function [data_out,b_plane,b_list_volume] = cutVolume(data,cut_grid,M,d)
% cutVolume cut a slice from the data.
% data_out.x, data_out.y, data_out.value is the output slice.
% b_list is the boundary polygon (5-by-2) in the plane coordination
% b_list1 is the boundary polygon (5-by-3) in the volume coordination
% The slice is specified by the initial grid (cut_grid.x,cut_grid.y) which
% is transformed by M.
% M is the tranform from slice coordination to volume coordination.
% The initial plane is assumed to be perpendicular to y axis.

if nargin==3 % Integration disenable
    %% Rotate slice
    % Get the size of the targeting image (it will be cropped afterwards)
    nx = length(cut_grid.x);
    ny = length(cut_grid.y);
    % Create the grid
    [xx0,zz0] = ndgrid(cut_grid.x,cut_grid.y);
    yy0 = zeros(size(xx0));
    pts0 = [xx0(:),yy0(:),zz0(:),ones(size(yy0(:)))]';
    % Transformation
    pts1 = M*pts0;
    xx1 = pts1(1,:);
    yy1 = pts1(2,:);
    zz1 = pts1(3,:);
    % Interp and get the data
    data_out = cut_grid;
    v = interp3(data.y,data.x,data.z,data.value,yy1,xx1,zz1);
    data_out.value = reshape(v,[nx,ny]);
    % Crop the image
    [data_out,xb,zb] = cropNaN(data_out);
    xb_list = [xb(1),xb(1),xb(2),xb(2),xb(1)];
    yb_list = [0,0,0,0,0];
    zb_list = [zb(1),zb(2),zb(2),zb(1),zb(1)];
    b_plane = [xb_list;zb_list];
    b_list_volume = M * [xb_list;yb_list;zb_list;ones(size(xb_list))];
elseif nargin>=4 % Integration enable
    %% Gather information
    stepX = mean(diff(data.x));
    stepSlice = stepX;
    num_slice = 2*round(d/2/stepSlice)+1;
    offset = (1:num_slice)-round(d/2/stepSlice)-1;
    offset = offset*stepSlice; % plus and minus d/2
    nx = length(cut_grid.x);
    ny = length(cut_grid.y);
    %% Rotate slice
    % Get the size of the targeting image (it will be cropped afterwards)
    % Create the grid
    [xx0,zz0] = ndgrid(cut_grid.x,cut_grid.y);
    yy0 = zeros(size(xx0));
    pts0 = [xx0(:),yy0(:),zz0(:),ones(size(yy0(:)))]';
    % Transformation
    pts1 = M*pts0;
    n = M*[0;1;0;0];
    data_out = cut_grid;
    % Interp and get the data
    v_sum = zeros(nx,ny);
    pts1_rep = repmat(pts1,[1,1,num_slice]);
    offset = reshape(offset,[1,1,num_slice]);
    n1 = reshape(n,[3,1,1]);
    pt_offset = offset.*n1;
    pts2 = pts1_rep + pt_offset;
    pts2 = reshape(pts2,[3,nx*ny*num_slice]);
    v = interp3(data.y,data.x,data.z,data.value,pts2(2,:),pts2(1,:),pts2(3,:));
    % I like vectorization...
    v_mat = reshape(v,[nx,ny,num_slice]);
%     nf = squeeze(isnan(v_mat(:,:,round(d/2/stepSlice)+1))); % not an NaN
%     nf = repmat(nf,[1,1,num_slice]);
    for i = 1:num_slice 
        v_mat_this = v_mat(:,:,i);
        nf = isnan(v_mat_this); % nf: is NaN
        v_mat_this = v_mat_this - mean(v_mat_this(~nf));        
        if i~=round(d/2/stepSlice)+1
            % Middle slice is the mask with nan
            v_mat_this(nf)=0;
        end
        v_mat(:,:,i) = v_mat_this;
    end
%     for i = 1:num_slice
%         pt_offset = offset(i)*n;
%         pts2 = pts1 + pt_offset;
%         v = interp3(data.y,data.x,data.z,data.value,pts2(2,:),pts2(1,:),pts2(3,:));
%         v_mat = reshape(v,[nx,ny]);
%         nf = isnan(v_mat); % nf: is NaN
%         v_mat = v_mat - mean(v_mat(~nf));
%         if i~=round(d/2/stepSlice)+1
%             % Middle slice is the mask with nan
%             v_mat(~nf)=0;
%         end
%         v_sum = v_sum+v_mat;
%     end
%     v_sum = v_sum/num_slice;
    v_sum = sum(v_mat,3);
    data_out.value = v_sum;
    % Crop the image
    [data_out,xb,zb] = cropNaN(data_out);
    xb_list = [xb(1),xb(1),xb(2),xb(2),xb(1)];
    yb_list = [0,0,0,0,0];
    zb_list = [zb(1),zb(2),zb(2),zb(1),zb(1)];
    b_plane = [xb_list;zb_list];
    b_list_volume = M * [xb_list;yb_list;zb_list;ones(size(xb_list))];
end


end
function [data_out,xb,zb] = cropNaN(data)
v = ~isnan(data.value);
vx = sum(v,2);
vy = sum(v,1);
ix1 = find(vx,1,'first');
ix2 = find(vx,1,'last');
iy1 = find(vy,1,'first');
iy2 = find(vy,1,'last');
data_out.x = data.x(ix1:ix2);
data_out.y = data.y(iy1:iy2);
data_out.value = data.value(ix1:ix2,iy1:iy2);
xb = data.x([ix1,ix2]);
zb = data.y([iy1,iy2]);
end
