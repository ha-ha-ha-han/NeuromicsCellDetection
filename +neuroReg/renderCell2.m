function data = renderCell2(x_lim,y_lim,step,pt_list,pt_area,Option)
% xlim and ylim are of the form [lb ub] to specify the render range.
% step specifies the resolution.
% pt_list is a 2-by-n matrix record the cell position.
% pt_area is a 1-by-n array specifies each cell's area
% data is returned with data.x, data.y, data.value
% size(data.value) == [length(data.x),length(data.y)]
% If you use a for-loop, it is (nearly) equivalent to renderCell3

if isfield(Option,'FSize')
    fSize = Option.FSize;
else
    fSize = 80;
end
if isfield(Option,'CellRadius')
    cellR = Option.CellRadius;
else
    cellR = 10;
end
data.x = x_lim(1):step:x_lim(2);
data.y = y_lim(1):step:y_lim(2);
nx = length(data.x);
ny = length(data.y);
data.value = zeros(nx,ny);
v_b = data.value;
if isempty(pt_area)||isempty(pt_list)
   return; 
end


fpx = fSize/step;
rpx = floor(fpx/2);
cpx = cellR/step;
f = fspecial('gaussian',[2*rpx+1,2*rpx+1],cpx);
rpx_b = rpx * 2;
f_b = ones([2*rpx_b+1,2*rpx_b+1]);
pt_area = pt_area/max(pt_area);
pt_list = (pt_list-[x_lim(1);y_lim(1)])/step;
% Now all the units are pixel
for i = 1:length(pt_area)
    % Boundary check
    pc = round(pt_list(:,i)');
    if pc(1)>1 && pc(1)<nx && pc(2)>1 && pc(2)<ny
        %% Cell
        % The center should be inside the region of interet
        x1 = pc(1) - rpx;
        x2 = pc(1) + rpx;
        y1 = pc(2) - rpx;
        y2 = pc(2) + rpx;
        xf1 = 1;
        xf2 = 2*rpx+1;
        yf1 = 1;
        yf2 = 2*rpx+1;
        % left
        if x1 < 1
            x1 = 1;
            xf1 = 2+rpx-pc(1);
        end
        % right
        if x2 > nx
            x2 = nx;
            xf2 = rpx+1-pc(1)+nx;
        end
        % down
        if y1 < 1
            y1 = 1;
            yf1 = 2+rpx-pc(2);
        end
        % up
        if y2 > ny
            y2 = ny;
            yf2 = rpx+1-pc(2)+ny;
        end
        % Render
        data.value(x1:x2,y1:y2) = data.value(x1:x2,y1:y2) + ...
            f(xf1:xf2,yf1:yf2)*pt_area(i);
        %% Background
                % The center should be inside the region of interet
        x1 = pc(1) - rpx_b;
        x2 = pc(1) + rpx_b;
        y1 = pc(2) - rpx_b;
        y2 = pc(2) + rpx_b;
        xf1 = 1;
        xf2 = 2*rpx_b+1;
        yf1 = 1;
        yf2 = 2*rpx_b+1;
        % left
        if x1 < 1
            x1 = 1;
            xf1 = 2+rpx_b-pc(1);
        end
        % right
        if x2 > nx
            x2 = nx;
            xf2 = rpx_b+1-pc(1)+nx;
        end
        % down
        if y1 < 1
            y1 = 1;
            yf1 = 2+rpx_b-pc(2);
        end
        % up
        if y2 > ny
            y2 = ny;
            yf2 = rpx_b+1-pc(2)+ny;
        end
        % Render
        v_b(x1:x2,y1:y2) = data.value(x1:x2,y1:y2) + ...
            f_b(xf1:xf2,yf1:yf2)*pt_area(i);
    end
end
%% finalize
data.value(v_b==0)=nan;       
% Visualizae
% figure(985211);
% plotData2(data);


end

