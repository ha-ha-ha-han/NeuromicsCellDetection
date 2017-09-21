function getROI
%% Choose three points
% p1 = (x1,y1) and p2 = (x2, y2) is the upper boarder of the rectangle
% p3 = (x3,y3) specify the height of the rectangle
nAr = 5;
data = evalin('base','dataS_gen');
disp('choose five regions in the current figure...');
disp('eave region is pecified by three points');
[x,y] = getpts;
if length(x)~=15
    disp('choose three points!');
    return;
end
disp('working...');
pts_all = [x,y];
for i = 1:nAr
    pts = pts_all((i-1)*3+1:i*3,:);
    v = pts(2,:) - pts(1,:);
    vd = pts(3,:) - pts(1,:);
    a = dot(vd,v)/norm(v);
    u = vd - a*v/norm(v);
    pts_roi = zeros(4,2);
    pts_roi(1,:) = pts(1,:) - 100*u/norm(u);
    pts_roi(2,:) = pts(2,:) - 100*u/norm(u);
    pts_roi(3,:) = pts(2,:) + u;
    pts_roi(4,:) = pts(1,:) + u;
    pts_roi(5,:) = pts_roi(1,:);
    hold on;
    plot(pts_roi(:,1),pts_roi(:,2),'r');
    hold off;
    %% Segment image
    [xx,yy] = ndgrid(data.x,data.z);
    vv = squeeze(data.value(:,4,:));
    vv = vv(:);
    in = inpolygon(xx(:),yy(:),pts_roi(:,1),pts_roi(:,2));
    xx_in = xx(in) - pts_roi(4,1); % Set lower left conner as the origin
    yy_in = yy(in) - pts_roi(4,2); % Set lower left conner as the origin
    p = [xx_in,yy_in];
    vv_in = double(vv(in));
    
    %% Coordination translation
    x_range = norm(v);
    y_range = norm(u)+100;
    M = [v/x_range;-u/y_range];
    M1 = M^-1;
    p1 = p*M1;
    xx_new = p1(:,1);
    yy_new = p1(:,2);
    data_out{i}.x = 0:1:x_range;
    data_out{i}.y = 0:1:y_range;
    [yy_mesh,xx_mesh] = meshgrid(data_out{i}.y,data_out{i}.x);
    data_out{i}.value = single(griddata(xx_new,yy_new,vv_in,xx_mesh,yy_mesh));
end
assignin('base','data_roi',data_out);
disp('done...');
end

