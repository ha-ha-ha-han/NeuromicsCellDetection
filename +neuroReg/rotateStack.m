function data_out = rotateStack(data,alpha,beta,gamma)
% rotateStack is used to rotate a ZStack data 
% data.value)
% data_out = rotateStack(data,alpha,beta,gamma)
% INPUT:
% data: ZStack data. data.x, data.y, data.z and data.value.
% alpha: rotation around y-axis
% beta: rotation around z-axis
% gamma: rotation around x-axis
% OUTPUT
% data_out: ZStack data. Rotated.
% WARNING:
% It is slow.
% Try to use rotateCells and renderCell3

%% Make rotation matrix
Ry = [cosd(alpha),0,sind(alpha);0,1,0;-sind(alpha),0,cosd(alpha)];
Rz = [cosd(beta),-sind(beta),0;sind(beta),cosd(beta),0;0,0,1];
Rx = [1,0,0;0,cosd(gamma),-sind(gamma);0,sind(gamma),cosd(gamma)];
R = Ry*Rz*Rx;
[xxx,yyy,zzz] = ndgrid(data.x,data.y,data.z);
%% Rotate
p = [xxx(:),yyy(:),zzz(:)]';
v = double(data.value(:));
p1 = R*p;
%% Interpolation
stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
stepZ = mean(diff(data.z));
x1 = min(p1(1,:));
x2 = max(p1(1,:));
y1 = min(p1(2,:));
y2 = max(p1(2,:));
z1 = min(p1(3,:));
z2 = max(p1(3,:));
data_out.x = x1:stepX:x2;
data_out.y = y1:stepY:y2;
data_out.z = z1:stepZ:z2;
[xxx_new,yyy_new,zzz_new] = ndgrid(data_out.x,data_out.y,data_out.z);
tic
data_out.value = single(griddata(p1(1,:),p1(2,:),p1(3,:),v,xxx_new,yyy_new,zzz_new,'nearest'));
toc

end

