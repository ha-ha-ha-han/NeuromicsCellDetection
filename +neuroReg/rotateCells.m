function [pt_list_rotated,R,x_lim,y_lim,d_lim] = rotateCells(pt_list,alpha,beta,gamma)
% rotateCells is very useful for:
% 1. Rotate pt_list
% >>[pt_list_rotated,R,x_lim,y_lim,d_lim] = rotateCells(pt_list,alpha,beta,gamma)
% 2. Make a rotation matrix from alpha beta and gamma
% >> [~,R,~,~,~] = rotateCells(~,alpha,beta,gamma)
% INPUT:
% pt_list is a point coordination list. 3-by-N
% alpha: rotation around y-axis
% beta: rotation around z-axis
% gamma: rotation around x-axis
% OUTPUT
% pt_list_rotated: rotated pt_list (3-by-N)
% R: roation matrix
% x_lim y_lim and d_lim: boundaries of the rotated pt_list.


%% Make rotation matrix
Ry = [cosd(alpha),0,sind(alpha);0,1,0;-sind(alpha),0,cosd(alpha)];
Rz = [cosd(beta),-sind(beta),0;sind(beta),cosd(beta),0;0,0,1];
Rx = [1,0,0;0,cosd(gamma),-sind(gamma);0,sind(gamma),cosd(gamma)];
R = Rx*Rz*Ry;
%% Rotate
pt_list_rotated = R*pt_list;

x_lim(1) = min(pt_list_rotated(1,:));
x_lim(2) = max(pt_list_rotated(1,:));
y_lim(1) = min(pt_list_rotated(3,:));
y_lim(2) = max(pt_list_rotated(3,:));
d_lim(1) = min(pt_list_rotated(2,:));
d_lim(2) = max(pt_list_rotated(2,:));
end