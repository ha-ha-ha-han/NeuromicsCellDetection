function M_1 = registerSixPoints(pts_1,pts_2)
% This function calculate the M that let the equation 
% point_set_2 = M^(-1) * point_set_1 best holds.
% x and y should be 3-by-6 arrays.
% in the neuroReg case, x(2,:) is assumed as zero.
% Because it is 1am and this is the 4th day that I work until 1am during the
% week :(
% 
% Han

M=[];
M_1 = [];
%% Match center of the mass
pts_2_center = mean(pts_2,2);
pts_1_center = mean(pts_1,2);
pts_1c = pts_1 - pts_1_center;
pts_2c = pts_2 - pts_2_center;
% figure(10099);hold off;
% scatter3(pts_1c(1,:),pts_1c(2,:),pts_1c(3,:));
% hold on;
% scatter3(pts_2c(1,:),pts_2c(2,:),pts_2c(3,:));
% hold off;
% xlabel x;
% ylabel y;
% zlabel z;
% axis equal;axis vis3d;
%% Get the average surface normal for pts_1
% pts_2 surface normal is [0,1,0]
% It is actually a linear regression
dist_plane = @(angles)dist2plane(angles,pts_1c);
angles = fminsearch(dist_plane,[0;0]);
d2 = dist2plane(angles,pts_1c);
sn(3,1) = cosd(angles(1));
sn(1,1) = sind(angles(1))*cosd(angles(2));
sn(2,1) = sind(angles(1))*sind(angles(2));

axis1 = cross(sn,[0;1;0]);
axis1 = axis1./norm(axis1);
axis3 = cross(axis1,sn);

RotMat_1 = [axis1';sn';axis3'];
pts_1r_step1 = RotMat_1*pts_1c;
%% Visualization for testing
% scatter3(pts_1r_step1(1,:),pts_1r_step1(2,:),pts_1r_step1(3,:));
% hold on;
% plot3([pts_1r_step1(1,1:3),pts_1r_step1(1,1)],[pts_1r_step1(2,1:3),pts_1r_step1(2,1)],[pts_1r_step1(3,1:3),pts_1r_step1(3,1)]);
% plot3([pts_1r_step1(1,4:6),pts_1r_step1(1,4)],[pts_1r_step1(2,4:6),pts_1r_step1(2,4)],[pts_1r_step1(3,4:6),pts_1r_step1(3,4)]);
% scatter3(pts_2c(1,:),pts_2c(2,:),pts_2c(3,:));
% 
% hold off;
% axis equal;axis vis3d;
% xlabel x;
% ylabel y;
% zlabel z;

%% Final rotation

dist_rotation = @(angle)distRotation(angle,pts_1r_step1,pts_2c);
rotation_angle = fminsearch(dist_rotation,0);
RotMat_2 = [cosd(rotation_angle),0,sind(rotation_angle);0,1,0;-sind(rotation_angle),0,cosd(rotation_angle)];
pts_1r_step2 = RotMat_2*pts_1r_step1;
pts_1r_step2_center = mean(RotMat_2*RotMat_1*pts_1,2);
pts_1r_step2_translate = pts_1r_step2 + pts_2_center;
% Note: for a Four-Element Number, the translation is after rotation.
M = RotMat_2*RotMat_1;
M_1 = M';
pts_2r = M'*pts_2;
pts_2r_center = mean(pts_2r,2);
M_1 = cat(2,M_1,-pts_2r_center+pts_1_center);

pts_2_translate = M_1*cat(1,pts_2,ones(size(pts_2(1,:))));
dist10 = norm(pts_1(:,1:3)-pts_1(:,4:6));
dist20 = norm(pts_2(:,1:3)-pts_2(:,4:6));
dist21 = norm(pts_2_translate(:,1:3)-pts_2_translate(:,4:6));
% fprintf('dist10 = %6.1f, dist20 = %6.1f, dist21 = %6.1f\n',dist10,dist20,dist21);
%% Visualization for testing
% figure(10099);hold off;
% scatter3(pts_1(1,:),pts_1(2,:),pts_1(3,:),'b');
% hold on;
% plot3([pts_1(1,1:3),pts_1(1,1)],[pts_1(2,1:3),pts_1(2,1)],[pts_1(3,1:3),pts_1(3,1)],'b');
% plot3([pts_1(1,4:6),pts_1(1,4)],[pts_1(2,4:6),pts_1(2,4)],[pts_1(3,4:6),pts_1(3,4)],'b');
% % Test M
% scatter3(pts_2_translate(1,:),pts_2_translate(2,:),pts_2_translate(3,:),'d');
% hold on;

% scatter3(pts_2(1,:),pts_2(2,:),pts_2(3,:),'r');
% plot3([pts_2_translate(1,1:3),pts_2_translate(1,1)],[pts_2_translate(2,1:3),pts_2_translate(2,1)],[pts_2_translate(3,1:3),pts_2_translate(3,1)],'r');
% plot3([pts_2_translate(1,4:6),pts_2_translate(1,4)],[pts_2_translate(2,4:6),pts_2_translate(2,4)],[pts_2_translate(3,4:6),pts_2_translate(3,4)],'r');
% % original distance dist0
% dist0 = norm(pts_1(:,1:3)-pts_1(:,4:6));
% % After re-center
% distC = norm(pts_1c(:,1:3)-pts_1c(:,4:6));
% %dist1 step1
% dist1_step1 = norm(pts_1r_step1(:,1:3)-pts_1r_step1(:,4:6));
% %dist1 step2
% dist1_step2 = norm(pts_1r_step2_translate(:,1:3)-pts_1r_step2_translate(:,4:6));
% %dist2
% dist2 = norm(pts_2(:,1:3)-pts_2(:,4:6));
% %Title dist1 dist2 difference
% diff_dist = dist2 - dist1_step1;
% 
% title(['dist1=',num2str(dist1_step1),' dist2=',num2str(dist2),' diff=',num2str(diff_dist)]);
% fprintf([...
%     'dist0=%6.1f, distC=%6.1f, dist1_step1=%6.1f,\n',...
%     'dist1_step2=%6.1f, dist2=%6.1f\n',...
%     'diff_12=%6.1f\n'], ...
%     dist0,distC,dist1_step1,dist1_step2,dist2,diff_dist);
% hold off;
% axis equal;axis vis3d;
% xlabel x;
% ylabel y;
% zlabel z;
% End of test

%% Do the point cloud registration
% ...............
%% Use this matrix to test all the data set. record the best 20 matches
%% Visualize the match
%% Fix coplanar and distance test!!!!


end

function dist = dist2plane(angles,pts)
% pts should be 3-by-n
% sn is 3-by-1
sn(3) = cosd(angles(1));
sn(1) = sind(angles(1))*cosd(angles(2));
sn(2) = sind(angles(1))*sind(angles(2));
x = pts(1,:);
y = pts(2,:);
z = pts(3,:);

t = x*sn(1)+y*sn(2)+z*sn(3);
dist = (sum(t.^2)).^0.5/6;
end

function dist = distRotation(angle,pts_1,pts_2)

RotMat = [cosd(angle),0,sind(angle);0,1,0;-sind(angle),0,cosd(angle)];

pts_1r = RotMat*pts_1;
diff_pts = pts_1r - pts_2;
dist = sum(diff_pts(:).^2)^0.5/6;



end

