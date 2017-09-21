function [pt_match_ind1,pt_match_ind2] = pointCloudTriangleCoplanarMatch(...
    pt_list1,pt_list2,pt_trial_match1,pt_trial_match2,Threshold,total_trial)
% pt_list1 (2 or 3-by-n) is used to pick points randomly.
% pt_list2 (2 or 3-by-n) will be used to generate the mutual distance list.
% pt_trial_match1 and 2 are two 1-by-total_trial cell arrays. Each cell
% in 1 is a 3-by-1 array and each cell in 2 is a 3-by-n array, which
% records the index of the triangle vertices that matches the triangle in
% 1.
% pt_trial_match1 and 2 can be the output from neuroReg.pointCloudTriangleMatch
% There should be no replica in pt_trial_match1.
% total_trial specifies number of triangles that are picked randomly.
% Threshold is the maximum mean distance that the points of one triangle
% are away from the plane, which is given by the other triangle.
%
% The out put are two 6-by-n array specifies the matching vertices.
%
% by Han Peng, penghan1992@gmail.com
% 2017.07.03
fprintf('\n****************************************************\n');
fprintf('*Point Cloud Registration. Coplanar Triangle Search*\n');
fprintf('****************************************************\n');

pt_match_ind1 = [];
pt_match_ind2 = [];

L = length(pt_trial_match1);
i_count = 1;
i_all_loop = 1;
[slect_mat_A,slect_mat_B] = ndgrid(1:L,1:L);
I = 1 - diag(ones(1,L));
slect_mat_A = slect_mat_A.*I;
slect_mat_B = slect_mat_B.*I;
max_loop = 10*total_trial;

while i_count < total_trial+1
    i_all_loop = i_all_loop+1;
    if i_all_loop >= max_loop
        break;
    end
    %% Select two triangles A and B from the point set 1 ramdomly.
    selected_ind = randi(L*L,1);
    selected(1) = slect_mat_A(selected_ind);
    selected(2) = slect_mat_B(selected_ind);
    if selected(1)==0
        continue;
    end
    fprintf('Trial #%d\t ',i_count);
    slect_mat_A(selected_ind)=0;
    slect_mat_B(selected_ind)=0;
    %% Get all the matching triangles from the point set 2 for A and B.
    ind_tri_A = pt_trial_match2{selected(1)};
    ind_tri_B = pt_trial_match2{selected(2)};
    pt_triA_1 = pt_list2(:,ind_tri_A(1,:));
    pt_triA_2 = pt_list2(:,ind_tri_A(2,:));
    pt_triA_3 = pt_list2(:,ind_tri_A(3,:));
    pt_triB_1 = pt_list2(:,ind_tri_B(1,:));
    pt_triB_2 = pt_list2(:,ind_tri_B(2,:));
    pt_triB_3 = pt_list2(:,ind_tri_B(3,:));
    %% Calculate all the surface normals for triangles A
    v1 = pt_triA_2 - pt_triA_1;
    v2 = pt_triA_3 - pt_triA_1;
    normal = cross(v1,v2,1);
    normal = normal./sum(normal.^2,1).^0.5;
    %% Calculate the distances from B triangles to the plane of A
    % triangles
    dist_mat_1 = normal'*pt_triB_1;
    dist_mat_2 = normal'*pt_triB_2;
    dist_mat_3 = normal'*pt_triB_3;
    dist_mean = (dist_mat_1.^2+dist_mat_2.^2+dist_mat_3.^2).^0.5/3;
    % The coplanar triangles (one from A, the other from B) are recorded
    flag_coplanar = (dist_mean<Threshold);
    %% The distance between the two triangles should match the distance
    % from A.
    pos_mean_A_2 = (pt_triA_1 + pt_triA_2 + pt_triA_3)/3;
    pos_mean_B_2 = (pt_triB_1 + pt_triB_2 + pt_triB_3)/3;
    dist_AB_2 = pdist2(pos_mean_A_2',pos_mean_B_2');
    pos_mean_A_1 = sum(pt_list1(:,pt_trial_match1{selected(1)}),2)/3;
    pos_mean_B_1 = sum(pt_list1(:,pt_trial_match1{selected(2)}),2)/3;
    dist_AB_1 = norm(pos_mean_A_1 - pos_mean_B_1);
    flag_dist = abs(dist_AB_2 - dist_AB_1)<dist_AB_1*0.2;
    [match_i,match_j] = find(flag_coplanar.*flag_dist);
    this_match_1 = [pt_trial_match1{selected(1)};pt_trial_match1{selected(2)}];
    fprintf('%3d matches\n',length(match_i));
    for i = 1:length(match_i)
        this_match_2 = [pt_trial_match2{selected(1)}(:,match_i(i));...
            pt_trial_match2{selected(2)}(:,match_j(i))];
        pt_match_ind1 = cat(2,pt_match_ind1,this_match_1);
        pt_match_ind2 = cat(2,pt_match_ind2,this_match_2);
        
        %% Visualization for testing
%         figure(10099);hold off;
        pts_1r = pt_list1(:,pt_match_ind1);
        v1 = pts_1r(1,:);
        v3 = pts_1r(2,:);
        v2 = zeros(size(v1));
        pts_1r = [v1;v2;v3];
        pts_2c = pt_list2(:,pt_match_ind2);
        scatter3(pts_1r(1,:),pts_1r(2,:),pts_1r(3,:),'b');
        hold on;
        plot3([pts_1r(1,1:3),pts_1r(1,1)],[pts_1r(2,1:3),pts_1r(2,1)],[pts_1r(3,1:3),pts_1r(3,1)],'b');
        plot3([pts_1r(1,4:6),pts_1r(1,4)],[pts_1r(2,4:6),pts_1r(2,4)],[pts_1r(3,4:6),pts_1r(3,4)],'b');
        scatter3(pts_2c(1,:),pts_2c(2,:),pts_2c(3,:),'r');
        plot3([pts_2c(1,1:3),pts_2c(1,1)],[pts_2c(2,1:3),pts_2c(2,1)],[pts_2c(3,1:3),pts_2c(3,1)],'r');
        plot3([pts_2c(1,4:6),pts_2c(1,4)],[pts_2c(2,4:6),pts_2c(2,4)],[pts_2c(3,4:6),pts_2c(3,4)],'r');
        hold off;
        axis equal;axis vis3d;
        xlabel x;
        ylabel y;
        zlabel z;
        % End of test
    end
    i_count = i_count+1;
end
end
















