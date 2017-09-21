function [pt_trial_match1, pt_trial_match2, pt_match_ind1, pt_match_ind2] = pointCloudTriangleMatch(pt_list1,pt_list2,r1,r2,total_trial)
% pt_list1 (2 or 3-by-n) is used to pick points randomly.
% pt_list2 (2 or 3-by-n) will be used to generate the mutual distance list.
% r1 specifies the shortest edge of the triangle.
% r2 specifies the longest edge of the triangle.
% total_trial specifies number of triangles that are picked randomly.
% pointCloudRegHan returns [pt_match_ind1, pt_match_ind2, s], in which
% pt_trial_match1 and 2 are two 1-by-total_trial cell arrays. Each cell 
% in 1 is a 3-by-1 array and each cell in 2 is a 3-by-n array, which
% records the index of the triangle vertices that matches the triangle in
% 1.
% The output can be used in the input of
% neuroReg.pointCloudTriangleCoplanarMatch(pt_list1,pt_list2,...
% pt_trial_match1,pt_trial_match2,total_trial)
%
% by Han Peng, penghan1992@gmail.com
% 2017.07.03
fprintf('\n****************************************************\n');
fprintf('*Point Cloud Registration. Matching Triangle Search*\n');
fprintf('****************************************************\n');

THETA_Z = 15; % normal is within 90+-15 degrees from z axis
THETA_Y = 30; % normal is within 90+-15 degrees from y axis


tic;
dist1 = squareform(pdist(pt_list1'));
dist2 = squareform(pdist(pt_list2'));
dist_legal = (dist1<r2).*(dist1>r1);
[i_p1,j_p1] = find(dist_legal);
L = length(i_p1);
toc;
i_count = 1;
i_all_loop=1;
di = [0 0 0];
pt_trial_match1=cell(1,total_trial);
pt_trial_match2=cell(1,total_trial);
trial_history = zeros(3,total_trial);
max_loop = 10*total_trial;
while i_count < total_trial+1
    i_all_loop = i_all_loop+1;
    if i_all_loop >= max_loop
        break;
    end
    %% Pick three points from set 1
    this_trial_match_list2 = [];
    ti_1 = randi(L,1); % pick three points randomly
    pt1 = i_p1(ti_1);
    pt2 = j_p1(ti_1);
    di(1) = norm(pt_list1(:,pt1) - pt_list1(:,pt2));
    k_p1 = find(dist_legal(pt1,:));
    k_p1 = k_p1(k_p1~=pt2);
    if length(k_p1)<1
        continue;
    end
    pt3 = k_p1(randi(length(k_p1),1));
    di(2) = norm(pt_list1(:,pt2) - pt_list1(:,pt3));
    di(3) = norm(pt_list1(:,pt3) - pt_list1(:,pt1));
    ddi(1) = abs(di(1)-di(2));
    ddi_thr(1) = max([5,0.1 * di(1)]) + max([0.5,0.1 * di(2)]);
    ddi(2) = abs(di(2)-di(3));
    ddi_thr(2) = max([5,0.1 * di(2)]) + max([0.5,0.1 * di(3)]);
    ddi(3) = abs(di(3)-di(1));
    ddi_thr(3) = max([5,0.1 * di(3)]) + max([0.5,0.1 * di(1)]);
    % Get a unique triangle
    if sum(di<=r1)>=1  || sum(di>=r2)>=1 || sum(ddi<ddi_thr)>=1
        continue;   
    end
    % Get rid of replicas
    pts = sort([pt1,pt2,pt3]);
    flag1 = trial_history(1,:) - pts(1);
    flag2 = trial_history(2,:) - pts(2);
    flag3 = trial_history(3,:) - pts(3);
    flags = (flag1==0).*(flag2==0).*(flag3==0);
    if sum(flags==1)>=1 
        continue;
    end
    % if everything is normal, then go on
    dist_match = {[],[],[]};
    dist_match_matix = {[],[],[]};
    %% Matching point set 2
    for i_tri = 1:3
        error_thres = max([5,0.1 * di(i_tri)]); % Threshold
        dist_match_matix{i_tri} = (abs(dist2 - di(i_tri))<error_thres);
        % Get the points with at leaset one neighbor with match distance
        dist_match{i_tri} = (sum(dist_match_matix{i_tri},2)>0);
    end
    dist_match_all = dist_match{1}+dist_match{2}+dist_match{3};
    % Find the point with at least two connections that matches the
    % triangle edges
    pt_two_matches = find(dist_match_all>=2);
    L_two_matches = length(pt_two_matches);
    if L_two_matches<3
        % In case there are too few matching points, jump this trial.
%         fprintf('%2d promising points. Less than 3. Initializing next trial.\n',L_two_matches);
        continue;
    end
    fprintf('Trial #%d\t ',i_count);
    % If there are more than 3 points, then test each of them.
%     fprintf('Testing %5d points...\t',L_two_matches);
    n_matches = 0;
    for k_matched_point = 1:L_two_matches
        % Searching the matched triangle
        this_point = pt_two_matches(k_matched_point);
        L_dist_ind = [0 0 0];
        dist_match_ind = {[],[],[]};
        for i_tri = 1:3
            dist_match_ind{i_tri} = find(dist_match_matix{i_tri}(this_point,:));
            L_dist_ind(i_tri) = length(dist_match_ind{i_tri});
        end
        i_test_tri = [1 2 3];
        % d1 d2 match, test d3
        for i = 1:3
            if (L_dist_ind(i_test_tri(1))>=1)&&(L_dist_ind(i_test_tri(2))>=1)
                %                 error_thres = max([10,0.1 * di(i_test_tri(3))]);
                %                 [ii,jj] = ndgrid(dist_match_ind{i_test_tri(1)},dist_match_ind{i_test_tri(2)});
                %                 i_pt_list = ii(:);
                %                 j_pt_list = jj(:);
                dist_d3_match = dist_match_matix{i_test_tri(3)}...
                    (dist_match_ind{i_test_tri(1)},dist_match_ind{i_test_tri(2)});
                [i_match_list_temp,j_match_list_temp] = find(dist_d3_match);
                i_match_list = dist_match_ind{i_test_tri(1)}(i_match_list_temp);
                j_match_list = dist_match_ind{i_test_tri(2)}(j_match_list_temp);
                %                 dist_d3_list = norm(pt_list2(:,i_pt_list) - pt_list2(:,j_pt_list));
                %                 dist_d3_match = find((dist_d3_list - di(i_test_tri(3)))<error_thres);
                if ~isempty(i_match_list)
                    %                     i_match_list = i_pt_list(dist_d3_match(:));
                    %                     j_match_list = j_pt_list(dist_d3_match(:));
                    this_point_list = this_point .* ones(1,length(i_match_list));
                    % THINK AGAIN HOW TO ARRANGE THE POINTS!!!
                    current_match_list2 = sort([i_match_list(:)';j_match_list(:)';this_point_list(:)'],1);
                    this_trial_match_list2 = cat(2,this_trial_match_list2,current_match_list2);
                    this_trial_match_list1 = [pt1;pt2;pt3];
                    n_matches = n_matches+length(i_match_list);
                end
            end
            i_test_tri = circshift(i_test_tri,1); % I love group theory!
        end
        
    end
    fprintf('%4d matching triangles. \t',n_matches);
    %% Verify the triangle surface normal
    pt1_position = pt_list2(:,this_trial_match_list2(1,:));
    pt2_position = pt_list2(:,this_trial_match_list2(2,:));
    pt3_position = pt_list2(:,this_trial_match_list2(3,:));
    v1_list = pt2_position - pt1_position;
    v2_list = pt3_position - pt1_position;
    normal_list = cross(v1_list,v2_list,1);
    normal_list = normal_list./sqrt(sum(normal_list.^2,1));
    cos_z = abs(dot(normal_list,repmat([0;0;1],[1,size(normal_list,2)]),1));
    cos_y = abs(dot(normal_list,repmat([0;1;0],[1,size(normal_list,2)]),1));
    cos_z_threshold = cosd(90-THETA_Z);
    cos_y_threshold = cosd(90-THETA_Y);
    ind_legal = (cos_z < cos_z_threshold).*(cos_y < cos_y_threshold);
    this_trial_match_list2 = this_trial_match_list2(:,ind_legal==1);
    %% Get rid of the repeat columns
    this_trial_match_list2 = unique(this_trial_match_list2','rows')'; 
    %% Rearrange each columns so that each points match each other.
    % d1 = |p1p2|, d2 = |p2p3|, d3 = |p3p1|
    % p1 - d1&d3, p2 - d1&d2, p3 - d2&d3
    di2pi_table = [0 2 1; 2 0 3; 1 3 0];
    
    pts = this_trial_match_list1;
    pts_new = [];
    d_1 = [];
    d_1(1,1) = norm(pt_list1(:,pts(1))-pt_list1(:,pts(2)));
    d_1(2,1) = norm(pt_list1(:,pts(2))-pt_list1(:,pts(3)));
    d_1(3,1) = norm(pt_list1(:,pts(3))-pt_list1(:,pts(1)));
    [~,I] = sort(d_1);
    pts_new(1,1) = pts(di2pi_table(I(1),I(3)));
    pts_new(2,1) = pts(di2pi_table(I(1),I(2)));
    pts_new(3,1) = pts(di2pi_table(I(2),I(3)));
    this_trial_match_list1 = pts_new;
    
    pts = this_trial_match_list2;
    pts_new = zeros(size(pts));
    d_2 = [];
    v1 = pt_list2(:,pts(1,:));
    v2 = pt_list2(:,pts(2,:));
    v3 = pt_list2(:,pts(3,:));
    d_2(1,:) = sum((v1 - v2).^2,1).^0.5;
    d_2(2,:) = sum((v2 - v3).^2,1).^0.5;
    d_2(3,:) = sum((v3 - v1).^2,1).^0.5;
    [~,I] = sort(d_2,1);
    ind_table1 = (I(1,:) + I(3,:)).^2;
    ind_table1(ind_table1==9)=2;
    ind_table1(ind_table1==16)=1;
    ind_table1(ind_table1==25)=3;
    ind_table2 = (I(1,:) + I(2,:)).^2;
    ind_table2(ind_table2==9)=2;
    ind_table2(ind_table2==16)=1;
    ind_table2(ind_table2==25)=3;
    ind_table3 = (I(2,:) + I(3,:)).^2;
    ind_table3(ind_table3==9)=2;
    ind_table3(ind_table3==16)=1;
    ind_table3(ind_table3==25)=3;
    for i = 1:size(pts_new,2)
        pts_new(:,i) = pts([ind_table1(i),ind_table2(i),ind_table3(i)],i);
    end
    this_trial_match_list2 = pts_new;
    
    % I hate matlab matrix dimension conventions. It is not homogeneus.
    n_matches = size(this_trial_match_list2,2);
    fprintf('%3d coplanar. \n',n_matches);
    %% Record the triangles
    pt_trial_match1{i_count} = this_trial_match_list1;
    pt_trial_match2{i_count} = this_trial_match_list2;
    trial_history(:,i_count) = sort([pt1;pt2;pt3]);
    i_count = i_count+1;
end

%% Output
% Build the cell array
pt_trial_match1 = pt_trial_match1(~cellfun('isempty',pt_trial_match1));
pt_trial_match2 = pt_trial_match2(~cellfun('isempty',pt_trial_match2));
% Build a matching list for point index
pt_match_ind1 = [];
pt_match_ind2 = [];
for i = 1:length(pt_trial_match1)
    pts_1 = pt_trial_match1{i};
    pts_2 = pt_trial_match2{i};    
    pts_1_rep = repmat(pts_1,[1,size(pts_2,2)]);    
    pt_match_ind1 = [pt_match_ind1,pts_1_rep];
    pt_match_ind2 = [pt_match_ind2,pts_2];
end

end

















