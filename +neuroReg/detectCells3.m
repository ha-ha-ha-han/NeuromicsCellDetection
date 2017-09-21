function [pt_list,pt_area] = detectCells3(data,Option)
% detectCells3 detects cell positions from a ZStack data.
% [pt_list,pt_area] = detectCells3(data,Option)
% ---------
% OUTPUT:
% pt_list: positions of the detected cells (3-by-N array)
% pt_area: volume of each detected cells (1-by-N array)
% ---------
% INPUT
% data: data.x, data.y, data.z, data.value (see doc neuroReg)
% Option: 
%   Note: pass [] to the function to use default settings.
%   -----------------------------
%   Options and default settings:
%   Option.Detect3Mode = 'Red';
%       >> Specify the detection mode. Can be 'Red' or 'Green'
%   Option.Sigma = [1 1 1]*4;
%       >> Size for difference of Gaussian filter. Unit: um
%   Option.Res0 = 0.0035;   
%       >> Threshold for difference of Gaussian filter
%   Option.SizeLimit = [100,4500];    
%       >> Estimization of the cell volume. Unite: um^3
%   Option.Sensitivity = 1.3;    
%       >> Used in Green Mode. The larger this number is, the less cells 
%          will be detected.
%   Option.MedianFilterSize = [7,7,7];   
%       >> Size for median filter
%   Option.NeighborSize = [30 30 50];
%       >> Neighbor size in adaptive thresholding. Green Mode only. Unit:
%       um.
    

%% Preset options
if ~isfield(Option,'Detect3Mode')
    % Specify the detection mode
    Option.Detect3Mode = 'Red';
end
if ~isfield(Option,'Sigma')  
    % Size for difference of Gaussian filter. Unit: um
    Option.Sigma = [1 1 1]*4;
end
if ~isfield(Option,'Res0')
    % Threshold for difference of Gaussian filter
    Option.Res0 = 0.0035;
end
if ~isfield(Option,'SizeLimit')
    % Estimization of the cell volume. Unite: um^3
    Option.SizeLimit = [100,4500];
end
if ~isfield(Option,'Sensitivity')
    % Used in Green Mode. The larger this number is, the less cells will be
    % detected.
    Option.Sensitivity = 1.3;
end
if ~isfield(Option,'MedianFilterSize')
    % Size for median filter
    Option.MedianFilterSize = [7,7,7];
end
if ~isfield(Option,'NeighborSize')
    % Used in Green Mode. Neighbor size in adaptive thresholding. 
    Option.NeighborSize = [50 50 50];
end
disp(Option);



%% Load data & get slice
% data = evalin('base','dataZ');
%% Preprocessing
% Use difference of Gaussian detector
Res0 = Option.Res0;
Sensitivity = Option.Sensitivity;
SizeLimit = Option.SizeLimit;
Sigma = Option.Sigma;


stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
stepZ = mean(diff(data.z));
SizeLimitPx = SizeLimit/stepX/stepY/stepZ;
MedianFilterSizePx = round(Option.MedianFilterSize./[stepX,stepY,stepZ]/2)*2+1;
v = data.value / max(data.value(:));
%% Green Channel Mode
if strcmp(Option.Detect3Mode,'Green')
    NbSize = round(Option.NeighborSize./[stepX,stepY,stepZ]);
    % Normalization
    v = (data.value - min(data.value(:))) / (max(data.value(:))-min(data.value(:)));
    tic
    % Local thresholding. Time consuming (around 200s).
    disp('Local thresholding... May take about 200 seconds.');
    v_local_thres = Sensitivity * convn(v,1/(NbSize(1)*NbSize(2)*NbSize(3))*ones(NbSize),'same');
    
    v = v - v_local_thres;
%     v = medfilt3(v,MedianFilterSizePx);
    %     data1 = data;
    %     data1.value = v;
    vf = v>0;
    v(~vf) = 0;
    data1 = data;
    data1.value = v;
    %     data3 = data;
    %     data3.value = v_local_thres;
    v = v./v_local_thres;
    v = (v - min(v(:)))/(max(v(:))-min(v(:)));
%     data4 = data;
%     data4.value = v;
    assignin('base','data1',data1)
    %     assignin('base','data2',data2)
    %     assignin('base','data3',data3)
    %     assignin('base','data4',data4)
    disp('Local threshold calculated.');
    toc
end

%% Gaussian
% um2vx
sigmaInd = [1 1 1];
sigmaInd(2) = Sigma(1)/stepX; % in imgaussfilt, sigmaInd(2) is for 1st dimension
sigmaInd(1) = Sigma(2)/stepY; % sigmaInd(1) is for 2nd dimension
sigmaInd(3) = Sigma(3)/stepZ;
tic
fprintf('Applying difference of gaussian detector...\n');
Ic = imgaussfilt3(v,sigmaInd);
Is = imgaussfilt3(v,sigmaInd*5);
res = Ic - Is;
fprintf('Difference of gaussian detector done...\n');
toc
% illumination correction
%  Res = adaptthresh3(res,Sensitivity);
%  im_bw = single((res-Res)>0);

tic
fprintf('Finalizing result...\n');
im_bw = single(res>Res0);
im_bw = imclearborder(im_bw);
CC = bwconncomp(im_bw,26);
stats = regionprops(CC,'Area','Centroid','FilledImage','BoundingBox');
pt_list = [nan nan nan]';
pt_area = [];
k=1;
for i = 1:length(stats)
    if stats(i).Area>SizeLimitPx(1)
        if stats(i).Area<SizeLimitPx(2)
            pt_list(:,k) = stats(i).Centroid';
            pt_area(k) = stats(i).Area;
            k=k+1;
        else
            bw = stats(i).FilledImage;
            [pt_this,pt_area_this] = splitCells(bw,stepX);
            n = length(pt_area_this);
            if n>=1
                for j = 1:n
                    pt_list(:,k) = pt_this(:,j) + stats(i).BoundingBox(1:3)';
                    pt_area(k) = pt_area_this(j);
                    k = k+1;
                end
            end
        end
        
    end
end
pt_area = pt_area*(stepX*stepY*stepZ);
if strcmp(Option.Detect3Mode,'Green')
    % If it Green mode, the area will not be accurate.
    % To fix this, a unified cell volume is assigned.
    pt_area = ones(size(pt_area))*mean(Option.SizeLimit);
end
temp = pt_list(1,:);
pt_list(1,:) = pt_list(2,:);
pt_list(2,:) = temp;
pt_list(1,:) = pt_list(1,:)*stepX + min(data.x) - stepX;
pt_list(2,:) = pt_list(2,:) * stepY + min(data.y) - stepY;
pt_list(3,:) = pt_list(3,:) * stepZ + min(data.z) - stepZ;
fprintf('Result done.\n Number of Cells = %d\n',length(pt_area));
toc
%% Visualized feedback
% ^^^^^^^^^^^^^^^^ TO CONTINUE
% Sub plot detected cells. Label axis. Set proper vis angle
% Sub plot origin data.
% Enable Adjustment.
% figure(10086);
% subplot(1,3,1);
% scatter3(pt_list(1,:),pt_list(2,:),pt_list(3,:));
% title(['detected cells, number = ',num2str(length(pt_list(1,:)))]);
% xlabel x;
% ylabel y;
% zlabel z;
% axis equal
% axis vis3d
% view([0,90])
% subplot(1,3,2);
% data_res = data;
% data_res.value = squeeze(sum(im_bw,3));
% plotSlice(data_res);
% title('binarized data');
% xlabel x;
% ylabel y;
% subplot(1,3,3);
% data_z = data;
% data_z.value = squeeze(sum(data.value,3));
% plotSlice(data_z);
% title('Original data');
% xlabel x;
% ylabel y;
assignin('base','data_temp_749',data);
neuroReg.PlotSlices('data_temp_749','y',pt_list,pt_area,25);
end
function [pt,pt_area] = splitCells(bw,step)
min_dist = 15/step;
D = bwdist(~bw);
D = -D;
D(~bw) = Inf;
L = watershed(D);
L(~bw) = 0;
stats = regionprops(L,'Area','Centroid');
n = length(stats);
SizeLimitPx = 0.8*sum(bw(:))/n;
pt_area = [];
pt = [nan,nan,nan]';
k=1;
for i = 1:length(stats)
    if stats(i).Area>SizeLimitPx % Leave out small parts
        pt(:,k) = stats(i).Centroid';
        pt_area(k) = stats(i).Area;
        k = k+1;
    end    
end
if k==3 % Check the case when one cell is falsely identified as two.
    d = norm(pt(:,1)-pt(:,2));
    if d < min_dist
        clear pt;
        clear pt_area;
        CC = bwconncomp(bw);
        stats = regionprops(CC,'Area','Centroid');
        pt(:,1) = stats(1).Centroid';
        pt_area(1) = stats(1).Area;
        k=2;
        
        % Visualization for testing
        figure(10099),cla, isosurface(bw,0.5), axis equal, title('BW')
        xlabel x, ylabel y, zlabel z
        view(3), camlight, lighting gouraud
        figure(10100),cla
        isosurface(L==1,0.5)
        isosurface(L==2,0.5)
        axis equal
        title(['Segmented objects, ',num2str(k-1)])
        xlabel x, ylabel y, zlabel z
        view(3), camlight, lighting gouraud
%         test ended
    end
end
% Visualization for testing
% figure(10099),cla, isosurface(bw,0.5), axis equal, title('BW')
% xlabel x, ylabel y, zlabel z
% view(3), camlight, lighting gouraud
% figure(10100),cla
% isosurface(L==1,0.5)
% isosurface(L==2,0.5)
% axis equal
% title(['Segmented objects, ',num2str(k-1)])
% xlabel x, ylabel y, zlabel z
% view(3), camlight, lighting gouraud
% test ended
end
