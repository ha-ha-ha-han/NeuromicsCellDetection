function [pt_list,pt_area] = detectCells2(data,Option,ex_list)
% detectCells2 detects cells from a slice data.
% [data_out,pt_list] = detectCells2(data,Option,ex_list)
% INPUT:
% data: a slice data with data.x, data.y, data.value (2D only).
% For 3D data, use bwCell3.
% Option: see neuroReg.setOption
% ex_list: points to exclude. Generated in the preprocess stage. 2-by-M
% OUTPUT:
% pt_list: detected cells. It should be the same with the result from
% bwCells2 with the same input. 2-by-N
% pt_area: detected cells' area

if nargin==2
    exclude_flag=0;
elseif nargin==3&&isempty(ex_list)
    exclude_flag=0;
elseif nargin==3&&~isempty(ex_list)
    exclude_flag=1;
else    
    error(['Check input.','neuroReg.detectCells2(data,Option)',...
        ' or neuroReg.detectCells2(data,Option, ex_list)']);
end
%% Preset options
if ~isfield(Option,'Sigma')    
    Option.Sigma = [1 1]*5;
end
if ~isfield(Option,'Res0')
    Option.Res0 = 0.03;
end
if ~isfield(Option,'SizeLimit')
    Option.SizeLimit = [100,2000];
end
if ~isfield(Option,'MedianFilterSize')
    Option.MedianFilterSize = [15,15];
end
if ~isfield(Option,'Threshold')
    Option.Threshold = 0.5;
end





%% Get slice
stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
Threshold = Option.Threshold;
v_thres = max(data.value(:))*Threshold;
SizeLimit = Option.SizeLimit/(stepX*stepY);
Res0 = Option.Res0;
%% Thresholding Median filter
data.value(data.value<v_thres)=0;
MedianFilterSize = Option.MedianFilterSize;
MedianFilterSizePx = ceil(MedianFilterSize./[stepX,stepY]);
MedianFilterSizePx = floor(MedianFilterSizePx/2)*2+1;
v = data.value / max(data.value(:));
v = medfilt2(v,MedianFilterSizePx);
%% Use difference of Gaussian detector
sigmaInd = [1 1];
sigmaInd(1) = Option.Sigma(1)/stepY; % in imgaussfilt, sigma(1) is for 2nd dimension
sigmaInd(2) = Option.Sigma(1)/stepX; % sigma(2) is for 1st dimension
tic;
fprintf('Applying difference of gaussian detector...\n');
Ic = imgaussfilt(v,sigmaInd);
Is = imgaussfilt(v,sigmaInd*5);
res = Ic - Is;
im_bw = imbinarize(res,Res0)*1.0;
fprintf('Difference of gaussian detector done...\n');
toc
tic
%% Connectivity detection

    fprintf('Finalizing result...\n');
    CC = bwconncomp(im_bw);
    stats = regionprops(CC,'Area','Centroid');
    pt_list = [nan nan]';
    pt_area = [];
    k=1;
    if exclude_flag==1 % exclude some false positive
        ex_list(1,:) = (ex_list(1,:) - data.x(1))/stepX;
        ex_list(2,:) = (ex_list(2,:) - data.y(1))/stepY;
        ex_range2 = 25/stepX/stepY;
        for i = 1:length(stats)
            if stats(i).Area>SizeLimit(1)&&stats(i).Area<SizeLimit(2)
                pt_this = stats(i).Centroid';
                d_this = ex_list-[pt_this(2);pt_this(1)];
                d_this = d_this(1,:).^2 + d_this(2,:).^2;
                if sum(d_this < ex_range2)==0 % Not near to any exclude list.
                    pt_list(:,k) = pt_this;
                    pt_area(k) = stats(i).Area;
                    k=k+1;
                end
            end
        end
    else
        
        for i = 1:length(stats)
            if stats(i).Area>SizeLimit(1)&&stats(i).Area<SizeLimit(2)
                pt_list(:,k) = stats(i).Centroid';
                pt_area(k) = stats(i).Area;
                k=k+1;
            end
        end
    end
%% Output
pt_area = pt_area*(stepX*stepY);
v = pt_list(1,:);
pt_list(1,:) = pt_list(2,:);
pt_list(2,:) = v;
pt_list(1,:) = pt_list(1,:)*stepX + min(data.x) - stepX;
pt_list(2,:) = pt_list(2,:)*stepY + min(data.y) - stepY;
fprintf('Result done.\n Number of Cells = %d\n',length(pt_area));
toc
%% Visualization for testing
figure(1065);cla;
plotSlice(data);
hold on;
scatter(pt_list(1,:),pt_list(2,:),64*pt_area(:)./max(pt_area(:)),'r');
title(['detected cells, number = ',num2str(length(pt_list(1,:)))]);
hold off;
end
function plotSlice(data)
% plotSlice(data) visualize a 2-D data set
pcolor(data.x,data.y,data.value');
axis equal;
axis tight;
shading flat;
end
function data_out = downSample(data,stepX,stepY,stepZ)
% downSample put data into a new grid
if ndims(data.value) == 2 %#ok<ISMAT>
    % Get range
    x1 = min(data.x);
    x2 = max(data.x);
    y1 = min(data.y);
    y2 = max(data.y);
    data_out.x = x1:stepX:x2;
    data_out.y = y1:stepZ:y2;
    [yy,xx] = meshgrid(data.y,data.x);
    [yyq,xxq] = meshgrid(data_out.y,data_out.x);
    data_out.value = interp2(yy,xx,data.value,yyq,xxq);
elseif ndims(data.value) == 3
    x1 = min(data.x);
    x2 = max(data.x);
    y1 = min(data.y);
    y2 = max(data.y);
    z1 = min(data.z);
    z2 = max(data.z);
    data_out.x = x1:stepX:x2;
    data_out.y = y1:stepY:y2;
    data_out.z = z1:stepZ:z2;
    [yyy,xxx,zzz] = meshgrid(data.y,data.x,data.z);
    [yyyq,xxxq,zzzq] = meshgrid(data_out.y,data_out.x,data_out.z);
    data_out.value = interp3(yyy,xxx,zzz,data.value,yyyq,xxxq,zzzq);
end
end