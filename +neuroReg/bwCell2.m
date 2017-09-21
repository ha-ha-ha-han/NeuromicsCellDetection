function [data_out,pt_list] = bwCell2(data,Option,ex_list)
% bwCell2 binarize a slice data.
% [data_out,pt_list] = bwCell2(data,Option,ex_list)
% INPUT:
% data: a slice data with data.x, data.y, data.value (2D only).
% For 3D data, use bwCell3.
% Option: see neuroReg.setOption
% ex_list: points to exclude. Generated in the preprocess stage. 2-by-M
% OUTPUT:
% data_out: binarized and then gaussian-blured slice data
% pt_list: detected cells. It should be the same with the result from
% detectCells3 with the same input. 2-by-N


if nargin==2
    exclude_flag=0;
elseif nargin==3&&isempty(ex_list)
    exclude_flag=0;
elseif nargin==3&&~isempty(ex_list)
    exclude_flag=1;
else    
    error(['Check input.','neuroReg.bwCell2(data,Option)',...
        ' or neuroReg.bwCell2(data,Option, ex_list)']);
end
if isfield(Option,'Threshold')
    Threshold = Option.Threshold;
else
    Threshold = 0;
end
if isfield(Option,'SigmaRender')
    SigmaRender = Option.SigmaRender;
else
    SigmaRender = 1.5*Option.Sigma;
end
%% Get slice
stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
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
sigmaRenderInd(1) = SigmaRender(1)/stepY;
sigmaRenderInd(2) = SigmaRender(1)/stepX;
Ic = imgaussfilt(v,sigmaInd);
Is = imgaussfilt(v,sigmaInd*5);
res = Ic - Is;
im_bw = imbinarize(res,Res0)*1.0;
%% Connectivity detection
CC = bwconncomp(im_bw);
L = labelmatrix(CC);
stats = regionprops(CC,'Area','Centroid');
k=1;
pt_list = [];
if exclude_flag==1 % exclude some false positive
    ex_list(1,:) = (ex_list(1,:) - data.x(1))/stepX;
    ex_list(2,:) = (ex_list(2,:) - data.y(1))/stepY;
    ex_range2 = 25/stepX/stepY;
    for i = 1:length(stats)
        if stats(i).Area<SizeLimit(1)||stats(i).Area>SizeLimit(2)
            L(L==i)=0;
        else
            pt_this = stats(i).Centroid';
            d_this = ex_list-[pt_this(2);pt_this(1)];
            d_this = d_this(1,:).^2 + d_this(2,:).^2;
            if sum(d_this < ex_range2)==0 
                % if there is no near exclude point, then keep the cell
                pt_list(:,k) = stats(i).Centroid';
                k = k+1;
            else
                L(L==i)=0;
            end
        end
    end
else
    for i = 1:length(stats)
        if stats(i).Area<SizeLimit(1)||stats(i).Area>SizeLimit(2)
            L(L==i)=0;
        else
            pt_list(:,k) = stats(i).Centroid';
            k = k+1;
        end
    end
end
if ~isempty(pt_list)
    v = pt_list(1,:);
    pt_list(1,:) = pt_list(2,:);
    pt_list(2,:) = v;
    pt_list(1,:) = pt_list(1,:)*stepX + min(data.x) - stepX;
    pt_list(2,:) = pt_list(2,:)*stepY + min(data.y) - stepY;
end
L(L>=1) = 1;
%% Output
data_out = data;
data_out.value=double(L);
data_out.value = imgaussfilt(data_out.value,sigmaRenderInd);
% figure(10087);pcolor(Ic');shading flat;axis equal;
% figure(10088);pcolor(Is');shading flat;axis equal;
% figure(10089);surf(res');shading flat;axis vis3d;
% figure(10090);surf(v');shading flat;axis vis3d;
end