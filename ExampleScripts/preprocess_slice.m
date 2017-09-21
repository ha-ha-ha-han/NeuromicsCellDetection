%% Load files
clear
load run_info.mat
%% Set filepath

Channel = 'Green';
Position = 13;
pf = (slice_files.Position == Position);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(pf.*cf==1,:);

%%
load reference_slice.mat
SliceFileName = slice_files_selected.FileName{1,1};
SliceDataName = slice_files_selected.DataName{1,1}; % without .tif
filepath = fullfile(slice_files_selected.RunPath{1,1});
source_filePathName = fullfile(slice_file_source,SliceFileName);
disp([SliceFileName, ' selected.']);
XStep = 4/3;
ZStep = 4/3;
%% Correction. TEMPERORY CODE
% MatFilePath = fullfile(filepath,[SliceDataName,'.mat']);
% load(MatFilePath);
% y1=data_slice.y(1);
% y2 = data_slice.y(end);
% data_slice.y = data_slice.y-y2-y1;
% ex_list(2,:) = ex_list(2,:) - y2 - y1;
% save(MatFilePath,'SliceInfo','data_slice','ex_list',...
%     'pt_list_slice','pt_area_slice','Option','Option_detect2');
%% Load Slice
dataS = neuroReg.loadTiff(XStep,ZStep,source_filePathName);
data_slice.value = squeeze(dataS.value(:,4,:));
data_slice.x = dataS.x;
data_slice.y = dataS.z;
% data_slice = neuroReg.downSample(data_slice,1.3,[],1.3);
data_slice.value = flip(flip(data_slice.value,1),2);
disp([SliceFileName, ' Loaded.']);

%% Rough registration for the slice
THRESHOLD = 0.03;
[I_fix,R_fix] = neuroReg.data2img(data_slice_ref);
[I_mov_0,R_mov] = neuroReg.data2img(data_slice);
% Thresholding
vm_fix = THRESHOLD * max(I_fix(:));
I_fix(I_fix > vm_fix)=vm_fix;
I_fix = I_fix/max(I_fix(:));
vm_mov = THRESHOLD * max(I_mov_0(:));
I_mov = I_mov_0;
I_mov(I_mov_0 > vm_mov)=vm_mov;
I_mov = I_mov/max(I_mov(:));
% Visualization
figure(79814);
imshowpair(I_fix,R_fix,I_mov,R_mov,'falsecolor');
x_mov = R_mov.XWorldLimits;
y_mov = R_mov.YWorldLimits;
%% Transform and Visualization
theta = 0;
t = [700 -220];
% theta = 0;
% t = [0 0];
M = [cosd(theta) -sind(theta) 0;...
    sind(theta) cosd(theta) 0;...
        t(1)          t(2)   1];
tform = affine2d(M);
[I_mov_2,R_mov_2] = imwarp(I_mov,R_mov,tform,'cubic');
figure(79814);
imshowpair(I_fix,R_fix,I_mov_2,R_mov_2,'falsecolor');

%% Truncate image
figure(79814);
[~,rect] = imcrop;
%%
rectangle('Position',rect,'EdgeColor','r', 'LineWidth',3);
[I_mov_0_2,R_mov_2] = imwarp(I_mov_0,R_mov,tform,'linear');
[indx1,indy1] = R_mov_2.worldToSubscript(rect(1),rect(2));
[indx2,indy2] = R_mov_2.worldToSubscript(rect(1)+rect(3),rect(2)+rect(4));
XStep = R_mov_2.PixelExtentInWorldX;
YStep = R_mov_2.PixelExtentInWorldY;
[x1,y1] = R_mov_2.intrinsicToWorld(indy1,indx1);
[x2,y2] = R_mov_2.intrinsicToWorld(indy2,indx2);
v = I_mov_0_2(indx1:indx2,indy1:indy2);
data_slice.value = flipud(v)';
data_slice.x = x1:XStep:x2;
data_slice.y = -y2:YStep:-y1;
%% Detect cells 2D
Option_detect2.Sigma = [1 1]*5;
Option_detect2.Res0 = 0.033;
Option_detect2.Threshold = 0;
Option_detect2.SizeLimit = [150,2000];
Option_detect2.MedianFilterSize = [10,10];

[pt_list_slice, pt_area_slice] = neuroReg.detectCells2(data_slice,Option_detect2);
%%
figure(666);
neuroReg.detectCells2(data_slice_adp,Option_detect2);
%%
ex_list = [];

%% Modify cell list for 2D
h = gca;
% [pt_list_slice,pt_area_slice] = neuroReg.addDelCells2(h,pt_list_slice,pt_area_slice);
[x,y] = getpts;
ex_list = [ex_list,[x,y]'];
[pt_list_slice, pt_area_slice] = neuroReg.detectCells2(data_slice,Option_detect2,ex_list);
%% Set Options
Option = neuroReg.setOption(Option_detect2);
Option.StepX = 5;
Option.StepD = 5;
Option.Integ = 15;
Option.MagicNumber = 2;
Option.CellRadius = 8;
Option.TransTol = 50;
Option.AngleTol = 4;
Option.MaxPeakNum = 200;
Option.SigmaRender = 15;

%% Test BW
disp(datestr(now));
disp('Option');
disp(Option);
[data_slice_bw,pt_list_slice] = neuroReg.bwCell2(data_slice,Option,ex_list);
figure(10089);
subplot(2,1,1);
neuroReg.plotData2(data_slice);
hold on
scatter(pt_list_slice(1,:),pt_list_slice(2,:),'r');
hold off;
subplot(2,1,2);
neuroReg.plotData2(data_slice_bw);
%% Save Slice
% Write image
[~] = mkdir(filepath);
img_path = fullfile(filepath,'SlicePosition.png');
img_path1 = fullfile(filepath,'SlicePosition.fig');
saveas(figure(79814),img_path);
saveas(figure(79814),img_path1);
SliceInfo.FileName = SliceDataName;
SliceInfo.RunDateTime = datestr(now);
disp(Option);
MatFilePath = fullfile(filepath,[SliceDataName,'.mat']);
save(MatFilePath,'SliceInfo','data_slice','ex_list',...
    'pt_list_slice','pt_area_slice','Option','Option_detect2');
fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
fprintf(['Slice_',SliceDataName,' saved at the folder\n']);
fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
disp(filepath);
disp(datestr(now));
