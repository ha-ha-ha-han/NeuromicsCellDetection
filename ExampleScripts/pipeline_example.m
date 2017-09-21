% THIS IS THE PIPELINE FOR THE NEUROMICS PROJECT
% Version 1.0  23-AUG-2017
% 
% This file contains all the pipeline steps. You can run this file section
% by section. Also, you can run separate files.
% 
% RUNNING SEQUENCE:
% 0. Set the current fold as 'RUN_THIS_FOLDER' (or the folder you put all 
%    the data)
% 1. run_config.m
%       Configure the environment
% 2. save_ref_slice.m
%       Set a reference slice for the rough registration of all slices
%       It often just needs to be done once and then you can copy and paste
%       your result for each run.
% 3. preprocess_zstack.m
%       Preprocess the ZStack data.
%       - Load data
%       - Detect cells
%       It often just needs to be done once and then you can copy and paste
%       your result for each run.
% 4. preprocess_slice.m
%       Preprocess the slice data.
%       - Load data
%       - Rough registration
%       - Detect cells
% 5. run_overnight_and_go_to_sleep.m
%       When everything is ready, run the actual registration algorithm.
% 6. sleep
%       Regular sleep hours help the project.
% 7. check the result in the data storage folder. (.tiff files)
% 8. visualization.m
%       Visualize the result.

% ===============================================================
%% 1. run_config.m
% Run this section to config the environment.
% In the example folder RUN_THIS_FOLDER\data\MD038\16182, only scene 08, 10 
% and 13 are available. However, you can still keep the following code
% unchanged, since no data is actually loaded in this section. Only the
% information of the slices are recorded.
% ===============================================================

% --------- Set the data file folder -------------
clear
reference_position = 9;
slice_file_source = ... 
    fullfile(pwd,'data','MD038','16186');   % Specify the path for slice file storage
slice_positions = [1 2 3 4 5 6 12 11 10 9 8 7 13 14 15 16 17 18];   % Scene 1 is at the position 1. Scene 7 is at the position 12.
slice_positions = [slice_positions;slice_positions];
Position = slice_positions(:);                                      % [Scene1_Green Scene1_Red,Scene2_Green,Scene2_Red,...]
fmt = ['2017_03_22__16186-Scene-%02d-ScanRegion%d.czi ',...
    '- 2017_03_22__16186-Scene-%02d-ScanRegion%d.czi #3 - C=%d.tif']; % The format of the tif files.
% --------- Create the information table ------------
for i = 1:18
        FileName{2*i-1,1} = sprintf(fmt,i,i-1,i,i-1,1);
        Channel{2*i-1,1} = 'Green';
        [~,name_temp,~] = fileparts(FileName{2*i-1,1});
        DataName{2*i-1,1} = ['Slice_',name_temp];
        RunPath{2*i-1,1} = fullfile(pwd,DataName{2*i-1,1});
        
        FileName{2*i,1} = sprintf(fmt,i,i-1,i,i-1,2);
        Channel{2*i,1} = 'Red';
        [~,name_temp,~] = fileparts(FileName{2*i,1});
        DataName{2*i,1} = ['Slice_',name_temp];
        RunPath{2*i,1} = fullfile(pwd,DataName{2*i,1});
end
Channel = categorical(Channel);
slice_files = table(DataName,Channel,Position,FileName,RunPath);
save('run_info.mat','slice_files','slice_file_source','reference_position');
fprintf('[run_config] Environment configuration done at %s\n',...
    datestr(now))
fprintf('[run_config] run_info.mat saved at %s\n',pwd)

% ===============================================================
%% 2. save_ref_slice.m
% Run this section to set a reference slice for the rough registration 
% of all slices
% ===============================================================

% ----------- Load files -----------
clear
load run_info.mat
% ----------- Set filepath -----------
Channel = 'Red';
Position = reference_position; % A slice for reference. It is set in run_config.m
pf = (slice_files.Position == Position);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(pf.*cf==1,:);
% ----------- Save reference slice -----------
slice_files_ref = slice_files(pf.*cf==1,:);
SliceFileName = slice_files_ref.FileName{1,1};
SliceDataName = slice_files_selected.DataName{1,1};
source_filePathName = fullfile(slice_file_source,SliceFileName);
XStep = 4/3;
ZStep = 4/3;
dataS_ref = neuroReg.loadTiff(XStep,ZStep,source_filePathName);
data_slice_ref.value = squeeze(dataS_ref.value(:,4,:));
data_slice_ref.x = dataS_ref.x;
data_slice_ref.y = dataS_ref.z;
data_slice_ref.value = flip(flip(data_slice_ref.value,1),2);
figure;
neuroReg.plotData2(data_slice_ref);
save reference_slice.mat data_slice_ref
fprintf('[save_ref_slice] Reference slice saved at %s\n',...
    datestr(now))
fprintf('[save_ref_slice] reference_slice.mat saved at %s\n',pwd)

% ===============================================================
%% 3. preprocess_zstack.m
%       Preprocess the ZStack data.
%       - Detect cells
%       It often just needs to be done once.
%       Note: it assumes a mat file has been saved. To convert a raw data
%       to .mat format, use neuroReg.loadZStack(iSession) for .bin files
%       and neuroReg.loadZStackTiff(iSession) for tiff file.
%       The two loaders are written only for specific data sets so it may
%       be updated in future release.
%       Anyway, a user is advised to develop his own loader for the ZStack
%       data. In this project, a standard ZStack 'data' should be a 
%       structure with fields:
%       dataZ = 
% 
%            struct with fields:
% 
%               value: [512×512×267 single]
%                   x: [1×512 double]
%                   y: [1×512 double]
%                   z: [1×267 double]
    
% ===============================================================
% ------- Load Z data---------
clear;
iSession = 4;
ZStackFileName = {'12_13','18_19','12_13_green','12_13_green_test'};
ZStackFileName = ZStackFileName{iSession};
filepath = fullfile(pwd,'data','MD038','PreprocessedZStacks');
filename = ['ZStack_',ZStackFileName];
load(fullfile(filepath,filename));
% --------- Set options for neuroReg.detectCells3 ---------
Option_detect3.Sigma = [1 1 1.2]*5;
Option_detect3.Res0 = 0.0030;
Option_detect3.SizeLimit = [100,3000];
Option_detect3.Sensitivity = 1.25;
Option_detect3.MedianFilterSize = [7,7,7];
Option_detect3.Detect3Mode = 'Green';
Option_detect3.NeighborSize = [30 30 50];
% --------- Apply median filter -------------
% dataZ_mid = neuroReg.medfilt3(dataZ,Option_detect3);

% --------- Detect cells ---------------
% Core function: neuroReg.detectCells3
[pt_list_vol, pt_area_vol] = neuroReg.detectCells3(dataZ,Option_detect3);

% --------- Save data ----------------
ZStackInfo.FileName = ZStackFileName;
ZStackInfo.RunDateTime = datestr(now);
save(fullfile(filepath,filename),'ZStackInfo','dataZ','Option_detect3','pt_list_vol','pt_area_vol');
fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
fprintf('[preprocess_zstack] ZStack processed at %s\n',datestr(now));
fprintf(['[preprocess_zstack] ',filename,' saved at the folder ']);
disp(filepath);
fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
% ===============================================================
%% 4. preprocess_slice.m
%       Preprocess the slice data.
%       - Load data
%       - Rough registration
%       - Detect cells
% This is the boring part of all the steps. Because this step should be run
% for each slice. Because you need to hand pick some false-positive
% detections.
% ================================================================
% ------------ Load environment information -------------
clear
load run_info.mat
% ------------ Settings, filepath and options -----------
Channel = 'Green'; % Select the signal channel you want to preprocess.
Position = 13;  % Select the position you want to preprocess.
XStep = 4/3;    % pixel to um along X and Y axis. For down sampling.
ZStep = 4/3;    % pixel to um along Z axis. For down sampling.
THRESHOLD = 0.03; % Thresholding in the rouph registration

% ------------ Set current slice --------------------
pf = (slice_files.Position == Position);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(pf.*cf==1,:);
% ------------ Load refence data ----------------------
load reference_slice.mat
SliceFileName = slice_files_selected.FileName{1,1};
SliceDataName = slice_files_selected.DataName{1,1}; % DataName is the same as FileName except without .tif
filepath = fullfile(slice_files_selected.RunPath{1,1});
source_filePathName = fullfile(slice_file_source,SliceFileName);
% -------------- Load slice ----------------------
disp([SliceFileName, ' selected.']);
dataS = neuroReg.loadTiff(XStep,ZStep,source_filePathName);
data_slice.value = squeeze(dataS.value(:,4,:));
data_slice.x = dataS.x;
data_slice.y = dataS.z;
data_slice.value = flip(flip(data_slice.value,1),2);
disp([SliceFileName, ' Loaded.']);
% -------------- Rough registration for the slice ------------
[I_fix,R_fix] = neuroReg.data2img(data_slice_ref);
[I_mov_0,R_mov] = neuroReg.data2img(data_slice);
vm_fix = THRESHOLD * max(I_fix(:)); % Thresholding for the reference image
I_fix(I_fix > vm_fix)=vm_fix;
I_fix = I_fix/max(I_fix(:));
vm_mov = THRESHOLD * max(I_mov_0(:)); % Thresholding for the moving image
I_mov = I_mov_0;
I_mov(I_mov_0 > vm_mov)=vm_mov;
I_mov = I_mov/max(I_mov(:));
% -------------- Visualization -----------------------
figure(79814);
imshowpair(I_fix,R_fix,I_mov,R_mov,'falsecolor');
x_mov = R_mov.XWorldLimits;
y_mov = R_mov.YWorldLimits;
%% preprocess_slice: Transform and Visualization
% You need to run this section several times to achieve a good registration
% You can change translation vector t and rotation angle theta.
theta = 0; % Rotation angle
t = [700 -220]; % Translation vector
% theta = 0; 
% t = [0 0];
M = [cosd(theta) -sind(theta) 0;...
    sind(theta) cosd(theta) 0;...
        t(1)          t(2)   1];
tform = affine2d(M); % tform will be saved in the mat file. This variable can be used to roughly register other slices.
[I_mov_2,R_mov_2] = imwarp(I_mov,R_mov,tform,'cubic');
figure(79814);
imshowpair(I_fix,R_fix,I_mov_2,R_mov_2,'falsecolor');
%% preprocess_slice: continue
% ---------- Crop the image ------------
% Note: click and drag to create a rectangle that covers the region of
% interest. Then right click the rectangle and select 'Crop the image'.
figure(79814);
[~,rect] = imcrop;
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
%% preprocess_slice: Detect cells
% Change parameters. Repeat until satisfied.
% When the result is nearly okay, you can run next section to hand pick
% points to be excluded.

% ---------- Detect cells 2D ----------------
Option_detect2.Sigma = [1 1]*5;
Option_detect2.Res0 = 0.032;
Option_detect2.Threshold = 0;
Option_detect2.SizeLimit = [150,2000];
Option_detect2.MedianFilterSize = [10,10];
% ---------- Green Channel Only --------------
if strcmp(Channel,'Green')
    fprintf('[preprocess_slice] Green channel mode.\n')
    data_slice_temp = data_slice;
    data_slice_temp.value = data_slice_temp.value/max(data_slice_temp.value(:));
    v = adaptthresh(data_slice_temp.value,'NeighborhoodSize',[51 51]);
    data_slice_adp = data_slice;
    data_slice_adp.value = data_slice_temp.value./v;
    [pt_list_slice, pt_area_slice] = neuroReg.detectCells2(data_slice_adp,Option_detect2);

else
    [pt_list_slice, pt_area_slice] = neuroReg.detectCells2(data_slice,Option_detect2);
end
% ---------- Detect cells 2D: core function ----------------

ex_list = [];
%% preprocess_slice: Detect cells
% Steps:
%   - Click the picture to select points to be excluded.
%   - Press enter to apply.
%   - Click on points and press enter again until satisfied.
%   - Select nothing and press enter to finish.

% Backup the original data.
if strcmp(Channel,'Green')
    data_slice_orginal = data_slice;
    data_slice = data_slice_adp;
end
% Modify cell list for 2D
while 1
    h = gca;
    [x,y] = getpts; % Select nothing to exist
    if isempty(x)
        break
    end
    ex_list = [ex_list,[x,y]'];
    [pt_list_slice, pt_area_slice] = neuroReg.detectCells2(data_slice,Option_detect2,ex_list);
end
% ---------- Set default options --------------
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
% ---------- Test and visualize binarization function -------------------
% This step is to confirm the settings are okay.
fprintf('[preprocess_slice] Slice processed at %s\n',datestr(now));
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
% ---------- Save Slice ------------------
% Write image to record the relative position of the slice.
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
    'pt_list_slice','pt_area_slice','Option','Option_detect2','tform');
if isvarname('data_slice_orginal')
    % In case there is a backup version of the original data, save it as
    % well.
    save(MatFilePath,'data_slice_orginal','-append');
end
fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
fprintf(['[preprocess_slice] Slice_',SliceDataName,' saved at the folder\n']);
fprintf(['Folder path:\n']);
disp(filepath);
disp(datestr(now));
fprintf('||||||||||||||||||||||||||||||||||||||||||||\n');
%% 5. run_overnight_and_go_to_sleep.m
%       Everything is ready. Run the actual registration algorithm.
%       Check run_overnight_and_go_to_sleep.m for more information.
%       This section is listed here only for clarification purpose.
run_overnight_and_go_to_sleep
%% 6. sleep
%       Regular sleep hours help the project.
%% 7. visualization.m
%       Visualize the result after the run.
%       Check visualization.m for more information.
%       This section is listed here only for clarification purpose.
visualization
