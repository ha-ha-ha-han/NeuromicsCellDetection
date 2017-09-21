%% Load Z data
% This pipeline needs to be modified
iSession = 3;
ZStackFileName = {'12_13','18_19','12_13_green'};
ZStackFileName = ZStackFileName{iSession};
filepath = fullfile(pwd,'\ZStack\');
filename = ['ZStack_',ZStackFileName];
load(fullfile(filepath,filename));
%% Detect cells 3D, preparation
Option_detect3.Sigma = [1 1 1]*4;
Option_detect3.Res0 = 0.0031;
Option_detect3.SizeLimit = [300,3000];
Option_detect3.Sensitivity = 0.13;
Option_detect3.MedianFilterSize = [7,7,7];

%% Apply median filter
dataZ_mid = neuroReg.medfilt3(dataZ,Option_detect3);

%% Detect cells 3D
[pt_list_vol, pt_area_vol] = neuroReg.detectCells3(dataZ_mid,Option_detect3);

%% Save data
ZStackInfo.FileName = ZStackFileName;
ZStackInfo.RunDateTime = datestr(now);
save(fullfile(filepath,filename),'ZStackInfo','dataZ','Option_detect3','pt_list_vol','pt_area_vol');
disp(datestr(now));
fprintf([filename,' saved at the folder ']);
disp(filepath);