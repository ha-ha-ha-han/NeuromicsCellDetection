%% Load files
clear
load run_info.mat
%% Set filepath

Channel = 'Red';
Position = reference_position;
pf = (slice_files.Position == Position);
cf = (slice_files.Channel == Channel);
slice_files_selected = slice_files(pf.*cf==1,:);

%% Save reference slice

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
% data_slice = neuroReg.downSample(data_slice,1.3,[],1.3);
data_slice_ref.value = flip(flip(data_slice_ref.value,1),2);
figure;
neuroReg.plotData2(data_slice_ref);
save reference_slice.mat data_slice_ref