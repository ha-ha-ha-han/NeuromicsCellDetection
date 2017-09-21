clear;
load run_info.mat
%%
Position = [9 11 13];
Channel = 'Green';
% Specify the parent folder for this run
FolderPath = pwd;
% Specify the storage path for preprocessed ZStack files
FolderPathForZStack = fullfile(FolderPath,'data','MD038','PreprocessedZStacks');
cf = (slice_files.Channel == Channel);
slice_files_cf = slice_files(cf,:);
pf=zeros(size(slice_files_cf.Position)); % Position Flag
for k = 1:length(Position)
    pf = pf+(slice_files_cf.Position == Position(k));
end
slice_files_selected=slice_files_cf(pf>0,:);
slice_files_selected = sortrows(slice_files_selected,'Position');
SlicesName = slice_files_selected.DataName; % i
ZStacksName = {'12_13','18_19','12_13_green','12_13_green_test'};   
AngleRange = [-10 40 26;-15 15 15;-15 15 15];
disp('Slices: ');
disp(slice_files_selected);
disp('ZStacks: ');
disp(ZStacksName');
fprintf('Preparation ready for visualization transforms.\n')
i = 1;
j = 1;
%%
for i = 3:3
    for j = 4:4
        for HF = 2:2
            this_slice = SlicesName{i};
            this_ZStack = ['ZStack_',ZStacksName{j}];
            this_slice_path = fullfile(FolderPath,this_slice);
            this_ZStack_file = fullfile(FolderPathForZStack,[this_ZStack,'.mat']);
            if HF == 1
                Option.Hist3Flag = 1;
                this_result_path = fullfile(this_slice_path,[this_ZStack,'_hist\']);
            else
                Option.Hist3Flag = 0;
                this_result_path = fullfile(this_slice_path,[this_ZStack,'_xcc\']);
            end
            fprintf('Result in %s',this_result_path);
            filepathname = fullfile(this_result_path,'result.mat');
%             [~] = mkdir(this_result_path);
            load(this_ZStack_file); % load zstack
            fprintf('%s Loaded.\n',this_ZStack);
            load(fullfile(this_slice_path,[this_slice,'.mat'])); % load slice
            fprintf('%s Loaded.\n',this_slice);
            load(filepathname); % load result
            fprintf('Result Loaded.\n');
            %% neuroReg.VisTransform
            % Set DataSets structure to record all the raw data (and binarized slice)
            DataSets.dataZ = dataZ;
            DataSets.data_slice = data_slice;
            data_slice_bw = neuroReg.bwCell2(data_slice,Option,ex_list);
            DataSets.data_slice_bw_low = neuroReg.downSample(data_slice_bw,Option.StepX,[],Option.StepX);
            obj = neuroReg.VisTransform(TransTable,DataSets,pt_list_vol,pt_list_slice,ex_list,Option,this_result_path);
            fprintf('VisTransform\n');
%% ----- Visualization
% 
% fprintf('*******************\n');
% disp(this_result_path);
% disp(RunInfo);
% h = figure(10086);
% fprintf('Saving figures...\n')
% for k = 1:size(TransTable,1)
%     if mod(k,10)==0
%         fprintf('%3d Figures Saved...\n',k);
%     end
%     h.Name = ['Peak_',num2str(k,'%04d')];
%     neuroReg.plotTransform(h,TransTable(k,:),DataSets,pt_list_vol,pt_list_slice,Option,ex_list)
%     filename = ['Peak_',num2str(k,'%04d'),'.tif'];
%     saveas(h, [this_result_path,filename])
% end
% disp(datestr(now));
% fprintf('Mission Completed!\n*******************\n');
        end
    end
end