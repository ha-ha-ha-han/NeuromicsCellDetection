% Before running the code, make sure you have already gone through all the
% preprocess.
% ***** TIPS *********
% You don't need to run the preprocess everytime. If you have some
% preprocess result, you can copy and paste the resulting .mat files to the
% designated folder.
% Where is the designated folder?
% Check the variable 'this_slice_path' below:
% --- Code relating to the file path ----------
% FolderPath = pwd;
% SlicesName = slice_files_selected.DataName;
% this_slice = SlicesName{i};
% this_ZStack = ['ZStack_',ZStacksName{j}];
% this_slice_path = fullfile(FolderPath,this_slice);
% this_ZStack_file = fullfile(FolderPath,'ZStack',[this_ZStack,'.mat']);
% ------------------------------------------------

%% Load Information
clear;
load run_info.mat
Position = [9 11 13];
Channel = 'Green'; % 'Red' for red channel, 'Green' for green channel
% Specify the parent folder for this run
FolderPath = pwd;
% Specify the storage path for preprocessed ZStack files
FolderPathForZStack = fullfile(FolderPath,'data','MD038','PreprocessedZStacks');
cf = (slice_files.Channel == Channel);
slice_files_cf = slice_files(cf,:);
pf=zeros(size(slice_files_cf.Position));
for k = 1:length(Position)
    pf = pf+(slice_files_cf.Position == Position(k));
end
slice_files_selected=slice_files_cf(pf>0,:);
slice_files_selected = sortrows(slice_files_selected,'Position');
SlicesName = slice_files_selected.DataName; % i
ZStacksName = {'12_13','18_19','12_13_green','12_13_green_test'};  % j
% AngleRange = [Alpha_start Alpha_end Alpha_points; Beta_start...; Gamma_...]
AngleRange = [-10 10 21;-10 10 21;-10 10 21];
% AngleRange = [0 0 1; 6.43 6.43 1; -4.2857 -4.2857 1];
disp('Slices: ');
disp(slice_files_selected);
disp('ZStacks: ');
disp(ZStacksName');
fprintf('Preparation ready.\n')
fprintf('Go home. Sleep. Run overnight.\n')
%% Go home. Sleep. Run overnight.
% AngleRange = [0 4 5; 4 8 5; 0 5 5];
% i for slice, j for ZStack
for i=3:3
    for j = 4:4        
        for HF = 2:2 % HF == 1: Hist. HF == 2: XCC.
        %% ----- Set Path -------
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
        [~] = mkdir(this_result_path);
        load(this_ZStack_file);
        load(fullfile(this_slice_path,[this_slice,'.mat']));
        RunInfo.SliceName = this_slice;
        RunInfo.ZStackName = this_ZStack;
        RunInfo.AngleRange = AngleRange;
        Option.StepX = 8;
        Option.StepD = 12;
        Option.Hist3Flag = 1;
        Option.Hist3Smooth = 1; % 1 = Box, 2 = Gaussian
        Option.MaxPeakNum = 50;
        %% ----- Calculation ----------        
        fprintf('*******************\n');
        disp(this_result_path);        
        RunInfo.StartTime = datestr(now);     
        disp(datestr(now));
        fprintf('Calculation started...\n') 
        data_slice_bw = neuroReg.bwCell2(data_slice,Option,ex_list);
        h = figure(5556233);
        subplot(3,1,1);
        neuroReg.plotData2(data_slice);
        subplot(3,1,2);
        neuroReg.plotData2(data_slice);
        hold on;
        scatter(pt_list_slice(1,:),pt_list_slice(2,:),'r');
        hold off;
        subplot(3,1,3);
        neuroReg.plotData2(data_slice_bw);
        saveas(h, fullfile(this_result_path,'Slice.tif'))
        tic
        if HF == 1
            TransTable = neuroReg.rotationHist3(pt_list_slice,pt_list_vol,data_slice,AngleRange,Option,ex_list);
        else
            TransTable = neuroReg.rotationCorr3(pt_list_slice,pt_list_vol,data_slice,AngleRange,Option,ex_list);
        end
        toc
        RunInfo.EndTime = datestr(now);
        disp(RunInfo.EndTime);
        fprintf('Calculation completed.\n')
        filepathname = fullfile(this_result_path,'result.mat');
        save(filepathname,'RunInfo','TransTable',...
            'Option','Option_detect2','Option_detect3',...
            'pt_list_vol','pt_area_vol','pt_list_slice','pt_area_slice','ex_list');
        disp(RunInfo);
        t1 = datetime(RunInfo.StartTime);
        t2 = datetime(RunInfo.EndTime) ;
        t = t2-t1;
        fprintf('Duration ');
        disp(t)
        fprintf('Option\n ');
        disp(Option);
        fprintf('Result saved!\n---------------\n');
        %% ----- Visualization
        h = figure(10086);
        fprintf('Saving figures...\n')
        DataSets.dataZ = dataZ;
        DataSets.data_slice = data_slice;
        data_slice_bw = neuroReg.bwCell2(data_slice,Option,ex_list);
        DataSets.data_slice_bw_low = neuroReg.downSample(data_slice_bw,Option.StepX,[],Option.StepX);
        for k = 1:size(TransTable,1)
            if mod(k,10)==0
                fprintf('%3d Figures Saved...\n',k);
            end
            h.Name = ['Peak_',num2str(k,'%04d')];
            neuroReg.plotTransform(h,TransTable(k,:),DataSets,pt_list_vol,pt_list_slice,Option,ex_list)
            filename = ['Peak_',num2str(k,'%04d'),'.tif'];
            saveas(h, fullfile(this_result_path,filename))
        end
        disp(datestr(now));
        fprintf('Mission Completed!\n*******************\n');
        end
    end
end