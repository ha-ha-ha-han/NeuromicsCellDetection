clear
reference_position = 9;
slice_file_source = ...
    fullfile('D:\ProjectLocalDataBase\Neuron\data\MD038\16186');
slice_positions = [1 2 3 4 5 6 12 11 10 9 8 7 13 14 15 16 17 18];
slice_positions = [slice_positions;slice_positions];
Position = slice_positions(:);
fmt = ['2017_03_22__16186-Scene-%02d-ScanRegion%d.czi ',...
    '- 2017_03_22__16186-Scene-%02d-ScanRegion%d.czi #3 - C=%d.tif'];
%%
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