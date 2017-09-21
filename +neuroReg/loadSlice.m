function data = loadSlice(XStep,ZStep,fileName)
%loadSlice loads raw slice data. It requires no input for now.
%
%loadSlice(FileName, XStep, YStep) is a function that loads raw slice data
%to work space.
% FileName is the file name of the data file.
% XStep, YStep and ZStep are the size for each pixel.
% The output data is a structure with data.value, data.x, data.y. 
% data.value is a data cube with 
% size(data.value) == [length(data.x),length(data.y)]
% dim1:X -> lateral to medial 
% dim3:Z -> superficial to deep
% dim2:Y -> different focus
% XStep = 1/3;
% ZStep = 1/3;
% fileName='16186_10_46_red';
fprintf(['    Loading file ',fileName, '...']); tic;
data.value = single(loadArr(['Slice Subsection\roi_' fileName]));
tt= toc; fprintf('%.2f s\n',tt)
data.value = permute(data.value,[2 3 1]);
data.value = flip(data.value,1);
data.value = flip(data.value,3);
nsize = size(data.value);
data.x = XStep*(1:nsize(1));
data.y = 1:nsize(2);
data.z = ZStep*(1:nsize(3));

