function data = loadZStack(iSession)
%loadZStack loads raw slice data. It requires no input for now.
%
%Ideally,loadZStackData(FileName, XStep, YStep, ZStep) is a function that loads 
% raw ZStack data to work space
% FileName is the file name of the data file.
% XStep, YStep and ZStep are the size for each pixel.
% The output data is a structure with data.value, data.x, data.y and
% data.z. data.value is a data cube with 
% size(data.value) == [length(data.x),length(data.y),length(data.z)]
% dim1:X -> medial (center of the brain) to lateral (near the ears)
% dim2:Y -> anterior (face) to posterior 
% dim3:Z -> superficial to deep

% data.value = single(loadArr(FileName));
% data.value = permute(data.value,[2 3 1]);
% nsize = size(data.value);
% data.x = XStep*(1:nsize(1));
% data.y = YStep*(1:nsize(1));
% data.z = ZStep*(1:nsize(1));

% iSession=1; % = 1 or 2: choose among one of the two 2p zstacks
tiffName = 'coronalMov_green';
fileIDX = [12 13; 18 19];
%every zstack is composed by two overlapping volumes
% movies at different depth Z and same X,Y
XStep=680/512;
YStep=665/512;
ZStep=2;

fileName = cell(1,2);
for iExp = 1:2
    fileName{iExp}=[tiffName num2str(fileIDX(iSession,iExp))];
end

%% load the zstacks
imageP=cell(1,2);
for iExp=1:2
    fprintf(['Loading file ',fileName{iExp}, '...'])
    tic;
%     imageP{iExp}=single(loadArr(['CoronalZstacks\' fileName{iExp}]));
    source_filePathName = fullfile(pwd,'CoronalZstacks\',[fileName{iExp},'.tif']);
    data = neuroReg.loadTiff(680/512,2,source_filePathName,665/512);
    imageP{iExp} = single(permute(data.value,[3,1,2])); % ugly code
    tt= toc; fprintf('%.2f s\n',tt)
end
%% information about real distance of cells from surface
if iSession==1
    dZ=1630-1487;
    dDu=1630-1580;
elseif iSession==2
    dZ=1986-1828;
    dDu=1986-1900;
end
%% compute the correctin on average intensity in order to combine the insitensity of the two z-stacks
dZp=round(dZ/ZStep);
maxZ=195;
temp1=mean(mean(imageP{2}(dZp:maxZ,:,:),3),2); %superficial
temp2=mean(mean(imageP{1}(1:(maxZ-dZp+1),:,:),3),2); %deeper
dFf=mean(temp1)-mean(temp2);


imageN=nan(maxZ+dZp,512,512);
imageN(1:dZp,:,:)=imageP{2}(1:dZp,:,:);
imageN((dZp+1):maxZ,:,:)=(imageP{2}((dZp+1):maxZ,:,:)+imageP{1}(1:(maxZ-dZp),:,:)+dFf)./2;
imageN(maxZ+1:(maxZ+dZp),:,:)=imageP{1}((maxZ-dZp+1):maxZ,:,:)+dFf;
imageN=single(imageN);




data.value = imageN;
data.value = permute(data.value,[2 3 1]);
data.value = flip(data.value,3);
nsize = size(data.value);
data.x = XStep*(1:nsize(1));
data.y = YStep*(1:nsize(2));
data.z = ZStep*(1:nsize(3));

end

