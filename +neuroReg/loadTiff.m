function data = loadTiff( stepX,stepZ,source_filePathName,stepY )
% loadTiff loads data from a TIff file.
% It works with both 2D and 3D Tiff.
% data = loadTiff( stepX,stepZ,source_filePathName)
% data = loadTiff( stepX,stepZ,source_filePathName,stepY)
% StepX and StepZ specifies the step of each pixel.
% source_filePathName is the file path and file name. (use fullfile!)
% stepY is optional. When missing, it assumes stepY = StepX.


finfo = imfinfo(source_filePathName);
for i = 1:length(finfo)
    if i == 1
        v = imread(source_filePathName,i);
        value = zeros([size(v),length(finfo)]);
        value(:,:,1) = v;
    else
        value(:,:,i) = imread(source_filePathName,i);
    end
end
%% Convert from the 'Image Convention' to the 'Data Convention'.
[z,x,y] = size(value);
data.value = permute(value,[2,3,1]);
data.x = (1:x)*stepX;
if nargin == 5
    data.y = (1:y)*stepY;
else
data.y = 1:y;
end
data.z = (1:z)*stepZ;


end

