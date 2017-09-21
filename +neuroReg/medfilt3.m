function data_out = medfilt3(data,Option)
% neuroReg.medfilt3 is simply a replacement of the default medfilt3.
% data_out = medfilt3(data,Option)
% data: structure. data.x, data.y, data.z, data.value
% Option: Option.MedianFilterSize, [SizeX, SizeY, SizeZ] in um

stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
stepZ = mean(diff(data.z));
MedianFilterSize = Option.MedianFilterSize;
% Convert from um to pixel
MedianFilterSizePx = ceil(MedianFilterSize./[stepX,stepY,stepZ]);
MedianFilterSizePx = floor(MedianFilterSizePx/2)*2+1;
v = data.value / max(data.value(:));

tic;
fprintf('Applying median filter...\n');
v = medfilt3(v,MedianFilterSizePx);
fprintf('Median filter done.\n');
toc

data_out = data;
data_out.value = v;
end

