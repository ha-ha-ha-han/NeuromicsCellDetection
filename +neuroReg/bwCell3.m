function data_out = bwCell3(data,Option)
% Use difference of Gaussian detector
Res0 = Option.Res0;
v = data.value / max(data.value(:));
stepX = mean(diff(data.x));
stepY = mean(diff(data.y));
stepZ = mean(diff(data.z));
sigmaInd = [1 1 1];
sigmaInd(2) = Option.Sigma(1)/stepX; % in imgaussfilt, sigmaInd(2) is for 1st dimension
sigmaInd(1) = Option.Sigma(2)/stepY; % sigmaInd(1) is for 2nd dimension
sigmaInd(3) = Option.Sigma(3)/stepZ;
Ic = imgaussfilt3(v,sigmaInd);
Is = imgaussfilt3(v,sigmaInd*5);
res = Ic - Is;
data_out = data;
data_out.value = (res>Res0)*1.0;
data_out.value = imgaussfilt3(data_out.value,sigmaInd*1.5);
end
