function Option = setOption(Option)
%% For image reconstruction from cell list
if ~isfield(Option,'Integ') % Integ: slice thickness
    Option.Integ = 10;
end
if ~isfield(Option,'CellRadius') % CellRadius: radius of a typical neuron
    Option.CellRadius = 5;
end
if ~isfield(Option,'StepD') % stepD: slice step (recommend: Integ/2)
    Option.StepD = 5;
end
if ~isfield(Option,'StepX') % stepX: resolution of the slice
    Option.StepX = 7;
end
%% For correlation function calculation
if ~isfield(Option,'MagicNumber')  % Intensity ratio for cell and background
    Option.MagicNumber = 2;
end
%% For cell detection in the slice data
if ~isfield(Option,'Sigma')  % For cell detection in the slice data
    Option.Sigma = [1 1]*8;
end
if ~isfield(Option,'SigmaRender')  % For cell detection in the slice data
    Option.SigmaRender = [1 1]*12;
end
if ~isfield(Option,'Res0')  % For cell detection.
    Option.Res0 = 0.03;
end
if ~isfield(Option,'SizeLimit')  % For cell detection. Cell Size Limit
    Option.SizeLimit = [100,2000];
end
if ~isfield(Option,'Threshold')  % For cell detection. Cell Size Limit
    Option.Threshold = 0;
end
if ~isfield(Option,'MedianFilterSize')  % For cell detection. Cell Size Limit
    Option.MedianFilterSize = [10,10];
end
%% For rotateCorr3, peak selection
if ~isfield(Option,'TransTol')
    Option.TransTol = 25^2;
end
if ~isfield(Option,'AngleTol')
    Option.AngleTol = 4^2;
end
if ~isfield(Option,'MaxPeakNum')
    Option.MaxPeakNum = 100;
end
if ~isfield(Option,'Hist3Flag')
    Option.Hist3Flag = 0;
end
if ~isfield(Option,'Hist3Smooth')
    Option.Hist3Smooth = 0;
end
if ~isfield(Option,'Hist3SmoothSize')
    Option.Hist3SmoothSize = [3,3,3];
end
FieldNames = {'Integ','CellRadius','StepD','StepX',...
    'MagicNumber',...
    'Sigma','SigmaRender','Res0','SizeLimit','Threshold','MedianFilterSize',...
    'TransTol','AngleTol','MaxPeakNum',...
    'Hist3Flag','Hist3Smooth','Hist3SmoothSize'};
Option = orderfields(Option,FieldNames);
end

