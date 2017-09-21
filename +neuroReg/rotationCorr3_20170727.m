function PeakTable = rotationCorr3(pt_list,data_slice,AngleRange,Option,ex_list)
% rotationCorr3 returns the transformation parameters and the corresponding
% correlation functions.
% The transfermation matrix is from volume to slice.
% OUTPUT:
% PeakTable record the possible transformation parameters and the
% coorelation intensity.
% n-by-(1+6) table with column names:
% Intensity, Angle, Translation
% [{intensity},{alpha,beta,gamma},{tx,ty,yz}]
% The transformation is from volume to cut.
% Corr record the corresponding correlation functions.
% 1-by-n
% INPUT:
% pt_list is 3xN array. Coordinations of the cells in 3D volume
% data_slice is the *2D* slice (with dataS.x, dataS.y and dataS.value)
% AngleRange should be a 3x3 matrix with the form
% [Alpha_Start, Alpha_End, Number; Beta_Start, Beta_End, Beta_Number; Gamma...]
% alpha: rotation around y-axis
% beta: rotation around z-axis
% gamma: rotation around x-axis
% Convention:
% R = Rx*Rz*Ry
% Y = R*X
% **********
% * OPTION *
% **********
% Option is a structure specifying the options used in calculating
% correlation function.
% Fields not specified will be set as default.
% Below is the description of each field:
% (I write it with enormous patience because tomorrow morning I will forget
% them all)
% ---------
% For peak record
% 'TransTol'    >> Translation tolerance to divide two peaks in XCC
% 'AngleTol'    >> Angle tolerance to divide two peaks in XCC
% 'MaxPeakNum'  >> Number of peaks in the output table
% ---------
% For image reconstruction from cell list
% 'Integ' >> Integration range of the reconstructed plane. Default: 10 (um)
% 'StepD' >> Off-plane translational resolution. Default: 5 (um)
% 'StepX' >> In-plane translational resolution. Smaller resolution leads to
% better result but also with more computation time. Default: 7 (um)
% 'CellRadius' >> The radius for each neuron, used in the render process.
% Default: 3 (um)
% ---------
% For cell detection in the slice data
% 'Sigma'   >> [sx,sy],Size of the gaussian filter. Default: [8,8] (um)
% 'Res0'    >> Intensity threshold for the cell. Change with the resolution
% of the input dataS. Trial with neuroReg.bwCell2 is recommended.
% Default: 0.03
% 'SizeLimit' >> Size limit of the cell. Default: [100,2000] (um2)
% ---------
% For correlation function calculation
% MagicNumber >> I = I - MageicNumber * Mean(I(:)) is calculated before
% going to correlation function Default:2
% PeakThreshold >> from 0 to 1. The peaks larger than
% PeakThreshold * PeakIntensity will be regarded as a possible match.
% Defualt: 0.8
% ---------
% By Han Peng, penghan1992@gmail.com
% 2017.07.18
%
%
% **********
% * TEST *
% **********
% pt_list = evalin('base','pt_list2');
% data = evalin('base','dataS_resam_gen');
% N = 21;                             % N: Number of angles
% alpha = linspace(-5,5,N);           % alpha: rotation around y-axis
% beta = linspace(-5,5,N);            % beta: rotation around z-axis
% gamma = linspace(-5,5,N);           % gamma: rotation around x-axis
% Integ = 10;                         % Integ: slice thickness
% stepD = 5;                          % stepD: slice step (recommend: Integ/2)
% stepX = 7;                          % stepX: resolution of the slice
% CellRadius = 3;
% MagicNumber = 2;                    % Intensity ratio for cell and background



%% Rotation Angles
alpha1 = AngleRange(1,1);
alpha2 = AngleRange(1,2);
n_alpha = AngleRange(1,3);
beta1 = AngleRange(2,1);
beta2 = AngleRange(2,2);
n_beta = AngleRange(2,3);
gamma1 = AngleRange(3,1);
gamma2 = AngleRange(3,2);
n_gamma = AngleRange(3,3);
alpha = linspace(alpha1,alpha2,n_alpha);           % alpha: rotation around y-axis
beta = linspace(beta1,beta2,n_beta);            % beta: rotation around z-axis
gamma = linspace(gamma1,gamma2,n_gamma);           % gamma: rotation around x-axis

%% Option
% Peak Records
Option = neuroReg.setOption(Option); % Check input, set default value.
TRANSLATION_TOLERANCE2 = Option.TransTol^2;
ANGLE_TOLERANCE2 = Option.AngleTol^2;
MAX_PEAK_RECORD_NUM = Option.MaxPeakNum;
% For image reconstruction from cell list
Integ = Option.Integ; % Integ: slice thickness
stepD = Option.StepD; % stepD: slice step (recommend: Integ/2)
stepX = Option.StepX; % stepX: resolution of the slice
CellRadius = Option.CellRadius; % CellRadius: radius of a typical neuron
% For cell detection in the slice data
option2.Sigma = Option.Sigma; % For cell detection in the slice data
option2.SigmaRender = Option.SigmaRender; % For cell detection in the slice data
option2.Res0 = Option.Res0; % For cell detection.
option2.SizeLimit = Option.SizeLimit; % For cell detection. Cell Size Limit
option2.Threshold = Option.Threshold; % For cell detection. Cell Size Limit
option2.MedianFilterSize = Option.MedianFilterSize; % For cell detection. Cell Size Limit
% For correlation function calculation
MagicNumber = Option.MagicNumber; % Intensity ratio for cell and background
%% Binarize the slice image
data_slice_bw = neuroReg.bwCell2(data_slice,option2,ex_list);
data_slice_bw_low = neuroReg.downSample(data_slice_bw, stepX, [], stepX);
data_slice_bw_low.value = data_slice_bw_low.value - mean(data_slice_bw_low.value(:))*MagicNumber;
figure(100860); cla;
subplot(2,1,1)
neuroReg.plotData2(data_slice);
subplot(2,1,2)
neuroReg.plotData2(data_slice_bw_low);
title('Binarization and render result from the input slice');

%% Rotation and Calculation of Correlation Function
h1 = waitbar(0,'Calculating correlation function... Have a coffee...',...
    'Name','Calculating XCC... ');
data_corr.x = data_slice_bw_low.x;
data_corr.z = data_slice_bw_low.y;
[nx,nz] = size(data_slice_bw_low.value);
corr_max = 0;
para_max = [nan;nan;nan;nan;nan;nan];
N = n_gamma * n_beta * n_alpha;
% Record the peak positions
% I love tables
PeakTable = table;
PeakTable.Intensity = zeros(0);
PeakTable.Angles = zeros(0,3);
PeakTable.Translation = zeros(0,3);
tic;
for i = 1:n_gamma
    for j = 1:n_beta
        for k = 1:n_alpha
            %% Rotate the cells from the volume to slice coordination
            [pt_list_rotated,R] = neuroReg.rotateCells(pt_list,alpha(k),beta(j),gamma(i));
            v1 = pt_list_rotated(1,:); % x-direction of the slice
            v2 = pt_list_rotated(3,:); % up-down-direction of the slice
            v3 = pt_list_rotated(2,:); % normal direction of the slice
            % Get the size of the rotated cube
            x_lim(1) = min(v1);
            x_lim(2) = max(v1);
            y_lim(1) = min(v2);
            y_lim(2) = max(v2);
            d_lim(1) = min(v3);
            d_lim(2) = max(v3);
            % Get the planes to calculate correlation function
            d_list = d_lim(1):stepD:d_lim(2);
            x_c(1,1) = mean(x_lim); % center position (x) of the 2d plane
            x_c(3,1) = mean(y_lim); % center position (z) of the 2d plane
            value_temp = zeros(nx,nz,length(d_list));
            % Preparation for FFT
            Im1 = data_slice_bw_low.value;
            Im2 = zeros(size(Im1));
            [Lx1,Ly1]=size(Im1);
            Lx2 = length(x_lim(1):stepX:x_lim(2));
            Ly2 = length(y_lim(1):stepX:y_lim(2));
            Mx1 = round(Lx1/2);
            My1 = round(Ly1/2);
            Mx2 = round(Lx2/2);
            My2 = round(Ly2/2);
            x_start = Mx1-Mx2+1;
            x_end = Mx1-Mx2+Lx2;
            y_start = My1-My2+1;
            y_end = My1-My2+Ly2;
            f1 = fft2(Im1);
            for m = 1:length(d_list)
                %% Calculate correlation function in 3D
                % x-z: in-plane directions. y: normal direction of the plane.
                % Reconstruct the cell image
                vd = v3 - d_list(m);
                selected = (abs(vd)<Integ);
                px = v1(selected==1);
                py = v2(selected==1);
                pt_now_dist = vd(selected==1);
                pt_now_area = ones(1,sum(selected==1)).*exp(-(pt_now_dist./CellRadius).^2/2);
                data_now = neuroReg.renderCell2(x_lim,y_lim,stepX,[px;py],pt_now_area,Option);
                data_now.value = data_now.value - mean(data_now.value(~isnan(data_now.value(:))))*MagicNumber;
                data_now.value(isnan(data_now.value))=0;
                % Calculate the correlation function for the m-th plane
                % value_temp_slice = filter2(data_now.value,data_slice_bw_low.value);
                Im2_temp = data_now.value;
                Im2(x_start:x_end,y_start:y_end)=Im2_temp;
                f2 = fft2(Im2);
                value_temp_slice = fftshift(ifft2(f1 .* conj(f2)));   
                value_temp(:,:,m) = value_temp_slice;
            end
            %% After the distance loop, record the maxima position
            % THIS PART WILL CHANGE TO RECORD ALL POSSIBLE PEAKS
            % y_c: coordination after transformation (Slice coordination)
            % x_c: coordiniation before transformation (Volume coordination)
            % M: the transform from Volume coordinate system to Slice
            % coordinate sytem.
            value_temp_bw = value_temp;
%             value_temp_bw = imgaussfilt3(value_temp,[10,10,10]);
            value_temp_bw(value_temp<0.8*max(value_temp(:))) = 0;
            CC = bwconncomp(value_temp_bw);
            PeakNum = CC.NumObjects;
            Angles = [alpha(k),beta(j),gamma(i)];            
            for i_pk = 1:PeakNum
                % Get the Maximum intensity from each object
                [Intensity,i_pos] = max(value_temp_bw(CC.PixelIdxList{i_pk}));
                % Get the Maximum position from each object
                [y_c_temp_x,y_c_temp_z,y_c_temp_d] = ind2sub(CC.ImageSize,CC.PixelIdxList{i_pk}(i_pos));
                x_c(2,1) = d_list(y_c_temp_d);
                y_c = [y_c_temp_x;0;y_c_temp_z];
                y_c(1) = y_c(1)*stepX+data_slice_bw_low.x(1)-stepX;
                y_c(3) = y_c(3)*stepX+data_slice_bw_low.y(1)-stepX;
                % Get the translation vector (after rotation. from volune
                % to slice)
                Translation = (-x_c+y_c)';
                Angles_diff = PeakTable.Angles - Angles;
                af = (sum(Angles_diff.^2,2)<ANGLE_TOLERANCE2);
                Translation_diff = PeakTable.Translation(af,:) - Translation;
                tf = (sum(Translation_diff.^2,2) < TRANSLATION_TOLERANCE2);
                PeakTableLength = size(PeakTable,1);
                % Modify the PeakTable
                if sum(tf)==0
                    % If this peak does not exist, then add it to the table
                    PeakTable(PeakTableLength+1,:) = {Intensity,Angles,Translation};
                else
                    % If this peak already exists, then compare the intensities
                    index_list = 1:PeakTableLength;
                    index_tf = index_list(af);
                    index_tf = index_tf(tf);
                    Intensity_diff = Intensity-PeakTable.Intensity(index_tf);
                    intensf = (Intensity_diff>0);
                    index_change = index_tf(intensf);
                    % If the intensity is larger, then change it
                    if ~isempty(index_change)
                        PeakTable(index_change(1),:) = {Intensity,Angles,Translation};
                        if length(index_change)>1
                            % delete duplicated peaks
                            PeakTable(index_change(2:end),:) = [];
                        end
                    end
                end
                % Peaks.Corr, Peaks.M
            end
            % Rank the Peaks with intensity
            PeakTable = sortrows(PeakTable,'Intensity','descend');
            PeakTableLength = size(PeakTable,1);
            if PeakTableLength>MAX_PEAK_RECORD_NUM
                % Only keep the peaks with large intensities.
                PeakTable((MAX_PEAK_RECORD_NUM+1):end,:)=[];
            end
            %%
            [corr_now_max,ind_max] = max(value_temp(:));
            [y_c_temp_x,y_c_temp_z,y_c_temp_d] = ind2sub(size(value_temp),ind_max);
            % y_c_temp(x_dim,z_dim,y_dim)
            x_c(2,1) = d_list(y_c_temp_d);
            y_c = [y_c_temp_x;0;y_c_temp_z];
            y_c(1) = y_c(1)*stepX+data_slice_bw_low.x(1)-stepX;
            y_c(3) = y_c(3)*stepX+data_slice_bw_low.y(1)-stepX;
            t = -x_c+y_c;
            if corr_now_max > corr_max
                corr_max = corr_now_max;
                para_max = [alpha(k);beta(j);gamma(i);t];
            end
            M = cat(2,R,t); % Transformation from Volume Coodination to Slice Coodination
            % Visualization for test
            data_corr.value(:,:,i) = value_temp(:,:,y_c_temp_d);
            vd = v3 - d_list(y_c_temp_d);
            selected = (abs(vd)<Integ);
            px = v1(selected==1);
            py = v2(selected==1);
            pt_now_dist = vd(selected==1);
            pt_now_area = ones(1,sum(selected==1)).*exp(-(pt_now_dist./10).^2/2);
            data_now = neuroReg.renderCell2(x_lim,y_lim,stepX,[px;py],pt_now_area,Option);
            data_now.value = data_now.value - mean(data_now.value(~isnan(data_now.value(:))))*MagicNumber;
            data_now.value(isnan(data_now.value))=0;
            data_now_corr = data_slice_bw_low;
            data_now_corr.value = filter2(data_now.value,data_slice_bw_low.value);
            % Visualize current result
            figure(100861);
            subplot(1,2,1);
            neuroReg.plotData2(data_now);
            subplot(1,2,2);
            neuroReg.plotData2(data_now_corr);
            title(['d=',num2str(d_list(y_c_temp_d)),' Max\_XCC=',num2str(corr_now_max)]);
            %%
            i_tot = sub2ind([n_alpha,n_beta,n_gamma],k,j,i);
            a = toc;
            waitbar(i_tot/N,h1,[...
                'Current angle: ',num2str([alpha(k),beta(j),gamma(i)]),...
                '  Time remaining: ',datestr(a/i_tot*(N-i_tot)/24/3600,'DD HH:MM:SS')]);
        end
    end
end
TransParameters = para_max;
Corr = corr_max;
close(h1);


end



function vts = getCube ( origin, size )
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
vts(1,:,:) = x;
vts(2,:,:) = y;
vts(3,:,:) = z;
end
function drawCube(vts,c)
for i=1:6
    hold on
    plot3(vts(1,:,i),vts(2,:,i),vts(3,:,i),c);
    hold off
end
h=patch(vts(1,:,6),vts(2,:,6),vts(3,:,6),'y');
alpha(h,0.3);
end
