% +NEUROREG Registration Toolbox for Neuromics project.
% Version 0.1 30-Aug-2017
%
% SET OPTIONS
%   setOption                       - For image reconstruction from cell list
% REGISTRATION ALGORITHMS
%   rotationCorr3                   - rotationCorr3 returns the transformation parameters and the corresponding
%   rotationHist3                   - rotationHist3 returns the transformation parameters and the corresponding
% HELPER FUNCTIONS FOR REGISTRATION
%   histcn                          - function [count edges mid loc] = histcn(X, edge1, edge2, ..., edgeN)
%   xcorrFFT                        - xcorrFFT provides a fast cross correlation calculation.
%   posHist3                        - posHist3 creates a correlated histogram from pt_list_slice3 to
% VISUALIZATION TOOLS:
%   plotData2                       - plotSlice(data) visualize a 2-D data set
%   plotTransform                   - plotTransform visualize the registration from dataZ to data_slice.
%   VisTransform                    - VisualizeTransform creates a GUI that display the transform from registration.
%   PlotSlices                      - Visualizes 3D data stack (and optionally marks the detected cells.
% ROTATION, TRANSLATION & GEOMETRY TOOLS
%   cutVolume                       - cutVolume cut a slice from the data.
%   renderCell2                     - xlim and ylim are of the form [lb ub] to specify the render range.
%   renderCell3                     - renderCell3 renders a volume from the pt_list (3-by-N)
%   rotateCells                     - rotateCells is very useful for:
%   rotateStack                     - rotateStack is used to rotate a ZStack data 
%   getROI                          - Choose three points
% CELL DETECTION TOOLS 
%   addDelCells2                    - addDelCells2 modify cell list interactively
%   detectCells2                    - detectCells2 detects cells from a slice data.
%   detectCells3                    - detectCells3 detects cell positions from a ZStack data.
% IMAGE PROCESSING TOOLS
%   bwCell2                         - bwCell2 binarize a slice data.
%   bwCell3                         - Use difference of Gaussian detector
%   downSample                      - downSample put data into a new grid. downSample(data,stepX,stepY,stepZ)
%   medfilt3                        - neuroReg.medfilt3 is simply a replacement of the default medfilt3.
%   data2img                        - This function deals with the stupid image storage convention.
% LOADING TOOLS
%   combineZStack                   - combineZStack combines data1 and data2 by matching data points
%   loadSlice                       - loadSlice loads raw slice data. It requires no input for now.
%   loadTiff                        - loadTiff loads data from a TIff file.
%   loadZStack                      - loadZStack loads raw slice data. It requires no input for now.
%   loadZStackTiff                  - loadZStack loads raw slice data. It requires no input for now.
% GEOMETRIC REGISTRATION METHOD (NOT SECCESSFUL YET)
%   PointCloudRegister2             - [M, Error] = PointCloudRegister(y, x, M0, DistScale, Options)
%   pointCloudTriangleCoplanarMatch - pt_list1 (2 or 3-by-n) is used to pick points randomly.
%   pointCloudTriangleMatch         - pt_list1 (2 or 3-by-n) is used to pick points randomly.
%   registerSixPoints               - This function calculate the M that let the equation 
% BACKUP FILES
%   rotationCorr3_20170727          - rotationCorr3 returns the transformation parameters and the corresponding





