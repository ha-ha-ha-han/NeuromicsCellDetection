function data_out = renderCell3(x_lim,y_lim,d_lim,stepX,stepD,pt_list,Option)
% renderCell3 renders a volume from the pt_list (3-by-N)
% It tries to simulate the ZStack of the neuron image.
% x_lim ([x1, x2]),y_lim ([y1,y2]),d_lim ([d1,d2]) specifies the cube boundaries.
% stepD specifies the pixel step of d-direction.
% pt_list: 3-by-N, position of the cells.
% Option: 
% Option.Integ: thickness of each X-Y slice.
% Option.CellRadius: radius of a cell
% Option.MagicNumber: value_temp = value_temp - mean(value_temp(~nf))*MagicNumber;
% To create a punishment for mismatching (This may not make sense)
% If you use a for-loop, renderCell2 is (nearly) equivalent to renderCell3


if isfield(Option,'Integ')
    Integ = Option.Integ; % Integ: slice thickness
else
    Integ = 10;
end
if isfield(Option,'CellRadius')
    CellRadius = Option.CellRadius; % Integ: slice thickness
else
    CellRadius = 6;
end
if isfield(Option,'MagicNumber')
    MagicNumber = Option.MagicNumber; % Integ: slice thickness
else
    MagicNumber = 2;
end
% renderCell3 renders a data cube from the point list
v1 = pt_list(1,:); % x-direction of the slice
v2 = pt_list(2,:); % up-down-direction of the slice
v3 = pt_list(3,:); % normal direction of the slice
% Get the size of the rotated cube
data_out.x = x_lim(1):stepX:x_lim(2);
data_out.y = y_lim(1):stepX:y_lim(2);
data_out.z = d_lim(1):stepD:d_lim(2);
nx = length(data_out.x);
ny = length(data_out.y);
nz = length(data_out.z);
% Get the distances from planes to cells
vd = abs(v3 - data_out.z'); % This result in a matrix. Column: each cell. Row: distance to d-plane
value = zeros(nx,ny,nz);
for i = 1:nz
    selected = (abs(vd(i,:))<Integ);
    px = v1(selected==1);
    py = v2(selected==1);
    pt_now_dist = vd(i,selected==1);
    pt_now_area = ones(1,sum(selected==1)).*exp(-(pt_now_dist./CellRadius).^2/2);
    data_now = neuroReg.renderCell2(x_lim,y_lim,stepX,[px;py],pt_now_area,Option);
    value_temp = data_now.value;
    nf = isnan(value_temp);
    value_temp = value_temp - mean(value_temp(~nf))*MagicNumber;
    value_temp(nf)=0;
    if sum(isnan(value_temp))>0
        1;
    end
    value(:,:,i) = value_temp;
end
data_out.value = value;
end

