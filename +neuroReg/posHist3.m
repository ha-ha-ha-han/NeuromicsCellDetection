function value_temp = posHist3(pt_list_slice3,pt_list_vol,Grid,Option)
% posHist3 creates a correlated histogram from pt_list_slice3 to
% pt_list_vol.
% INPUT
% pt_list_slice3: point list 1.
% pt_list_vol: point list 2.
% Grid: specify the grid for value_temp. Grid.x, Grid.y, Grid.z
% Option: Option.Hist3Smooth. Smoothing method. 
% 'box', 'gaussian' or 'none'
% OUTPUT
% value_temp: 3D histogram. Nx-by-Ny-by-Nz
p1 = pt_list_slice3;
n1 = size(p1,2);
p1 = reshape(p1,[3,n1,1]);
p2 = pt_list_vol;
n2 = size(p2,2);
p2 = reshape(p2,[3,1,n2]);
pd = p1 - p2; % The translation is from volume to slice
pd = reshape(pd,[3,n1*n2]);

stepX = mean(diff(Grid.x));
stepY = mean(diff(Grid.y));
stepZ = mean(diff(Grid.z));

binX = [Grid.x(:)-stepX/2;Grid.x(end)+stepX/2];
binY = [Grid.y(:)-stepY/2;Grid.y(end)+stepY/2];
binZ = [Grid.z(:)-stepZ/2;Grid.z(end)+stepZ/2];
count = neuroReg.histcn(pd',binX,binY,binZ);

switch Option.Hist3Smooth
    case 1 % box
        value_temp = smooth3(count,'box',Option.Hist3SmoothSize);
    case 2 % gaussian
        value_temp = smooth3(count,'gaussian',Option.Hist3SmoothSize);
    otherwise
        value_temp = count;
end



end

