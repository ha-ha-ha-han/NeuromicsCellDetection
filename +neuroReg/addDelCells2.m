function [pt_list,pt_area] = addDelCells2(h,pt_list,pt_area)
% addDelCells2 modify cell list interactively
%% Delete Cells
ht = get(h,'Title');
titleString = ht.String;
title(h,[titleString,' || Deleting Cells...']);

[pts_x,pts_y] = getpts(h);
h1 = findobj(h);
for i = 1:length(h1)
    if isa(h1(i),'matlab.graphics.chart.primitive.Scatter')
        delete(h1(i));
    end
end
for i = 1:length(pts_x)
    xd = pt_list(1,:)-pts_x(i);
    yd = pt_list(2,:)-pts_y(i);
    d = xd.^2 + yd.^2;
    [~,i_mind]=min(d);
    pt_list = cat(2,pt_list(:,1:i_mind-1),pt_list(:,(i_mind+1):end));
    pt_area = cat(2,pt_area(:,1:i_mind-1),pt_area(:,(i_mind+1):end));
end
hold on
h2 = scatter(pt_list(1,:),pt_list(2,:),pt_area./max(pt_area)*50,'r');

%% Adding Cells
title(h,[titleString,' || Adding Cells...']);
[pts_x,pts_y] = getpts(h);
n = length(pts_x);
pts = [pts_x,pts_y]';
pt_list = cat(2,pt_list,pts);
pt_area = [pt_area(:);300*ones(n,1)];
title(h,titleString);
scatter(pts(1,:),pts(2,:),50*300/max(pt_area),'g');
hold off
