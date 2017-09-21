function plotData2(data)
% plotSlice(data) visualize a 2-D data set
pcolor(data.x,data.y,data.value');
axis equal;
axis tight;
shading flat;
end

