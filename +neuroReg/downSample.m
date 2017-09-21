function data_out = downSample(data,stepX,stepY,stepZ)
% downSample put data into a new grid. downSample(data,stepX,stepY,stepZ)
if ndims(data.value) == 2 %#ok<ISMAT>
    % Get range
    x1 = min(data.x);
    x2 = max(data.x);
    y1 = min(data.y);
    y2 = max(data.y);
    data_out.x = x1:stepX:x2;
    data_out.y = y1:stepZ:y2;
    [yy,xx] = meshgrid(data.y,data.x);
    [yyq,xxq] = meshgrid(data_out.y,data_out.x);
    data_out.value = interp2(yy,xx,data.value,yyq,xxq);
elseif ndims(data.value) == 3
    x1 = min(data.x);
    x2 = max(data.x);
    y1 = min(data.y);
    y2 = max(data.y);
    z1 = min(data.z);
    z2 = max(data.z);
    data_out.x = x1:stepX:x2;
    data_out.y = y1:stepY:y2;
    data_out.z = z1:stepZ:z2;
    [yyy,xxx,zzz] = meshgrid(data.y,data.x,data.z);
    [yyyq,xxxq,zzzq] = meshgrid(data_out.y,data_out.x,data_out.z);
    data_out.value = interp3(yyy,xxx,zzz,data.value,yyyq,xxxq,zzzq);
end
