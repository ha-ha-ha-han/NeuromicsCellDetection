function [I,R] = data2img(data)
% This function deals with the stupid image storage convention.
% imshow(I) will result in the same image as pcolor(data.value')
% I == flipud(data.value')
% R specifies picture position in the real world.

% px2um_x = mean(diff(data.x));
% px2um_y = mean(diff(data.y));
% if (px2um_x - px2um_y)>1e-4
%     error('x and y pixels should be of the same scale');
% end
% px2um = px2um_x;
I = flipud(data.value');
R = imref2d(size(I),[data.x(1),data.x(end)],-[data.y(end),data.y(1)]);

end

