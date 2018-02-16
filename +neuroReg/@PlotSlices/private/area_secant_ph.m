function [x,y,r] = area_secant_ph( x_grid , y_grid, x0, y0, x1, y1 ) 
%    Area_secant_ph returns the coodinates of points of a secant of a
%rectangular area defined by x_grid and y_grid.

%% determine the line
k=(y1-y0)/(x1-x0);
reverse_flag=0;
if k>10e4   % in case the line is vertical to x-axis
    [x_grid,y_grid]=exchange_p(x_grid,y_grid);
    [x0,y0]=exchange_p(x0,y0);
    [x1,y1]=exchange_p(x1,y1);
    k=(y1-y0)/(x1-x0);
    reverse_flag=1;
end
b=y1-x1*k;

%% determine the range of output x array

y_max=max(y_grid);
y_min=min(y_grid);
x_max=max(x_grid);
x_min=min(x_grid);

x1=(y_max-b)/k;
x2=(y_min-b)/k;
if x1>x2
    [x1,x2]=exchange_p(x1,x2);
end

if x2<x_min||x1>x_max
    errordlg('hehe, the line is out of the area');
    return;
end

if x1<x_min
    x1=x_min;
    y1=k*x1+b;
else
    y1=y_min;
end

if x2>x_max
    x2=x_max;
    y2=k*x2+b;
else
    y2=y_max;
end



% determine the number of sampling points
r=sqrt((x1-x2)^2+(y1-y2)^2);
n=floor(abs(r/(min(x_grid(1)-x_grid(2),y_grid(1)-y_grid(2)))))+1;

%% determine the output array x and array y
x=linspace(x1,x2,n);
y=k*x+b;
sgn=(x-x0)./abs(x-x0);
sgn(isnan(sgn))=0;
r=sgn.*sqrt((x-x0).^2+(y-y0).^2);
if reverse_flag==1
    [x,y]=exchange_p(x,y);
end


function [a,b] = exchange_p(c,d)  % just exchange
    a=d;
    b=c;
  
    


