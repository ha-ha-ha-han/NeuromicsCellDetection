function sqCorrection(obj,varargin)%Square Correction for correcting the 
%blank border 
switch obj.Direction
    case 'x'
        d=1;
        xData=obj.Data.y;
        yData=obj.Data.z;
    case 'y'
        d=2;
        xData=obj.Data.x;
        yData=obj.Data.z;
    otherwise
        d=3;
        xData=obj.Data.x;
        yData=obj.Data.y;
end
pos=get(obj.Slider,'Value');
[~,posInd]=min(abs(obj.AxisData-pos));
index=arraySliceIndex(3,posInd,d);
sliceData=squeeze(obj.Data.value(index{:}));
[~,sizeOfxData]=size(xData);
[~,sizeOfyData]=size(yData);
disp('Please select the points on the  left side.');
[m1,n1]=getpts;%To get the points on the left border
%disp([num2str(m1),   num2str(n1)]);
[numberOfPoints,~]=size(m1);
xP=NaN.*ones(numberOfPoints,1);
yP=NaN.*ones(numberOfPoints,1);
for k=1:numberOfPoints
temp1=(m1(k)-xData(1))^2;
temp2=(n1(k)-yData(1))^2;
xP(k)=1;
yP(k)=1;
    for i=1:sizeOfxData
         if ((m1(k)-xData(i))^2<temp1)
             temp1=(m1(k)-xData(i))^2;   
             xP(k)=i;
         end
    end
    for j=1:sizeOfyData
        if ((n1(k)-yData(j))^2<temp2)
            temp2=(n1(k)-yData(j))^2;
            yP(k)=j;
        end
    end
end%To get the position in the matrix of the left boundary points 
%disp([num2str(xP),   num2str(yP)]);
A1=polyfit(yP,xP,1);%Linear fit of the left boundary
disp('Please select the points on the right side.');
[m2,n2]=getpts;%To get the points on the right border
%disp([num2str(m2),   num2str(n2)]);
[numberOfPoints,~]=size(m2);
xP=NaN.*ones(numberOfPoints,1);
yP=NaN.*ones(numberOfPoints,1);
for k=1:numberOfPoints
temp1=(m2(k)-xData(1))^2;
temp2=(n2(k)-yData(1))^2;
xP(k)=1;
yP(k)=1;
    for i=1:sizeOfxData
         if ((m2(k)-xData(i))^2<temp1)
                temp1=(m2(k)-xData(i))^2;
                xP(k)=i;
         end
    end
    for j=1:sizeOfyData
        if ((n2(k)-yData(j))^2<temp2)
            temp2=(n2(k)-yData(j))^2;
            yP(k)=j;
        end
    end
end%To get the position in the matrix of the right boundary points
%disp([num2str(xP),   num2str(yP)]);
A2=polyfit(yP,xP,1);%Linear fit of the right boundary
newSliceData=NaN.*ones(sizeOfxData,sizeOfyData);
ll=NaN.*ones(sizeOfyData);
rl=NaN.*ones(sizeOfyData);
for j=1:sizeOfyData
    ll(j)=polyval(A1,j);%left limit of each line j
    rl(j)=polyval(A2,j);%right limit of each line j
    if ll(j)>=rl(j)
        errordlg('Left side value is bigger than right side.');
        disp(['ll(j)=' num2str(ll(j)), '  rl(j)=' num2str(rl(j))]);
        return;
    end
    if ll(j)<1
        ll(j)=1;
    end
    if rl(j)>sizeOfxData
        rl(j)=sizeOfxData;
    end
    ll(j)=round(ll(j));
    rl(j)=round(rl(j));
    xx(ll(j):rl(j),j)=((ll(j):rl(j))-ll(j))*(sizeOfxData-1)/(rl(j)-ll(j))+1;
    newSliceData(round(((ll(j):rl(j))-ll(j))*(sizeOfxData-1)/(rl(j)-ll(j))+1),j)=sliceData((ll(j):rl(j)),j);
end
%figure;
%pcolor(xData,yData,newSliceData');
%shading interp;
%% Correction for one cut
boolNewSliceData=~isnan(newSliceData);
count=sum(sum(boolNewSliceData));
%disp(num2str(count));
x=1:count;
y=1:count;
v=1:count;
count=0;
for j=1:sizeOfyData
    il=round(polyval(A1,j));
    if il<1
        il=1;
    end
    for i=1:sizeOfxData  
        if ~isnan(newSliceData(i,j))
            count=count+1;
            x(count)=(xx(il,j)-1)*(xData(sizeOfxData)-xData(1))/(sizeOfxData-1)+xData(1);
            y(count)=yData(j);
            v(count)=newSliceData(i,j);
            il=il+1;
        end
    end
end
[X,Y,V]=griddata(x,y,v,xData',yData);
figure;
pcolor(X,Y,V);
shading interp;
drawnow;
%% Correction for all cuts
data=obj.Data;
if (d==1)
    [~,columnsize]=size(obj.Data.x);
    for k=1:columnsize
        sliceData=squeeze(obj.Data.value(k,:,:));
        for j=1:sizeOfyData
            newSliceData(round(((ll(j):rl(j))-ll(j))*(sizeOfxData-1)/(rl(j)-ll(j))+1),j)=sliceData((ll(j):rl(j)),j);
        end
        boolNewSliceData=~isnan(newSliceData);
        count=sum(sum(boolNewSliceData));
        %disp(num2str(count));
        x=1:count;
        y=1:count;
        v=1:count;
        count=0;
        for j=1:sizeOfyData
            il=round(polyval(A1,j));
            if il<1
                il=1;
            end
            for i=1:sizeOfxData  
                if ~isnan(newSliceData(i,j))
                    count=count+1;
                    x(count)=(xx(il,j)-1)*(xData(sizeOfxData)-xData(1))/(sizeOfxData-1)+xData(1);
                    y(count)=yData(j);
                    v(count)=newSliceData(i,j);
                    il=il+1;
                end
            end
        end
        [~,~,T]=griddata(x,y,v,xData',yData);
        if (sum(size(T))==0)
            disp(k);
            Value1(k,:,:)=NaN.*ones(sizeOfxData,sizeOfyData);
            break;
        end
        Value1(k,:,:)=T';
    end
    data.value=Value1;
    data.x=obj.Data.x;
    data.y=X(1,:);
    data.z=Y(:,1)';
end
if (d==2)
    [~,columnsize]=size(obj.Data.y);
    for k=1:columnsize
        sliceData=squeeze(obj.Data.value(:,k,:));
        for j=1:sizeOfyData
            newSliceData(round(((ll(j):rl(j))-ll(j))*(sizeOfxData-1)/(rl(j)-ll(j))+1),j)=sliceData((ll(j):rl(j)),j);
        end
        boolNewSliceData=~isnan(newSliceData);
        count=sum(sum(boolNewSliceData));
        %disp(num2str(count));
        x=1:count;
        y=1:count;
        v=1:count;
        count=0;
        for j=1:sizeOfyData
            il=round(polyval(A1,j));
            if il<1
                il=1;
            end
            for i=1:sizeOfxData  
                if ~isnan(newSliceData(i,j))
                    count=count+1;
                    x(count)=(xx(il,j)-1)*(xData(sizeOfxData)-xData(1))/(sizeOfxData-1)+xData(1);
                    y(count)=yData(j);
                    v(count)=newSliceData(i,j);
                    il=il+1;
                end
            end
        end
        [~,~,T]=griddata(x,y,v,xData',yData);
 %       disp(num2str(size(T)));
        if (sum(size(T))==0)
            Value2(k,:,:)=NaN.*ones(sizeOfxData,sizeOfyData);
            disp(k);
            break;
        end
        Value2(k,:,:)=T';
    end
    data.value=permute(Value2,[2,1,3]);
    data.x=X(1,:);
    data.y=obj.Data.y;
    data.z=Y(:,1)';
end
%%
% data=obj.Data;
% boolNewDataValue=~isnan(newDataValue);
% count=sum(sum(sum(boolNewDataValue)));
% x=1:count;
% y=1:count;
% z=1:count;
% v=1:count;
% count=0;
% if(d==1)
%     [~,columnsize]=size(obj.Data.x);
%     for i=1:sizeOfxData
%         for j=1:sizeOfyData
%             for k=1:columnsize
%             if ~isnan(newDataValue(k,i,j))
%                 count=count+1;
%                 x(count)=obj.Data.x(k);
%                 y(count)=xData(i);
%                 z(count)=yData(j);
%                 v(count)=newDataValue(k,i,j);
%             end
%             end
%         end
%     end
%     F = scatteredInterpolant(x,y,z,v);
%     V = F(x,y,z);
%    % V=griddata(x,y,z,v,linspace(obj.Data.x(1),obj.Data.x(columnsize)),linspace(xData(1),xData(sizeOfxData)),linspace(yData(1),yData(sizeOfyData)));
% end
% if(d==2)
%     [~,columnsize]=size(obj.Data.y);
%     for i=1:sizeOfxData
%         for j=1:columnsize
%             for k=1:sizeOfyData
%             if ~isnan(newDataValue(i,j,k))
%                 count=count+1;
%                 x(count)=xData(i);
%                 y(count)=obj.Data.y(j);
%                 z(count)=yData(k);
%                 v(count)=newDataValue(i,j,k);
%             end
%             end
%         end
%     end
%  %   F = scatteredInterpolant(x,y,z,v);
%  %   V = F(x,y,z);
% %    V=griddata(x,y,z,v,linspace(xData(1),xData(sizeOfxData)),linspace(obj.Data.y(1),obj.Data.y(columnsize)),linspace(yData(1),yData(sizeOfyData)));
% end 
SaveName=inputdlg({'Data Name'},...
'Square Correction',...
1,...
{[obj.DataName,'_sq']});
if isempty(SaveName)
    return;
end
assignin('base',SaveName{1},data);
end