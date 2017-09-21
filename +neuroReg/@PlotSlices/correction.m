function correction(obj,varargin)
% Set the Upper and Lower limit. It is changeable, but now
% disabled.
%            tll=1;  %Tolerance lower limit.
%Max of tolerance is 1 and zero tolerance is forbidden.
%Choose this value in between.
%            tul=1;  %Tolerance upper limit.
%Max of tolerance is 1 and zero tolerance is forbidden.
%Choose this value in between.

%modified by sandy for photon dep data

%% Gather information
disp(['Please process Square correction first if needed.'...
    'Along x or y is permitted.'...
    'Please select at least three points.',...
    'Data selection end with Enter.']);
switch obj.Direction
    case 'x'
        d=1;
        xData=obj.Data.y;
        yData=obj.Data.z;
        %sandy modified
        if isfield(obj.Data,'ztot')
            ytotData=obj.Data.ztot;
        end
        %sandy end modified
    case 'y'
        d=2;
        xData=obj.Data.x;
        yData=obj.Data.z;
        %sandy modified
        if isfield(obj.Data,'ztot')
            ytotData=obj.Data.ztot;
        end
        %sandy end modified
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
[m,n]=getpts;
hold on;
correctionPointHandle=plot(m,n,'k+');
% can give bigger order or not?
A=polyfit(m,n,2);
z=polyval(A,xData);
correctionLineHandle=plot(xData,z,'r--');
hold off;
rangeOfyData=max(yData)-min(yData);%range of yData

%% positive quadratic term coefficient
if (A(1,1)>0) %positive quadratic term coefficient
    xd=xData;
    pv=polyval(A,xd);
    fitmin1=min(pv);
    fitmax1=max(pv);
    newyDataDimension=ceil((fitmax1-fitmin1)*sizeOfyData/rangeOfyData)+sizeOfyData;%new yData dimension
    nsliceData=NaN.*ones(sizeOfxData,newyDataDimension);
    for i=1:sizeOfxData
        dif=fix((polyval(A,xd(1,i))-fitmin1)*sizeOfyData/rangeOfyData);%difference
        nsliceData(i,newyDataDimension-dif-sizeOfyData+(1:sizeOfyData))=sliceData(i,1:sizeOfyData);
    end
    if (d==1)
        [~,columnsize]=size(obj.Data.x);
        newDataValue=NaN.*ones(columnsize,sizeOfxData,newyDataDimension);
        for i=1:sizeOfxData
            dif=fix((polyval(A,xd(1,i))-fitmin1)*sizeOfyData/rangeOfyData);%difference
            newDataValue(1:columnsize,i,newyDataDimension-dif-sizeOfyData+(1:sizeOfyData))=obj.Data.value(1:columnsize,i,1:sizeOfyData);
        end
    end
    if (d==2)
        [~,columnsize]=size(obj.Data.y);
        newDataValue=NaN.*ones(sizeOfxData,columnsize,newyDataDimension);
        for i=1:sizeOfxData
            dif=fix((polyval(A,xd(1,i))-fitmin1)*sizeOfyData/rangeOfyData);%difference
            newDataValue(i,1:columnsize,newyDataDimension-dif-sizeOfyData+(1:sizeOfyData))=obj.Data.value(i,1:columnsize,1:sizeOfyData);
        end
    end
end %end of positive quadratic term coefficient case

%% negative quadratic term coefficient
if (A(1,1)<0)%negative quadratic term coefficient
    xd=xData;
    pv=polyval(A,xd);
    fitmin2=min(pv);
    fitmax2=max(pv);
    newyDataDimension=ceil((fitmax2-fitmin2)*sizeOfyData/rangeOfyData)+sizeOfyData;
    nsliceData=NaN.*ones(sizeOfxData,newyDataDimension);
    for i=1:sizeOfxData
        dif=fix((fitmax2-polyval(A,xd(1,i)))*sizeOfyData/rangeOfyData);%difference
        nsliceData(i,dif+(1:sizeOfyData))=sliceData(i,1:sizeOfyData);
    end
    if (d==1)
        [~,columnsize]=size(obj.Data.x);
        newDataValue=NaN.*ones(columnsize,sizeOfxData,newyDataDimension);
        for i=1:sizeOfxData
            dif=fix((fitmax2-polyval(A,xd(1,i)))*sizeOfyData/rangeOfyData);%difference
            newDataValue(1:columnsize,i,dif+(1:sizeOfyData))=obj.Data.value(1:columnsize,i,1:sizeOfyData);
        end
    end
    if (d==2)
        [~,columnsize]=size(obj.Data.y);
        newDataValue=NaN.*ones(sizeOfxData,columnsize,newyDataDimension);
        for i=1:sizeOfxData
            dif=fix((fitmax2-polyval(A,xd(1,i)))*sizeOfyData/rangeOfyData);%difference
            newDataValue(i,1:columnsize,dif+(1:sizeOfyData))=obj.Data.value(i,1:columnsize,1:sizeOfyData);
        end
    end
end %end of negative quadratic term coefficient case
boolNewDataValue=~isnan(newDataValue);
flagU=1;
Testu=NaN.*ones(1,newyDataDimension);
for ciru=newyDataDimension:-1:1
    Testu(1,ciru)=sum(sum(boolNewDataValue(:,:,ciru)));
    if sum(sum(boolNewDataValue(:,:,ciru)))~=0
        newUpperLimit=ciru;
        flagU=0;
        break;
    end
end
%         disp(Testu(1,:))
flagL=1;
for cirl=1:newyDataDimension
    %             disp(sum(sum(boolNewDataValue(:,:,cirl))))
    if sum(sum(boolNewDataValue(:,:,cirl)))~=0
        newLowerLimit=cirl;
        flagL=0;
        break;
    end
end
obj.CorrData.x=obj.Data.x;
obj.CorrData.y=obj.Data.y;
if (flagL+flagU~=0)
    obj.CorrData.value=newDataValue;
    nyData=NaN.*ones(1,newyDataDimension);
    for i=1:newyDataDimension
        nyData(1,i)=min(yData)+(i-1)*rangeOfyData/sizeOfyData;%new yData
        %sandy modified
        if isfield(obj.Data,'ztot')
            nytotData(i,:) = min(ytotData)+(i-1)*rangeOfyData/sizeOfyData;%#ok<AGROW> %new ytotData
        end
        %sandy end modified
    end
    obj.CorrData.z=nyData;
    %sandy modified
    if isfield(obj.Data,'ztot')
        obj.CorrData.ztot=nytotData;
    end
    %sandy end modified
end
if (flagL+flagU==0)
    if(d==1)
        [~,columnsize]=size(obj.Data.x);
        obj.CorrData.value(1:columnsize,1:sizeOfxData,1:(newUpperLimit-newLowerLimit+1))=newDataValue(1:columnsize,1:sizeOfxData,newLowerLimit:newUpperLimit);
    end
    if(d==2)
        [~,columnsize]=size(obj.Data.y);
        obj.CorrData.value(1:sizeOfxData,1:columnsize,1:(newUpperLimit-newLowerLimit+1))=newDataValue(1:sizeOfxData,1:columnsize,newLowerLimit:newUpperLimit);
    end
    nyData=NaN.*ones(1,(newUpperLimit-newLowerLimit+1));
    for i=1:(newUpperLimit-newLowerLimit+1)
        nyData(1,i)=min(yData)+(i-1)*rangeOfyData/sizeOfyData;%new yData
        %sandy modified
        if isfield(obj.Data,'ztot')
            nytotData(i,:) = min(ytotData)+(i-1)*rangeOfyData/sizeOfyData;%new ytotData
        end
        %sandy end modified
    end
    obj.CorrData.z=nyData;
    %sandy modified
    if isfield(obj.Data,'ztot')
        obj.CorrData.ztot=nytotData;
    end
    %sandy end modified
end

[flag,DataName] = saveCorrection(obj);
if flag
    if(d==1)
        PlotSlices(DataName,'x');
    end
    if(d==2)
        PlotSlices(DataName,'y');
    end
end
delete(correctionLineHandle);
delete(correctionPointHandle);

end

function [flag,DataName]=saveCorrection(obj,varargin)
flag=0;
newData=obj.Data;
newData.x=obj.CorrData.x;
newData.y=obj.CorrData.y;
newData.z=obj.CorrData.z;
newData.value=obj.CorrData.value;
%sandy modified
if isfield(obj.Data,'ztot')
    newData.ztot=obj.CorrData.ztot;
end
%sandy end modified
if length(newData.x)==1
    newData.x=newData.y;
    newData.y=newData.z;
    newData=rmfield(newData,'z');
    newData.value=squeeze(newData.value);
end
DataName=inputdlg({'Data Name'},...
    'Save correction',...
    1,...
    {[obj.DataName,'_corr']});
if isempty(DataName)
    return;
end
assignin('base',DataName{1},newData);
flag=1;
end