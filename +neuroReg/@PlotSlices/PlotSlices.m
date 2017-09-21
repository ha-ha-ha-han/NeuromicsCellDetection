classdef PlotSlices < handle
    % PlotSlices visualizes 3D data stack.
    %   obj = neuroReg.PlotSlices(DataName,Direction)
    %   obj = neuroReg.PlotSlices(DataName,Direction,pt_list,pt_area_list,integ_range)
    %   - DataName (string): name of the variable. See the note below. 
    %   - Direction (string):'x', 'y' or 'z', specifying the normal
    %     direction of the cut.
    %   For the second usage, the visualization overlaps the stack with the
    %   detected cell positions, which are described by the adiitional
    %   input.
    %   - pt_list (3-by-n matrix): position of the detected cells.
    %   - pt_area_list (1-by-n array): size of the cells. volume in 3D 
    %     case.
    %   - integ_range (scalar): the range within which the point_list
    %     will be projected to the cut.
    %   NOTE 1: 
    %       The variable specified by DataName should be a structure. 
    %       Type 'help +neuroReg' for more information about  the format of  
    %       the data structure. 
    %   NOTE 2:
    %       The variable should be stored in the 'base' workspace. If 
    %       calling neuroReg.PlotSlices inside a function, use the example
    %       below.
    %   EXAMPLE:
    %       % Visualize data
    %       assignin('base','data_temp_749',data);
    %       neuroReg.PlotSlices('data_temp_749','z',pt_list,pt_area,25);    
    
    properties
        Direction
        Data
        DataName
        AxisData % Store the axis that the plot is along with
        CorrData
        DetectedPoints % 3-by-n matrix
        DetectedPointsArea % 1-by-n array
        IntegralInterval % scalar, unit is the same with obj.X/Y/Z
        
        Figure
        Axis
        SurfaceObj
        Checkbox
        Slider
        SliderMin
        SliderMax
        SliderGamma
        FigureEDC
        FigureMDC
        MenuMDCEDC
        MenuSaveMDCEDC
        ColorMap
        FlipCheckbox
        Scatter % Cell detection plot
        
        Position
        PositionIndex
        MinValue
        MaxValue
        GammaValue
        PosX
        PosY
        VLine
        HLine
        Markers=[]
    end
    methods
        function obj=PlotSlices(DataName,Direction,pt_list,pt_area,integral)
            if iscell(DataName)
                DataName=DataName{1};
            end
            obj.DataName=DataName;
            obj.Direction = Direction;
            obj.Data=evalin('base',DataName);
            obj.DetectedPoints = pt_list;
            obj.DetectedPointsArea = pt_area;
            obj.IntegralInterval = integral;
            if ~isfield(obj.Data,'value')
                errordlg('Check data.value');
                return
            end
            switch  ndims(obj.Data.value)
                case 2
                    obj.data2Dto3D;
                case 4
                    obj.data4Dto3D;
                case 3
                    0;
                otherwise
                    errordlg('Check data dimension');
                    return;
            end
            obj.AxisData=obj.Data.(obj.Direction);
            obj.figureConstruct(DataName);
            obj.plot;
            obj.setUIMenu;
            obj.PosX=64;
            obj.PosY=0;
        end
        function data2Dto3D(obj)
            obj.Direction='x';
            [sizeY,sizeZ]=size(obj.Data.value);
            valueNew=zeros(1,sizeY,sizeZ);
            valueNew(1,:,:)=obj.Data.value;
            
            data=obj.Data;
            data.x=1;
            data.y=obj.Data.x;
            data.z=obj.Data.y;
            data.value=valueNew;
            obj.Data=data;
            
        end
        function data4Dto3D(obj)
            return; % to be added
        end
        function figureConstruct(obj,DataName)
            %% Frames
            obj.Figure = figure('Name',DataName,'MenuBar','none');
            homeVBox=uiextras.VBox('Parent',obj.Figure,...
                'BackgroundColor','w');
            infoHBox=uiextras.HBox('Parent',homeVBox);
            colorHBox=uiextras.HBox('Parent',homeVBox);
            colorleftVBox=uiextras.VBox('Parent',colorHBox);
            colorrightVBox=uiextras.VBox('Parent',colorHBox);
            colorMinHBox=uiextras.HBox('Parent',colorleftVBox);
            colorMaxHBox=uiextras.HBox('Parent',colorleftVBox);
            colorGammaHBox=uiextras.HBox('Parent',colorrightVBox);
            colormapHBox=uiextras.HBox('Parent',colorrightVBox);
            obj.Axis=axes('Parent',homeVBox);
            set(homeVBox,'Sizes',[25,50,-1]);
            set(colorHBox,'Sizes',[-1,-1]);
            set(colorleftVBox,'Sizes',[-1,-1]);
            set(colorrightVBox,'Sizes',[-1,-1]);
            axisMin=min(obj.AxisData);
            axisMax=max(obj.AxisData);
            step=1/length(obj.AxisData);
            
            %% uicontrols for slice control
            uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Text',...
                'String',       ['Along ',obj.Direction],...
                'BackgroundColor','w');
            obj.Checkbox=uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Checkbox',...
                'String',       'Ax Equal',...
                'Callback',     @obj.plot);
            obj.Slider=uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Slider',...
                'Min',          axisMin,...
                'Max',          axisMax,...
                'Value',        mean([axisMin,axisMax]),...
                'SliderStep',   [step,5*step]);
            addlistener(obj.Slider,'Value','PostSet',@obj.plot);
            uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Text',...
                'String',       'Pos');
            obj.Position=uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Edit',...
                'String',       'N/A',...
                'BackgroundColor','w',...
                'Callback',     @obj.setPos);
            uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Text',...
                'String',       'Ind');
            uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'PushButton',...
                'String',       '<',...
                'Callback',     {@obj.setInd,'-'});
            obj.PositionIndex=uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'Edit',...
                'String',       'N/A',...
                'BackgroundColor','w',...
                'Callback',     {@obj.setInd,'='});
            uicontrol(...
                'Parent',       infoHBox,...
                'Style',        'PushButton',...
                'String',       '>',...
                'Callback',     {@obj.setInd,'+'});
            set(infoHBox,'Sizes',[50,70,-1,20,50,20,20,40,20]);
            
            %% uicontrols for color correction
            uicontrol(...
                'Parent',       colorMinHBox,...
                'Style',        'Text',...
                'String',       'Min(%)');
            obj.MinValue=uicontrol(...
                'Parent',       colorMinHBox,...
                'Style',        'Edit',...
                'String',       '0',...
                'BackgroundColor','w',...
                'Callback',     @obj.sliderinput);
            obj.SliderMin=uicontrol(...
                'Parent',       colorMinHBox,...
                'Style',        'Slider',...
                'Min',          -1,...
                'Max',          1,...
                'Value',        -1);
            addlistener(obj.SliderMin,'Value','PostSet',@obj.plot);
            set(colorMinHBox,'Sizes',[40,55,-1]);
            
            uicontrol(...
                'Parent',       colorMaxHBox,...
                'Style',        'Text',...
                'String',       'Max(%)');
            obj.MaxValue=uicontrol(...
                'Parent',       colorMaxHBox,...
                'Style',        'Edit',...
                'String',       '100',...
                'BackgroundColor','w',...
                'Callback',     @obj.sliderinput);
            obj.SliderMax=uicontrol(...
                'Parent',       colorMaxHBox,...
                'Style',        'Slider',...
                'Min',          -1,...
                'Max',          1,...
                'Value',        0);
            addlistener(obj.SliderMax,'Value','PostSet',@obj.plot);
            set(colorMaxHBox,'Sizes',[40,55,-1]);
            
            uicontrol(...
                'Parent',       colorGammaHBox,...
                'Style',        'Text',...
                'String',       'Gamma');
            obj.GammaValue=uicontrol(...
                'Parent',       colorGammaHBox,...
                'Style',        'Edit',...
                'String',       '1',...
                'BackgroundColor','w',...
                'Callback',     @obj.sliderinput);
            obj.SliderGamma=uicontrol(...
                'Parent',       colorGammaHBox,...
                'Style',        'Slider',...
                'Min',          0.01,...
                'Max',          5,...
                'Value',        1);
            addlistener(obj.SliderGamma,'Value','PostSet',@obj.plot);
            set(colorGammaHBox,'Sizes',[40,55,-1]);
            
            obj.FlipCheckbox=uicontrol(...
                'Parent',       colormapHBox,...
                'Style',        'Checkbox',...
                'String',       'flip',...
                'Callback',     @obj.plot);
            obj.ColorMap=uicontrol(...
                'Parent',       colormapHBox,...
                'Style',        'popup',...
                'String',       {...
                'gray','jet','hsv','hot','copper','cool',...
                'spring','summer','autumn','winter','cm_blackbody',...
                'cm_blue','cm_bluehot','cm_blueredgreen','cm_copper',...
                'cm_cyan','cm_cyanmagenta','cm_geo','cm_gold','cm_green',...
                'cm_landandsea','cm_magenta','cm_pastels','cm_pastelsmap',...
                'cm_planetearth','cm_rainbow','red','cm_redwhiteblue',...
                'cm_redwhitegreen','cm_relief','cm_spectrum','cm_terrain',...
                'cm_yellow','cm_yellowhot'},...
                'Callback',     @obj.plot);
            set(colormapHBox,'Sizes',[50,-1]);
        end
        function setUIMenu(obj)
            hmenu=uimenu(obj.Figure,'Label','ARPES Tools');
            cutMenu=uimenu(hmenu,'Label','Cut','Callback',@obj.cut);
            if ~strcmp(obj.Direction,'z')
                set(cutMenu,'Enable','off');
            end
            uimenu(hmenu,'Label','Save Slice','Callback',@obj.save);
            plotMenu=uimenu(hmenu,'Label','Plot in new figure');
            uimenu(plotMenu,'Label','2D Plot',...
                'Callback',{@obj.newFigure,2});
            uimenu(plotMenu,'Label','3D Plot',...
                'Callback',{@obj.newFigure,3});
            uimenu(hmenu,'Label','Sub plot in new figure',...
                'Callback',@obj.newFigureSubPlot);
            obj.MenuMDCEDC=uimenu(hmenu,'Label','MDC/EDC',...
                'Callback',@obj.mdcedc);
            obj.MenuSaveMDCEDC=uimenu(hmenu,'Label','Save MDC/EDC',...
                'Callback',@obj.saveMdcEdc);
            set(obj.MenuSaveMDCEDC,'Enable','off');
            corMenu=uimenu(hmenu,'Label','Square correction',...
                'Callback',@obj.sqCorrection);
            if strcmp(obj.Direction,'z')
                set(corMenu,'Enable','off');
            end
            corrMenu=uimenu(hmenu,'Label','Fermi Surface correction','Callback',@obj.correction);
            %Along x or y is permitted
            if strcmp(obj.Direction,'z')
                set(corrMenu,'Enable','off');
            end
        end
        
        function save(obj,varargin)
            pos=get(obj.Slider,'Value');
            [~,posInd]=min(abs(obj.AxisData-pos));
            set(obj.Position,'String',num2str(pos));
            set(obj.PositionIndex,'String',num2str(posInd));
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
            index=arraySliceIndex(3,posInd,d);
            sliceData=squeeze(obj.Data.value(index{:}));
            data=obj.Data;
            data=rmfield(data,'z');
            data.x=xData;
            data.y=yData;
            data.value=sliceData;
            data.SliceInfo.Direction=obj.Direction;
            data.SliceInfo.pos=pos;
            data.SliceInfo.parentData=obj.DataName;
            SaveName=inputdlg({'Data Name'},...
                'Save Slice',...
                1,...
                {[obj.DataName,'_sa']});
            if isempty(SaveName)
                return;
            end
            assignin('base',SaveName{1},data);
        end
        function cut(obj,varargin)
            [xlist,ylist]=getpts(obj.Axis);
            x0=xlist(1);
            y0=ylist(1);
            x1=xlist(2);
            y1=ylist(2);
            [x,y,r]=area_secant_ph(obj.Data.x,obj.Data.y,x0,y0,x1,y1);
            [x_grid,y_grid,z_grid]=meshgrid(obj.Data.x,obj.Data.y,obj.Data.z);
            xx=repmat(x,max(size(obj.Data.z)),1);
            yy=repmat(y,max(size(obj.Data.z)),1);
            zz=repmat(obj.Data.z',1,max(size(x)));
            data=obj.Data;
            data=rmfield(data,'z');
            data.value=interp3(x_grid,y_grid,z_grid, permute(obj.Data.value,[2,1,3]),xx,yy,zz);
            %             % Temp for 3D visualization
            %             figure;
            %             surf(xx,yy,zz,data.value);
            %             %% Temp End
            data.value=data.value';
            data.x=r;
            data.y=obj.Data.z;
            figure;
            pcolor(data.x,data.y,data.value');
            shading interp;
            % save
            data.SliceInfo.Direction='nv';
            data.SliceInfo.pos=[xlist,ylist];
            data.SliceInfo.parentData=obj.DataName;
            SaveName=inputdlg({'Data Name'},...
                'Save Slice',...
                1,...
                {[obj.DataName,'_sa']});
            if isempty(SaveName)
                return;
            end
            assignin('base',SaveName{1},data);
        end
        function sliderinput(obj,varargin)
            minInput=str2double(get(obj.MinValue,'String'));
            if minInput<0
                minInput=0;
            end
            if minInput>1000
                minInput=1000;
            end
            if minInput<100
                minInput=minInput/100-1;
            else
                minInput=minInput/100;
                minInput=log10(minInput);
            end
            set(obj.SliderMin,'Value',minInput);
            
            maxInput=str2double(get(obj.MaxValue,'String'));
            if maxInput<0
                maxInput=0;
            end
            if maxInput>1000
                maxInput=1000;
            end
            if maxInput<100
                maxInput=maxInput/100-1;
            else
                maxInput=maxInput/100;
                maxInput=log10(maxInput);
            end
            set(obj.SliderMax,'Value',maxInput);
            
            gammaInput=str2double(get(obj.GammaValue,'String'));
            if 0.01>gammaInput
                gammaInput=0.01;
            end
            if 5<gammaInput
                gammaInput=5;
            end
            set(obj.SliderGamma,'Value',gammaInput);
            plot(obj,obj.Axis);
        end
        function plot(obj,varargin)
            %% Plot
            axPlot=obj.Axis;
            axFlag=0;
            if ~isempty(varargin)
                if isa(varargin{1},'char')
                    if strcmp(varargin{1},'new')
                        axPlot=varargin{2};
                        axFlag=1;
                    end
                end
            end
            minValue=min(obj.Data.value(:));
            maxValue=max(obj.Data.value(:));
            minSlider=get(obj.SliderMin,'Value');
            if minSlider<0
                minSlider=minSlider+1;
            else
                minSlider=10^minSlider;
            end
            minSlider=minSlider*(maxValue-minValue)+minValue;
            
            maxSlider=get(obj.SliderMax,'Value');
            if maxSlider<0
                maxSlider=maxSlider+1;
            else
                maxSlider=10^maxSlider;
            end
            maxSlider=maxSlider*(maxValue-minValue)+minValue;
            
            gamma=get(obj.SliderGamma,'Value');
            set(obj.GammaValue,'String',num2str(gamma));
            
            if minSlider<maxSlider
                cmin=minSlider;
                cmax=maxSlider;
            else
                cmin=maxSlider;
                cmax=minSlider;
            end
            if cmax==cmin
                cmax=cmin+1;
            end
            minPercent=100*(cmin-minValue)/(maxValue-minValue);
            maxPercent=100*(cmax-minValue)/(maxValue-minValue);
            set(obj.MinValue,'String',num2str(minPercent));
            set(obj.MaxValue,'String',num2str(maxPercent));
            pos=get(obj.Slider,'Value');
            [~,posInd]=min(abs(obj.AxisData-pos));
            set(obj.Position,'String',num2str(pos));
            set(obj.PositionIndex,'String',num2str(posInd));
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
            index=arraySliceIndex(3,posInd,d);
            sliceData=squeeze(obj.Data.value(index{:}));
            %% Color Mapping
            sliceData(sliceData<cmin)=cmin;
            sliceData=(sliceData-cmin)/(cmax-cmin);
            sliceData=sliceData.^gamma;
            sliceData=sliceData*(cmax-cmin)+cmin;
            if axFlag==1
                pcolor(axPlot,xData,yData,sliceData');
            else
                if ishandle(obj.SurfaceObj)
                    set(obj.SurfaceObj,...
                        'XData',xData,'YData',yData,...
                        'ZData',zeros(size(sliceData')),'CData',sliceData')
                else
                    obj.SurfaceObj=pcolor(obj.Axis,xData,yData,sliceData');
                end
            end
            shading interp;
            %use caxis to change the Min and Max of color
            cmapvar=get(obj.ColorMap,'Value');
            cmap=get(obj.ColorMap,'String');
            load('Colormap.mat');
            cmapData=eval(cmap{cmapvar});
            colormap(cmapData);
            if get(obj.FlipCheckbox,'Value')
                colormap(flipud(cmapData));
            end
            caxis([cmin,cmax]);
            if get(obj.Checkbox,'Value')
                XLim=xlim(obj.Axis);
                YLim=ylim(obj.Axis);
                axis equal;
                xlim(obj.Axis,XLim);
                ylim(obj.Axis,YLim);
            else
                axis normal;
            end
            flagMDCEDC = get(obj.MenuMDCEDC,'Checked');
            if strcmp(flagMDCEDC,'on')
                obj.plotCrossHair;
            end
            %% Plot Detected Cells
            if ~isempty(obj.DetectedPoints)
                % Permute coordinations
                switch obj.Direction
                    case 'x'
                        d=1;
                        xData=obj.Data.y;
                        yData=obj.Data.z;
                        zData=obj.Data.x;
                        pt_list(1,:) = obj.DetectedPoints(2,:);
                        pt_list(2,:) = obj.DetectedPoints(3,:);
                        pt_list(3,:) = obj.DetectedPoints(1,:);
                    case 'y'
                        d=2;
                        xData=obj.Data.x;
                        yData=obj.Data.z;
                        zData=obj.Data.y;
                        pt_list(1,:) = obj.DetectedPoints(1,:);
                        pt_list(2,:) = obj.DetectedPoints(3,:);
                        pt_list(3,:) = obj.DetectedPoints(2,:);
                    otherwise
                        d=3;
                        xData=obj.Data.x;
                        yData=obj.Data.y;
                        zData=obj.Data.z;
                        pt_list = obj.DetectedPoints;
                end
                % Find cells within integral interval
                lb = pos - obj.IntegralInterval/2;
                ub = pos + obj.IntegralInterval/2;
                v = (pt_list(3,:)>lb) .* (pt_list(3,:)<ub);
                pt_show = pt_list(:,v==1);
                ptArea = 36*obj.DetectedPointsArea/mean(obj.DetectedPointsArea);
                ptArea_show = ptArea(v==1);
                hold(axPlot,'on')
                if axFlag==1
                    scatter(pt_show(1,:),pt_show(2,:),ptArea_show);
                else
                    delete(obj.Scatter);
                    obj.Scatter = scatter(pt_show(1,:),pt_show(2,:),ptArea_show);
                end
                hold(axPlot,'off')
                
                
                
            end
        end
        function plot3D(obj,varargin)
            axPlot=obj.Axis;
            axFlag=0;
            if ~isempty(varargin)
                if isa(varargin{1},'char')
                    if strcmp(varargin{1},'new')
                        axPlot=varargin{2};
                        axFlag=1;
                    end
                end
            end
            minValue=min(obj.Data.value(:));
            maxValue=max(obj.Data.value(:));
            minSlider=get(obj.SliderMin,'Value');
            if minSlider<0
                minSlider=minSlider+1;
            else
                minSlider=10^minSlider;
            end
            minSlider=minSlider*(maxValue-minValue)+minValue;
            
            maxSlider=get(obj.SliderMax,'Value');
            if maxSlider<0
                maxSlider=maxSlider+1;
            else
                maxSlider=10^maxSlider;
            end
            maxSlider=maxSlider*(maxValue-minValue)+minValue;
            
            gamma=get(obj.SliderGamma,'Value');
            set(obj.GammaValue,'String',num2str(gamma));
            
            if minSlider<maxSlider
                cmin=minSlider;
                cmax=maxSlider;
            else
                cmin=maxSlider;
                cmax=minSlider;
            end
            if cmax==cmin
                cmax=cmin+1;
            end
            minPercent=100*(cmin-minValue)/(maxValue-minValue);
            maxPercent=100*(cmax-minValue)/(maxValue-minValue);
            set(obj.MinValue,'String',num2str(minPercent));
            set(obj.MaxValue,'String',num2str(maxPercent));
            pos=get(obj.Slider,'Value');
            [~,posInd]=min(abs(obj.AxisData-pos));
            set(obj.Position,'String',num2str(pos));
            set(obj.PositionIndex,'String',num2str(posInd));
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
            index=arraySliceIndex(3,posInd,d);
            sliceData=squeeze(obj.Data.value(index{:}));
            %normalization to use gamma index
            sliceData=(sliceData-minValue)/(maxValue-minValue);
            sliceData=sliceData.^gamma;
            sliceData=sliceData*(maxValue-minValue)+minValue;
            if axFlag==1
                surface(xData,yData,sliceData',sliceData','Parent',axPlot);
            else
                if ishandle(obj.SurfaceObj)
                    set(obj.SurfaceObj,...
                        'XData',xData,'YData',yData,...
                        'ZData',zeros(size(sliceData')),'CData',sliceData')
                else
                    obj.SurfaceObj=surface(xData,yData,sliceData',sliceData');
                end
            end
            shading interp;
            %use caxis to change the Min and Max of color
            cmapvar=get(obj.ColorMap,'Value');
            cmap=get(obj.ColorMap,'String');
            load('Colormap.mat');
            cmapData=eval(cmap{cmapvar});
            colormap(cmapData);
            if get(obj.FlipCheckbox,'Value')
                colormap(flipud(cmapData));
            end
            caxis([cmin,cmax]);
        end
        function newFigure(obj,varargin)
            plotType=varargin{3};
            pos=get(obj.Position,'String');
            posInd=get(obj.PositionIndex,'String');
            hFig=figure('Name',[obj.DataName,...
                ' ',obj.Direction,'=',pos,' ',...
                'index=',posInd],...
                'Color','w');
            h=axes('Parent',hFig);
            if plotType==3
                obj.plot3D('new',h);
            else
                obj.plot('new',h);
            end
        end
        function newFigureSubPlot(obj,varargin)
            inputText=inputdlg('Input Positions');
            plotArray=inputText{1};
            plotArray=str2num(plotArray);
            sliderMin=get(obj.Slider,'Min');
            sliderMax=get(obj.Slider,'Max');
            if isnan(plotArray)
                errordlg('Input Array should not be NaN');
                return;
            end
            if sum(plotArray<sliderMin)
                errordlg('Exeeds Min');
                return;
            end
            if sum(plotArray>sliderMax)
                errordlg('Exeeds Max');
                return;
            end
            n=length(plotArray);
            for i = 1:ceil(n/9)
                figure('Name',[obj.DataName,' ',obj.Direction],...
                    'Color','w');
                for j=1:9
                    k=(i-1)*9+j;
                    if k<=n
                        set(obj.Slider,'Value',plotArray(k));
                        hAx=subplot(3,3,j);
                        obj.plot('new',hAx);
                        title([obj.Direction,'=',num2str(plotArray(k))]);
                    end
                end
            end
        end
        function ButtonMotionFcn(obj,~,~)
            obj.plotCrossHair;
        end
        function ButtonUpFcn(obj,~,~)
            set(obj.Figure,...
                'WindowButtonMotionFcn',[],...
                'WindowButtonUpFcn',[]);
        end
        function setInd(obj,~,~,SetMethod)
            switch SetMethod
                case '+'
                    posInd = str2double(get(obj.PositionIndex,'String'))+1;
                case '-'
                    posInd = str2double(get(obj.PositionIndex,'String'))-1;
                case '='
                    posInd = str2double(get(obj.PositionIndex,'String'));
            end
            if posInd > length(obj.AxisData)
                posInd = length(obj.AxisData);
            end
            if posInd < 1
                posInd = 1;
            end
            pos = obj.AxisData(posInd);
            set(obj.Slider,'Value',pos);
        end
        function setPos(obj,~,~)
            pos = str2double(get(obj.Position,'String'));
            if pos > max(obj.AxisData)
                pos = max(obj.AxisData);
            end
            if pos < min(obj.AxisData)
                pos = min(obj.AxisData);
            end
            set(obj.Slider,'Value',pos);
        end
    end
    methods % separate files
        correction(obj,varargin) % Fermi surface correction
        plotCrossHair(obj,~,~) % plotCrossHair and EDC/MDC plot drawing
        ButtonDownFcn(obj,~,~) % mouse click action
        mdcedc(obj,varargin) % EDC/MDC function
        saveMdcEdc(obj,varargin) % save data of EDC/MDC plots
        sqCorrection(obj,varargin) % Square Correction for data set
    end
end


