classdef VisTransform < handle
    %VisualizeTransform creates a GUI that display the transform from registration.
    %   obj = VisTransform(TransTable,DataSets,pt_list_vol,...
    %            pt_list_slice,ex_list,Option,filepath)
    %   The transform is from volume to slice in the neuron registration.
    %   To create a GUI, use VisualizeTransform(TransTable, ...
    %   DataSets, pt_list_vol, pt_list_slice, Option, ex_list,filepath). 
    %   TransformTable: can have multiple observants.
    %   DataSets: 
    %       DataSets.dataZ = dataZ;
    %       DataSets.data_slice = data_slice;
    %       data_slice_bw = neuroReg.bwCell2(data_slice,Option,ex_list);
    %       DataSets.data_slice_bw_low = neuroReg.downSample(...
    %           data_slice_bw,Option.StepX,[],Option.StepX);
    %   Option: see neuroReg.setOption.
    %   All these input vaiables can get from the results of
    %   neuroReg.rotationCorr3 pipeline.
    % by Han Peng, 2017.07.29.
  
    
    
    properties
        Gui
        Data
        Option
    end    
    methods
        function obj = VisTransform(TransTable,DataSets,pt_list_vol,...
                pt_list_slice,ex_list,Option,filepath)
            obj.Data.TransTable = TransTable;
            obj.Data.TransTableNow = obj.Data.TransTable(1,:);
            obj.Data.pt_list_vol = pt_list_vol;
            obj.Data.pt_list_slice = pt_list_slice;
            obj.Data.ex_list = ex_list;
            obj.Data.DataSets = DataSets;
            obj.Option = neuroReg.setOption(Option);
            data_slice = obj.Data.DataSets.data_slice;
            obj.Data.DataSets.data_slice_low = ...
                neuroReg.downSample(data_slice,4,[],4);
            % obj.Data.DataSets.data_slice_low = data_slice;
            obj.Gui.Path = filepath;
            constructFigure(obj);
        end
        function constructFigure(obj)
            [N,~] = size(obj.Data.TransTable);
            obj.Gui.Figure = figure('Position',[200 600 1500 150],...
                'Name',['VisTransform - ',obj.Gui.Path],...
                'NumberTitle','off',...
                'MenuBar','none',...
                'IntegerHandle','off',...
                'HandleVisibility','off',...
                'Tag','VisTransform',...
                'WindowKeyPressFcn',@obj.pressKey);
            HomeVBox = uiextras.VBox('Parent',obj.Gui.Figure,...
                'BackgroundColor','w');
            %% HBox1: for VBox11 and VBox12
            HBox1 = uiextras.HBox('Parent',HomeVBox,...
                'BackgroundColor','w'); % HBox for VBox1 and Image            
            VBox1 = uiextras.VBox('Parent',HBox1,...
                'BackgroundColor',[0.94 0.94 0.94]); % VBox for data and control                     
            HBoxControl = uiextras.HBox('Parent',VBox1,...
                'BackgroundColor',[0.94 0.94 0.94]); % HBox for uicontrols
            HBoxTextInput = uiextras.HBox('Parent',VBox1,...
                'BackgroundColor',[0.94 0.94 0.94]); % HBox for current table#
            %% HBoxControl: for controllers
            obj.Gui.Save = uicontrol(...
                'Parent',       HBoxControl,...
                'Style',        'pushbutton',...
                'String',       'Save',...
                'Callback',     @obj.save);
            obj.Gui.Save = uicontrol(...
                'Parent',       HBoxControl,...
                'Style',        'pushbutton',...
                'String',       'Fiugure Output',...
                'Callback',     @obj.plotOut);
            obj.Gui.Help = uicontrol(...
                'Parent',       HBoxControl,...
                'Style',        'pushbutton',...
                'String',       'Help',...
                'Callback',     @obj.help);
            set(HBoxControl,'Sizes',[100 100 100]);
            %% HBoxTextInput: Current Table
            uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Text',...
                'String',       'Table#');
            obj.Gui.EditTableNum = uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Edit',...
                'String',       '1',...
                'Callback',     @obj.setTransTable);
            obj.Gui.TotNum = uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Text',...
                'String',       ['/ ',num2str(N)]);
            uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Text',...
                'String',       '   Current Figure#');
            obj.Gui.EditFigNum1 = uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Edit',...
                'String',       '12580');
            uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Text',...
                'String',       '   Output Figure#');
            obj.Gui.FigNum2 = uicontrol(...
                'Parent',       HBoxTextInput,...
                'Style',        'Edit',...
                'String',       '12588');
            set(HBoxTextInput,'Sizes',[40 30 40 60 40 60 40])
            %% TransTable
            CT1 = obj.Data.TransTableNow;
            cn1 = {'Intense','Alpha','Beta','Gamma',...
                'X','Y','Z'};
            obj.Gui.TransTable = uitable(...
                'Parent',       HBox1,...
                'Data',         CT1{:,:},...
                'ColumnName',   cn1);
            cn2 = {'Integ'    'CellRadius'    'StepD'    'StepX'  ...
                'MagicNumber'    'SigmaX'  'SigmaY'  'SigmaRender'  ...
                'Res0'    'SizeLimit1' 'SizeLimit2'   'Threshold'  ...
                'MedianFilterSizeX' 'MedianFilterSizeY'   ...
                'TransTol'    'AngleTol' 'MaxPeakNum'};
            CT2 = struct2table(obj.Option);
            obj.Gui.OptionTable = uitable(...
                'Parent',       HomeVBox,...
                'Data',         CT2{:,:},...
                'ColumnName',   cn2);
            %% Set Sizes for VBox1
            set(VBox1,'Sizes',[50 25]);
            set(HBox1,'Sizes',[320 -1]);
            set(HomeVBox,'Sizes',[80,70]);

            %% Plot
%             HBoxAxis = uiextras.HBox('Parent',VBox2,...
%                 'BackgroundColor','w'); 
%             obj.Gui.Axis1 = axes(HBoxAxis);
%             obj.Gui.Axis2 = axes(HBoxAxis);
            %% Finalize
            obj.Gui.FigNum1 = 12580;
            obj.Gui.FigNum2 = 12588;
            tic
            obj.plot;
            toc
        end

        function plot(obj,varargin)
            % Update Tables. Update Figure.
            icp_flag = 0;
            if nargin>1
                plot_method = varargin{1};
                if strcmp(plot_method,'icp')
                    icp_flag = 1;
                    M_icp_s2v = varargin{2};
                    M_icp_v2s = varargin{3};
                    UseMe = varargin{4};
                    % Matched cell numbers:
                    n_UseMe = sum(UseMe);
                    % Print information
                    disp('Plot ICP result.')
                    fprintf('Number of matched cells: %d\n',n_UseMe);
                    disp('slice2volume');
                    disp(M_icp_s2v);
                    disp('volume2slice');
                    disp(M_icp_v2s);
                end
            end
            %% Update Tables
            CT1 = obj.Data.TransTableNow;
            set(obj.Gui.TransTable,'Data',CT1{:,:});
            CT2 = struct2table(obj.Option);
            set(obj.Gui.OptionTable,'Data',CT2{:,:});

            %% Update Figure
            Integ = obj.Option.Integ;
            CellRadius = obj.Option.CellRadius;
            dataZ = obj.Data.DataSets.dataZ;
            data_slice = obj.Data.DataSets.data_slice_low;
            pt_list_vol = obj.Data.pt_list_vol;
            pt_list_slice = obj.Data.pt_list_slice;
            TransTable = table2array(obj.Data.TransTableNow);
            TransParameters = TransTable(1,2:end);
            [~,R,x_lim,y_lim,d_lim] = ...
                neuroReg.rotateCells(pt_list_vol,...
                TransParameters(1),TransParameters(2),TransParameters(3));
            t = TransParameters(4:6)';
            if icp_flag == 0
                M = [R',-R'*t]; % M: slice to volume. Default.
                M1 = [R,t]; % Volume to Slice
            elseif icp_flag == 1
                M = M_icp_s2v ; %  M: slice to volume. Default.
                M1 = M_icp_v2s; % Volume to Slice
            end
                
            pt_list_vol_new = M1*[pt_list_vol;ones(size(pt_list_vol(1,:)))];
            % One day I will optimize this line:
            [~,b_plane,b_list_volume] = neuroReg.cutVolume(dataZ,data_slice,M,0);
            % It is stupid to use cutVolume to calculate the intersection
            % of the slice and the cube. Sometimes I am astonished by my
            % stupidity and laziness. Well, at least it works for now.
            %% Plot: Axis 1
            figure(obj.Gui.FigNum1)
            h1 = subplot(1,2,1);
            [xx0,zz0] = ndgrid(data_slice.x,data_slice.y);
            yy0 = zeros(size(xx0));
            surf(h1,xx0,yy0,zz0,double(data_slice.value));
            shading(h1,'flat');
            axis(h1,'equal');
            axis(h1,'vis3d');
            hold(h1,'on');
            vts = getCube([dataZ.x(1),dataZ.y(1),dataZ.z(1)],...
                [dataZ.x(end)-dataZ.x(1),dataZ.y(end)-dataZ.y(1),dataZ.z(end)-dataZ.z(1)]);
            vts1 = M1*[reshape(vts,[3,24]);ones(1,24)];
            vts1 = reshape(vts1,[3 4 6]);
            plot3(h1,b_plane(1,:),zeros(size(b_plane(1,:))),b_plane(2,:),'r');
            drawCube(h1,vts1,'r');
            hold(h1,'off');
            title('Geometry');
            %% Plot: Axis 2
            figure(obj.Gui.FigNum1)
            h2 = subplot(1,2,2);
            xs_full = pt_list_slice(1,:);
            ys_full = pt_list_slice(2,:);
            xb = b_plane(1,:);
            yb = b_plane(2,:);
            pt_list_slice_now = pt_list_slice(:,inpolygon(xs_full,ys_full,xb,yb));
            pt_list_vol_rotated = M1*[pt_list_vol;ones(size(pt_list_vol(1,:)))];
            pf = abs(pt_list_vol_rotated(2,:))<Integ;
            pt_list_vol_now = pt_list_vol_rotated(:,pf);
            xs = pt_list_slice_now(1,:);
            ys = pt_list_slice_now(2,:);
            xv = pt_list_vol_now(1,:);
            yv = pt_list_vol_now(3,:);
            dv = pt_list_vol_now(2,:);
            pt_now_area = ones(1,sum(pf==1)).*exp(-(dv/CellRadius).^2/2);            
            hold(h2,'off');
            scatter(h2,xs,ys,72,'k+');
            hold(h2,'on');
            scatter(h2,xv,yv,72*pt_now_area,'ro');
            axis(h2,'equal');
            xlim(h2,[b_plane(1,1),b_plane(1,3)]);
            ylim(h2,[b_plane(2,1),b_plane(2,2)]);
            if icp_flag == 1
                scatter(h2,pt_list_slice(1,UseMe),pt_list_slice(2,UseMe),300,'sb');
            end
            hold(h2,'off');
            title(h2,'+ = Slice, O = Volume');
            box(h2,'on');
            hold(h1,'on')
            h=scatter3(h1,xv,pt_list_vol_now(2,:),yv,'r');
            hold(h1,'off')
            figure(obj.Gui.Figure);
            function vts = getCube ( origin, size )
                x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
                y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
                z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
                vts(1,:,:) = x;
                vts(2,:,:) = y;
                vts(3,:,:) = z;
            end
            
            function drawCube(h1,vts,c)
                for i=1:6
                    hold on
                    plot3(h1,vts(1,:,i),vts(2,:,i),vts(3,:,i),c);
                    hold off
                end
                h1=patch(h1,vts(1,:,6),vts(2,:,6),vts(3,:,6),'y');
                alpha(h1,0.3);
            end
            
        end
        function plotOut(obj,varargin)    
            % Call neuroReg.plotTransform.
            % plotOut(obj,plot_method)  
            %   if plot_method is 'icp', then plot out the icp transform
            % plotOut(obj) 
            %   Plot the transform specified by the TransTable
            icp_flag = 0;
            if nargin>1
                plot_method = varargin{1};
                if strcmp(plot_method,'icp')
                    icp_flag = 1;
                    M_icp_s2v = obj.Data.ICP.M_icp_s2v;
                    M_icp_v2s = obj.Data.ICP.M_icp_v2s;
                    UseMe = obj.Data.ICP.UseMe;
                    % Matched cell numbers:
                    n_UseMe = sum(UseMe);
                    % Print information
                    disp('Plot ICP result.')
                    fprintf('Number of matched cells: %d\n',n_UseMe);
                    disp('slice2volume');
                    disp(M_icp_s2v);
                    disp('volume2slice');
                    disp(M_icp_v2s);
                end
            end
            h = figure(obj.Gui.FigNum2);
            if icp_flag == 1
            neuroReg.plotTransform(h,...
                obj.Data.TransTableNow,...
                obj.Data.DataSets,...
                obj.Data.pt_list_vol,...
                obj.Data.pt_list_slice,...
                obj.Option,...
                obj.Data.ex_list,...
                M_icp_s2v,...
                M_icp_v2s,...
                UseMe); 
            else
                neuroReg.plotTransform(h,...
                obj.Data.TransTableNow,...
                obj.Data.DataSets,...
                obj.Data.pt_list_vol,...
                obj.Data.pt_list_slice,...
                obj.Option,...
                obj.Data.ex_list);
            end
            figure(obj.Gui.Figure);
        end
        function save(obj,varargin)
            % Save current transformation to the TransTable.
            % Write the TransTable to obj.Gui.Path
            
            [N,~] = size(obj.Data.TransTable);
            obj.Data.TransTable(N+1,:) = obj.Data.TransTableNow;            
            set(obj.Gui.TotNum,'String',['/ ',num2str(N+1)]);   
            set(obj.Gui.EditTableNum,'String',num2str(N+1));
            filepathname = fullfile(obj.Gui.Path,'TransTable_Save.mat');
            TransTable = obj.Data.TransTable;
            save(filepathname,'TransTable');
            disp([filepathname,' saved at ',datestr(now)]);
        end
        function setTransTable(obj,varargin)
            % When obj.Gui.EditTableNum is changed, set the current
            % transtable (obj.Data.TransTableNow) to the designated index.
            
            ind = get(obj.Gui.EditTableNum,'String');
            ind = str2double(ind);
            obj.Data.TransTableNow = obj.Data.TransTable(ind,:);
            obj.plot;
        end
        function pressKey(obj,~,keydata)
            % Pree 'Help' to check the behaviour of WindowKeyPressFcn
            STEP_BIG = 50;
            STEP_SMALL = 2;
            cf = 0; % If anything changed ==> cf=1
            table_now = table2array(obj.Data.TransTableNow);
            TransParameters = table_now(2:7);
            if isempty(keydata.Modifier)
                switch keydata.Key
                    %% Rotation
                    case 'i'
                        cf = 1;
                        TransParameters(3) = TransParameters(3) + 1;
                    case 'k'
                        cf = 1;
                        TransParameters(3) = TransParameters(3) - 1;
                    case 'u'
                        cf = 1;
                        TransParameters(2) = TransParameters(2) + 1;
                    case 'o'
                        cf = 1;
                        TransParameters(2) = TransParameters(2) - 1;
                    case 'j'
                        cf = 1;
                        TransParameters(1) = TransParameters(1) - 1;
                    case 'l'
                        cf = 1;
                        TransParameters(1) = TransParameters(1) + 1;
                        %% Translation
                    case 'q' % Y+
                        cf = 1;
                        TransParameters(5) = TransParameters(5) + STEP_SMALL;
                    case 'e' % Y-
                        cf = 1;
                        TransParameters(5) = TransParameters(5) - STEP_SMALL;
                    case 'w' % Z+
                        cf = 1;
                        TransParameters(6) = TransParameters(6) + STEP_SMALL;
                    case 's' % Z-
                        cf = 1;
                        TransParameters(6) = TransParameters(6) - STEP_SMALL;
                    case 'a' % X-
                        cf = 1;
                        TransParameters(4) = TransParameters(4) - STEP_SMALL;
                    case 'd' % X+
                        cf = 1;
                        TransParameters(4) = TransParameters(4) + STEP_SMALL;
                        %% Options
                    case 'rightbracket'
                        cf = 1;
                        obj.Option.Integ = obj.Option.Integ + 2;
                    case 'leftbracket'
                        cf = 1;
                        obj.Option.Integ = obj.Option.Integ - 2;
                    case 'p'
                        cf = 0;
                        obj.plotOut;
                    otherwise
                end
            elseif strcmp(keydata.Modifier{1},'control')
                switch keydata.Key
                    case 'p'
                        obj.plotOut('icp');
                    case 'i'
                        obj.icp;
                    otherwise
                end
            elseif strcmp(keydata.Modifier{1},'shift')
                switch keydata.Key
                    case 'q' % Y+
                        cf = 1;
                        TransParameters(5) = TransParameters(5) + STEP_BIG;
                    case 'e' % Y-
                        cf = 1;
                        TransParameters(5) = TransParameters(5) - STEP_BIG;
                    case 'w' % Z+
                        cf = 1;
                        TransParameters(6) = TransParameters(6) + STEP_BIG;
                    case 's' % Z-
                        cf = 1;
                        TransParameters(6) = TransParameters(6) - STEP_BIG;
                    case 'a' % X-
                        cf = 1;
                        TransParameters(4) = TransParameters(4) - STEP_BIG;
                    case 'd' % X+
                        cf = 1;
                        TransParameters(4) = TransParameters(4) + STEP_BIG;
                    otherwise
                end
            end
            if cf==1 % In case anything changed
                obj.Data.TransTableNow.Angles = TransParameters(1:3);
                obj.Data.TransTableNow.Translation = TransParameters(4:6);
                obj.plot;
            end
        end
        function help(obj,varargin)
            % Show help.
            
            h = figure('Name','Help',...
                'NumberTitle','off',...
                'MenuBar','none',...
                'IntegerHandle','off',...
                'HandleVisibility','off',...
                'Color','w');
            h1 = axes(h);
            fp = which('neuroReg.VisTransform');
            fp = fileparts(fileparts(fp));
            fp = fullfile(fp,'+neuroReg','VisTransform.png');
            Img = imread(fp);
            imshow(Img,'Parent',h1);            
        end
        function [M_icp_s2v,M_icp_v2s]=icp(obj,varargin)
            % ICP based on the current transform matrix.
            % obj.icp plots and print the icp result.
            % obj.icp(d) plots and print the icp result, with the maximum 
            % neigbhor range d.
            % THIS PART IS FOR TEST AND HAS NO GUI AT THE MOMENT.            
            if nargin>=2
                d = varargin{1};
            else
                d = 20;
            end
            pt_list_slice = obj.Data.pt_list_slice;
            pt_list_vol = obj.Data.pt_list_vol;    
            TransParameters = table2array(obj.Data.TransTableNow(:,2:end));
            %% Get M2
            [~,R,x_lim,y_lim,d_lim] = ...
                neuroReg.rotateCells(pt_list_vol,...
                TransParameters(1),TransParameters(2),TransParameters(3));
            t = TransParameters(4:6)';
            M = [R',-R'*t]; % M: slice to volume. Default.
            M_23 = M(:,[1 3 4]);

            [M2,Error,MyNeighb,UseMe] = neuroReg.PointCloudRegister2(pt_list_vol',pt_list_slice',M_23',d,[]);
            M2 = M2';
            y3 = M2*cat(1,pt_list_slice,ones(size(pt_list_slice(1,:))));
            figure(10099);cla;
            % Plot cells in the volume
            scatter3(pt_list_vol(1,:),pt_list_vol(2,:),pt_list_vol(3,:),'k');
            ha=gca;
            hold on;
            % Plot cells in the slice
            scatter3(y3(1,:),y3(2,:),y3(3,:),'c');    
            scatter3(y3(1,UseMe),y3(2,UseMe),y3(3,UseMe),300,'rs');  
            % Plot Boundary
            x_axis = obj.Data.DataSets.data_slice.x;
            y_axis = obj.Data.DataSets.data_slice.y;
            xb1 = x_axis(1);
            xb2 = x_axis(end);
            yb1 = y_axis(1);
            yb2 = y_axis(end);
            vts_slice = [xb1 yb1;xb2 yb1; xb2 yb2;xb1 yb2; xb1 yb1]';
            vts_slice_trans = M2*cat(1,vts_slice,ones(size(vts_slice(1,:))));
            hp = patch(vts_slice_trans(1,:),vts_slice_trans(2,:),vts_slice_trans(3,:),'k');
            alpha(hp,0.3);
            x_axis = obj.Data.DataSets.dataZ.x;
            y_axis = obj.Data.DataSets.dataZ.y;
            z_axis = obj.Data.DataSets.dataZ.z;
            xrange = range(x_axis);
            yrange = range(y_axis);
            zrange = range(z_axis);
            vts = getCube([0 0 0],[xrange yrange zrange]);
            drawCube(ha,vts,'b');
            disp(sum(UseMe==1))
            hold off;
            axis equal;
            axis vis3d;
            %% Get the affine matrix (3-by-4) from the (3-by-3) matrix
            n1 = M2*[1;0;0];
            n3 = M2*[0;1;0];
            n2 = cross(n3,n1);
            n2 = n2/norm(n2);
            M3(:,1) = M2(:,1);
            M3(:,2) = n2;
            M3(:,3) = M2(:,2);
            M3(:,4) = M2(:,3);
            M_icp_s2v = M3;
            M4 = [M3;[0,0,0,1]];
            M_icp_v2s = M4^(-1);
            M_icp_v2s = M_icp_v2s(1:3,:);
            obj.Data.ICP.M_icp_v2s = M_icp_v2s;
            obj.Data.ICP.M_icp_s2v = M_icp_s2v;
            obj.Data.ICP.UseMe = UseMe;
            obj.plot('icp',M_icp_s2v,M_icp_v2s,UseMe);
            %% Visualization
            
            %%
            function vts = getCube ( origin, size )
                x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
                y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
                z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
                vts(1,:,:) = x;
                vts(2,:,:) = y;
                vts(3,:,:) = z;
            end
            
            function drawCube(h1,vts,c)
                for i=1:6
                    hold on
                    plot3(h1,vts(1,:,i),vts(2,:,i),vts(3,:,i),c);
                    hold off
                end
                h=patch(h1,vts(1,:,6),vts(2,:,6),vts(3,:,6),'y');
                alpha(h,0.3);
            end
            
        end
    end
            
    
end

