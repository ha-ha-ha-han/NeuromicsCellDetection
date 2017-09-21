function saveMdcEdc(obj,varargin)
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
            x0=obj.PosX;
            y0=obj.PosY;
            [~,Xind]=min(abs(x0-xData));
            [~,Yind]=min(abs(y0-yData));
            
            pos=get(obj.Slider,'Value');
            [~,posInd]=min(abs(obj.AxisData-pos));
            set(obj.Position,'String',num2str(pos));
            set(obj.PositionIndex,'String',num2str(posInd));
            
            index=arraySliceIndex(3,posInd,d);
            sliceData=squeeze(obj.Data.value(index{:}));

            xMDC=xData;
            yMDC=sliceData(:,Yind);
            xEDC=sliceData(Xind,:);
            yEDC=yData;
            saveEDCdata=struct('x',yEDC,'value',xEDC);
            saveMDCdata=struct('x',xMDC,'value',yMDC);
            assignin('base','MapData_MDC',saveMDCdata);
            assignin('base','MapData_EDC',saveEDCdata);
        end