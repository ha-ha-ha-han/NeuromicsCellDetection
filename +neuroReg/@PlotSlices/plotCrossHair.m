function plotCrossHair(obj,~,~)
            delete(obj.VLine);
            delete(obj.HLine);

            % Deal with Direction
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
            end % end of switch obj.Direction

            pt=get(obj.Axis,'CurrentPoint');
            obj.PosX=pt(1,1);
            obj.PosY=pt(1,2);
            x0=obj.PosX;
            y0=obj.PosY;
            % Get Index
            [~,Xind]=min(abs(x0-xData));
            [~,Yind]=min(abs(y0-yData));
            x0=xData(Xind);
            y0=yData(Yind);
            % Plot Cross Hair

            obj.VLine=line([x0,x0],ylim(obj.Axis),...
                'Color','red','Parent',obj.Axis);
            obj.HLine=line(xlim(obj.Axis),[y0,y0],...
                'Color','red','Parent',obj.Axis);

            % Plot Edc and Mdc

            pos=get(obj.Slider,'Value');
            [~,posInd]=min(abs(obj.AxisData-pos));
            set(obj.Position,'String',num2str(pos));
            set(obj.PositionIndex,'String',num2str(posInd));

            index=arraySliceIndex(3,posInd,d);
            sliceData=squeeze(obj.Data.value(index{:}));
            figure(obj.FigureEDC);
            plot(sliceData(Xind,:),yData);
            hold on
            line(xlim,[y0,y0])
            hold off
            figure(obj.FigureMDC);
            plot(xData,sliceData(:,Yind));
            hold on
            line([x0,x0],ylim);
            hold off
        end