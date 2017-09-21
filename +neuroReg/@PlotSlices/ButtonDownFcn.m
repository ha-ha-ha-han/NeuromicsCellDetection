function ButtonDownFcn(obj,~,~)
            EventFlag=get(obj.Figure,'SelectionType');
            FigPos=get(obj.Figure,'Position');
            switch EventFlag
                case 'normal' 
                    % mouse left click ==> plot EDC and MDC
                    if ~isvalid(obj.FigureEDC)
                        obj.FigureEDC=figure('Color','w',...
                            'Position',[FigPos(1)+FigPos(3)+15,FigPos(2),FigPos(4),FigPos(4)]);
                        set(obj.FigureEDC,'Name','MDC');
                    end
                    if ~isvalid(obj.FigureMDC)
                        obj.FigureMDC=figure('Color','w',...
                            'Position',[FigPos(1),FigPos(2)-FigPos(4)-80,FigPos(3),FigPos(4)]);
                        set(obj.FigureMDC,'Name','EDC');
                    end
                    obj.plotCrossHair;
                    set(obj.Figure,...
                        'WindowButtonMotionFcn',@obj.ButtonMotionFcn,...
                        'WindowButtonUpFcn',    @obj.ButtonUpFcn);
                case 'alt'
                    % mouse right click ==> mark
                    if ~(ishandle(obj.VLine)||ishandle(obj.HLine))
                        return;
                    end
                    ptx=get(obj.VLine,'XData');
                    pty=get(obj.HLine,'YData');
                    x0=ptx(1);
                    y0=pty(1);
                    hold(obj.Axis,'on');
                    len=length(obj.Markers);
                    obj.Markers(len+1)=scatter(obj.Axis,x0,y0,100,'+','g');
                    hold off;
                otherwise
                    % double click or shift+left click ==> delete last marker
                    len=length(obj.Markers);
                    if len>=1
                        if ishandle(obj.Markers(len))
                            delete(obj.Markers(len));
                        end
                        NewMarkersList=obj.Markers(1:len-1);
                        obj.Markers=NewMarkersList;
                    end
            end
        end