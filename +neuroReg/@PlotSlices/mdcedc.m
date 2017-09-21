function mdcedc(obj,varargin)
            checkedFlg = get(obj.MenuMDCEDC,'Checked');
            if strcmp(checkedFlg,'on')
                set(obj.MenuMDCEDC,'Checked','off');
                set(obj.MenuSaveMDCEDC,'Enable','off');
                if isvalid(obj.FigureEDC)
                    close(obj.FigureEDC);
                end
                if isvalid(obj.FigureMDC)
                    close(obj.FigureMDC);
                end
                delete(obj.VLine);
                delete(obj.HLine);
                for i = 1:length(obj.Markers)
                    delete(obj.Markers(i));
                end
                obj.Markers=[];
                set(obj.Figure,'WindowButtonDownFcn','');
            else % When menu is unchecked
                set(obj.MenuMDCEDC,'Checked','on');
                set(obj.MenuSaveMDCEDC,'Enable','on');
                FigPos=get(obj.Figure,'Position');
                obj.FigureEDC=figure('Color','w','Position',...
                    [FigPos(1)+FigPos(3)+15,FigPos(2),FigPos(4),FigPos(4)]);
                obj.FigureMDC=figure('Color','w','Position',...
                    [FigPos(1),FigPos(2)-FigPos(4)-80,FigPos(3),FigPos(4)]);
                set(obj.FigureEDC,'Name','EDC');
                set(obj.FigureMDC,'Name','MDC');
                set(obj.Figure,'WindowButtonDownFcn',@obj.ButtonDownFcn);
            end
        end