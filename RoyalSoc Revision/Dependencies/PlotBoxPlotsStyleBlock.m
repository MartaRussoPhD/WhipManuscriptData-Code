function PlotBoxPlotsStyleBlock(sp,YBPplot,rgbmystyle,varargin)
% Plot a boxplot with applying formatting (Works with 16x5x2 data Subject x Block x Style)


flag_ThinnerLines = 0;
if any(strcmpi(varargin,'ThinnerLines'))
    flag_ThinnerLines = 1;   
end

YBPplotStyle = [];
if any(strcmpi(varargin,'AddAllBlocks'))
    YBPplotStyle = varargin{find(strcmpi(varargin,'AddAllBlocks'))+1};
end

flag_BarNotBox = 0;
if any(strcmpi(varargin,'BarNotBox'))
    flag_BarNotBox = 1;   
end

if ~isempty(YBPplotStyle) && size(YBPplot,2) == 5
    YBPplot(:,6,:) = nan(16,2);
    YBPplot(:,7,:) = nan(16,2);
    YBPplot(:,7,:) = YBPplotStyle;
end

if ~flag_BarNotBox
%Plot boxes
styles2plot = 2;
if size(YBPplot,3) == 1
    styles2plot = 1;
end
for istyle = 1:styles2plot
	Bp = boxplot(YBPplot(:,:,istyle),'BoxStyle','filled','Colors',rgbmystyle(istyle,:),'Widths',0.2);
	NChild = length(sp.Children(istyle).Children);
	
	if istyle == 1
		sp.Children(1).Tag = 'BoxPlotDiscrete';
	else
		sp.Children(1).Tag = 'BoxPlotRhythmic';
	end
	for ichild = NChild:(-1):1
		hC = sp.Children(1).Children(ichild);
		hC.XData = hC.XData - 0.2 * (istyle == 1) + 0.2 * (istyle == 2);
		if strcmpi(hC.Tag,'Whisker')
			colval = rgbmystyle(istyle,:);
			hC.Color = [0.5 0.5 0.5];
			hC.LineStyle = '-';
			hC.LineWidth = .5;
			WhiskEdges(1:2,(ichild-15),istyle) = hC.YData;
		elseif strcmpi(hC.Tag,'Box')
			%hC.LineWidth = 20;
			hC.Visible = 'off';
			%Create a rectangle instead
			ydat = hC.YData;
			xdat = hC.XData;
			colval = hC.Color;
            if any(isnan(ydat)), continue; end
			
            if flag_ThinnerLines
                NewRect = rectangle('Position',[xdat(1) ydat(1) diff(xdat) diff(ydat)],'EdgeColor','k','FaceColor',colval,'LineWidth',1);
            else
                NewRect = rectangle('Position',[xdat(1) ydat(1) diff(xdat) diff(ydat)],'EdgeColor','k','FaceColor',colval,'LineWidth',2);
            end
			NewRect.Position = NewRect.Position + [-0.13 0 0.26 0];
			NewRect.Tag = 'Box';
			NewRect.Parent = hC.Parent;
			delete(hC);
			ChildrenTemp = sp.Children(1).Children(2:ichild);
			ChildrenTemp(ichild) = sp.Children(1).Children(1);
			ChildrenTemp(ichild+1:NChild) = sp.Children(1).Children(ichild+1:end);
			sp.Children(1).Children = ChildrenTemp;
		elseif strcmpi(hC.Tag,'Median')
			hC.Color = [0 0 0];
			hC.LineWidth = 3;
            if any(isnan(ydat)), continue; end
            if flag_ThinnerLines
                NewMed = line(hC.XData + [-0.03 0.03],hC.YData,'LineWidth',2,'Color','k');
            else
                NewMed = line(hC.XData + [-0.03 0.03],hC.YData,'LineWidth',4,'Color','k');
            end
			NewMed.Tag = 'Median';
			NewMed.Parent = hC.Parent;
			delete(hC);
			ChildrenTemp = sp.Children(1).Children(2:ichild);
			ChildrenTemp(ichild) = sp.Children(1).Children(1);
			ChildrenTemp(ichild+1:NChild) = sp.Children(1).Children(ichild+1:end);
			sp.Children(1).Children = ChildrenTemp;
		elseif strcmpi(hC.Tag,'Outliers')
			hC.Marker = 'o';
			hC.MarkerEdgeColor = 'none';
			hC.MarkerFaceColor = 'none';%[0.7 0.7 0.7];
			hC.Visible = 'off';
		end
	end
end
%Add points.
for istyle = 1:2
	for iblock = 1:size(YBPplot,2)
		xdat = iblock - 0.2 * (istyle == 1) + 0.2 * (istyle == 2);
        if flag_ThinnerLines
            Pts = line(xdat .* ones(size(YBPplot,1),1),YBPplot(:,iblock,istyle),'LineStyle','none','Marker','o','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor','k');
        else
            Pts = line(xdat .* ones(size(YBPplot,1),1),YBPplot(:,iblock,istyle),'LineStyle','none','Marker','o','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor','k');
        end
		
	end
end
ylim('auto');
Y = YBPplot(:);
Ymed = median(Y,'omitnan');
Yiqr = iqr(Y);
% Adjust for optimal display.
indOutlier = find(Y > Ymed + 3 * Yiqr | Y < Ymed - 3 * Yiqr);
Y(indOutlier) = [];
PlotMax = max(Y,[],'all','omitnan');
PlotMin = min(Y,[],'all','omitnan');
if PlotMin ~= PlotMax
    ylim(AutoLimsStdBP(PlotMin,PlotMax,sp));
end



else
    % Plot Bars
    BarData = squeeze(mean(YBPplot,1,'omitnan'));
    ErrorbarData = squeeze(std(YBPplot,[],1,'omitnan'));
    Bp = bar(BarData,.6);
    Eb = errorbar(BarData,ErrorbarData);
    for istyle = 1:2
        Bp(istyle).FaceColor = rgbmystyle(istyle,:);
        Eb(istyle).LineStyle = 'none';
        Eb(istyle).Color = 'k';
        XTrue(:,istyle) = Bp(istyle).XEndPoints;
        Eb(istyle).XData = XTrue(:,istyle);
        Eb(istyle).CapSize = 3;
    end

    % Adjust limits
    ylim('auto');
    
end


end
