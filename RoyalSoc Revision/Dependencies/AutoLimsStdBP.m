function lims = AutoLimsStdBP(Ymin,Ymax,ax)
% Uses input to set nice automatic Axis limits and ticklabels - for BoxPlots

%     Xmin = min(Y-Ystd,[],'all');
%     Xmax = max(Y+Ystd,[],'all');

Xmin = Ymin;
Xmax = Ymax;

Xrange = Xmax - Xmin;

%Same order of magnitude?
if Xmin > 0.1 * Xmax
else %or different
end

limL = Xmin - 0.05 * Xrange;
limR = Xmax + 0.05 * Xrange;

if Xmin == 0, limL = 0; end
if Xmax == 0, limR = 0; end

lims = [limL limR];
ylim(ax,lims);
%If only 2 ticks present, try keep them the same
if length(ax.YTick) < 3
	tickrange = ax.YTick(2) - ax.YTick(1);
	if ax.YTick(2) + tickrange - limR > limL - (ax.YTick(1) - tickrange)
		%try lower the bottom limit
		if Xmin - (ax.YTick(1) - tickrange) < 0.15 * Xrange
			limL = ax.YTick(1) - tickrange;
		end
	else
		if ax.YTick(2) + tickrange - Xmax < 0.15 * Xrange
			limR = ax.YTick(2) + tickrange;
		end
		%try raise the upper limit
	end
	lims = [limL limR];
end

end
