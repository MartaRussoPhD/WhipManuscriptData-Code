function lims = AutoLims(Y,ax,Dir,varargin)
% Uses input to set nice automatic Axis limits and ticklabels - AllTrials
% Dir = 'Y' (default) or 'X'
if nargin < 3
    Dir = 'Y';
end
tickobj = strcat(Dir,'Tick');
limobj = strcat(Dir,'Lim');

flag_topbottomalign = 0;
if any(strcmpi(varargin,'TopBottomAlign'))
    flag_topbottomalign = 1;
end

if ~ishandle(ax)
    disp('Second argument must be the axes, skipping AutoLims');
    return
end


if isempty(Y)
    % if histogram
    for iax = 1:length(ax)
        obj(iax,:) = findobj('Parent',ax(iax),'Type','Histogram');
    end
    objlist = obj(:);
    for iobj = 1:numel(obj)
        Ydat = objlist(iobj).Values;
        Y(end+1:end+length(Ydat)) = Ydat;
    end
end
Xmin = min(Y,[],'all');
Xmax = max(Y,[],'all');

Xrange = Xmax - Xmin;

%Same order of magnitude?
if Xmin > 0.1 * Xmax
else %or different
end

%Extend the limits slightly beyond the range of the data
limL = Xmin - 0.05 * Xrange;
limR = Xmax + 0.05 * Xrange;

%But keep it at zero when applicable
if limL < 0 && Xmin > 0, limL = 0; end
if limR > 0 && Xmax < 0, limR = 0; end

lims = [limL limR];

% assuming axes are already linked
ax = ax(1);
if lims(1) == lims(2)
    lims(1) = lims(1) - .5 * lims(1);
    lims(2) = lims(2) + .5 * lims(2);
end
ax.(limobj) = lims;

%Make sure there're enough ticks
if length(ax.(tickobj)) < 3
	tickrange = ax.(tickobj)(2) - ax.(tickobj)(1);
	if ax.(tickobj)(2) + tickrange - limR > limL - (ax.(tickobj)(1) - tickrange)
		%try lower the bottom limit
		if Xmin - (ax.(tickobj)(1) - tickrange) < 0.15 * Xrange
			limL = ax.(tickobj)(1) - tickrange;
		end
    else
        %try raise the upper limit
		if ax.(tickobj)(2) + tickrange - Xmax < 0.15 * Xrange
			limR = ax.(tickobj)(2) + tickrange;
        end
	end
	lims = [limL limR];
end
end