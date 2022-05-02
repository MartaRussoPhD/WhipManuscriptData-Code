function [T,PercRemoved,ind] = RemoveOutliersMy(T,thrMAD_STD,varargin)
% Fix Scale - remove outliers for Displaying AllTrials


% When T is double, it can be a 2-D value and will be treated from plain geometry
% perspective

if strcmpi(thrMAD_STD,'Agressive4Plot') || any(strcmpi(varargin,'Agressive4Plot'))
	AggressiveMargin = 1;
else
	AggressiveMargin = 0;
end


if istable(T)
	mea = mean(T.Y,'omitnan');
	med = median(T.Y,'omitnan');
	madsc = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(T.Y - med),'omitnan');
	stdd = std(T.Y,[],'omitnan');
else
	mea = mean(T,'all','omitnan');
	med = median(T,'all','omitnan');
	madsc = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(T - med),'all','omitnan');
	stdd = std(T,[],'all','omitnan');
end



% if binary, do nothing.
if islogical(T) || (istable(T) && islogical(T.Y))
    PercRemoved = 0;
    return
end

% %Outlier: outside of 3 MAD OR outside of 3 STD
% ind = find(((T.Y > med + 3 * madsc) & (T.Y > mea + 3 * stdd)));

%Outlier: outside of 2.5 Std OR outside of 2.5 MAD

if AggressiveMargin
	thresholdMAD = 1.5;
	thresholdSTD = 1.5;
elseif nargin < 2
	thresholdMAD = 3;
	thresholdSTD = 3;
else
	thresholdMAD = thrMAD_STD(1);
	thresholdSTD = thrMAD_STD(2);
end


if istable(T)
	ind = find(((T.Y > med + thresholdMAD * madsc) & (T.Y > mea + thresholdSTD * stdd)) |...
		((T.Y < mea - thresholdSTD * stdd) & (T.Y < med - thresholdMAD * madsc)));
elseif size(T,2) == 1
	ind = find(((T > med + thresholdMAD * madsc) & (T > mea + thresholdSTD * stdd)) |...
		((T < mea - thresholdSTD * stdd) & (T < med - thresholdMAD * madsc)));
    
    
elseif size(T,2) == 2
    % Planar. 
    % For parametric distribution, find cov.ellipsoid centered at mean, enlarge its
    % semiaxes by thresholdSTD, find Focal Pts, and find the points that are outside of
    % that ellipsoid.
    [largest_eigenvec,smallest_eigenvec,a,b,R] = error_ellipse(T,'NoChi');
    EllC = mean(T,'omitnan');
    Tc = CoordTransform(T,R,EllC'); % Now my X is major semiaxis
    aout = a * thresholdSTD;
    bout = b * thresholdSTD;
    foc1 = [-sqrt(aout.^2-bout.^2) 0];
    foc2 = [sqrt(aout.^2-bout.^2) 0];
    dist2foc1 = sqrt(sum((Tc - foc1).^2,2));
    dist2foc2 = sqrt(sum((Tc - foc2).^2,2));
    PtOutsideBin = dist2foc1 + dist2foc2 - aout * 2 > 0;
    ind = find(PtOutsideBin);

    % For nonparametric distribution, use the same ellipsoid, but centered at median.
    % Assuming proper orientation, convert(rotate) points to its r.f., re-measure semiaxes
    % with MAD somehow........
    
    
%     % DEBUG
if  1 == 0
    figure, scatter(Tc(:,1),Tc(:,2),5,'k','filled','MarkerEdgeColor','none'); hold on;
    % Ellipse
    Tell = linspace(0,2*pi,50)';
    Xell = aout * cos(Tell);
    Yell = bout * sin(Tell);
    EllPlot = [Xell Yell] * inv(R) + EllC;
    EllPlot = [Xell Yell];
    line(EllPlot(:,1),EllPlot(:,2),'Color',[0.5 0.5 0.5]);
    % Discarded pts
    scatter(Tc(ind,1),Tc(ind,2),25,'k');
end
    
end







if istable(T)
	for i = 1:length(ind)
		T.Y(ind(i)) = nan;
	end
	PercRemoved = round(length(ind) / height(T) .* 100,2);
else
	T(ind,:) = nan(length(ind),size(T,2));
	PercRemoved = round(length(ind) / length(T) .* 100,2);
end


end
