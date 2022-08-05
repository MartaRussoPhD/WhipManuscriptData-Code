function WhipTask_HSAveragePlot(srcTimeSeries,rgbmystyle)

% Using previously averaged data, this function plots HandSpeed for each participant and
% style, aligning the time-series by the end of the longest trial average. To set the
% XAxis leftwards from -2 to 0 seconds
% Average time is obtained from averaging durations
% Average speed profile is obtained from time-normalizing and averaging trial time-series
% SD patch is obtained the same way.

Nresamp = 200;
SampleRate = 500;
% Ranking according to grand median Error of each subject across both styles.
subjlistMed = [14 10 13 16 6 5 11 8 4 1 3 7 2 15 9 12]; 

% Align by the end.
for irow = 1:length(srcTimeSeries)
    Tmax(irow,1) = srcTimeSeries(irow).HSTimeMean(end);
end
for irow = 1:length(srcTimeSeries)
    srcTimeSeries(irow).HSTimeMean = srcTimeSeries(irow).HSTimeMean + max(Tmax) - Tmax(irow);
end

% Plot HandSpeed
curFig = Create_Reuse_Figure([],'Hand Speed in a Trial',[1921 52 850 240]);

% Create tile of plots. I don't like tiledLayout for its incompatibility, so here is a
% more customizable version.
optionsIn.SpacingsOut = [10 30 40 40]; % [Right Top Left Bottom] in pixels. Default ~60
%os.SpacingsOutStr = 'Right Top Left Bottom in px'; %. Default ~60
optionsIn.SpacingsBetween = [25 10];
%os.SpacingsBetweenStr = 'Vert Horiz';
optionsIn.YLabels = 'BotLeftCorn'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn"
optionsIn.XLabels = 'BotLeftCorn'; 
optionsIn.XTickMarks = 'BotLeftCorn'; % Same
optionsIn.YTickMarks = 'BotLeftCorn'; 
[sp, optionsOut] = TiledAxesMy(curFig,2,8,[0 0 1 1],optionsIn);
sp(9).XLabel.String = 'Time (s)';
sp(9).YLabel.String = 'Hand Speed (m/s)';

for iisubj = 1:16
    title(sp(iisubj),sprintf('P%d',iisubj));
    isubj = subjlistMed(iisubj);
    
    % Plot SD patches
    for istyle = 1:2
        irow = istyle + (isubj-1)*2;
        X = srcTimeSeries(irow).HSTimeMean;
        X = [X; flip(X)];
        
        Ym = srcTimeSeries(irow).HSmean;
        Ys = srcTimeSeries(irow).HSstd;
        Y = [Ym + Ys; flip(Ym - Ys)];
        patch(sp(iisubj),'XData',X,'YData',Y,'FaceColor',rgbmystyle(istyle,:),'FaceAlpha', 0.5,'EdgeColor','none');
    end
    
    % Plot mean lines
    for istyle = 1:2
        irow = istyle + (isubj-1)*2;
        X = srcTimeSeries(irow).HSTimeMean;
        Y = srcTimeSeries(irow).HSmean;
        plot(sp(iisubj),X,Y,'Color',rgbmystyle(istyle,:),'LineWidth',2);
    end
end
linkaxes(sp,'xy');
xlim(sp,[0 2]);
ylim(sp,[0 8]);
end





