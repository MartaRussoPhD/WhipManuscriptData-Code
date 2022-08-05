function WhipTask_WSAveragePlot(srcTimeSeries,rgbmystyle)

% Using previously averaged data, this function plots Whip Marker Speed (Including Hand) for each participant and
% style.
% Data were sampled from time interval ThrowOnset--MinDist and XAxis is time-normalized.
% Average speed profile is obtained from time-normalizing and averaging trial time-series
% SD patch is obtained the same way.

Nresamp = 200;
SampleRate = 500;
StyleNames = {'Discrete','Rhythmic'};
% Ranking according to grand median Error of each subject across both styles.
subjlistMed = [14 10 13 16 6 5 11 8 4 1 3 7 2 15 9 12]; 


%Color for markers - 5 distant same as in QTM, then faded (HSV: -S -V)
colWhip = [1 0 1; 0.85 0.325 0.098; 0.929 0.694 0.125; 0 1 0; 0 0.447 0.7410];
tempCol = rgb2hsv(colWhip);
tempCol = tempCol .* [1,0.5,0.8];
tempCol(tempCol>1) = 1;
colWhip(6:10,:) = hsv2rgb(tempCol);
% Hand color dark green
colWhip(11,:) = [0.2 0.4 0];


% Normalize time.
for irow = 1:length(srcTimeSeries)
    srcTimeSeries(irow).WSTimeMean = srcTimeSeries(irow).WSTimeMean ./ srcTimeSeries(irow).WSTimeMean(end);
end




for istyle = 1:2
% Plot WhipSpeed
curFig(istyle) = Create_Reuse_Figure([],sprintf('Whip Speed during Throw %s',StyleNames{istyle}),...
    [1921 52 + 400 * (istyle-1) 1000 240]);

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
[sp, optionsOut] = TiledAxesMy(curFig(istyle),2,8,[0 0 1 1],optionsIn);
sp(9).XLabel.String = 'Throw time (a.u.)';
sp(9).YLabel.String = 'Marker Speed (m/s)';

for iisubj = 1:16
    title(sp(iisubj),sprintf('P%d',iisubj));
    isubj = subjlistMed(iisubj);
    
    irow = istyle + (isubj-1)*2;
    X = srcTimeSeries(irow).WSTimeMean;
    X = [X; flip(X)];
    % Plot SD patches
    for imark = 1:11
        Ym = srcTimeSeries(irow).WSmean(:,imark);
        Ys = srcTimeSeries(irow).WSstd(:,imark);
        Y = [Ym + Ys; flip(Ym - Ys)];
        patch(sp(iisubj),'XData',X,'YData',Y,'FaceColor',colWhip(imark,:),'FaceAlpha', 0.25,'EdgeColor','none');
    end
    
    % Plot mean lines
    X = srcTimeSeries(irow).WSTimeMean;
    for imark = 1:11
        Y = srcTimeSeries(irow).WSmean(:,imark);
        plot(sp(iisubj),X,Y,'Color',colWhip(imark,:),'LineWidth',1);
    end
end
linkaxes(sp,'xy');
%xlim(sp,[0 2]);
ylim(sp,[0 50]);
end





