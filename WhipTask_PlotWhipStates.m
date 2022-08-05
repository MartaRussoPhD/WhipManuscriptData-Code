function WhipTask_PlotWhipStates(WhipTask_TimeSeries)

% Using whip marker and handle positions measured at the Peak Hand Speed landmark, 
% plots the whip snapshots from every trial, where handle positions are translated to the same coordinate for display.
% Each whip snapshot is colored with respect to Error acquired in the respective trial.

StyleNames = {'Discrete','Rhythmic'};

% Prepare color palette
CAx = [0 0.30]; % Colors for lines wrt MinDist.
jetf = flipud(jet); % Consider making Log, or sqrt, or .^2
% Make it more contrast
NewVec = linspace(0,1,256).^2.5;
jetf = flipud(interp1(linspace(0,1,256),jet,NewVec));


for istyle = 1:2

Pt = WhipTask_TimeSeries(32+istyle).WhipPos_PeakHand;
MinDist = WhipTask_TimeSeries(32+istyle).MinDistList;

% Make color for lines.
md = MinDist;
md(md > CAx(2)) = CAx(2);
md(md < CAx(1)) = CAx(1);
YColRGB = interp1(linspace(CAx(1),CAx(2),256),jetf,md);
% make bad trials more transparent?..
AlphaLine = 0.13 - 0.09 * ((md - min(md)) ./ (max(md) - min(md)));

curFig = Create_Reuse_Figure([],sprintf('Whip_AtHSmax %d',StyleNames{istyle}),[1950 20 + 200*(istyle-1) 900 600]);
sp = subplot(1,1,1); 
xlabel('X (m)','fontweight','bold');
ylabel('Y (m)','fontweight','bold');
zlabel('Z (m)','fontweight','bold');


% Plot......
for iitrial = 1:length(MinDist)
    xplot = Pt(:,1,iitrial);
    yplot = Pt(:,2,iitrial);
    zplot = Pt(:,3,iitrial);
    WhipP(iitrial) = line(xplot,yplot,zplot,'Color',YColRGB(iitrial,:),'LineWidth',2);
    
%     WhipP(iitrial).Color(4) = 0.1;
    WhipP(iitrial).Color(4) = AlphaLine(iitrial);
end
title(StyleNames{istyle});
cb = colorbar;
colormap(jetf);
cb.Ticks = linspace(0,1,4);
cb.TickLabels = linspace(CAx(1),CAx(2),4);
cb.Label.String = 'MinDist (m)';

axis square;
xlim([1 4]);
ylim([-1 2]);
view(0,90)
end

end
