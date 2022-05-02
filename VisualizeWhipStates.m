function VisualizeWhipStates(subjparamsDerived,subjparamsMarkers,InterpRes)



styleNames = {'Discrete','Rhythmic'};
% Plot the base for illustration
plotSilhoette_Trial = 1;
SampleRate = 500;
isubj = 14;
istyle = 2;
iblock = 3;
NThrowFrames = 0;
[PointTraj, HDirH, HDirIF, HDirP, iMinDistList, iThrowOnsetList, iHSmaxList, iThrowList, curFig, sp3d, ~, ~] = ...
                ExtractTrajOriTimes_PlotSilhoette(subjparamsDerived,subjparamsMarkers,isubj,istyle,iblock,...
                NThrowFrames,SampleRate,plotSilhoette_Trial,'NoWhipMarkers','NoWhipConn');

% Remove body connectors?
alphaobjects = findobj('Color','k');
for iobj = length(alphaobjects):(-1):1
    delete(alphaobjects(iobj));
end
alphaobjects = [findobj('FaceColor',[0.1 0.4 0.1]); findobj('FaceColor',[0.6 0.1 0.1])];
for iobj = length(alphaobjects):(-1):1
    delete(alphaobjects(iobj));
end

            
            
%%% Main part

CAx = [0 0.12]; % Colors for lines wrt MinDist.
jetf = flipud(jet); % Consider making Log, or sqrt, or .^2

% Make it more contrast
NewVec = linspace(0,1,256).^1.5;
jetf = flipud(interp1(linspace(0,1,256),jet,NewVec));



% Extract 10 whip markers at Throw Onset transform into the Target R.F. (for each trial) in a list of trials
isubj = 1:16;
%istyle = 1;
iblock = 1:5;
[ind, Nind] = gr_ind_friendly([],subjparamsDerived,'SubjNum',isubj,'StyleDR',istyle,'Block',iblock,'Discard',0,'OHeadPSag',0:1);

% Extract landmarks
for iitrial = 1:Nind
    itrial = ind(iitrial);
    iOnset(iitrial,1) = round(subjparamsDerived.AllTrials(itrial).Times.tOnsetThrow * SampleRate);
    iHSmax(iitrial,1) = round(subjparamsDerived.AllTrials(itrial).Times.tHSmax2 * SampleRate);
    iHSmax1(iitrial,1) = round(subjparamsDerived.AllTrials(itrial).Times.tHSmax1 * SampleRate);
    iMinDist(iitrial,1) = round(subjparamsDerived.AllTrials(itrial).DUR);
    MinDist(iitrial,1) = subjparamsDerived.AllTrials(itrial).MinDist;
end

% Make color for lines.
md = MinDist;
md(md > CAx(2)) = CAx(2);
md(md < CAx(1)) = CAx(1);
YColRGB = interp1(linspace(CAx(1),CAx(2),256),jetf,md);
        
LandMarksInds = [iOnset iHSmax iMinDist];
LandmarkInterest = iHSmax;

% Extract Whip markers and target rotation/transform values
for iitrial = 1:Nind
    itrial = ind(iitrial);
    W(:,:,iitrial) = squeeze((subjparamsMarkers.AllTrials(itrial).Markers.Whip(:,:,LandmarkInterest(iitrial)) ./ 1000));
    [Rot(:,:,iitrial),Tr(:,iitrial),~,~,~] = GetStaticMarkerSilhoette_ThrowOnset(subjparamsMarkers.AllTrials(itrial));
end

% Coord transform to the target r.f.
for iitrial = 1:Nind %parfor?
    Pt(:,:,iitrial) = CoordTransform(squeeze(W(:,:,iitrial)),Rot(:,:,iitrial),Tr(:,iitrial)); %Into the target ref frame
    
    dX = Pt(10,:,iitrial) - PointTraj{1}(1,:);
    Pt(:,:,iitrial) = Pt(:,:,iitrial) - dX;
end

% Coord transform to put the whip at the handle?...
iThrowOnsetList(plotSilhoette_Trial);

% make bad trials more transparent?..
AlphaLine = 0.13 - 0.09 * ((md - min(md)) ./ (max(md) - min(md)));

% Plot......
for iitrial = 1:Nind
    xplot = Pt(:,1,iitrial);
    yplot = Pt(:,2,iitrial);
    zplot = Pt(:,3,iitrial);
    WhipP(iitrial) = line(xplot,yplot,zplot,'Color',YColRGB(iitrial,:),'LineWidth',2);
    
%     WhipP(iitrial).Color(4) = 0.1;
    WhipP(iitrial).Color(4) = AlphaLine(iitrial);
end
title(styleNames{istyle});
cb = colorbar;
colormap(jetf);

view(0,90)
% Record a video since alpha channel doesn't really save well
NRevs = 2;
TRev = 10;
FrameRate = 20;
v = VideoWriter('WhipStatesAtHSmax_Discrete','MPEG-4');
v.FrameRate = FrameRate;
open(v);
NFrames = NRevs * TRev * FrameRate;
AngleList = linspace(0,NRevs * 360,NFrames);
View0 = sp3d.View(1);
ElAmp = 10;
ElFreqPerRev = 0.7;
El0 = sp3d.View(2);
ElList = ElAmp * sin(2*pi*ElFreqPerRev/FrameRate/TRev * linspace(0,NFrames,NFrames));
for iframe = 1:NFrames
    sp3d.View(1) = View0 + AngleList(iframe);
    sp3d.View(2) = El0 + ElList(iframe);%pause(0.03);
    drawnow limitrate
    writeVideo(v, getframe(curFig));
end
close(v);



% Have the the ~2000 lines each 10 3D points. Plot with a proper transparency
            
end