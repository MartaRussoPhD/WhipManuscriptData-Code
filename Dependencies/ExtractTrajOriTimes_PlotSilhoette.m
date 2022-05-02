%% Obtain handle trajectories and orientations, transformed into the Target ref.frame
function [PointTraj, HDirH, HDirIF, HDirP, iMinDistList, iThrowOnsetList, iHSmaxList, iThrowList, curFig, sp3d, ind, Nind, iindHits] = ...
    ExtractTrajOriTimes_PlotSilhoette(subjparamsDerivedOrig,subjparamsMarkers,isubj,istyle,iblock,NThrowFrames,SampleRate,plotSilhoette_Trial,varargin)

FigNameSpec = [];
if any(strcmpi(varargin,'FigNameSpec'))
    FigNameSpec = varargin{find(strcmpi(varargin,'FigNameSpec'))+1};
end


flag_interval = 0;
if any(strcmpi(varargin,'Prep'))
    flag_interval = 1;
elseif any(strcmpi(varargin,'WholeTrial'))
    flag_interval = 2;
end


flag_noWhipMarkers = 0;
if any(strcmpi(varargin,'NoWhipMarkers'))
    flag_noWhipMarkers = 1;
end

flag_noWhipConn = 0;
if any(strcmpi(varargin,'NoWhipConn'))
    flag_noWhipConn = 1;
end



% None means Handle (w10--Hd0); 
% Available Hip (hipR), Sternum, Shoulder, Elbow, Wrist, Hand, Handle, COM, w10....w1
flag_obj = 0; % meaning Handle
if any(strcmpi(varargin,'Object'))
    indvarargin = find(strcmpi(varargin,'Object'));
    if strcmpi(varargin{indvarargin+1},'Wrist')
        flag_obj = 1;
    elseif strcmpi(varargin{indvarargin+1},'Elbow')
        flag_obj = 2;
    elseif strcmpi(varargin{indvarargin+1},'Shoulder')
        flag_obj = 3;
    elseif strcmpi(varargin{indvarargin+1},'w8')
        flag_obj = 4;
    elseif strcmpi(varargin{indvarargin+1},'COM')
        flag_obj = 5;
    elseif strcmpi(varargin{indvarargin+1},'Tip') || strcmpi(varargin{indvarargin+1},'w2')
        flag_obj = 6;
    elseif strcmpi(varargin{indvarargin+1},'w4')
        flag_obj = 7;
    elseif strcmpi(varargin{indvarargin+1},'w3')
        flag_obj = 8;
    elseif strcmpi(varargin{indvarargin+1},'Hand')
        flag_obj = 9;
    elseif strcmpi(varargin{indvarargin+1},'w6')
        flag_obj = 10;    
	 elseif strcmpi(varargin{indvarargin+1},'w9')
        flag_obj = 11;   
	 elseif strcmpi(varargin{indvarargin+1},'w10')
		 flag_obj = 12;
	 elseif strcmpi(varargin{indvarargin+1},'w7')
		 flag_obj = 13;
	 elseif strcmpi(varargin{indvarargin+1},'w5')
		 flag_obj = 14;
	 elseif strcmpi(varargin{indvarargin+1},'w1')
		 flag_obj = 15;
	 elseif strcmpi(varargin{indvarargin+1},'Sternum')
		 flag_obj = 16;
	 elseif strcmpi(varargin{indvarargin+1},'Hip')
		 flag_obj = 17;
     elseif strcmpi(varargin{indvarargin+1},'Hd1')
		 flag_obj = 18;
	 elseif strcmpi(varargin{indvarargin+1},'Handle')
		 flag_obj = 0;
	 else
		 flag_obj = 0;
		 disp(sprintf('<strong>No valid object specified, proceeding with Handle</strong>'));
	 end
end

[ind, Nind] = gr_ind_friendly([],subjparamsDerivedOrig,'SubjNum',isubj,'StyleDR',istyle,'Block',iblock,'Discard',0);
[indHit,~] = gr_ind_friendly([],subjparamsDerivedOrig,'SubjNum',isubj,'StyleDR',istyle,'Block',iblock,'Discard',0,'Hit',1);

iindHits = find(ismember(ind,indHit));


subjparamsDerived = struct;
curFig = [];
sp3d = [];

% This step takes too long, given many Blocks and Ps. On average 5 s for Discrete block,
% and 2 s for Rhythmic. 
% Interpolation?...

% First take all relevant parts of subjparams to sCur in a normal loop, and then run Coord
% transformations in a parfor loop. See if it boosts the speed.
% Yes it did, by about 8 times.
% Using parfor for the second loop, reduced the speed by about 20% more, to 1/10 of
% original time.
ibeg = 1;
iend = 0;
iMinDistList = nan(Nind,1);
iThrowOnsetList = nan(Nind,1);
iHSmaxList = nan(Nind,1);
iThrowList = nan(Nind,1);
for iitrial = 1:Nind
    itrial = ind(iitrial);
    subjparamsDerived(itrial).SZ = subjparamsDerivedOrig.AllTrials(itrial).SZ;
    subjparamsDerived(itrial).Times = subjparamsDerivedOrig.AllTrials(itrial).Times;
    subjparamsDerived(itrial).iMinDist = subjparamsDerivedOrig.AllTrials(itrial).iMinDist;
    subjparamsDerived(itrial).MinDist = subjparamsDerivedOrig.AllTrials(itrial).MinDist;
    subjparamsDerived(itrial).W11COM = subjparamsDerivedOrig.AllTrials(itrial).W11COM;
    subjparamsDerived(itrial).HDIR_Handle = subjparamsDerivedOrig.AllTrials(itrial).HDIR_Handle;
    subjparamsDerived(itrial).HDIR_IndF = subjparamsDerivedOrig.AllTrials(itrial).HDIR_IndF;
    subjparamsDerived(itrial).HDIR_Palm = subjparamsDerivedOrig.AllTrials(itrial).HDIR_Palm;
    subjparamsDerived(itrial).Markers.Hd0 = subjparamsMarkers.AllTrials(itrial).Markers.Hd0;
    subjparamsDerived(itrial).Markers.Hd1 = subjparamsMarkers.AllTrials(itrial).Markers.Hd1;
    subjparamsDerived(itrial).Markers.WrL = subjparamsMarkers.AllTrials(itrial).Markers.WrL;
    subjparamsDerived(itrial).Markers.WrM = subjparamsMarkers.AllTrials(itrial).Markers.WrM;
    subjparamsDerived(itrial).Markers.HandL = subjparamsMarkers.AllTrials(itrial).Markers.HandL;
    subjparamsDerived(itrial).Markers.HandM = subjparamsMarkers.AllTrials(itrial).Markers.HandM;
    subjparamsDerived(itrial).Markers.Elb = subjparamsMarkers.AllTrials(itrial).Markers.Elb;
    subjparamsDerived(itrial).Markers.ShR = subjparamsMarkers.AllTrials(itrial).Markers.ShR;
	 subjparamsDerived(itrial).Markers.hipR = subjparamsMarkers.AllTrials(itrial).Markers.hipR;
	 subjparamsDerived(itrial).Markers.Sternum = subjparamsMarkers.AllTrials(itrial).Markers.Sternum;
	 
    subjparamsDerived(itrial).Markers.Whip = subjparamsMarkers.AllTrials(itrial).Markers.Whip;
    subjparamsDerived(itrial).Markers.Target = subjparamsMarkers.AllTrials(itrial).Markers.Target;
    
    
    iMinDist = subjparamsDerived(itrial).iMinDist;
    iThrowOnset = round(SampleRate * subjparamsDerived(itrial).Times.tOnsetThrow);
    iHSmax = round(SampleRate * subjparamsDerived(itrial).Times.tHSmax2);
    
    
    % If I agree to ignore time duration... I would just rescale this trajectory to.. 100
    % frames, save real duration in a separate list.
    iPrep = 1:iThrowOnset;
    iThrow = iThrowOnset:iMinDist;
    iWholeTrial = 1:iMinDist;
    switch flag_interval
        case 0
            iRange = iThrow;
        case 1
            iRange = iPrep;
        case 2
            iRange = iWholeTrial;
    end
    
    
    % Keep indices
    iMinDistList(iitrial,1) = iMinDist + ibeg;
    iThrowOnsetList(iitrial,1) = iThrowOnset + ibeg;
    iHSmaxList(iitrial,1) = iHSmax + ibeg;
    iend = iend + iMinDist;
    ibeg = ibeg + iMinDist;
    iThrowList(iitrial,1) = length(iRange);
    
    
    plotSilhoetteAtLandmark = iThrowOnset;
    
    
    % Plot if needed
    if iitrial == plotSilhoette_Trial
        sCur = subjparamsDerivedOrig.AllTrials(itrial);
        sCur.Markers = subjparamsMarkers.AllTrials(itrial).Markers;
        [~,~,obj,objConn,objLabels] = GetStaticMarkerSilhoette_ThrowOnset(sCur,plotSilhoetteAtLandmark,iMinDist,iHSmax,flag_noWhipMarkers,flag_noWhipConn);
        [curFig, sp3d] = PlotStaticMarkersConnectors3D(obj,objConn,'TargetRF','CustomView',[-0.5 0.5; ],'FullScreen','FigNameSpec',FigNameSpec);
    end

end





    
    
    
%     sCur.Markers.Hd0
% sCur.Markers.Hd1
% sCur.Markers.WrL
% sCur.Markers.WrM
% sCur.Markers.Elb
% sCur.Markers.ShR
% sCur.Markers.Whip
% sCur.Markers.Target
% 
% sCur.SZ;
% sCur.Times;
% sCur.iMinDist
% 
% sCur.HDIR_Handle
% sCur.HDIR_IndF
% sCur.HDIR_Palm





PointTraj = nan(NThrowFrames,3,Nind);
HDirH = nan(NThrowFrames,3,Nind);
HDirIF = nan(NThrowFrames,3,Nind);
HDirP = nan(NThrowFrames,3,Nind);

if NThrowFrames == 0
    PointTrajC = cell(Nind,1);
    HDirHC = cell(Nind,1);
    HDirIFC = cell(Nind,1);
    HDirPC = cell(Nind,1);
end
%RotTarg = nan(3,3,Nind); Unused for now.
%TrTarg = nan(3,Nind);




parfor iitrial = 1:Nind
    itrial = ind(iitrial);
    sCur = subjparamsDerived(itrial);
    %sCur.Markers = subjparamsMarkers.AllTrials(itrial).Markers;

    
    iMinDist = sCur.iMinDist;
    iThrowOnset = round(SampleRate * sCur.Times.tOnsetThrow);
    iHSmax = round(SampleRate * sCur.Times.tHSmax2);
    
    iRange = nan;
    Pt = nan;
    iPrep = 1:iThrowOnset;
    iThrow = iThrowOnset:iMinDist;
    iWholeTrial = 1:iMinDist;
    switch flag_interval
        case 0
            iRange = iThrow;
        case 1
            iRange = iPrep;
        case 2
            iRange = iWholeTrial;
    end
    
    [Rot,Tr,obj,objConn,objLabels] = GetStaticMarkerSilhoette_ThrowOnset(sCur);
    
    %RotTarg(:,:,iitrial) = Rot;
    %TrTarg(:,iitrial) = Tr;
    
    Hd0 = squeeze(sCur.Markers.Hd0 ./ 1000)';
    w10 = squeeze(sCur.Markers.Whip(10,:,:) ./ 1000)';
    Hd = .5 * (Hd0 + w10);
    if nnz(isnan(Hd0)) > 0.5 * numel(Hd0)
        Hd = w10;
    end
    if nnz(isnan(w10)) > 0.5 * numel(w10)
        Hd = Hd0;
    end
    
    Hd1 = squeeze(sCur.Markers.Hd1 ./ 1000)';
    
    % Extract other landmarks too
    ShR = squeeze(sCur.Markers.ShR ./ 1000)';
    Elb = squeeze(sCur.Markers.Elb ./ 1000)';
    WrL = squeeze(sCur.Markers.WrL ./ 1000)';
    WrM = squeeze(sCur.Markers.WrM ./ 1000)';
    HandL = squeeze(sCur.Markers.HandL ./ 1000)';
    HandM = squeeze(sCur.Markers.HandM ./ 1000)';
	 Sternum = squeeze(sCur.Markers.Sternum ./ 1000)';
	 hipR = squeeze(sCur.Markers.hipR ./ 1000)';
    
    Wr = .5 * (WrM + WrL);
    if nnz(isnan(WrL)) > 0.5 * numel(WrL)
        Wr = WrM;
    end
    if nnz(isnan(WrM)) > 0.5 * numel(WrM)
        Wr = WrL;
    end
    Hand = .5 * (HandM + HandL)
    if nnz(isnan(HandL)) > 0.5 * numel(HandL)
        Hand = HandM;
    end
    if nnz(isnan(HandM)) > 0.5 * numel(HandM)
        Hand = HandL;
	 end
	 
	 w10 = squeeze(sCur.Markers.Whip(10,:,:) ./ 1000)';
	 w9 = squeeze(sCur.Markers.Whip(9,:,:) ./ 1000)';
    w8 = squeeze(sCur.Markers.Whip(8,:,:) ./ 1000)';
	 w7 = squeeze(sCur.Markers.Whip(7,:,:) ./ 1000)';
	 w6 = squeeze(sCur.Markers.Whip(6,:,:) ./ 1000)';
	 w5 = squeeze(sCur.Markers.Whip(5,:,:) ./ 1000)';
    w4 = squeeze(sCur.Markers.Whip(4,:,:) ./ 1000)';
    w3 = squeeze(sCur.Markers.Whip(3,:,:) ./ 1000)';
	 w2 = squeeze(sCur.Markers.Whip(2,:,:) ./ 1000)';
	 w1 = squeeze(sCur.Markers.Whip(1,:,:) ./ 1000)';
    WCOM = sCur.W11COM;
    
    
	 
    switch flag_obj
        case 0
            Pt = Hd;
        case 1
            Pt = Wr;
        case 2
            Pt = Elb;
        case 3
            Pt = ShR;
        case 4
            Pt = w8;
        case 5
            Pt = WCOM;
        case 6
            Pt = w2;
        case 7
            Pt = w4;
        case 8
            Pt = w3;
        case 9
			  Pt = Hand;
		 case 10
			 Pt = w6
		 case 11
			 Pt = w9;
		 case 12
			 Pt = w10;
		 case 13
			 Pt = w7;
		 case 14
			 Pt = w5;
		 case 15
			 Pt = w1;
		 case 16
			 Pt = Sternum;
		 case 17
			 Pt = hipR;
        case 18
            Pt = Hd1;
    end
    
    Pt = CoordTransform(Pt,Rot,Tr); %Into the target ref frame
    
    
    
    
    % Add variation: if NThrowFrames == 0, do not respline and put stuff in Cell array
    % instead.
    
    
    iThrowRes = linspace(iRange(1),iRange(end),NThrowFrames);
    
    if NThrowFrames == 0 % Do not resample/interpolate
        PointTrajC{iitrial,1} = Pt(iRange,:);
        HDirHC{iitrial,1} = CoordTransform(sCur.HDIR_Handle(iRange,:),Rot,[0 0 0]');
        HDirIFC{iitrial,1} = CoordTransform(sCur.HDIR_IndF(iRange,:),Rot,[0 0 0]');
        HDirPC{iitrial,1} = CoordTransform(sCur.HDIR_Palm(iRange,:),Rot,[0 0 0]');
    else
        PointTraj(:,:,iitrial) = interp1My(iRange,Pt(iRange,:),iThrowRes,'spline');
        HDirH(:,:,iitrial) = CoordTransform(interp1My(iRange,sCur.HDIR_Handle(iRange,:),iThrowRes,'spline'),Rot,[0 0 0]');
        HDirIF(:,:,iitrial) = CoordTransform(interp1My(iRange,sCur.HDIR_IndF(iRange,:),iThrowRes,'spline'),Rot,[0 0 0]');
        HDirP(:,:,iitrial) = CoordTransform(interp1My(iRange,sCur.HDIR_Palm(iRange,:),iThrowRes,'spline'),Rot,[0 0 0]');
    end
    

end

if NThrowFrames == 0
    PointTraj = PointTrajC;
    HDirH = HDirHC;
    HDirIF = HDirIFC;
    HDirP = HDirPC;
end