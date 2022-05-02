function [ind, iind] = gr_ind_friendly(subjparamsMarkers,subjparamsDerived,varargin)
%[ind, iind] = gr_ind_friendly(subjparamsMarkers,subjparamsDerived,'SubjNum',[1:16],'StyleDR',1:2,'Block',1:5,'Discard',0:1,'Hit',0:1)    
%For more parameters, look inside;
    
%trialstructure, subjnum, Subjnum, SexFM, DR/RD, Exp-Arm, Exp-Throw, Exp-list, style, block, suspicious, edited, duration, hit
%indEXP = [{soccer},{tennis},{zumba},{karate}...]

%This function still does a loop, not a vector... I could not find at time
%the option to do vector with a structure. Should be doable now, but won't
%affect performance a lot.



%Varargin assignments
for i = 1:1
    
    
    if any(strcmpi(varargin,'IndStartSearch')) %UNLIKELY TO BE USED
        indstart = varargin{find(strcmpi(varargin,'IndStartSearch'))+1};
    else
        indstart = [];
    end
    if any(strcmpi(varargin,'SubjNum'))
        indSubj = varargin{find(strcmpi(varargin,'SubjNum'))+1};
    else
        indSubj = [];
    end
    if any(strcmpi(varargin,'SexFM'))
        indSEXFM = varargin{find(strcmpi(varargin,'SexFM'))+1};
    else
        indSEXFM = [];
    end
    if any(strcmpi(varargin,'OrderDRRD'))
        indDRRD = varargin{find(strcmpi(varargin,'OrderDRRD'))+1};
    else
        indDRRD = [];
    end
    if any(strcmpi(varargin,'ExperArm'))
        indEXPARM = varargin{find(strcmpi(varargin,'ExperArm'))+1};
    else
        indEXPARM = [];
    end
    if any(strcmpi(varargin,'ExperThrow'))
        indEXPTHR = varargin{find(strcmpi(varargin,'ExperThrow'))+1};
    else
        indEXPTHR = [];
    end
    if any(strcmpi(varargin,'ExperList')) % should be a cell list {'tennis','volley'}....
        indEXP = varargin{find(strcmpi(varargin,'ExperList'))+1};
    else
        indEXP = [];
    end
    if any(strcmpi(varargin,'StyleDR'))
        indStyle = varargin{find(strcmpi(varargin,'StyleDR'))+1};
    else
        indStyle = [];
    end
    if any(strcmpi(varargin,'Block'))
        indBlock = varargin{find(strcmpi(varargin,'Block'))+1};
    else
        indBlock = [];
    end
    if any(strcmpi(varargin,'Suspicious')) %UNLIKELY TO BE USED
        indSusp = varargin{find(strcmpi(varargin,'Suspicious'))+1};
    else
        indSusp = [];
    end
    if any(strcmpi(varargin,'Discard'))
        indDiscard = varargin{find(strcmpi(varargin,'Discard'))+1};
    else
        indDiscard = [];
    end
    if any(strcmpi(varargin,'Edit')) %UNLIKELY TO BE USED
        indEdit = varargin{find(strcmpi(varargin,'Edit'))+1};
    else
        indEdit = [];
    end
    if any(strcmpi(varargin,'Duration')) %UNLIKELY TO BE USED
        indDur = varargin{find(strcmpi(varargin,'Duration'))+1};
    else
        indDur = [];
    end
    if any(strcmpi(varargin,'Hit'))
        indHit = varargin{find(strcmpi(varargin,'Hit'))+1};
    else
        indHit = [];
    end
    if any(strcmpi(varargin,'OHeadPSag'))
        indOHeadPSag = varargin{find(strcmpi(varargin,'OHeadPSag'))+1};
    else
        indOHeadPSag = [];
    end
    
end

%% List of parameters
% indstart - not used, leave empty. It is the starting set of indices. Just
% in case you already have a subset of indices and you want to further look
% within that subset and care about that little performance increase...
% Subject number
% SEX (0F 1M)
% DRRD order (1 DR 2RD)
% EXPARM experience with torso (0 none, 1 some not throwing, 2 throwing or racket or bat or whip), don't change
% EXPTHR experience with throwing (0 none, 1 throwing or racket or bat, 2 whip), Look at 0:1 (not whip) and 2 separately.
% EXP list of sport experiences, don't change
% BLK block number
% STLDR style (1 D 2 R)
% SUSP suspicious (0/1), don't change
% EDIT begin/end manually edited (0/1), don't change
% SZ size of data saved - with extra ~150 frames.
% TTD - tip to target distance. It is the minimal one of 12 distances (w1..w2 with one virtual ball between them) to (Targ1..Targ2 with two virtual balls between them)
%       won't work as a parameter here, since it's a vector
% WTDmin - minimal TTD during the trial. I subtracted 10 mm of whip tip thickness and 20 mm of target thickness. These are not (always) exact, so negative values might be there.
% iWTDmin - frame when mTTD was registered.
% wWTDmin and TmWTDmin - between which of 3 whip and 4 target markers was the minDist registered. WmTTD 1 tip w1, 2 virtual, 3 w2;   TmTTD 1 Targ-1, 2 virtual-1, 3 virtual-2, 4 Target-2.
% shouldn't matter now.
% 
% iHit - frame of hit. Contains error and might even have a user error of +-10 frame. Use only as (isnan) check for Hit. 
% DUR - the smallest iWTDmin
%% 
% 
% Markers - substructure with all the markers but Target. Its length is not trimmed to DUR and is roughly DUR+150, just in case.  When you extract those for speed, take mark(:,1:DUR). We don't want the data for Rhythmic trials to overlap between consequent trials.


% Not all the fields from the table could be indexed with this function.
% But it's unlikely that we start pulling out, say, only those where the
% whip hit the Target-2 marker....  To count them (maybe), a simple loop
% would be enough.



%The function takes output to find all row-numbers (Trials) and their total
%number, that corresponds to parameters in brackets. [] means any.

%e.g. [ind,iind] = groupind(subjparams,2:3,[],[],[],[],[],2,2:4,[],[],[],1)
%will give you a list of trial numbers and its size for these conditions:
%subjects 2 and 3 (numeration as in the large table); Rhythmic style 2,3,4 blocks; Only when a hit was recorded via oscillations



%subjpar = subjparams.subjparams;
subjparams = subjparamsDerived;
%subjparams.AllTrials.Markers = subjparamsMarkers.AllTrials.Markers;

if ~isempty(subjparamsMarkers) && isfield(subjparamsMarkers.AllTrials,'Markers')
for i = 1:length(subjparamsMarkers.AllTrials)
    subjparams.AllTrials(i).Markers = subjparamsMarkers.AllTrials(i).Markers; 
end
end

ntrials = size(subjparams.AllTrials,2);
if isempty(indstart)
    indstart = 1:ntrials; 
end
ind = [];
iind = 0;
for itrial = indstart
      
    if (all(isempty(indSubj)) || all(ismember(subjparams.AllTrials(itrial).S,indSubj))) &...
            (all(isempty(indSEXFM)) || all(ismember(subjparams.AllTrials(itrial).SEXFM,indSEXFM))) &...
            (all(isempty(indDRRD)) || all(ismember(subjparams.AllTrials(itrial).DRRD,indDRRD))) &...
            (all(isempty(indEXPARM)) || all(ismember(subjparams.AllTrials(itrial).EXPARM,indEXPARM))) &...
            (all(isempty(indEXPTHR)) || all(ismember(subjparams.AllTrials(itrial).EXPTHR,indEXPTHR))) &...
            (all(isempty(indEXP)) || all(ismember(subjparams.AllTrials(itrial).EXP,indEXP))) &...   % MAYBE LOOK FOR INTERSECTIONS?
            (all(isempty(indStyle)) || all(ismember(subjparams.AllTrials(itrial).STLDR,indStyle))) &...
            (all(isempty(indBlock)) || all(ismember(subjparams.AllTrials(itrial).BLK,indBlock))) &...
            (all(isempty(indSusp)) || all(ismember(subjparams.AllTrials(itrial).SUSP,indSusp))) &...
            (all(isempty(indEdit)) || all(ismember(subjparams.AllTrials(itrial).EDIT,indEdit))) &...
            (all(isempty(indDiscard)) || all(ismember(subjparams.AllTrials(itrial).Discard,indDiscard))) &...
            (all(isempty(indDur)) || all(ismember(subjparams.AllTrials(itrial).DUR,indDur))) &...
            (all(isempty(indHit)) || indHit == ~isnan(subjparams.AllTrials(itrial).iHit)) &...
            (all(isempty(indOHeadPSag)) || all(ismember(subjparams.AllTrials(itrial).Ohead_Psag_01,indOHeadPSag)))
   
        %disp(itrial);
        
        iind = iind + 1;
        ind(iind) = itrial;
    end
end
ind = ind';

end

