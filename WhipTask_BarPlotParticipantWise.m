function WhipTask_BarPlotParticipantWise(ParTab,rgbmystyle,ParamName,varargin)

% Using the per-trial data, this function creates a barplot for 16 subjects in 2 styles (Discrete and Rhythmic)
% Barplots (and error-bars) are computed as mean and SD values across 5 block means (medians with flag) for each subject-style
% Participants are ranked according to grand-median error of each participant, from best to worst.


indPar = find(strcmpi(ParTab.Properties.VariableNames,ParamName));
if isempty(indPar)
    disp('Parameter not found in the table. Check name spelling');
    return
end


ParName = 'Par name (unit)';
if any(strcmpi(varargin,'ParameterNameUnit'))
    ParName = varargin{find(strcmpi(varargin,'ParameterNameUnit'))+1};
end

flag_median = 0;
if any(strcmpi(varargin,'NotNormal')) || any(strcmpi(varargin,'NonNormal'))
    flag_median = 1;
end


flag_Ranking = 0;
string_Ranking = '';
if any(strcmpi(varargin,'Ranking'))
    string_Ranking = varargin{find(strcmpi(varargin,'Ranking'))+1};
end
if ~isempty(string_Ranking)
    if any(strcmpi(string_Ranking,{'D','R'}))
        flag_Ranking = 1 * strcmpi(string_Ranking,'D') + 2 * strcmpi(string_Ranking,'R');
    else
        disp('Ranking string can be ''R'' or ''D''');
        return
    end
end
subjlistMedDR = [14 10 13 16 6 5 11 8 4 1 3 7 2 15 9 12]; %ranking accuracy-based
subjlistMedR = [10 13 14 5 6 7 16 15 3 1 11 4 8 9 2 12];
subjlistMedD = [16 11 14 2 10 13 8 4 6 3 1 5 7 15 9 12];
subjlistOld = [14 15 19 20 24 21 22 23 16 17 11 12 13 110 25 26]; %recording numbers

switch flag_Ranking
    case 0
        subjlistMed = subjlistMedDR;
    case 1
        subjlistMed = subjlistMedD;
    case 2
        subjlistMed = subjlistMedR;
end


if flag_median
    BlockMeanName = 'Block Median';
    BlockVarName = 'Block IQR';
    BinMeanName = 'Bin Median';
    BinVarName = 'Bin IQR';
    MeanName = 'Median';
    SDName = 'IQR';
else
    BlockMeanName = 'Block Mean';
    BlockVarName = 'Block Std';
    BinMeanName = 'Bin Mean';
    BinVarName = 'Bin Std';
    MeanName = 'Mean';
    SDName = 'SD';
end



ic(1) = find(strcmpi(ParTab.Properties.VariableNames,'Subj'));
ic(2) = find(strcmpi(ParTab.Properties.VariableNames,'Subj_Ranked'));
ic(3) = find(strcmpi(ParTab.Properties.VariableNames,'StyleDR'));
ic(4) = find(strcmpi(ParTab.Properties.VariableNames,'Block'));

T = ParTab(:,ic); % Subj, Subj_Ranked,StyleDR,Block

T.Y = ParTab.(ParamName);
T.Error = ParTab.MinDist;
T.StyleDR = T.StyleDR - 1; % To make consistent with R.

for iiS = 1:16
    iS = subjlistMedDR(iiS);
    for istyle = 1:2
        for iblock = 1:5
            T1 = T(T.StyleDR == (istyle-1) & T.Subj == iS & T.Block == iblock,:);
            ParBlockMea(iS,istyle,iblock) = mean(T1.Y,'omitnan');
            if flag_median
                ParBlockMea(iS,istyle,iblock) = median(T1.Y,'omitnan');
            end
        end
        T1 = T(T.StyleDR == (istyle-1) & T.Subj == iS,:);
        ParBPTrialM(iS,istyle) = mean(T1.Y,'omitnan');
        ParBPTrialS(iS,istyle) = std(T1.Y,[],'omitnan');
        if flag_median
            ParBPTrialM(iS,istyle) = median(T1.Y,'omitnan');
            ParBPTrialS(iS,istyle) = iqr(T1.Y);
        end
    end
end

%% Plot parameter for each participant, barplot, sorting by D ranking
% 5 block-med or block-mean values used for each bar
%%%%%%%%%%% OR ALL TRIAL VALUES USED FOR EACH BAR?


FigBar = figure('Position',[300 200 551 209]);
FigBar.Name = 'Barplot ranked by accuracy';
spBPRanked = axes(FigBar,'position',[0.1 0.25 0.8 0.7]); hold on;
axes(spBPRanked);
ylabel(sprintf('%s Style Average',ParamName),'fontweight','bold','FontSize',10);
XL = xlabel('Participant ranked','fontweight','bold');
XL.Units = 'normalized';
XL.Position(1:2) = [0.3 -0.15];

% RANK THEM FIRST
%Prepare data for barplot. Have 16 points per block-style
YBPplot = permute(ParBlockMea,[3 1 2]);
YBPplot = YBPplot(:,subjlistMed,:);


bP = bar(spBPRanked,squeeze(mean(YBPplot,1,'omitnan')),'FaceColor','flat');
% OR USE TRIAL VALUES?
%bP = bar(spBPRanked,ParBPTrialM,'FaceColor','flat');

bP(1).FaceColor = rgbmystyle(1,:);
bP(2).FaceColor = rgbmystyle(2,:);
lg = legend({'Discrete','Rhythmic'});

% Add errorbars
ydat(:,1) = bP(1).YEndPoints;
ydat(:,2) = bP(2).YEndPoints;
xdat(:,1) = bP(1).XEndPoints;
xdat(:,2) = bP(2).XEndPoints;

eB = errorbar(xdat,ydat,squeeze(std(YBPplot,[],1,'omitnan')),'LineStyle','none','HandleVisibility','off');
% OR USE TRIAL VALUES?
%eB = errorbar(xdat,ydat,ParBPTrialS,'LineStyle','none','HandleVisibility','off');

eB(1).Color = 'k';
eB(2).Color = 'k';

ylim(AutoLims(YBPplot,spBPRanked,'Y'));

lg.Orientation = 'horizontal';
lg.Position(1:2) = [0.59 0.02];


end


