function CorrModelOut = CorrParameters_Inter_Intra(ParTab,rgbmystyle,ParamName,varargin)
% 1. Runs a simple linear model on the parameter mean and the error median (one value per
% subject per block). Specify 'NonNormal' to use median for the parameter. Chooses the
% model with or without interaction, by ML-criterion (chi-sq statistic with the F-test).
% 2. Runs a mixed-effect linear model on all (non-discarded) trials. The models are fitted
% iteratively. a) testing whether random effects (intercepts) are needed. b) testing
% whether fixed interaction is needed. c) if b) confirms, continuing with more random
% effect terms, up to random interaction. Iteratively compares via ML-criterion.
% The model handles binary parameters well.
%%%% ADD support for only-one-style parameters (the rest nan).
%%%%% ADD EMBEDDED SIMPLE TEST FOR NORMALITY? To justify using mean/median for simple LM?..
% ______ Returns a structure with an LM object, an LMER object, LMER coefficient test table (t-test, using
% satterthwaite method for dofs, 95% CI included), random coefficients for participants
% sorted by performance (specify 'RankingType', [default 'DR', 'D', 'R']), covariance
% coefficients for the random terms.

warning('off','MATLAB:table:RowsAddedExistingVars');

flag_output_graph = 1;
if any(strcmpi(varargin,'SupressGraphics'))
    flag_output_graph = 0;
end
flag_output_text = 1;
if any(strcmpi(varargin,'SupressText'))
    flag_output_text = 0;
end

flag_median = 0;
if any(strcmpi(varargin,'NonNormal'))
    flag_median = 1;
end

RankingType = 'DR';
if any(strcmpi(varargin,'RankingType'))
    RankingType = varargin{find(any(strcmpi(varargin,'RankingType')))+1};
end

ParNameY = 'MinDist';
if any(strcmpi(varargin,'DependentParName'))
    ParNameY = varargin{find(any(strcmpi(varargin,'DependentParName')))+1};
end


T = ParTab(:,[2 3 8 9]);
T.X = ParTab.(ParamName);
T.Y = ParTab.(ParNameY);
T.Error = ParTab.MinDist;


flag_OnlyOneStyle = 0;
T.StyleDR = T.StyleDR - 1; % To make consistent with R.
%%%% IS THE ENTRY BINARY?
ParamsNames = {ParamName,ParNameY,T};

CorrModelOut.T = T;

if flag_output_text
    disp(sprintf('Testing Correlation of %s with %s, Between- and Within-participant',ParamName,ParNameY));
end

%% Inter-participant (averaged parameter and Error) via Linear Model
TaggMea = varfun(@(x)mean(x,'omitnan'),T,'InputVariables',{'X','Y'},...
    'GroupingVariables',{'Subj','StyleDR'});
TaggMed = varfun(@(x)median(x,'omitnan'),T,'InputVariables',{'X','Y'},...
    'GroupingVariables',{'Subj','StyleDR'});
TaggMed.X = TaggMed.Fun_X;
TaggMed.Y = TaggMed.Fun_Y;

curFig = [];
if flag_output_graph
    figname = sprintf('Corr %s vs %s',ParamsNames{1},ParamsNames{2});
    figpos = [2000 100 1280 360];
    curFig = Create_Reuse_Figure([],figname,figpos);
end

%% Histogram
if flag_output_graph
    BinMethodStr = 'fd'; % 'scott','fd','intergers','sturges',sqrt'
    NormalizationStr = 'count'; %count, countdensity,cumcount,probability,pdf,cdf

    
    spHistY = subplot('Position',[0.06 0.6 0.18 0.3]); hold on;
    xlabel(ParNameY,'fontweight','bold','Interpreter','none');
    ylabel('Count','fontweight','bold');
    
    for istyle = 1:2
        XPlot = T.Y(T.StyleDR == istyle-1,:);
        [XPlotNoOut, PercRemoved] = RemoveOutliersMy(XPlot,'Agressive4Plot');
        hisY = histogram(XPlotNoOut,'BinMethod',BinMethodStr,'EdgeColor','none','FaceColor',rgbmystyle(istyle,:),...
            'FaceAlpha',0.65,'Normalization',NormalizationStr);
    end
    xline(0,'Color','k','LineStyle','--','LineWidth',.5);
    if strcmpi(NormalizationStr,'pdf')
        spHistY.YLabel.String = 'Probability';
    end
    
    
    spHistX = subplot('Position',[0.06 0.15 0.18 0.3]); hold on;
    xlabel(ParamName,'fontweight','bold','Interpreter','none');
    ylabel('Count','fontweight','bold');
    %spHist2.Position(1) = spHist2.Position(1) - 0.07;
    
    for istyle = 1:2
        XPlot = T.X(T.StyleDR == istyle-1,:);
        [XPlotNoOut, PercRemoved] = RemoveOutliersMy(XPlot,'Agressive4Plot');
        hisX = histogram(XPlotNoOut,'BinMethod',BinMethodStr,'EdgeColor','none','FaceColor',rgbmystyle(istyle,:),...
            'FaceAlpha',0.65,'Normalization',NormalizationStr);
    end
    xline(0,'Color','k','LineStyle','--','LineWidth',.5);
    if strcmpi(NormalizationStr,'pdf')
        spHistX.YLabel.String = 'Probability';
    end
end



% Using median for error averaging by default.
% Use mean for the parameter?
%flag_useParMean = 1;

if ~flag_median
    TaggMed.X = TaggMea.Fun_X;
    if flag_output_text
        disp(sprintf('<strong>Using parameter MEAN for between-participant correlation</strong>'));
    end
else
    if flag_output_text
        disp(sprintf('<strong>Using parameter MEDIAN for between-participant correlation</strong>'));
    end
end
if flag_output_text
    disp(sprintf('\n'));
end

indD = T.StyleDR == 0;
indR = T.StyleDR == 1;

%%%%% Use stepwiselm to do the job for me?
LMFormulas = {'Y ~ X + StyleDR';...
    'Y ~ X * StyleDR'};
if nnz(isnan(T.X(indD))) > 0.8 * nnz(indD) ||...
        nnz(isnan(T.X(indR))) > 0.8 * nnz(indR)
    flag_OnlyOneStyle = 1;
    disp('Only one style present');
    LMFormulas = {'Y ~ X'};
end


H = length(LMFormulas);
dLogLikX2 = nan(H,1);
LRTest_p = nan(H,1);
LMMdof = nan(H,1);
LogLik = nan(H,1);
SSR = nan(H,1);
Selected = nan(H,1);
LMs(H).LM = [];

for ilm = 1:H
    LMs(ilm).LM = fitlm(TaggMed,LMFormulas{ilm});
    LMMdof(ilm) = LMs(ilm).LM.NumObservations - LMs(ilm).LM.DFE;
    LogLik(ilm) = LMs(ilm).LM.LogLikelihood;
    SSR(ilm) = LMs(ilm).LM.SSR;
end
% Compare
for ilm = 2:length(LMFormulas)
    % Detailed comparison
    % LR = 2*(LM1.LogLikelihood - LM0.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
    % LMCompar_pval = 1 - chi2cdf(LR, LM0.DFE - LM1.DFE);
    
    % More compact from toolbox.
    [~,LRTest_p(ilm,1),dLogLikX2(ilm,1),~] = lratiotest(LMs(ilm).LM.LogLikelihood,LMs(ilm-1).LM.LogLikelihood,LMMdof(ilm) - LMMdof(ilm-1));
end

% Pick the most parsimonious model
if ~flag_OnlyOneStyle
    nmin = 2;
    if any(LRTest_p > 0.05)
        nmin = find(LRTest_p > 0.05,1,'first') - 1;
    end
else
    nmin = 1;
end

% DEBUG - FORCE NO INTERACTION
if 1 == 0
    nmin = 1;
    disp('Debug! Forcing no interaction');
end
% DEBUG - FORCE INTERACTION
if 1 == 0
    nmin = 2;
    disp('Debug! Forcing interaction');
end

LM = LMs(nmin).LM;
Selected(nmin) = 1;
LMSelection = table(LMFormulas,LMMdof,SSR,LogLik,dLogLikX2,LRTest_p,Selected);
LMAnova = anova(LM);



CorrModelOut.InterLM = LM;
if ~flag_OnlyOneStyle
    CorrModelOut.InterLM_Interaction_Chi2_dof_P = [dLogLikX2(2) LMMdof(2)-LMMdof(1) LRTest_p(2)];
else
    CorrModelOut.InterLM_Interaction_Chi2_dof_P = [nan nan nan];
end

if flag_output_text
    disp('*****Linear Models attempted*****');
    disp(LMSelection);
    disp('*****Selected model*****');
    disp(LM);
    disp('*****95% CI of coefficients*****');
    disp(coefCI(LM));
    disp('*****ANOVA on Fixed Effects*****');
    disp(LMAnova);
    disp(sprintf('\n\n\n'));
    
end

%%%%% PLOT?
if flag_output_graph
    spCorrInter = subplot('Position',[0.35 0.15 0.2 0.65]); hold on
    ylabel(sprintf('%s',ParNameY),'fontweight','bold','FontSize',14,'Interpreter','none');
    xlabel(sprintf('%s - Participant average',ParamName),'fontweight','bold','Interpreter','none');
    
    % Plot scatters
    for istyle = 1:2
        T1 = TaggMed(TaggMed.StyleDR == istyle-1,:);
        if istyle == 1, MarkStyle = '^'; MarkEdgeStyle = rgbmystyle(istyle,:); MarkFaceStyle = 'none';
        else, MarkStyle = 'o'; MarkEdgeStyle = 'none'; MarkFaceStyle = rgbmystyle(istyle,:);
        end
        sc = scatter(T1.X,T1.Y,20,'MarkerFaceColor',rgbmystyle(istyle,:),'Marker',MarkStyle,...
            'MarkerEdgeColor','none');
        % For special display
        %             sc = scatter(T.mean_Y,T.mean_MinDist,20,rgbmystyle,'filled','Marker',MarkStyle);
    end
    
    % Plot lines
    for istyle = 1:2
        T1 = TaggMed(TaggMed.StyleDR == istyle-1,:);
        if height(LM.Coefficients) == 2
            Rftcoeffs(1) = LM.Coefficients{1,1};
            Rftcoeffs(2) = LM.Coefficients{2,1};
        elseif height(LM.Coefficients) == 3
            Rftcoeffs(1) = LM.Coefficients{1,1} + (istyle-1) * LM.Coefficients{2,1};
            Rftcoeffs(2) = LM.Coefficients{3,1};
        else
            Rftcoeffs(1) = LM.Coefficients{1,1} + (istyle-1) * LM.Coefficients{2,1};
            Rftcoeffs(2) = LM.Coefficients{3,1} + (istyle-1) * LM.Coefficients{4,1};
        end
        XplotLM = [min(T1.Fun_X) - 0.05 * (max(T1.Fun_X) - min(T1.Fun_X))  max(T1.Fun_X) + 0.05 * (max(T1.Fun_X) - min(T1.Fun_X))];
        YplotLM = Rftcoeffs(1) + Rftcoeffs(2) * XplotLM;
        plot(XplotLM,YplotLM,'LineWidth',1,'LineStyle','--','Color',rgbmystyle(istyle,:));
        %[ypred,yci] = predict(LM,TaggMed);
    end
    
    % Adjust XLim properly!
    TT = TaggMed;
    TT.Y = TaggMed.X;
    [TT1,~] = RemoveOutliersMy(TT,[3 3]);
    xlim(AutoLims(TT1.Y,gca,'X'));
    
    
    % Add text
    Intercept = LM.Coefficients.Estimate(1);
    InterceptP = LM.Coefficients.pValue(1);
    iStyleDR = find(strcmpi(LM.Coefficients.Properties.RowNames,'StyleDR'));
    if ~isempty(iStyleDR)
        Style = LM.Coefficients.Estimate(iStyleDR);
        StyleP = LM.Coefficients.pValue(iStyleDR);
    end
    iX = find(strcmpi(LM.Coefficients.Properties.RowNames,'X'));
    if ~isempty(iX)
        Xest = LM.Coefficients.Estimate(iX);
        XP = LM.Coefficients.pValue(iX);
    end
    iInteraction = find(strcmpi(LM.Coefficients.Properties.RowNames,'StyleDR:X'));
    if ~isempty(iInteraction)
        Interaction = LM.Coefficients.Estimate(iInteraction);
        InteractionP = LM.Coefficients.pValue(iInteraction);
    end
    
    TextFont = 'Normal';
    if XP < 0.05
        TextFont = 'bold';
    end
    MainStr = sprintf('X effect: b = %.3f, p=%.4f',Xest,XP);
    TxMain = text(-0.25,1.23,MainStr,'HorizontalAlignment','Left','FontSize',10,'Units','normalized','Color','k','fontweight',TextFont);
    
    TextFont = 'Normal';
    if ~isempty(iInteraction) && InteractionP < 0.05
        TextFont = 'bold';
    end
    if ~isempty(iInteraction)
        InteractStr = sprintf('Interaction: b = %.3f, p=%.4f',Interaction,InteractionP);
    else
        InteractStr = sprintf('No interaction, p(chi2)=%.4f',CorrModelOut.InterLM_Interaction_Chi2_dof_P(3));
    end
    TxInteract = text(-0.25,1.16,InteractStr,'HorizontalAlignment','Left','FontSize',10,'Units','normalized','Color','k','fontweight',TextFont);
    
    InterceptStr = sprintf('b_0=%.3f (p=%.4f)',Intercept,InterceptP);
    TxInter = text(0.52,1.22,InterceptStr,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
    
    if ~isempty(iStyleDR)
        StyleStr = sprintf('b_s=%.3f (p=%.4f)',Style,StyleP);
    else
        StyleStr = sprintf('b_s=0 (p=1)');
    end
    TxStyle = text(0.52,1.15,StyleStr,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
    
    R2Str = sprintf('R^2=%.2f',LM.Rsquared.Adjusted);
    TxR2 = text(0.92,1.16,R2Str,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
    
end


%% Within-participant global correlation, via Mixed Linear Model

% Consider reading about iterative fitting of LMEs. Do we have the right sequence?..
% Barr suggested adding random slopes for highest-degree-combination of within-factors.
% But those would be too perplexing to analyze and would divert focus from the main
% results, so unless I'm doing pure statistics here, nope.

if flag_output_text
    disp(sprintf('<strong>Testing within-participant correlation</strong>'));
    disp(sprintf('\n'));
end
%%%% Test whether random effects need to be included at all.
LMEFormulas = {'Y ~ X + StyleDR';... % checked: fitlme generalizes well even if the formula does not contain RND factors
    'Y ~ X + StyleDR + (1|Subj)'}; % 1 comparison
LMEparN = [4 5];
if flag_OnlyOneStyle
    disp('CorrError_Within::Only one style present');
    LMEFormulas = {'Y ~ X';...
        'Y ~ X + (1|Subj)'};
    LMEparN = [3 4];
end


Msgs1 = {'Random effects will not be included';...
    'Random effects will be included'};


H = length(LMEFormulas);
dLogLikX2 = nan(H,1);
LRTest_p = nan(H,1);
npar = nan(H,1);
LogLik = nan(H,1);
SSR = nan(H,1);
Selected = nan(H,1);
LMMs = struct;
LMMs(H).LME = [];
LMEComparison = table(LMEFormulas,npar,SSR,LogLik,dLogLikX2,LRTest_p,Selected);
for ilme = 1:H
    LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
    LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
    LMEComparison.npar(ilme,1) = LMEparN(ilme);
    LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
end
[~,LMEComparison.LRTest_p(2,1),LMEComparison.dLogLikX2(2,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(ilme-1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1));

nmin = 1;
if LMEComparison.LRTest_p(2,1) < 0.05, nmin = 2; end

%%% DEBUG ONLY:
if 1 == 0
    disp('Debug!!! Force no random effects');
    nmin = 1;
end

LME = LMMs(nmin).LME;
LMEComparison.Selected(nmin) = 1;
Msg = Msgs1{nmin};

if flag_output_text
    disp(LMEComparison);
    disp(sprintf('<strong>%c</strong>',Msg));
    disp(sprintf('\n'));
end

LMEComparison = LMEComparison(nmin,:);
LMEComparison{1,5:7} = nan;

CorrModelOut.IntraLM_Interaction_Chi2_dof_P = [nan nan nan];
% If one style only, and no difference ther, summarize and return
if flag_OnlyOneStyle
    if nmin == 1
        CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
        return
    else % Otherwise, prepare the table for one next step
        LMMs(1) = LMMs(2);
        LMMs(2) = [];
    end
end



%%%% Test whether fixed interaction needs to be included at all
if flag_OnlyOneStyle
    disp('CorrError_Within::Only one style present');
else
    
    LMEFormulasNew = {'Y ~ X * StyleDR';...
        'Y ~ X * StyleDR + (1|Subj)'};
    LMEparNNew = [5 6];
    % 1 comparison, with respective previous!
    % after this comparison, if the choices are 1 and 1, finish.
    Msgs2 = {'Fixed interaction will not be included';...
        'Fixed interaction will be included'};
    
    if nmin == 1, LMEComparison.LMEFormulas(2) = LMEFormulasNew(1); LMEComparison.npar(2) = LMEparNNew(1); LMMs(2).LME = [];
    else
        
        LMEComparison.LMEFormulas(2) = LMEFormulasNew(2); LMEComparison.npar(2) = LMEparNNew(2); LMMs(1).LME = LMMs(2).LME; LMMs(2).LME = [];
    end
    H = height(LMEComparison);
    
    for ilme = 2:H
        LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
        LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
        LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
    end
    [~,LMEComparison.LRTest_p(2,1),LMEComparison.dLogLikX2(2,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(ilme-1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1));
    
    nmin = 1;
    CorrModelOut.IntraLM_Interaction_Chi2_dof_P = [LMEComparison.dLogLikX2(2,1) LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1) LMEComparison.LRTest_p(2,1)];
    if LMEComparison.LRTest_p(2,1) < 0.05, nmin = 2; end
    %%% DEBUG ONLY:
    if 1 == 0
        disp('Debug!!! Force no interaction');
        nmin = 1;
    end
    
    
    LME = LMMs(nmin).LME;
    LMEComparison.Selected(nmin) = 1;
    Msg = Msgs2{nmin};
    
    if flag_output_text
        disp(LMEComparison);
        disp(sprintf('<strong>%c</strong>',Msg));
        disp(sprintf('\n'));
    end
    
    flag_stop = strcmpi(LMEComparison.LMEFormulas{1},LMEFormulasNew{1}) || strcmpi(LMEComparison.LMEFormulas{1},LMEFormulas{1});
    
    if flag_stop
        LMEAnovaFix = anova(LME,'DFMethod','satterthwaite');
        CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
        return
    end
    
    LMEComparison = LMEComparison(nmin,:);
    LMEComparison{1,5:7} = nan;
    
end



if nnz(isnan(T.X(indD))) > 0.8 * nnz(indD) ||...
        nnz(isnan(T.X(indR))) > 0.8 * nnz(indR)
    disp('CorrError_Within::Only one style present');
    LMEFormulasNoInter = {'Y ~ X + (1|Subj)';...
        'Y ~ X + (1 + X|Subj)'};...
        LMEparRandNoInter = [4 6];
    ilme = 2;
    LMEComparison.LMEFormulas(ilme) = LMEFormulasNoInter(2);
    LMEComparison.npar(ilme) = LMEparRandNoInter(2);
    LMEComparison{ilme,3:7} = nan;
    LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
    LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
    LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
    [~,LMEComparison.LRTest_p(ilme,1),LMEComparison.dLogLikX2(ilme,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(1,1));
    nmin = 1;
    if LMEComparison.LRTest_p(ilme,1) < 0.05
        LMEComparison.Selected(ilme,1)  = 1;
        nmin = 2;
    end
    if flag_output_text
        disp(LMEComparison);
        if nmin == 1
            disp(sprintf('<strong>Random slope will not be included</strong>'));
        else
            disp(sprintf('<strong>Random slope will be included</strong>'));
        end
    end
    LME = LMMs(nmin).LME;
    CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
    return
else
    
    % if no interaction, but random good
    LMEFormulasNoInter = {'Y ~ X + StyleDR + (1|Subj)';...
        'Y ~ X + StyleDR + (1 + StyleDR|Subj)';...
        'Y ~ X + StyleDR + (1 + X|Subj)';...
        'Y ~ X + StyleDR + (1 + StyleDR + X|Subj)'};
    LMEparRandNoInter = [5 7 7 10];
    
    % if interaction, and random good
    LMEFormulasInter = {'Y ~ X * StyleDR + (1|Subj)';...
        'Y ~ X * StyleDR + (1 + StyleDR|Subj)';...
        'Y ~ X * StyleDR + (1 + X|Subj)';...
        'Y ~ X * StyleDR + (1 + StyleDR + X|Subj)'};
    LMEparRandInter = [6 8 8 11];
    
    
    
    if strcmpi(LMEComparison.LMEFormulas{1},LMEFormulasNew{2})
        LMEFormulasNew = LMEFormulasInter(2:end); LMEparNNew = LMEparRandInter(2:end);
    else % meaning 'X + StyleDR + (1|Subj)'
        LMEFormulasNew = LMEFormulasNoInter(2:end); LMEparNNew = LMEparRandNoInter(2:end);
    end
    LMMs(1).LME = LME;
    
    % Pre-populate the table
    for ilm = 1:2%length(LMEFormulasNew)
        LMEComparison.LMEFormulas(ilm + 1) = LMEFormulasNew(ilm);
        LMEComparison.npar(ilm + 1) = LMEparNNew(ilm);
        LMEComparison{ilm+1,3:7} = nan;
    end
    
    if flag_output_text
        disp(sprintf('\n'));
        disp('Comparing 2 with 1 and 3 with 1');
    end
    
    % Compute the first two new models
    for ilme = 2:3
        LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
        LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
        LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
        [~,LMEComparison.LRTest_p(ilme,1),LMEComparison.dLogLikX2(ilme,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(1,1));
        if LMEComparison.LRTest_p(ilme,1) < 0.05
            LMEComparison.Selected(ilme,1)  = 1;
        end
    end
    if flag_output_text
        disp(LMEComparison);
    end
    %%%% DEBUG ONLY
    if 1 == 0
        disp('Debug!!! Force no rnd except intercepts');
        CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
        return
    end
    %%%% DEBUG ONLY
    if 1 == 0
        disp('Debug!!! Force only random X or Random StyleDR'); % nmin = 2 or 3.
        nmin = 3;
        LMEComparison(2,:) = LMEComparison(nmin,:);
        LMEComparison(3,:) = [];
        LME = LMMs(nmin).LME;
        CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
        return
    end
    
    
    
    if any(LMEComparison.LRTest_p(2:3,1) < 0.05)
        if flag_output_text
            disp(sprintf('<strong>Random 1-st order terms will be included</strong>'));
        end
        [~,nmin] = max(LMEComparison.dLogLikX2);
        LMEComparison(2,:) = LMEComparison(nmin,:);
        LMEComparison(3,:) = [];
        LME = LMMs(nmin).LME;
    else
        if flag_output_text
            disp(sprintf('<strong>Random 1-st order terms will not be included</strong>'));
        end
        CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
        return
    end
    
    
    % Compute the remaining two models. Prepopulate the table
    for ilm = 3:4
        LMEComparison.LMEFormulas(ilm) = LMEFormulasNew(ilm);
        LMEComparison.npar(ilm) = LMEparNNew(ilm);
        LMEComparison{ilm,3:7} = nan;
    end
    % Compute the remaining two models.
    for ilme = 3:4
        LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
        LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
        LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
        [~,LMEComparison.LRTest_p(ilme,1),LMEComparison.dLogLikX2(ilme,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(ilme-1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1));
        if LMEComparison.LRTest_p(ilme,1) < 0.05
            LMEComparison.Selected(ilme,1)  = 1;
        end
    end
    if flag_output_text
        disp(LMEComparison);
    end
    
    
    
    % Make the final decision
    if LMEComparison.LRTest_p(3,1) >= 0.05 % the first one is not good
        if flag_output_text
            disp(sprintf('<strong>Combination of random slopes will not be included</strong>'));
        end
        LME = LMMs(2).LME;
        %CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType);
    else
        if flag_output_text
            disp(sprintf('<strong>Sum, but not interaction, of random slopes will be included</strong>'));
        end
        LME = LMMs(3).LME;
        %CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType);
    end

    CorrModelOut.IntraLMM = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text);
    
    
    
    
    
    
    if 1 == 0 %flag_output
        % Isn't working yet. I'm not sure how to plot individual participant
        % slopes given that random ones do not always exist.
        
        figure
        subplot(1,2,1); hold on
        ylabel('MinDist (m)','fontweight','bold','FontSize',14);
        xlabel(sprintf('%s',ParamName),'fontweight','bold');
        title('Discrete');
        
        % Plot lines
        istyle = 1;
        T1 = T(T.StyleDR == istyle-1,:);
        CorrModelOut.IntraLMM
        
        if height(LM.Coefficients) == 3
            Rftcoeffs(1) = LM.Coefficients{1,1} + (istyle-1) * LM.Coefficients{2,1};
            Rftcoeffs(2) = LM.Coefficients{3,1};
        else
            Rftcoeffs(1) = LM.Coefficients{1,1} + (istyle-1) * LM.Coefficients{2,1};
            Rftcoeffs(2) = LM.Coefficients{3,1} + (istyle-1) * LM.Coefficients{4,1};
        end
        XplotLM = [min(T1.Fun_X) - 0.05 * (max(T1.Fun_X) - min(T1.Fun_X))  max(T1.Fun_X) + 0.05 * (max(T1.Fun_X) - min(T1.Fun_X))];
        YplotLM = Rftcoeffs(1) + Rftcoeffs(2) * XplotLM;
        plot(XplotLM,YplotLM,'LineWidth',1,'LineStyle','--','Color',rgbmystyle(istyle,:));
        %[ypred,yci] = predict(LM,TaggMed);
        
        subplot(1,2,2); hold on
        xlabel(sprintf('%s',ParamName),'fontweight','bold');
        title('Discrete');
        istyle = 2;
        T1 = TaggMed(TaggMed.StyleDR == istyle-1,:);
        if height(LM.Coefficients) == 3
            Rftcoeffs(1) = LM.Coefficients{1,1} + (istyle-1) * LM.Coefficients{2,1};
            Rftcoeffs(2) = LM.Coefficients{3,1};
        else
            Rftcoeffs(1) = LM.Coefficients{1,1} + (istyle-1) * LM.Coefficients{2,1};
            Rftcoeffs(2) = LM.Coefficients{3,1} + (istyle-1) * LM.Coefficients{4,1};
        end
        XplotLM = [min(T1.Fun_X) - 0.05 * (max(T1.Fun_X) - min(T1.Fun_X))  max(T1.Fun_X) + 0.05 * (max(T1.Fun_X) - min(T1.Fun_X))];
        YplotLM = Rftcoeffs(1) + Rftcoeffs(2) * XplotLM;
        plot(XplotLM,YplotLM,'LineWidth',1,'LineStyle','--','Color',rgbmystyle(istyle,:));
        %[ypred,yci] = predict(LM,TaggMed);
    end
    
end






% My version
% LMEFormulas = {'Y ~ X + StyleDR';...
%     'Y ~ X + StyleDR + (1|Subj)';...
%     'Y ~ X * StyleDR + (1|Subj)';... % If significant, use the rest with *, otherwise with +
%     'Y ~ X * StyleDR + (1 + StyleDR|Subj)';...
%     'Y ~ X * StyleDR + (1 + X|Subj)';...
%     'Y ~ X * StyleDR + (1 + StyleDR + X|Subj)';...
%     'Y ~ X * StyleDR + (1 + StyleDR * X|Subj)'};
%%%% But in this version be careful what is compared to what, it's not merely sequential.


warning('on','MATLAB:table:RowsAddedExistingVars');



end






function LMEOut = DisplayAndPackLMER(CorrModelOut,LME,RankingType,ParamsNames,rgbmystyle,curFig,flag_output_graph,flag_output_text)
% LME object. RankingType can be 'DR','D','R'

[B,Bnames] = randomEffects(LME); % Random effects
Bnames.Value = B;
L = height(B) ./ 16; % Can be 0,16,32,48,64;
RankTypes = {'DR','D','R'};
if nargin < 2
    RankingType = 'DR';
end
subjlistMedDR = [14 10 13 16 6 5 11 8 4 1 3 7 2 15 9 12]; %ranking accuracy-based
subjlistMedR = [10 13 14 5 6 7 16 15 3 1 11 4 8 9 2 12];
subjlistMedD = [16 11 14 2 10 13 8 4 6 3 1 5 7 15 9 12];
switch find(strcmpi(RankingType,RankTypes))
    case 1
        subjlist = subjlistMedDR;
    case 2
        subjlist = subjlistMedD;
    case 3
        subjlist = subjlistMedR;
    otherwise
        error('Wrong ranking type');
end



for isubj = 1:16
    Ranking(isubj,1) = find(isubj == subjlist);
end

Subj = (1:16)';

Torig = CorrModelOut.T;
TStyleMed = varfun(@(x)median(x,'omitnan'),Torig,'InputVariables','Error',...
    'GroupingVariables',{'Subj','StyleDR'});
TSubjMed = varfun(@(x)median(x,'omitnan'),Torig,'InputVariables','Error',...
    'GroupingVariables',{'Subj'});
Error = TSubjMed.Fun_Error;
ErrorD = TStyleMed.Fun_Error(TStyleMed.StyleDR == 0);
ErrorR = TStyleMed.Fun_Error(TStyleMed.StyleDR == 1);


RanCoefs = table(Ranking,Subj,Error,ErrorD,ErrorR);

CoefsErrorLM = [];
switch L
    case 0
        
    case 1 % Intercept only
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj);
        end
        CoefsErrorLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
    case 2 % Intercept and one slope
        SlopeNameRaw = Bnames.Name{2};
        if strcmpi(SlopeNameRaw,'X')
            SlopeName = 'X';
        else
            SlopeName = Bnames.Name{2}(2:end-1);
        end
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*2-1);
            RanCoefs.(SlopeName)(isubj) = Bnames.Value(isubj*2);
        end
        CoefsErrorLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
        CoefsErrorLM.(SlopeName) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName));
    case 3 %Intercept and two slopes
        SlopeName1 = Bnames.Name{2};
        SlopeName2 = Bnames.Name{3};
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*3-2);
            RanCoefs.(SlopeName1)(isubj) = Bnames.Value(isubj*3-1);
            RanCoefs.(SlopeName2)(isubj) = Bnames.Value(isubj*3);
        end
        CoefsErrorLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
        CoefsErrorLM.(SlopeName1) = fitlm(RanCoefs,sprintf('Error ~  %s',SlopeName1));
        CoefsErrorLM.(SlopeName2) = fitlm(RanCoefs,sprintf('Error ~  %s',SlopeName2));
    case 4 %Intercept, two slopes and interaction
        % FOR THE SAKE OF COMPLYING WITH MARTA'S APPROACH (and boosting fixed interaction
        % significance) I leave out the random interaction.
        
        SlopeName1 = Bnames.Name{2};
        SlopeName2 = Bnames.Name{3};
        %SlopeName3 = Bnames.Name{4};
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*4-3);
            RanCoefs.(SlopeName1)(isubj) = Bnames.Value(isubj*4-2);
            RanCoefs.(SlopeName2)(isubj) = Bnames.Value(isubj*4-1);
            RanCoefs.Interaction(isubj) = Bnames.Value(isubj*4);
        end
        CoefsErrorLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
        CoefsErrorLM.(SlopeName1) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName1));
        CoefsErrorLM.(SlopeName2) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName2));
        CoefsErrorLM.Interaction = fitlm(RanCoefs,'Error ~ Interaction');
end
RanCoefs = sortrows(RanCoefs,1);
[psi,mse,stats] = covarianceParameters(LME);
RanCov = stats{1,1};
[betasFixRaw,~,FixEst] = fixedEffects(LME,'DFMethod','satterthwaite');
LMEOut.LME = LME;
LMEOut.FixCoefs = FixEst;
LMEOut.RanCoefs = RanCoefs;
LMEOut.RanCov = RanCov;

if flag_output_text
    disp(sprintf('\n'));
    disp('*****Selected model*****');
    disp(sprintf('<strong>Linear mixed-effects model fit by %s\n</strong>',LME.FitMethod));
    
    % Observations = LME.NumObservations;
    % FixedEffects = LME.Formula.FELinearFormula.Nterms;
    % RandomEffects = LME.Formula.RELinearFormula.Nterms * 16;
    % CovariancePars = LME.Formula.RELinearFormula.Nterms + cumsum(1:LME.Formula.RELinearFormula.Nterms);
    
    disp(sprintf('<strong>Formula</strong>'));
    disp(LME.Formula);
    
    disp(LME.ModelCriterion)
    disp(sprintf('Rsq-adj    %.3f',LME.Rsquared.Adjusted));
    
    disp(FixEst); fprintf('\n');
    disp(sprintf('<strong>Random coefficients per participant (sorted)</strong>'));
    disp(RanCoefs);fprintf('\n');
end

if L > 0
    if flag_output_text
        disp(sprintf('\n'));
        disp(sprintf('<strong>Correlation of random coefficients with ranking</strong>'));
        disp(sprintf('\n'));
        fn = fieldnames(CoefsErrorLM);
        for ifld = 1:length(fn)
            disp(CoefsErrorLM.(fn{ifld}));
            disp(sprintf('\n'));
        end
    end
    LMEOut.CoefsRankingLM = CoefsErrorLM;
end


if flag_output_text
    disp(sprintf('<strong>Covariance of random coefficients</strong>'));
    disp(RanCov);
end




% Visualize coefficients? betawFixRaw and RanCoefs can have different elements..
CoefsTotal = RanCoefs(:,1:2);
CoefsTotal.Intercept = zeros(16,1);
CoefsTotal.StyleDR = zeros(16,1);
CoefsTotal.X = zeros(16,1);
CoefsTotal.Interaction = zeros(16,1);

% Make sure fix effects have standard structure and order. Zero if n/a
iFixIntercept = find(strcmpi(FixEst.Name,'(Intercept)'));
if ~isempty(iFixIntercept)
    CoefsTotal.Intercept = repmat(betasFixRaw(iFixIntercept),16,1);
end
iFixStyle = find(strcmpi(FixEst.Name,'StyleDR'));
if ~isempty(iFixStyle)
    CoefsTotal.StyleDR = repmat(betasFixRaw(iFixStyle),16,1);
end
iFixX = find(strcmpi(FixEst.Name,'X'));
if ~isempty(iFixX)
    CoefsTotal.X = repmat(betasFixRaw(iFixX),16,1);
end
iFixInteraction = find(strcmpi(FixEst.Name,'StyleDR:X'));
if ~isempty(iFixInteraction)
    CoefsTotal.Interaction = repmat(betasFixRaw(iFixInteraction),16,1);
end

% Make sure random effects have standard structure and order. Zero if n/a
iRandIntercept = find(strcmpi(RanCoefs.Properties.VariableNames,'Intercept'));
if ~isempty(iRandIntercept)
    CoefsTotal.Intercept = CoefsTotal.Intercept + RanCoefs.Intercept;
end
iRandStyle = find(strcmpi(RanCoefs.Properties.VariableNames,'StyleDR'));
if ~isempty(iRandStyle)
    CoefsTotal.StyleDR = CoefsTotal.StyleDR + RanCoefs.StyleDR;
end
iRandX = find(strcmpi(RanCoefs.Properties.VariableNames,'X'));
if ~isempty(iRandX)
    CoefsTotal.X = CoefsTotal.X + RanCoefs.X;
end

iRandInteraction = find(strcmpi(RanCoefs.Properties.VariableNames,'Interaction'));
if ~isempty(iRandInteraction)
    CoefsTotal.Interaction = CoefsTotal.Interaction + RanCoefs.Interaction;
end

LMEOut.CoefsForPlot = CoefsTotal;




% Alternative method, to prepare data for plotting -



% Test figure
if flag_output_graph
    %figure(curFig);
    T = ParamsNames{3};
    % Adjust XLim properly!
    X1 = RemoveOutliersMy([T.X],'Agressive4Plot');
    xmin = min(X1);
    xmax = max(X1);
    xplot = [xmin xmax];
    jj = flipud(jet);
    ColorsSubj = jj(round(linspace(1,256,16)),:);
    
    %sp1 = subplot('Position',[0.15 0.15 0.3 0.8]); hold on
    sp1 = subplot('Position',[0.6 0.15 0.15 0.65]); hold on;
    title('Discrete');  
    %xlabel(ParamsNames{1},'Interpreter','none','Fontweight','bold');
    % ylabel(ParamsNames{2},'Interpreter','none','Fontweight','bold');
    
    % Create lines.
    istyle = 0; % To keep it consistent with R
    for iisubj = 1:16
        yplot = CoefsTotal.Intercept(iisubj) + CoefsTotal.StyleDR(iisubj) * istyle + CoefsTotal.X(iisubj) * xplot + CoefsTotal.Interaction(iisubj) * istyle * xplot;
        plot(sp1,xplot,yplot,'Color',ColorsSubj(iisubj,:));
    end
    
    %sp2 = subplot('Position',[0.55 0.15 0.3 0.8]); hold on
    sp2 = subplot('Position',[0.78 0.15 0.15 0.65]); hold on;
    title('Rhythmic');  
    %xlabel(ParamsNames{1},'Interpreter','none','Fontweight','bold');
    yticklabels({});
    yticks('auto');
    % ylabel(ParamsNames{2},'Interpreter','none','Fontweight','bold');
    
    
    % Create lines.
    istyle = 1; % To keep it consistent with R
    for iisubj = 1:16
        yplot = CoefsTotal.Intercept(iisubj) + CoefsTotal.StyleDR(iisubj) * istyle + CoefsTotal.X(iisubj) * xplot + CoefsTotal.Interaction(iisubj) * istyle * xplot;
        plot(sp2,xplot,yplot,'Color',ColorsSubj(iisubj,:));
    end
    
    % Adjust YLim properly!
    Y1 = RemoveOutliersMy([T.Y],'Agressive4Plot');
    ymin = min(Y1);
    ymax = max(Y1);
    linkaxes([sp1 sp2],'xy');
    ylim(sp1,[0 ymax + 0.1 * (ymax-ymin)]);
    
    colormap(jet);
    cb = colorbar('Position',[0.94 0.15 0.01 0.65]);
    cb.Ticks = [0 1];
    cb.TickLabels = {'Worst','Best'};
    cLab = ylabel(cb,'Accuracy Ranked','Units','normalized','Position',[1 0.5 0],'Fontweight','bold');
    
    
    
    xplotD = prctile([LME.Variables.X(LME.Variables.StyleDR == 0)],[2.5 97.5]);
    xplotR = prctile([LME.Variables.X(LME.Variables.StyleDR == 1)],[2.5 97.5]);
    % Plot Fix-Effect Slopes
    % Add average within-participant slopes, via CI
    iIntercept = find(strcmpi(FixEst.Name,'(Intercept)'));
    iStyle = find(strcmpi(FixEst.Name,'StyleDR'));
    iX = find(strcmpi(FixEst.Name,'X'));
    iInteraction = find(strcmpi(FixEst.Name,'StyleDR:X'));
    
    b = zeros(4,1);
    bCI = zeros(4,2);
    if ~isempty(iIntercept)
        b(1) = FixEst.Estimate(iIntercept);
        bCI(1,1) = FixEst.Lower(iIntercept);
        bCI(1,2) = FixEst.Upper(iIntercept);
    end
    if ~isempty(iStyle)
        b(2) = FixEst.Estimate(iStyle);
        bCI(2,1) = FixEst.Lower(iStyle);
        bCI(2,2) = FixEst.Upper(iStyle);
    end
    if ~isempty(iX)
        b(3) = FixEst.Estimate(iX);
        bCI(3,1) = FixEst.Lower(iX);
        bCI(3,2) = FixEst.Upper(iX);
    end
    if ~isempty(iInteraction)
        b(4) = FixEst.Estimate(iInteraction);
        bCI(4,1) = FixEst.Lower(iInteraction);
        bCI(4,2) = FixEst.Upper(iInteraction);
    end
    
    
    % Get x
    %             T = TStyleMea(TStyleMea.StyleDR == 1,:);
    %             xplot(1,[1 3]) = [min(T.mean_Y) max(T.mean_Y)];
    %             T = TStyleMea(TStyleMea.StyleDR == 2,:);
    %             xplot(2,[1 3]) = [min(T.mean_Y) max(T.mean_Y)];
    
    xplot(1,[1 3]) = xplotD;
    xplot(2,[1 3]) = xplotR;
    xplot(:,2) = mean(xplot(:,[1 3]),2);
    
    XMeanByStyle = xplot(:,2);
    
    % Get lines
    yplot = b(1) + b(2) .* [0; 1] + (b(3) + b(4) * [0; 1]) .* xplot;
    
    % Coefficient uncertainties split onto the two styles
    dBH = (bCI(:,2) - b) ./ 2;
    dx = xplot - XMeanByStyle;
    
    ypatch = [yplot flip(yplot,2)];
    % Account for intercept and style effect
    ypatch(:,1:3) = ypatch(:,1:3) - (dBH(1));%; + dBH(2));
    ypatch(:,4:6) = ypatch(:,4:6) + (dBH(1));% + dBH(2));
    % Do Not Account for Style Effect, it can be misleading
    
    % Account for block and interaction effects
    ypatch(:,[1:2 4]) = ypatch(:,[1:2 4]) + (dBH(3) + dBH(4)) .* dx(:,1:3);
    ypatch(:,[3 5:6]) = ypatch(:,[3 5:6]) - (dBH(3) + dBH(4)) .* dx(:,[3 2 1]);
    
    % Plot stuff
    istyle = 1;
    patch(sp1,'XData',[xplot(istyle,:), flip(xplot(istyle,:))],'YData',ypatch(istyle,:),'FaceColor',rgbmystyle(istyle,:),'FaceAlpha',0.35,'EdgeColor','none');
    plot(sp1,xplot(istyle,:),yplot(istyle,:),'Color',rgbmystyle(istyle,:));
    
    istyle = 2;
    patch(sp2,'XData',[xplot(istyle,:), flip(xplot(istyle,:))],'YData',ypatch(istyle,:),'FaceColor',rgbmystyle(istyle,:),'FaceAlpha',0.35,'EdgeColor','none');
    plot(sp2,xplot(istyle,:),yplot(istyle,:),'Color',rgbmystyle(istyle,:));
    
    
    
    % Add text
    istyle = 1;
    Intercept = LMEOut.FixCoefs.Estimate(1);
    InterceptP = LMEOut.FixCoefs.pValue(1);
    iStyleDR = find(strcmpi(LMEOut.FixCoefs.Name,'StyleDR'));
    if ~isempty(iStyleDR)
        Style = LMEOut.FixCoefs.Estimate(iStyleDR);
        StyleP = LMEOut.FixCoefs.pValue(iStyleDR);
    end
    iX = find(strcmpi(LMEOut.FixCoefs.Name,'X'));
    if ~isempty(iX)
        Xest = LMEOut.FixCoefs.Estimate(iX);
        XP = LMEOut.FixCoefs.pValue(iX);
    end
    iInteraction = find(strcmpi(LMEOut.FixCoefs.Name,'StyleDR:X'));
    if ~isempty(iInteraction)
        Interaction = LMEOut.FixCoefs.Estimate(iInteraction);
        InteractionP = LMEOut.FixCoefs.pValue(iInteraction);
    end
    
    TextFont = 'Normal';
    if XP < 0.05
        TextFont = 'bold';
    end
    MainStr = sprintf('X effect: b = %.3f, p=%.4f',Xest,XP);
    TxMain = text(-0.8,1.24,MainStr,'HorizontalAlignment','Left','FontSize',10,'Units','normalized','Color','k','fontweight',TextFont);
    
    TextFont = 'Normal';
    if ~isempty(iInteraction) && InteractionP < 0.05
        TextFont = 'bold';
    end
    if ~isempty(iInteraction)
        InterStr = sprintf('Interaction: b = %.3f, p=%.4f',Interaction,InteractionP);
    else
        InterStr = sprintf('No interaction, p(chi2)=%.4f',CorrModelOut.IntraLM_Interaction_Chi2_dof_P(3));
    end
    TxInter = text(-.8,1.17,InterStr,'HorizontalAlignment','Left','FontSize',10,'Units','normalized','Color','k','fontweight',TextFont);
    
    InterceptStr = sprintf('b_0=%.3f (p=%.4f)',Intercept,InterceptP);
    TxInter = text(0.4,1.24,InterceptStr,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
    
    if ~isempty(iStyleDR)
        StyleStr = sprintf('b_s=%.3f (p=%.4f)',Style,StyleP);
    else
        StyleStr = sprintf('b_s=0 (p=1)');
    end
    TxStyle = text(0.4,1.16,StyleStr,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
    
    R2Str = sprintf('R^2=%.2f',LMEOut.LME.Rsquared.Adjusted);
    TxR2 = text(1.05,1.17,R2Str,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
    
end


end
