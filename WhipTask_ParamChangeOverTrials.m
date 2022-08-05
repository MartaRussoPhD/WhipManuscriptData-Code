function ParamChangeModelOut = ParamChangeOverTrials(ParTab,rgbmystyle,ParamName,varargin)

% Runs a mixed-effect linear model on all (non-discarded) trials Parameter ~ Style x
% Block, or a similar GLMM (binomial, logit, LogLik, Laplacian approx.) for a binary parameter.
% The models are fitted iteratively. a) testing whether random effects (intercepts) are needed.
% b) testing whether fixed interaction is needed. c) if b) confirms, continuing with more random
% effect terms, up to random interaction. Iteratively compares via ML-criterion.

% ______ Accepts a table with per-trial parameters (must have Style and Block
% columns), my color palette, parameter name. Optional flags are available (see below)

% ______ Returns a structure with an LM object, an LMER object, LMER coefficient test table (t-test, using
% satterthwaite method for dofs, 95% CI included), random coefficients for participants
% sorted by performance (specify 'RankingType', [default 'DR', 'D', 'R']), covariance
% coefficients for the random terms.



warning('off','MATLAB:table:RowsAddedExistingVars');


indPar = find(strcmpi(ParTab.Properties.VariableNames,ParamName));
if isempty(indPar)
    disp('Parameter not found in the table. Check name spelling');
    return
end



flag_output_graph = 1;
if any(strcmpi(varargin,'SupressGraphics'))
    flag_output_graph = 0;
end
flag_output_text = 1;
if any(strcmpi(varargin,'SupressText'))
    flag_output_text = 0;
end
if any(strcmpi(varargin,'OnlyMainText'))
    flag_output_text = 2;
end

% Used only for boxplot visualization (using either mean or median), not for stats.
flag_median = 0;
if any(strcmpi(varargin,'NotNormal')) || any(strcmpi(varargin,'NonNormal'))
    flag_median = 1;
end

RankingType = 'DR';
if any(strcmpi(varargin,'RankingType'))
    RankingType = varargin{find(any(strcmpi(varargin,'RankingType')))+1};
end

ic(1) = find(strcmpi(ParTab.Properties.VariableNames,'Subj'));
ic(2) = find(strcmpi(ParTab.Properties.VariableNames,'Subj_Ranked'));
ic(3) = find(strcmpi(ParTab.Properties.VariableNames,'StyleDR'));
ic(4) = find(strcmpi(ParTab.Properties.VariableNames,'Block'));

T = ParTab(:,ic); % Subj, Subj_Ranked,StyleDR,Block

T.Block = ParTab.Block;
T.Y = ParTab.(ParamName);
T.Error = ParTab.MinDist;

if flag_output_text
    disp(sprintf('Testing change of %s, between Styles and across Blocks',ParamName));
end

T.StyleDR = T.StyleDR - 1; % To cross-compare with our tests in R

indD = T.StyleDR == 0;
indR = T.StyleDR == 1;

flag_OnlyOneStyle = 0;
if nnz(isnan(T.Y(indD))) > 0.8 * nnz(indD) ||...
        nnz(isnan(T.Y(indR))) > 0.8 * nnz(indR)
    flag_OnlyOneStyle = 1;
end


%%%% IS THE ENTRY BINARY?
flag_bin = 0;
if islogical(T.Y)
    flag_bin = 1;
end
glme_dist = 'binomial';
glme_link = 'logit';

%% Within-participant global correlation, via Mixed Linear Model


% Iterative LMM models fitting Based on Cesqui 2016 and Barr 2013
% 1. Compute and compare
%     1 + X + C
%     1 + X + C + (1|Subj)

% 2. Compare the best(simplest) with interaction
%     1 + X * C
%     1 + X * C + (1|Subj)
%     if no significant, that's it.
%     if significant, continue

% 3. Continue either with or without interaction
%     1 + X *+ C + (1 + X | Subj)
%     1 + X *+ C + (1 + C | Subj)
%     % compare these with the base. If neither is significant, only leave the base
%     % if only one is, pick that one.
%     % if both are significant, pick the largest LogLik
%     1 + X *+ C + (1 + X + C | Subj)
%     % Compare with the previous best. If significantly different, accept.

%%%%%%% If computing effect sizes without random effects,
% Use eqns from Field & Myers "Discovering Statistics with R", p.542


%%%% Test whether random effects need to be included at all.

% Create a table of model formulas and fit metrics
LMEFormulas = {'Y ~ Block + StyleDR';... % checked: fitglme generalizes well even if the formula does not contain RND factors
    'Y ~ Block + StyleDR + (1|Subj)'}; % 1 comparison
LMEparN = [4 5];
if flag_OnlyOneStyle
    disp('Only one style present');
    LMEFormulas = {'Y ~ Block';...
        'Y ~ Block + (1|Subj)'};
    LMEparN = [3 4];
end

ParamsNames = {ParamName,'Error (m)',T};

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
% Fit models
for ilme = 1:H
    if flag_bin
        LMMs(ilme).LME = fitglme(T,LMEComparison.LMEFormulas{ilme},'Distribution',glme_dist,'Link',glme_link,'FitMethod','Laplace');
    else
        LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
    end
    LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
    LMEComparison.npar(ilme,1) = LMEparN(ilme);
    LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
end
% Compare the models
[~,LMEComparison.LRTest_p(2,1),LMEComparison.dLogLikX2(2,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(ilme-1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1));

% Choose the model based on LR-test
nmin = 1;
if LMEComparison.LRTest_p(2,1) < 0.05, nmin = 2; end
LME = LMMs(nmin).LME;
LMEComparison.Selected(nmin) = 1;
Msg = Msgs1{nmin};

if flag_output_text
    if flag_output_text == 1
        disp(LMEComparison);
    end
    disp(sprintf('<strong>%c</strong>',Msg));
    disp(sprintf('\n'));
end
LMEComparison = LMEComparison(nmin,:);
LMEComparison{1,5:7} = nan;

% Initialize the field to report insignificant interaction (if the respective model does
% not improve fit)
ParamChangeModelOut.IntraLM_Interaction_Chi2_dof_P = [nan nan nan];

% If one style only, and no difference there, summarize and return
if flag_OnlyOneStyle
    if nmin == 1
        ParamChangeModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,T,flag_median,flag_OnlyOneStyle,flag_output_graph,flag_output_text,rgbmystyle);
        return
    else % Otherwise, prepare the table for one next step
        LMMs(1) = LMMs(2);
        LMMs(2) = [];
    end
end

%%%% Test whether fixed interaction needs to be included at all
if flag_OnlyOneStyle
    disp('Only one style present');
else
    
    % Table, equations, metrics
    LMEFormulasNew = {'Y ~ Block * StyleDR';...
        'Y ~ Block * StyleDR + (1|Subj)'};
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
        if flag_bin
            LMMs(ilme).LME = fitglme(T,LMEComparison.LMEFormulas{ilme},'Distribution',glme_dist,'Link',glme_link,'FitMethod','Laplace');
        else
            LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
        end
        LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
        LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
    end
    [~,LMEComparison.LRTest_p(2,1),LMEComparison.dLogLikX2(2,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(ilme-1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1));
    
    % Pick the model
    nmin = 1;
    ParamChangeModelOut.IntraLM_Interaction_Chi2_dof_P = [LMEComparison.dLogLikX2(2,1) LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1) LMEComparison.LRTest_p(2,1)];
    if LMEComparison.LRTest_p(2,1) < 0.05, nmin = 2; end
    LME = LMMs(nmin).LME;
    LMEComparison.Selected(nmin) = 1;
    Msg = Msgs2{nmin};
    if flag_output_text
        if flag_output_text == 1
            disp(LMEComparison);
        end
        disp(sprintf('<strong>%c</strong>',Msg));
        disp(sprintf('\n'));
    end
    flag_stop = strcmpi(LMEComparison.LMEFormulas{1},LMEFormulasNew{1}) || strcmpi(LMEComparison.LMEFormulas{1},LMEFormulas{1});
    if flag_stop
        ParamChangeModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,T,flag_median,flag_OnlyOneStyle,flag_output_graph,flag_output_text,rgbmystyle);
        return
    end
    LMEComparison = LMEComparison(nmin,:);
    LMEComparison{1,5:7} = nan;
end


%%%% Test whether random Block and Style effects need to be included
if flag_OnlyOneStyle
    disp('Only one style present');
    LMEFormulasNoInter = {'Y ~ Block + (1|Subj)';...
        'Y ~ Block + (1 + Block|Subj)'};...
        LMEparRandNoInter = [4 6];
    ilme = 2;
    LMEComparison.LMEFormulas(ilme) = LMEFormulasNoInter(2);
    LMEComparison.npar(ilme) = LMEparRandNoInter(2);
    LMEComparison{ilme,3:7} = nan;
    
    if flag_bin
        LMMs(ilme).LME = fitglme(T,LMEComparison.LMEFormulas{ilme},'Distribution',glme_dist,'Link',glme_link,'FitMethod','Laplace');
    else
        LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
    end
    
    LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
    LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
    [~,LMEComparison.LRTest_p(ilme,1),LMEComparison.dLogLikX2(ilme,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(1,1));
    
    % Pick the model
    nmin = 1;
    if LMEComparison.LRTest_p(ilme,1) < 0.05
        LMEComparison.Selected(ilme,1)  = 1;
        nmin = 2;
    end
    if flag_output_text == 1
        disp(LMEComparison);
    end
    if flag_output_text
        if nmin == 1
            disp(sprintf('<strong>Random slope will not be included</strong>'));
        else
            disp(sprintf('<strong>Random slope will be included</strong>'));
        end
    end
    LME = LMMs(nmin).LME;
    ParamChangeModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,T,flag_median,flag_OnlyOneStyle,flag_output_graph,flag_output_text,rgbmystyle);
    return
else
    
    % if no interaction, but random good
    LMEFormulasNoInter = {'Y ~ Block + StyleDR + (1|Subj)';...
        'Y ~ Block + StyleDR + (1 + StyleDR|Subj)';...
        'Y ~ Block + StyleDR + (1 + Block|Subj)';...
        'Y ~ Block + StyleDR + (1 + StyleDR + Block|Subj)'};
    LMEparRandNoInter = [5 7 7 10];
    
    % if interaction, and random good
    LMEFormulasInter = {'Y ~ Block * StyleDR + (1|Subj)';...
        'Y ~ Block * StyleDR + (1 + StyleDR|Subj)';...
        'Y ~ Block * StyleDR + (1 + Block|Subj)';...
        'Y ~ Block * StyleDR + (1 + StyleDR + Block|Subj)'};
    LMEparRandInter = [6 8 8 11];
    
    if strcmpi(LMEComparison.LMEFormulas{1},LMEFormulasNew{2})
        LMEFormulasNew = LMEFormulasInter(2:end); LMEparNNew = LMEparRandInter(2:end);
    else
        LMEFormulasNew = LMEFormulasNoInter(2:end); LMEparNNew = LMEparRandNoInter(2:end);
    end
    LMMs(1).LME = LME;
    
    % Pre-populate the table
    for ilm = 1:2%length(LMEFormulasNew)
        LMEComparison.LMEFormulas(ilm + 1) = LMEFormulasNew(ilm);
        LMEComparison.npar(ilm + 1) = LMEparNNew(ilm);
        LMEComparison{ilm+1,3:7} = nan;
    end
    if flag_output_text == 1
        disp(sprintf('\n'));
        disp('Comparing 2 with 1 and 3 with 1');
    end
    
    % Compute the first two new models
    for ilme = 2:3
        if flag_bin
            LMMs(ilme).LME = fitglme(T,LMEComparison.LMEFormulas{ilme},'Distribution',glme_dist,'Link',glme_link,'FitMethod','Laplace');
        else
            LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
        end
        LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
        LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
        [~,LMEComparison.LRTest_p(ilme,1),LMEComparison.dLogLikX2(ilme,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(1,1));
        if LMEComparison.LRTest_p(ilme,1) < 0.05
            LMEComparison.Selected(ilme,1)  = 1;
        end
    end
    if flag_output_text == 1
        disp(LMEComparison);
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
        ParamChangeModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,T,flag_median,flag_OnlyOneStyle,flag_output_graph,flag_output_text,rgbmystyle);
        return
    end
    
    
    % Compute the remaining model
    ilm = 3;
    LMEComparison.LMEFormulas(ilm) = LMEFormulasNew(ilm);
    LMEComparison.npar(ilm) = LMEparNNew(ilm);
    LMEComparison{ilm,3:7} = nan;
    if flag_bin
        LMMs(ilme).LME = fitglme(T,LMEComparison.LMEFormulas{ilme},'Distribution',glme_dist,'Link',glme_link,'FitMethod','Laplace');
    else
        LMMs(ilme).LME = fitlme(T,LMEComparison.LMEFormulas{ilme});
    end
    
    LMEComparison.SSR(ilme,1) = LMMs(ilme).LME.SSR;
    LMEComparison.LogLik(ilme,1) = LMMs(ilme).LME.LogLikelihood;
    [~,LMEComparison.LRTest_p(ilme,1),LMEComparison.dLogLikX2(ilme,1),~] = lratiotest(LMMs(ilme).LME.LogLikelihood,LMMs(ilme-1).LME.LogLikelihood,LMEComparison.npar(ilme,1) - LMEComparison.npar(ilme-1,1));
    if LMEComparison.LRTest_p(ilme,1) < 0.05
        LMEComparison.Selected(ilme,1)  = 1;
    end
    if flag_output_text == 1
        disp(LMEComparison);
    end
    
    
    
    % Make the final decision
    if LMEComparison.LRTest_p(3,1) >= 0.05 % the most recent model is not good
        if flag_output_text
            disp(sprintf('<strong>Combination of random slopes will not be included</strong>'));
        end
        LME = LMMs(2).LME;
    else % or is good
        if flag_output_text
            disp(sprintf('<strong>Sum, but not interaction, of random slopes will be included</strong>'));
        end
        LME = LMMs(3).LME;
    end
    ParamChangeModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,T,flag_median,flag_OnlyOneStyle,flag_output_graph,flag_output_text,rgbmystyle);
    
    % Plot near-final figure with boxplot Style x Block (using median or mean across each
    % block of each style of each participant)
    if flag_output_graph
        ParamChangeSinglePlot(ParamsNames,ParamChangeModelOut,T,rgbmystyle,flag_OnlyOneStyle,flag_median);
    end
    
    
end
warning('on','MATLAB:table:RowsAddedExistingVars');



end




function LMEOut = DisplayAndPackLMER(LME,RankingType,ParamsNames,ParTab,flag_median,flag_OnlyOneStyle,flag_output_graph,flag_output_text,rgbmystyle)
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


Subj = (1:16)';
for isubj = 1:16
    Ranking(isubj,1) = find(isubj == subjlist);
end

Torig = ParTab;
TStyleMed = varfun(@(x)median(x,'omitnan'),Torig,'InputVariables','Error',...
    'GroupingVariables',{'Subj','StyleDR'});
TSubjMed = varfun(@(x)median(x,'omitnan'),Torig,'InputVariables','Error',...
    'GroupingVariables',{'Subj'});
Error = TSubjMed.Fun_Error;
ErrorD = TStyleMed.Fun_Error(TStyleMed.StyleDR == 0);
ErrorR = TStyleMed.Fun_Error(TStyleMed.StyleDR == 1);

RanCoefs = table(Ranking,Subj,Error,ErrorD,ErrorR);
CoefsRankLM = [];
switch L
    case 0
    case 1 % Intercept only
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
    case 2 % Intercept and one slope
        SlopeName = Bnames.Name{2}(2:end-1);
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*2-1);
            RanCoefs.(SlopeName)(isubj) = Bnames.Value(isubj*2);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
        CoefsRankLM.(SlopeName) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName));
    case 3 %Intercept and two slopes
        SlopeName1 = Bnames.Name{2};
        SlopeName2 = Bnames.Name{3};
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*3-2);
            RanCoefs.(SlopeName1)(isubj) = Bnames.Value(isubj*3-1);
            RanCoefs.(SlopeName2)(isubj) = Bnames.Value(isubj*3);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
        CoefsRankLM.(SlopeName1) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName1));
        CoefsRankLM.(SlopeName2) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName2));
    case 4 %Intercept, two slopes and interaction
        SlopeName1 = Bnames.Name{2};
        SlopeName2 = Bnames.Name{3};
        SlopeName3 = Bnames.Name{4};
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*4-3);
            RanCoefs.(SlopeName1)(isubj) = Bnames.Value(isubj*4-2);
            RanCoefs.(SlopeName2)(isubj) = Bnames.Value(isubj*4-1);
            RanCoefs.Interaction(isubj) = Bnames.Value(isubj*4);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Error ~ Intercept');
        CoefsRankLM.(SlopeName1) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName1));
        CoefsRankLM.(SlopeName2) = fitlm(RanCoefs,sprintf('Error ~ %s',SlopeName2));
        CoefsRankLM.Interaction = fitlm(RanCoefs,'Error ~ Interaction');
end


RanCoefs = sortrows(RanCoefs,1);
[psi,mse,stats] = covarianceParameters(LME);
RanCov = stats{1,1};

if isa(LME,'GeneralizedLinearMixedModel')
    [betasFixRaw,~,FixEst] = fixedEffects(LME,'DFMethod','residual');
else
    [betasFixRaw,~,FixEst] = fixedEffects(LME,'DFMethod','satterthwaite');
end
LMEOut.LME = LME;
LMEOut.FixCoefs = FixEst;
LMEOut.RanCoefs = RanCoefs;
LMEOut.RanCov = RanCov;


% Visualize coefficients? betawFixRaw and RanCoefs can have different elements..
CoefsTotal = RanCoefs(:,1:2);
CoefsTotal.Intercept = zeros(16,1);
CoefsTotal.StyleDR = zeros(16,1);
CoefsTotal.Block = zeros(16,1);
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
iFixBlock = find(strcmpi(FixEst.Name,'Block'));
if ~isempty(iFixBlock)
    CoefsTotal.Block = repmat(betasFixRaw(iFixBlock),16,1);
end
iFixInteraction = find(strcmpi(FixEst.Name,'StyleDR:Block'));
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
iRandBlock = find(strcmpi(RanCoefs.Properties.VariableNames,'Block'));
if ~isempty(iRandBlock)
    CoefsTotal.Block = CoefsTotal.Block + RanCoefs.Block;
end

iRandInteraction = find(strcmpi(RanCoefs.Properties.VariableNames,'Interaction'));
if ~isempty(iRandInteraction)
    CoefsTotal.Interaction = CoefsTotal.Interaction + RanCoefs.Interaction;
end

LMEOut.CoefsForPlot = CoefsTotal;


% Magic math: we only care about resudial variance and about how far x is from the mean.
% https://stats.stackexchange.com/questions/85560/shape-of-confidence-interval-for-predicted-values-in-linear-regression/85565#85565
% BUT NOOOOO, it must be done differently on the MLM: remember, there're so many DoFs here
% that the added variance would be tiny. Nonetheless, slope parameter estimates have
% finite (not too small) SDs.


if flag_output_text
    
    disp(sprintf('\n'));
    disp('*****Selected model*****');
    disp(sprintf('<strong>Linear mixed-effects model fit by %s\n</strong>',LME.FitMethod));
    disp(sprintf('<strong>Formula</strong>'));
    disp(LME.Formula);
    disp(LME.ModelCriterion)
    disp(sprintf('Rsq-adj    %.3f',LME.Rsquared.Adjusted));
    disp(FixEst); fprintf('\n');
    if flag_output_text == 1
        disp(sprintf('<strong>Random coefficients per participant (sorted)</strong>'));
        disp(RanCoefs);fprintf('\n');
    end
end

if L > 0
    if flag_output_text == 1
        disp(sprintf('\n'));
        disp(sprintf('<strong>Correlation of random coefficients with ranking</strong>'));
        disp(sprintf('\n'));
        fn = fieldnames(CoefsRankLM);
        for ifld = 1:length(fn)
            disp(CoefsRankLM.(fn{ifld}));
            disp(sprintf('\n'));
        end
    end
    LMEOut.CoefsRankingLM = CoefsRankLM;
end

if flag_output_text == 1
    disp(sprintf('<strong>Covariance of random coefficients</strong>'));
    disp(RanCov);
end





end




%% Figure with boxplot and summary of fixed effects
function ParamChangeSinglePlot(ParamsNames,ParamChangeModelOut,T,rgbmystyle,flag_OnlyOneStyle,flag_median)

% Version taken from ParameterOneWindow dashboard
flag_BlotBarNotBox = 0; % <--- Change for barplot instead of boxplot
% Prep data for boxplot.
subjlistMedDR = [14 10 13 16 6 5 11 8 4 1 3 7 2 15 9 12]; %ranking accuracy-based for both styles
subjlistMedR = [10 13 14 5 6 7 16 15 3 1 11 4 8 9 2 12]; % or only for R
subjlistMedD = [16 11 14 2 10 13 8 4 6 3 1 5 7 15 9 12]; % or only for D
styles2plot = 1:2;
if flag_OnlyOneStyle
    styles2plot = 1;
end
for iiS = 1:16
    iS = subjlistMedDR(iiS);
    for istyle = styles2plot
        T1 = T(T.StyleDR == (istyle-1) & T.Subj == iS,:);
        ParSSMea(iS,istyle) = mean(T1.Y,'omitnan');
        if flag_median
            ParSSMea(iS,istyle) = median(T1.Y,'omitnan');
        end
        for iblock = 1:5
            T1 = T(T.StyleDR == (istyle-1) & T.Subj == iS & T.Block == iblock,:);
            ParBlockMea(iS,istyle,iblock) = mean(T1.Y,'omitnan');
            
            if flag_median
                ParBlockMea(iS,istyle,iblock) = median(T1.Y,'omitnan');    
            end
        end
    end
end



FigBP = figure('Position',[300 200 361 209]);
FigBP.Name = 'Block mean(med) displayed. MLM-fit to all-trial values';
spBP = axes('Parent',FigBP,'position',[0.1 0.2 0.8 0.6]); hold on;
axes(spBP);
ylabel(ParamsNames{1},'fontweight','bold','FontSize',10,'Interpreter','none');
xlabel('Block','fontweight','bold');

%Prepare data for Boxplot. Have 16 points per block-style, Plot
YBPplot = permute(ParBlockMea,[1 3 2]);
if ~flag_BlotBarNotBox %Plot Boxes, apply formatting for further work
    PlotBoxPlotsStyleBlock(spBP,YBPplot,rgbmystyle,'ThinnerLines','AddAllBlocks',ParSSMea);
else % Plot Bars, apply formatting for further work
    PlotBoxPlotsStyleBlock(spBP,YBPplot,rgbmystyle,'ThinnerLines','BarNotBox','AddAllBlocks',ParSSMea);
end
xticks([1:5 7]);
xticklabels(cat(2,num2cell(1:5),{'All'}));
xlim([0 8]);


% Model Params - names changing because the script adopted from another fn.
% Extract slopes and CIs
ModelChange = ParamChangeModelOut;
LME = ModelChange.IntraLMM;
FixEst = LME.FixCoefs;
CoefsTotal = LME.CoefsForPlot;

iIntercept = find(strcmpi(FixEst.Name,'(Intercept)'));
iStyle = find(strcmpi(FixEst.Name,'StyleDR'));
iBlock = find(strcmpi(FixEst.Name,'Block'));
iInteraction = find(strcmpi(FixEst.Name,'StyleDR:Block'));

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
if ~isempty(iBlock)
    b(3) = FixEst.Estimate(iBlock);
    bCI(3,1) = FixEst.Lower(iBlock);
    bCI(3,2) = FixEst.Upper(iBlock);
end
if ~isempty(iInteraction)
    b(4) = FixEst.Estimate(iInteraction);
    bCI(4,1) = FixEst.Lower(iInteraction);
    bCI(4,2) = FixEst.Upper(iInteraction);
end


% Get x for plotting - generic approach despite exactly 5 blocks
xplot(1,[1 3]) = prctile(LME.LME.Variables.Block(LME.LME.Variables.StyleDR == 0),[2.5 97.5]);
xplot(2,[1 3]) = prctile(LME.LME.Variables.Block(LME.LME.Variables.StyleDR == 1),[2.5 97.5]);
xplot(:,2) = mean(xplot(:,[1 3]),2);
xplot(1,:) = xplot(1,:) - 0.17;
xplot(2,:) = xplot(2,:) + 0.17;
XMeanByStyle = xplot(:,2);

% Get lines
yplot = b(1) + b(2) .* [0; 1] + (b(3) + b(4) * [0; 1]) .* xplot;
% Coefficient CIs split onto the two styles
dBH = (bCI(:,2) - b) ./ 2;
dx = xplot - XMeanByStyle;
ypatch = [yplot flip(yplot,2)];
% Account for intercept and style effect
ypatch(:,1:3) = ypatch(:,1:3) - (dBH(1) + dBH(2));
ypatch(:,4:6) = ypatch(:,4:6) + (dBH(1) + dBH(2));
% Account for block and interaction effects
ypatch(:,[1:2 4]) = ypatch(:,[1:2 4]) + (dBH(3) + dBH(4)) .* dx(:,1:3);
ypatch(:,[3 5:6]) = ypatch(:,[3 5:6]) - (dBH(3) + dBH(4)) .* dx(:,[3 2 1]);

% Plot stuff
for istyle = 1:2
    patch('XData',[xplot(istyle,:), flip(xplot(istyle,:))],'YData',ypatch(istyle,:),'FaceColor',rgbmystyle(istyle,:),'FaceAlpha',0.35,'EdgeColor','none');
    plot(xplot(istyle,:),yplot(istyle,:),'Color',rgbmystyle(istyle,:));
end

% Display text summaries of the LMM
Intercept = LME.FixCoefs.Estimate(1);
InterceptP = LME.FixCoefs.pValue(1);
iStyleDR = find(strcmpi(LME.FixCoefs.Name,'StyleDR'));
if ~isempty(iStyleDR)
    StyleEst = LME.FixCoefs.Estimate(iStyleDR);
    StyleP = LME.FixCoefs.pValue(iStyleDR);
end
iBlock = find(strcmpi(LME.FixCoefs.Name,'Block'));
if ~isempty(iBlock)
    BlockEst = LME.FixCoefs.Estimate(iBlock);
    BlockP = LME.FixCoefs.pValue(iBlock);
end
iInteraction = find(strcmpi(LME.FixCoefs.Name,'StyleDR:Block'));
if ~isempty(iInteraction)
    Interaction = LME.FixCoefs.Estimate(iInteraction);
    InteractionP = LME.FixCoefs.pValue(iInteraction);
end

TextFont = 'Normal';
if ~isempty(iStyleDR) && StyleP < 0.05
    TextFont = 'bold';
end
if ~isempty(iStyleDR)
    StyleStr = sprintf('Style: b = %.4f, p=%.4f',StyleEst,StyleP);
else
    StyleStr = sprintf('No Style in the model');
end
TxStyle = text(0.105,1.3,StyleStr,'HorizontalAlignment','Left','FontSize',9,'Units','normalized','Color','k','fontweight',TextFont);


TextFont = 'Normal';
if BlockP < 0.05
    TextFont = 'bold';
end
MainStr = sprintf('Block: b = %.4f, p=%.4f',BlockEst,BlockP);
TxBlock = text(0.095,1.18,MainStr,'HorizontalAlignment','Left','FontSize',9,'Units','normalized','Color','k','fontweight',TextFont);

TextFont = 'Normal';
if ~isempty(iInteraction) && InteractionP < 0.05
    TextFont = 'bold';
end
if ~isempty(iInteraction)
    InteractionStr = sprintf('Interaction: b = %.4f, p=%.4f',Interaction,InteractionP);
else
    InteractionStr = sprintf('No interaction, p(chi2)=%.4f',ParamChangeModelOut.IntraLM_Interaction_Chi2_dof_P(3));
end
TxInteraction = text(0,1.065,InteractionStr,'HorizontalAlignment','Left','FontSize',9,'Units','normalized','Color','k','fontweight',TextFont);

InterceptStr = sprintf('b_0=%.4f (p=%.4f)',Intercept,InterceptP);
TxInter = text(0.73,1.24,InterceptStr,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');


R2Str = sprintf('R^2=%.2f',LME.LME.Rsquared.Adjusted);
TxR2 = text(0.73,1.08,R2Str,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');



%% Figure with individual slopes Style x Block from LMM
if 0
    FigBPInd = figure('Position',[700 200 361 209]);
    FigBPInd.Name = 'Individual slopes from the LMM';
    spMLMSlopesD = axes(FigBPInd,'position',[0.17 0.2 0.3 0.7]); hold on;
    tt = title('Discrete');
    tt.Units = 'normalized';
    tt.Position(2) = 1;
    ylabel(ParamsNames{1},'fontweight','bold','FontSize',12,'Interpreter','none');
    xlabel('Block','fontweight','bold');
    
    spMLMSlopesR = axes(FigBPInd,'position',[0.54 0.2 0.3 0.7]); hold on;
    tt = title('Rhythmic');
    tt.Units = 'normalized';
    tt.Position(2) = 1;
    yticklabels({}); yticks('auto');
    
    
    xmin = 1;
    xmax = 5;
    xplot = [xmin xmax];
    jj = flipud(jet);
    ColorsSubj = jj(round(linspace(1,256,16)),:);
    
    if ~isempty(ModelChange)
        % Create lines.
        istyle = 0; % To keep it consistent with R
        for iisubj = 1:16
            yplot = CoefsTotal.Intercept(iisubj) + CoefsTotal.StyleDR(iisubj) * istyle + CoefsTotal.Block(iisubj) * xplot + CoefsTotal.Interaction(iisubj) * istyle * xplot;
            plot(spMLMSlopesD,xplot,yplot,'Color',ColorsSubj(iisubj,:));
        end
        
        istyle = 1;
        for iisubj = 1:16
            yplot = CoefsTotal.Intercept(iisubj) + CoefsTotal.StyleDR(iisubj) * istyle + CoefsTotal.Block(iisubj) * xplot + CoefsTotal.Interaction(iisubj) * istyle * xplot;
            plot(spMLMSlopesR,xplot,yplot,'Color',ColorsSubj(iisubj,:));
        end
        
        xlim(spMLMSlopesD,[0 6]);
        xlim(spMLMSlopesR,[0 6]);
        linkaxes([spMLMSlopesR spMLMSlopesD],'x');
        xticks(spMLMSlopesR,[1 5]);
        xticks(spMLMSlopesD,[1 5]);
        
        YLim0 = spBP.YLim;
        linkaxes([spBP spMLMSlopesR spMLMSlopesD],'y');
        spBP.YLim = YLim0;
        
        colormap(jet);
        cb = colorbar('Position',[0.87 0.2 0.02 0.7]);
        cb.Ticks = [0 1];
        cb.TickLabels = {'Worst','Best'};
        cLab = ylabel(cb,'Accuracy Ranked','Units','normalized','Position',[1 0.5 0],'Fontweight','bold');
    end
end

end



