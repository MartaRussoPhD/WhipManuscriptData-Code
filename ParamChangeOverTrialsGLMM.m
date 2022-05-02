function CorrModelOut = ParamChangeOverTrialsGLMM(ParTab,rgbmystyle,ParamName,varargin)
% 1. Runs a mixed-effect linear model on all (non-discarded) trials. Or GLMM 
% (Laplacian approx of ML, binomial distribution, logit link). The models are fitted
% iteratively. a) testing whether random effects (intercepts) are needed. b) testing
% whether fixed interaction is needed. c) if b) confirms, continuing with more random
% effect terms, up to random interaction. Iteratively compares via ML-criterion.
% The model handles binary parameters well.
% ______ Returns a structure with an LM object, an LMER object, LMER coefficient test table (t-test, using
% satterthwaite method for dofs, 95% CI included; default DOFs estimation for GLMM), random coefficients for participants
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

T = ParTab(:,[2 3 8 9]);
T.Block = ParTab.Block;
T.Y = ParTab.(ParamName);

if flag_output_text
    disp(sprintf('Testing change of %s, between Styles and across Blocks',ParamName));
end

T.StyleDR = T.StyleDR - 1; % To make consistent with R.

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

%%%%%%% If computing effect sizes without random effects,
% Use eqns from Field & Myers "Discovering Statistics with R", p.542

%%%% Test whether random effects need to be included at all.
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
        CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
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
        if flag_bin
            LMEAnovaFix = anova(LME);
        else
            LMEAnovaFix = anova(LME,'DFMethod','satterthwaite');
        end
        CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
        return
    end
    
    LMEComparison = LMEComparison(nmin,:);
    LMEComparison{1,5:7} = nan;
    
end



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
    nmin = 1;
    if LMEComparison.LRTest_p(ilme,1) < 0.05
        LMEComparison.Selected(ilme,1)  = 1;
        nmin = 2;
    end
    if flag_output_text
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
    CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
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
    if flag_output_text
        disp(LMEComparison);
    end
    
    %%%% DEBUG ONLY
    if 1 == 0
        disp('Debug!!! Force no rnd except intercepts');
        CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
        return
    end
    %%%% DEBUG ONLY
    if 1 == 0
        disp('Debug!!! Force only random Block or Random StyleDR'); % nmin = 2 or 3.
        nmin = 3;
        LMEComparison(2,:) = LMEComparison(nmin,:);
        LMEComparison(3,:) = [];
        LME = LMMs(nmin).LME;
        CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
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
        CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
        return
    end
    
    
    % Compute the remaining model. Prepopulate the table
    ilm = 3;
        LMEComparison.LMEFormulas(ilm) = LMEFormulasNew(ilm);
        LMEComparison.npar(ilm) = LMEparNNew(ilm);
        LMEComparison{ilm,3:7} = nan;
    
    % Compute the remaining model.
	 ilme = 3;
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
	 
    if flag_output_text
        disp(LMEComparison);
    end
    
    
    
    % Make the final decision
    if LMEComparison.LRTest_p(3,1) >= 0.05 % the first one is not good
        if flag_output_text
            disp(sprintf('<strong>Combination of random slopes will not be included</strong>'));
        end
        LME = LMMs(2).LME;
    else
        if flag_output_text
            disp(sprintf('<strong>Sum, but not interaction, of random slopes will be included</strong>'));
        end
        LME = LMMs(3).LME;
    end
       
    CorrModelOut.IntraLMM = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text);
    
    
    
end
warning('on','MATLAB:table:RowsAddedExistingVars');


end




function LMEOut = DisplayAndPackLMER(LME,RankingType,ParamsNames,flag_OnlyOneStyle,flag_output_graph,flag_output_text)
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
RanCoefs = table(Ranking,Subj);
CoefsRankLM = [];
switch L
    case 0
        
    case 1 % Intercept only
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Intercept ~ Ranking');
    case 2 % Intercept and one slope
        SlopeName = Bnames.Name{2}(2:end-1);
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*2-1);
            RanCoefs.(SlopeName)(isubj) = Bnames.Value(isubj*2);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Intercept ~ Ranking');
        CoefsRankLM.(SlopeName) = fitlm(RanCoefs,sprintf('%s ~ Ranking',SlopeName));
    case 3 %Intercept and two slopes
        SlopeName1 = Bnames.Name{2};
        SlopeName2 = Bnames.Name{3};
        for isubj = 1:16
            RanCoefs.Intercept(isubj) = Bnames.Value(isubj*3-2);
            RanCoefs.(SlopeName1)(isubj) = Bnames.Value(isubj*3-1);
            RanCoefs.(SlopeName2)(isubj) = Bnames.Value(isubj*3);
        end
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Intercept ~ Ranking');
        CoefsRankLM.(SlopeName1) = fitlm(RanCoefs,sprintf('%s ~ Ranking',SlopeName1));
        CoefsRankLM.(SlopeName2) = fitlm(RanCoefs,sprintf('%s ~ Ranking',SlopeName2));
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
        CoefsRankLM.Intercept = fitlm(RanCoefs,'Intercept ~ Ranking');
        CoefsRankLM.(SlopeName1) = fitlm(RanCoefs,sprintf('%s ~ Ranking',SlopeName1));
        CoefsRankLM.(SlopeName2) = fitlm(RanCoefs,sprintf('%s ~ Ranking',SlopeName2));
        CoefsRankLM.Interaction = fitlm(RanCoefs,'Interaction ~ Ranking');    
end

%%%%% Prepare a better alternative for RanCoefs table. In case some participants are
%%%%% completely nan. Use a bunch of strcmpi, find participant number, each factor, be
%%%%% happy


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

xplotline = [1 5];

% Test figure
if flag_output_graph
    figure
    T = ParamsNames{3};
    xmin = 1;
    xmax = 5;
    
    xplot(1,[1 3]) = prctile(LME.Variables.Block(LME.Variables.StyleDR == 0),[2.5 97.5]);
    xplot(2,[1 3]) = prctile(LME.Variables.Block(LME.Variables.StyleDR == 1),[2.5 97.5]);
    xplot(:,2) = mean(xplot(:,[1 3]),2);
    
    
    jj = flipud(jet);
    ColorsSubj = jj(round(linspace(1,256,16)),:);
    
    sp1 = subplot('Position',[0.15 0.15 0.3 0.8]); hold on
    title('Discrete'); ylabel(ParamsNames{1},'Interpreter','none','Fontweight','bold'); xlabel(ParamsNames{1},'Interpreter','none','Fontweight','bold');
    
    
    % Create lines.
    
    if ~flag_OnlyOneStyle
    istyle = 0; % To keep it consistent with R
    for iisubj = 1:16
        yplot = CoefsTotal.Intercept(iisubj) + CoefsTotal.StyleDR(iisubj) * istyle + CoefsTotal.Block(iisubj) * xplot(istyle+1,:) + CoefsTotal.Interaction(iisubj) * istyle * xplot(istyle+1,:);
        plot(sp1,xplot(istyle+1,:),yplot,'Color',ColorsSubj(iisubj,:));
    end
    end
    

    sp2 = subplot('Position',[0.55 0.15 0.3 0.8]); hold on
    title('Rhythmic'); ylabel(ParamsNames{1},'Interpreter','none','Fontweight','bold'); xlabel(ParamsNames{1},'Interpreter','none','Fontweight','bold');

    
    % Create lines.
    istyle = 1; % To keep it consistent with R
    for iisubj = 1:16
        yplot = CoefsTotal.Intercept(iisubj) + CoefsTotal.StyleDR(iisubj) * istyle + CoefsTotal.Block(iisubj) * xplot(istyle+1,:) + CoefsTotal.Interaction(iisubj) * istyle * xplot(istyle+1,:);
        plot(sp2,xplot(istyle+1,:),yplot,'Color',ColorsSubj(iisubj,:));
    end
    
    % Adjust YLim properly!
    Y1 = RemoveOutliersMy([T.Y],[3 3]);
    ymin = min(Y1);
    ymax = max(Y1);
    ylim(sp1,[0 ymax + 0.1 * (ymax-ymin)]);
    ylim(sp2,[0 ymax + 0.1 * (ymax-ymin)]);
    xlim([1 5]);
    
    colormap(jet);
    cb = colorbar('Position',[0.87 0.15 0.02 0.8]);
    cb.Ticks = [0 1];
    cb.TickLabels = {'Worst','Best'};
    cLab = ylabel(cb,'Accuracy Ranked','Units','normalized','Position',[1 0.5 0],'Fontweight','bold');
    
    
    % Plot single line per style with CI?
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
    
    %XMeanByStyle = [mean(LME.Variables.Block(LME.Variables.StyleDR == 0),'omitnan'); mean(LME.Variables.Block(LME.Variables.StyleDR == 1),'omitnan')];
    XMeanByStyle = xplot(:,2);
    % Get lines
    yplot = b(1) + b(2) .* [0; 1] + (b(3) + b(4) * [0; 1]) .* xplot;
    
    % Coefficient uncertainties split onto the two styles
    dBH = (bCI(:,2) - b) ./ 2;
    dx = xplot - XMeanByStyle;
    
    ypatch = [yplot flip(yplot,2)];
    % Account for intercept and style effect
    ypatch(:,1:3) = ypatch(:,1:3) - (dBH(1) + dBH(2));
    ypatch(:,4:6) = ypatch(:,4:6) + (dBH(1) + dBH(2));
    
    % Account for block and interaction effects
    ypatch(:,[1:2 4]) = ypatch(:,[1:2 4]) + (dBH(3) + dBH(4)) .* dx(:,1:3);
    ypatch(:,[3 5:6]) = ypatch(:,[3 5:6]) - (dBH(3) + dBH(4)) .* dx(:,[3 2 1]);
    
    
    figure;
    subplot(1,1,1); hold on
    
    % Plot stuff
    istyle = 1;
    patch('XData',[xplot(istyle,:), flip(xplot(istyle,:))],'YData',ypatch(istyle,:),'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none');
    plot(xplot(istyle,:),yplot(istyle,:),'Color','r');
    istyle = 2;
    patch('XData',[xplot(istyle,:), flip(xplot(istyle,:))],'YData',ypatch(istyle,:),'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none');
    %yplot = b(1) + b(2) * (istyle-1) + (b(3) + b(4) * (istyle-1)) * xplot;
    plot(xplot(istyle,:),yplot(istyle,:),'Color','g');
    
end




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
        fn = fieldnames(CoefsRankLM);
        for ifld = 1:length(fn)
            disp(CoefsRankLM.(fn{ifld}));
            disp(sprintf('\n'));
        end
    end
    LMEOut.CoefsRankingLM = CoefsRankLM;
end

if flag_output_text
    disp(sprintf('<strong>Covariance of random coefficients</strong>'));
    disp(RanCov);
end

end




