function WhipTask_BoxPlotCoV(ParTab,rgbmystyle,ParamName,varargin)

% Computes coefficient of variation (SD/mean or IQR/median for non-parametric approximation) across ~30 trials 
% of each participant, each style, each block. 
% Plots a boxplot of the values, pooling 16 participant values.
% Runs a linear mixed model (fitted iteratively) and displays fixed effect estimates and p-values in the figure.

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

flag_type = 0; % CoV
if any(strcmpi(varargin,'Variability')) || any(strcmpi(varargin,'SD')) || any(strcmpi(varargin,'STD')) || any(strcmpi(varargin,'IQR'))
    flag_type = 1;
end
if any(strcmpi(varargin,'Mean')) || any(strcmpi(varargin,'Median'))
    flag_type = 2;
end





if flag_median
    switch flag_type
        case 2
            YName = 'Block Median';
        case 1
            YName = 'Block IQR';
        case 0
            YName = 'Block nonparam CoV';
    end
else
    switch flag_type
        case 2
            YName = 'Block Mean';
        case 1
            YName = 'Block STD';
        case 0
            YName = 'Block CoV';
    end
end

flag_BlotBarNotBox = 0; % <--- Change for barplot instead of boxplot



ic(1) = find(strcmpi(ParTab.Properties.VariableNames,'Subj'));
ic(2) = find(strcmpi(ParTab.Properties.VariableNames,'Subj_Ranked'));
ic(3) = find(strcmpi(ParTab.Properties.VariableNames,'StyleDR'));
ic(4) = find(strcmpi(ParTab.Properties.VariableNames,'Block'));

T = ParTab(:,ic); % Subj, Subj_Ranked,StyleDR,Block

T.Y = ParTab.(ParamName);
T.Error = ParTab.MinDist;
T.MinDist = ParTab.MinDist;
T.StyleDR = T.StyleDR - 1; % To make consistent with R.

% Prep data.
subjlistMedDR = [14 10 13 16 6 5 11 8 4 1 3 7 2 15 9 12]; %ranking accuracy-based for both styles
for iiS = 1:16
    iS = subjlistMedDR(iiS);
    for istyle = 1:2
        T1 = T(T.StyleDR == (istyle-1) & T.Subj == iS,:);
        ParSSMea(iS,istyle) = mean(T1.Y,'omitnan');
        ParSSStd(iS,istyle) = std(T1.Y,[],'omitnan');
        if flag_median
            ParSSMea(iS,istyle) = median(T1.Y,'omitnan');
            ParSSStd(iS,istyle) = iqr(T1.Y);
        end
        for iblock = 1:5
            T1 = T(T.StyleDR == (istyle-1) & T.Subj == iS & T.Block == iblock,:);
            ParBlockMea(iS,istyle,iblock) = mean(T1.Y,'omitnan');
            ParBlockS(iS,istyle,iblock) = std(T1.Y,[],'omitnan');
            if flag_median
                ParBlockMea(iS,istyle,iblock) = median(T1.Y,'omitnan');
                ParBlockS(iS,istyle,iblock) = iqr(T1.Y);
            end
        end
    end
end
ParBlockCOV = ParBlockS ./ ParBlockMea;
ParSSMeaCOV = ParSSStd ./ ParSSMea;

% Also extract as a table for LMM
% across a block
if flag_median
    TBlockMea = varfun(@(x)median(x,'omitnan'),T,'InputVariables','Y',...
        'GroupingVariables',{'Subj','Subj_Ranked','StyleDR','Block'});
    TBlockStd = varfun(@(x)iqr(x),T,'InputVariables','Y',...
        'GroupingVariables',{'Subj','Subj_Ranked','StyleDR','Block'});
else
    TBlockMea = varfun(@(x)mean(x,'omitnan'),T,'InputVariables','Y',...
        'GroupingVariables',{'Subj','Subj_Ranked','StyleDR','Block'});
    TBlockStd = varfun(@(x)std(x,[],'omitnan'),T,'InputVariables','Y',...
        'GroupingVariables',{'Subj','Subj_Ranked','StyleDR','Block'});
end
TBlockMea.mean_Y = TBlockMea.Fun_Y;
TBlockMea.std_Y = TBlockStd.Fun_Y;
TBlockMinDistMed = varfun(@(x)median(x,'omitnan'),T,'InputVariables','Error',...
    'GroupingVariables',{'Subj','Subj_Ranked','StyleDR','Block'});
TBlockMinDistIqr = varfun(@(x)iqr(x),T,'InputVariables','Error',...
    'GroupingVariables',{'Subj','Subj_Ranked','StyleDR','Block'});
TBlockMea.MinDist = TBlockMinDistMed.Fun_Error;



switch flag_type
    case 2
        TBlockMea.std_Y = TBlockMea.mean_Y;
    case 1
    case 0
    TBlockMea.std_Y = TBlockMea.std_Y ./ TBlockMea.mean_Y;
end









FigBP = figure('Position',[300 200 361 209]);
FigBP.Name = 'Block mean(med) displayed. MLM-fit to all-trial values';
spBP = axes('Parent',FigBP,'position',[0.15 0.2 0.8 0.6]); hold on;
axes(spBP);
ylabel(sprintf('%s of %s',YName,ParamName),'fontweight','bold','FontSize',8,'Interpreter','none');
xlabel('Block','fontweight','bold');

switch flag_type
    case 2
    case 1
        ParBlockMea = ParBlockS;
        ParSSMea = ParSSStd;
    case 0
        ParBlockMea = ParBlockCOV;
        ParSSMea = ParSSMeaCOV;
end


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


%%% LMM Test
TBlockMea.StyleDR = TBlockMea.StyleDR + 1;
ModelChange = WhipTask_ParamChangeOverTrials(TBlockMea,[],'std_Y','OnlyMainText','SupressGraphics');
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
            

            % Get x
            xplot(1,[1 3]) = prctile(LME.LME.Variables.Block(LME.LME.Variables.StyleDR == 0),[2.5 97.5]);
            xplot(2,[1 3]) = prctile(LME.LME.Variables.Block(LME.LME.Variables.StyleDR == 1),[2.5 97.5]);
            xplot(:,2) = mean(xplot(:,[1 3]),2);
            xplot(1,:) = xplot(1,:) - 0.17;
            xplot(2,:) = xplot(2,:) + 0.17;
            
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
            
            % Plot stuff
            for istyle = 1:2
                patch('XData',[xplot(istyle,:), flip(xplot(istyle,:))],'YData',ypatch(istyle,:),'FaceColor',rgbmystyle(istyle,:),'FaceAlpha',0.35,'EdgeColor','none');
                plot(xplot(istyle,:),yplot(istyle,:),'Color',rgbmystyle(istyle,:));
            end
            
            
            
            % Display text
            istyle = 1;
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
                InteractionStr = sprintf('No interaction, p(chi2)=%.4f',ModelChange.IntraLM_Interaction_Chi2_dof_P(3));
            end
            TxInteraction = text(0,1.065,InteractionStr,'HorizontalAlignment','Left','FontSize',9,'Units','normalized','Color','k','fontweight',TextFont);
            
            InterceptStr = sprintf('b_0=%.4f (p=%.4f)',Intercept,InterceptP);
            TxInter = text(0.73,1.24,InterceptStr,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
            
            
            R2Str = sprintf('R^2=%.2f',LME.LME.Rsquared.Adjusted);
            TxR2 = text(0.73,1.08,R2Str,'HorizontalAlignment','Left','FontSize',8,'Units','normalized','Color','k');
            


end
