%% Krotov, Russo, Nah, Hogan & Sternad (2022) Motor Control Beyond Reach: Hitting a Target with a Bullwhip.
% Royal Society Open Science

% The findings reported in the manuscript are supported with the data stored in the
% WhipTaskData variable (.csv or .mat). The following scripts perform the statistical
% tests and produce the figures that are identical (up to design) to those
% submitted in the manuscript.

% MATLAB R2021a was used with Statistic  and Curve Fitting Toolboxes
% See auxilliary dependent functions from the same root folder.


% Set nicer background and plot labels
WhipTask_setMyDefaults();

%% Generate qualitative figures
% Plot Error as barplot, using 5 block-median (or mean) values for each bar, sorting by DR
% median error
WhipTask_BarPlotParticipantWise(WhipTask_TrialStats,rgbmystyle2,'MinDist','NonNormal');
WhipTask_BarPlotParticipantWise(WhipTask_TrialStats,rgbmystyle2,'HitFlag'); % Median cannot handle binary values.

% Plot Whip Speed (and hand) time-profiles, taken on ThrowOnset-MinDist interval, averaged
% within subject-style.
WhipTask_WSAveragePlot(WhipTask_TimeSeries,rgbmystyle2);

% Plot Hand Speed time-profiles, averaged within subject-style, aligned by MinDist
WhipTask_HSAveragePlot(WhipTask_TimeSeries,rgbmystyle2);

% Plot whip states at PeakHandSpeed, top-view, colored wrt Error (nonlinear jet color scale)
WhipTask_PlotWhipStates(WhipTask_TimeSeries);

%% Generate statistical reports in command line and figures illustrating these effects (also with brief stats)
% Run tests Parameter ~ Style x Block, plot results and display stat summary
% Linear Mixed Model fia fitlme (fitglme for HitFlag), 
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'MinDist','OnlyMainText','NonNormal');
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'HitFlag','OnlyMainText'); % disregard LMM displayed predictions due to GLMM linking fn. 
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'t_InterTrial','OnlyMainText');
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'PeakHandSpeed','OnlyMainText');
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'PeakTipSpeed','OnlyMainText');
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'WhipExt_PeakHand','OnlyMainText','NonNormal'); % measured in a.u. on [0,1]. To match X-axis and slopes (betas), scale by whip length, L= 1.6m
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'WhipAzAng_PeakHand','OnlyMainText');
WhipTask_ParamChangeOverTrials(WhipTask_TrialStats,rgbmystyle2,'HandAzAng_PeakHand','OnlyMainText');

% Plot CoV of intertrial interval, run LMM test Parameter ~ Style x Block on these
% block-averaged values
WhipTask_BoxPlotCoV(WhipTask_TrialStats,rgbmystyle2,'t_InterTrial');

% Plot IQR of Minimal Distance, run LMM test Parameter ~ Style x Block on these
WhipTask_BoxPlotCoV(WhipTask_TrialStats,rgbmystyle2,'MinDist','NonNormal','IQR');

% Run correlation tests Error ~ Style x Parameter, plot results and display stat summary
% Linear model on participant-mean values via fitlm, interaction may be added iteratively, following the ML-test.
% Linear Mixed Model fia fitlme, Nparameters increased iteratively, evaluating by ML-test
WhipTask_CorrParameters_Inter_Intra(WhipTask_TrialStats,rgbmystyle2,'PeakHandSpeed','MinDist','OnlyMainText');
WhipTask_CorrParameters_Inter_Intra(WhipTask_TrialStats,rgbmystyle2,'PeakTipSpeed','MinDist','OnlyMainText');
WhipTask_CorrParameters_Inter_Intra(WhipTask_TrialStats,rgbmystyle2,'WhipExt_PeakHand','MinDist','OnlyMainText');
WhipTask_CorrParameters_Inter_Intra(WhipTask_TrialStats,rgbmystyle2,'WhipAzAng_PeakHand','MinDist','OnlyMainText');
WhipTask_CorrParameters_Inter_Intra(WhipTask_TrialStats,rgbmystyle2,'HandAzAng_PeakHand','MinDist','OnlyMainText');

% Run correlation test Parameter 2 ~ Style x Parameter 1, plot results and display stat summary
% Linear Mixed Model fia fitlme, Nparameters increased iteratively, evaluating by ML-test
WhipTask_CorrParameters_Inter_Intra(WhipTask_TrialStats,rgbmystyle2,'PeakHandSpeed','PeakTipSpeed','OnlyMainText');





