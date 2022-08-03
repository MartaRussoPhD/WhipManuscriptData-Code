% HEADER for data processing and plotting for
% Krotov, Russo, Nah, Hogan & Sternad (2022) Motor Control Beyond Reach: Hitting a Target with a Bullwhip.


%% DISCUSS

% Statistical findings reported in the manuscript are supported with the data stored in the
% WhipTaskData variable (.csv or .mat). The following scripts perform the statistical
% tests and produce the figures that are identical (up to prettifying design) to those
% submitted in the manuscript.

% !!! Complete processing (filtering, structuring, computing tangential velocity, computing error etc.)
% and complete dataset (all marker positions) are not included with this manuscript.

% !!! OR SHOULD WE INCLUDE SOME DATA? OTHERWISE PLOTS WITH HAND SPEEDS, WHIP MARKER SPEEDS
% CANNOT BE REPLICATED.




% We have whip state graphical plots at HSmax. To plot those, I need to include whip
% marker data at Peak Hand Speed landmark. Include?

% For that I need the hand speed data + a script how to compute Peak Hand Speed. Include?

% Rev 1 asked about definition of Error in the figure with the whip states. It's the same
% as Error in stats plots. Define? Should I include data (whip w1-w3 + targ1-2) + script
% to allow computing MinDist and tMinDist?

% We have whip marker speeds from all trials. Should I include original WhipSpeed data +
% script?


%% Contents: statistics figures
% Set nicer background and plot labels
setMyDefaults();

% Run tests Parameter ~ Style x Block, plot results and display stat summary
% Linear Mixed Model fia fitlme (fitglme for HitFlag), 
%   Nparameters increased iteratively, evaluating by ML-test
ParamChangeOverTrials(WhipTaskData,rgbmystyle2,'MinDist','OnlyMainText');
ParamChangeOverTrials(WhipTaskData,rgbmystyle2,'HitFlag','OnlyMainText'); % disregard LMM displayed predictions due to GLMM linking fn. 
ParamChangeOverTrials(WhipTaskData,rgbmystyle2,'PeakHandSpeed','OnlyMainText');
ParamChangeOverTrials(WhipTaskData,rgbmystyle2,'PeakTipSpeed','OnlyMainText');
ParamChangeOverTrials(WhipTaskData,rgbmystyle2,'WhipExt_PeakHand','OnlyMainText');
ParamChangeOverTrials(WhipTaskData,rgbmystyle2,'WhipAzAng_PeakHand','OnlyMainText');

% Run correlation tests Error ~ Style x Parameter, plot results and display stat summary
% Linear Mixed Model fia fitlme, Nparameters increased iteratively, evaluating by ML-test
CorrParameters_Inter_Intra(WhipTaskData,rgbmystyle2,'PeakHandSpeed','MinDist','OnlyMainText');
CorrParameters_Inter_Intra(WhipTaskData,rgbmystyle2,'PeakTipSpeed','MinDist','OnlyMainText');
CorrParameters_Inter_Intra(WhipTaskData,rgbmystyle2,'WhipExt_PeakHand','MinDist','OnlyMainText');
CorrParameters_Inter_Intra(WhipTaskData,rgbmystyle2,'WhipAzAng_PeakHand','MinDist','OnlyMainText');

% Run correlation test Parameter 2 ~ Style x Parameter 1, plot results and display stat summary
% Linear Mixed Model fia fitlme, Nparameters increased iteratively, evaluating by ML-test
CorrParameters_Inter_Intra(WhipTaskData,rgbmystyle2,'PeakHandSpeed','PeakTipSpeed','OnlyMainText');











%% Definition of Error and tMinDist
% Prepare an exemplary code that takes subjparamsDerived, interpolates w1-w3, t1-t3, 
% computes pairwise distances, and finds mindist and tmindist.
%% Definition of HSmax
% Using hand speed (include all HS data or 1 participant 1 block?), find HSmax (with 10 Hz
% filtering)

%% Definition of ThrowOnset
% Using HdL -- Targ1 distance, find its maximum.

%% Definition of Rhythmicity
% Using 

%% Definition of Whip Orientation and Extension
%% Definition of Hand/Handle Orientation

%% Boxplot Y ~ Block*Style
%% Correlations Error ~ Y*Style
%% LMM Y ~ Block*Style
%% LMM Error ~ Y*Style
%% Plot WhipStates
