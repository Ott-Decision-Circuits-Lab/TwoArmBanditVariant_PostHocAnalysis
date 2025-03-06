function AnalysisFigure = Matching_MS_B_ChoiceSymmetricQLearning_Visualisation(DataFolderPath, Prior)
% MS = MultiSession
% B = Bayesian <- Prior using simulation & MCMC (Hamiltonian MC) sampling
% from prior to get marginal posterior
% Matching Analysis Function
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ Jan 2025
% Model iteration see the end of script

%% load files
if nargin < 1
    DataFolderPath = uigetdir(OttLabDataServerFolderPath());
elseif ~ischar(DataFolderPath) && ~isstring(DataFolderPath)
    disp('Error: Unknown input format. No further analysis can be performed.')
    return
end

try
    load(fullfile(DataFolderPath, '\Selected_Data.mat'));
    load(fullfile(DataFolderPath, '\Concatenated_Data.mat'));
catch
    disp('Error: Selected DataFolderPath does not contain the required .mat for further steps.')
    return
end

SessionDateRange = DataFolderPath(end-16:end);
[~, RatName] = fileparts(fileparts(fileparts(DataFolderPath)));

RatID = str2double(RatName);
if isnan(RatID)
    RatID = -1;
end
RatName = num2str(RatID);

AnalysisName = 'Matching_MS_B_ChoiceSymmetricQLearning';

%% Hierarchaical Symmetric Q-Learning with Forgetting and Stickiness model
% Parametric Prior
% 20241220 tested with simulation that works well as initial parameters
if nargin < 2
    Prior.LearningRateAlpha = 4; % beta distribution beta(1, 1) = uniform distribution
    Prior.LearningRateBeta = 12; % beta(4, 12) gives mean ~ 0.25, SD = 0.11 <- from simulation
    
    Prior.InverseTemperatureMean = 8; % normal distribution (can be neg. but fully)
    Prior.InverseTemperatureSigma = 1; % N(8, 1) <- from simulation
    
    Prior.ForgettingRateAlpha = 2; % beta distribution beta(1, 1) = uniform distribution
    Prior.ForgettingRateBeta = 10; % beta(2, 10) gives mean ~ 0.167, SD ~ 0.10
    
    Prior.ChoiceStickinessMean = -1; % normal distribution (can be neg. but fully)
    Prior.ChoiceStickinessSigma = 0.4; % N(-1, 0.4) <- from simulation
    
    Prior.ChoiceForgettingRateAlpha = 10; % beta distribution beta(1, 1) = uniform distribution
    Prior.ChoiceForgettingRateBeta = 6; % beta(10, 6) gives mean ~0.63, SD ~0.12
    
    Prior.BiasMean = 0; % normal distribution
    Prior.BiasSigma = 1;

    Prior.BurnIn = 200;
    Prior.nSample = 1000;
    Prior.nChain = 12;

elseif fieldnames(Prior)
    disp('Error: Unknown input format. No further analysis can be performed.')
    return
end

Models = {};
try
    load(fullfile(DataFolderPath, '\Matching_MS_B_ChoiceSymmetriQLearning.mat'));
catch
    for iSession = 1:length(DataHolder)
        SessionData = DataHolder{iSession};
        
        try
            Models{iSession} = Matching_SS_B_ChoiceSymmetricQLearning_Model(SessionData, Prior);
        catch
            Models{iSession} = [];
            disp(strcat('Warning: Problem running model for Session', num2str(iSession)))
        end
    end

    save(fullfile(DataFolderPath, '\Matching_MS_B_ChoiceSymmetriQLearning.mat'), 'Models')
    disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')
end

if isempty(Models)
    disp('Error: No model is either computed or loaded. No further visualisation is possible.')
    return
end

%% create figure
% create figure
AnalysisFigure = figure('Position', [   0       0    1191     842],... % DIN A3, 72 ppi
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateRange, '_', AnalysisName),...
                        'MenuBar', 'none',...
                        'Resize', 'off');

% spacer for correct saving dimension
FrameAxes = axes(AnalysisFigure, 'Position', [0 0 1 1]);
set(FrameAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'w',...
    'YColor', 'w')

% Figure Info
FigureInfoAxes = axes(AnalysisFigure, 'Position', [0.01    0.98    0.48    0.01]);
set(FigureInfoAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'w',...
    'YColor', 'w')

FigureTitle = strcat(RatName, '_', SessionDateRange, '_', AnalysisName);

FigureTitleText = text(FigureInfoAxes, 0, 0,...
                       FigureTitle,...
                       'FontSize', 14,...
                       'FontWeight','bold',...
                       'Interpreter', 'none');

% colour palette
ColourPalette = CommonColourPalette();

%% Analysis across sessions
SessionDateLabel = [];

% Psychometric
PsychometricAxes = axes(AnalysisFigure, 'Position', [0.01    0.75    0.15    0.19]);
hold(PsychometricAxes, 'on')

AllPredictedProb = [];
AllLogOdds = [];
AllChoiceLeft = [];

set(PsychometricAxes,...
    'FontSize', 10,...
    'XLim', [-5 5],...
    'YLim', [0, 100],...
    'YAxisLocation', 'right')
title(PsychometricAxes, 'Psychometric')
xlabel(PsychometricAxes, 'log(odds)')
ylabel(PsychometricAxes, 'Left Choices (%)')

% Posterior mode (i.e. MAP, maximum a posteriori) estimate of ChoiceSymmetricQ
MAPAxes = axes(AnalysisFigure, 'Position', [0.22    0.75    0.15    0.19]);
hold(MAPAxes, 'on');

LearningRateMAPs = [];
InverseTemperatureMAPs = [];
ForgettingRateMAPS = [];
ChoiceStickinessMAPs = [];
ChoiceForgettingRateMAPs = [];
BiasMAPs = [];

set(MAPAxes,...
    'FontSize', 10,...
    'XLim', [0 7],...
    'XTick', 1:6,...
    'XTickLabel', {'\alpha', '\beta', '\gamma', '\phi', 'c_\gamma', 'bias'},...
    'YLim', [-2, 3],...
    'YAxisLocation', 'right')
xlabel(MAPAxes, 'Parameters');
ylabel(MAPAxes, 'Maximum a posteriori (a.u.)');
title(MAPAxes, 'Mode of posterior distribution')

% Vevaiometric      
VevaiometricAxes = axes(AnalysisFigure, 'Position', [0.01    0.48    0.15    0.19]);
hold(VevaiometricAxes, 'on')

AllExploringTI = [];
AllExploitingTI = [];

AllExploringLogOdds = [];
AllExploitingLogOdds = [];

set(VevaiometricAxes,...
    'FontSize', 10,...
    'XLim', [-5 5],...
    'YLim', [0 12],...
    'YAxisLocation', 'right')
title(VevaiometricAxes, 'Vevaiometric');
ylabel(VevaiometricAxes, 'Invested Time (s)');

% Vevaiometric z-score
VevaiometricSqrtZScoreAxes = axes(AnalysisFigure, 'Position', [0.01    0.25    0.15    0.19]);
hold(VevaiometricSqrtZScoreAxes, 'on')

AllExploringTISqrtZScore = [];
AllExploitingTISqrtZScore = [];

set(VevaiometricSqrtZScoreAxes,...
    'FontSize', 10,...
    'XLim', [-5 5],...
    'YLim', [-2 2],...
    'YAxisLocation', 'right')
xlabel(VevaiometricSqrtZScoreAxes, 'log(odds)');
ylabel(VevaiometricSqrtZScoreAxes, 'sqrt(Invested Time) (z-score)');

% Vevaiometric (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
LRVevaiometricAxes = axes(AnalysisFigure, 'Position', [0.33    0.48    0.15    0.19]);
hold(LRVevaiometricAxes, 'on')

AllLeftTI = [];
AllRightTI = [];
AllLeftAbsResidual = [];
AllRightAbsResidual = [];

set(LRVevaiometricAxes,...
    'FontSize', 10,...
    'XLim', [0 1],...
    'YLim', [0 12])
title(LRVevaiometricAxes, 'LRVevaiometric');

% Vevaiometric z-score (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
LRVevaiometricSqrtZScoreAxes = axes(AnalysisFigure, 'Position', [0.33    0.25    0.15    0.19]);
hold(LRVevaiometricSqrtZScoreAxes, 'on')

AllLeftTISqrtZScore = [];
AllRightTISqrtZScore = [];

set(LRVevaiometricSqrtZScoreAxes,...
    'FontSize', 10,...
    'XLim', [0 1],...
    'YLim', [-2 2])
xlabel(LRVevaiometricSqrtZScoreAxes, 'abs(Residuals)');

% TI distribution per left-right and explore/exploit
TIDistributionAxes = axes(AnalysisFigure, 'Position', [0.22    0.48    0.09    0.19]);
hold(TIDistributionAxes, 'on')

set(TIDistributionAxes,...
    'TickDir', 'out',...
    'XLim', [-0.5, 1.5],...
    'XTick', [0 1],...
    'XTickLabel', {},...
    'YLim', [0, 12],...
    'FontSize', 10);
xlabel(TIDistributionAxes, '')
title(TIDistributionAxes, 'TI Distribution');

% TI-ZScore distribution per left-right and explore/exploit
TISqrtZScoreDistributionAxes = axes(AnalysisFigure, 'Position', [0.22    0.25    0.09    0.19]);
hold(TISqrtZScoreDistributionAxes, 'on')

set(TISqrtZScoreDistributionAxes,...
    'TickDir', 'out',...
    'XLim', [-0.5, 1.5],...
    'XTick', [0 1],...
    'XTickLabel', {'Left', 'Right'},...
    'YLim', [-2, 2],...
    'FontSize', 10);
xlabel(TISqrtZScoreDistributionAxes, '')
title(TISqrtZScoreDistributionAxes, '');

%% Explore/exploit level around block switch
% Block transition
BlockTransitionAxes = axes(AnalysisFigure, 'Position', [0.01    0.06    0.15    0.12]);
hold(BlockTransitionAxes, 'on');

Block1TransitionYData = nan(100, 61);
Block2TransitionYData = nan(100, 61);
Block3TransitionYData = nan(100, 61);
Block4TransitionYData = nan(100, 61);

set(BlockTransitionAxes,...
    'TickDir', 'out',...
    'XLim', [-10 50],...
    'XTick', [0 20 40],...
    'XTickLabel', [1 21 41],...
    'YLim', [0 40],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(BlockTransitionAxes, 'Choice_{explore} (%)')
title(BlockTransitionAxes, 'Block Transition')

% Explore/exploit level against reward rate
RewardRateAxes = axes(AnalysisFigure, 'Position', [0.33    0.06    0.15    0.12]);
hold(RewardRateAxes, 'on');

AllAbsResidual = [];
AllRewardRate = [];

set(RewardRateAxes,...
    'TickDir', 'out',...
    'XLim', [0 1],...
    'YLim', [0 100],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(RewardRateAxes, 'abs(Residuals)')
ylabel(RewardRateAxes, 'Reward rate')

AllLeftTIRewardRate = [];
AllRightTIRewardRate = [];

%% Move time
% Vevaiometric MT
VevaiometricMTAxes = axes(AnalysisFigure, 'Position', [0.51    0.48    0.15    0.19]);
hold(VevaiometricMTAxes, 'on')

AllExploringMT = [];
AllExploitingMT = [];

AllExploringMTLogOdds = [];
AllExploitingMTLogOdds = [];

set(VevaiometricMTAxes,...
    'FontSize', 10,...
    'XLim', [-5 5],...
    'YLim', [0 0.5],...
    'YAxisLocation', 'right')
title(VevaiometricMTAxes, 'Vevaiometric MT');
ylabel(VevaiometricMTAxes, 'Move Time (s)');

% Vevaiometric MT z-score
%{
% NOT USE AS NOT NORMAL DISTRIBUTION
VevaiometricMTSqrtZScoreAxes = axes(AnalysisFigure, 'Position', [0.51    0.25    0.15    0.19]);
hold(VevaiometricMTSqrtZScoreAxes, 'on')

AllExploringMTSqrtZScore = [];
AllExploitingMTSqrtZScore = [];

set(VevaiometricMTSqrtZScoreAxes,...
    'FontSize', 10,...
    'XLim', [-5 5],...
    'YLim', [0 3],...
    'YAxisLocation', 'right')
xlabel(VevaiometricMTSqrtZScoreAxes, 'log(odds)');
ylabel(VevaiometricMTSqrtZScoreAxes, 'sqrt(Move Time) (z-score)');
%}

% Vevaiometric (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
LRVevaiometricMTAxes = axes(AnalysisFigure, 'Position', [0.83    0.48    0.15    0.19]);
hold(LRVevaiometricMTAxes, 'on')

AllLeftMT = [];
AllRightMT = [];
AllLeftMTAbsResidual = [];
AllRightMTAbsResidual = [];

set(LRVevaiometricMTAxes,...
    'FontSize', 10,...
    'XLim', [0 1],...
    'YLim', [0 0.5])
title(LRVevaiometricMTAxes, 'LRVevaiometric MT');

% Vevaiometric MT z-score (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
%{
% NOT USE AS NOT NORMAL DISTRIBUTION
LRVevaiometricMTSqrtZScoreAxes = axes(AnalysisFigure, 'Position', [0.83    0.25    0.15    0.19]);
hold(LRVevaiometricMTSqrtZScoreAxes, 'on')

AllLeftMTSqrtZScore = [];
AllRightMTSqrtZScore = [];

set(LRVevaiometricMTSqrtZScoreAxes,...
    'FontSize', 10,...
    'XLim', [0 1],...
    'YLim', [0 3])
xlabel(LRVevaiometricMTSqrtZScoreAxes, 'abs(Residuals)');
%}

% TI distribution per left-right and explore/exploit
MTDistributionAxes = axes(AnalysisFigure, 'Position', [0.72    0.48    0.09    0.19]);
hold(MTDistributionAxes, 'on')

set(MTDistributionAxes,...
    'TickDir', 'out',...
    'XLim', [-0.5, 1.5],...
    'XTick', [0 1],...
    'XTickLabel', {},...
    'YLim', [0, 0.5],...
    'FontSize', 10);
xlabel(MTDistributionAxes, '')
title(MTDistributionAxes, 'MT Distribution');

% TI-ZScore distribution per left-right and explore/exploit
%{
% NOT USE AS NOT NORMAL DISTRIBUTION
MTSqrtZScoreDistributionAxes = axes(AnalysisFigure, 'Position', [0.72    0.25    0.09    0.19]);
hold(MTSqrtZScoreDistributionAxes, 'on')

set(MTSqrtZScoreDistributionAxes,...
    'TickDir', 'out',...
    'XLim', [-0.5, 1.5],...
    'XTick', [0 1],...
    'XTickLabel', {'Left', 'Right'},...
    'YLim', [0, 3],...
    'FontSize', 10);
xlabel(MTSqrtZScoreDistributionAxes, '')
title(MTSqrtZScoreDistributionAxes, '');
%}

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')

%% saving
DataFolder = OttLabDataServerFolderPath;
RatName = SessionData.Info.Subject;
% %%The following lines doesn't not work, as the timestamp documented
% in the SessionData may not be the same as the one being used for saving
% SessionDate = string(datetime(SessionData.Info.SessionDate), 'yyyyMMdd')';
% SessionTime = string(datetime(SessionData.Info.SessionStartTime_UTC), 'HHmmSS')';
% SessionDateTime = strcat(SessionDate, '_', SessionTime);
DataPath = strcat(DataFolder, RatName, '\bpod_session\', SessionDateTime, '\',...
                  RatName, '_TwoArmBanditVariant_', SessionDateTime, '_Matching_ChoiceSymmetricQLearning.png');
exportgraphics(AnalysisFigure, DataPath);

DataPath = strcat(DataFolder, RatName, '\bpod_graph\',...
                  RatName, '_TwoArmBanditVariant_', SessionDateTime, '_Matching_ChoiceSymmetricQLearning.png');
exportgraphics(AnalysisFigure, DataPath);

close(AnalysisFigure)

end % function