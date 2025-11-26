function AnalysisFigure = Matching_SS_B_ChoiceSymmetricQLearning_Diagnosis(DataFile, Model)
% SS = SingleSession
% MLE = Maximum Log Likelihood
% Matching Analysis Function
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ July 2024
% Version 2.0 ~ Dec 2024 <- Convert to symmertric with forgetting
% Model iteration see the end of script

if nargin < 1
    global BpodSystem
    if isempty(BpodSystem) || isempty(BpodSystem.Data)
        [datafile, datapath] = uigetfile('\\ottlabfs.bccn-berlin.pri\ottlab\data\');
        load(fullfile(datapath, datafile));
        SessionDateTime = datapath(end-15:end-1);
    else
        SessionData = BpodSystem.Data;
        [~, name, ~] = fileparts(BpodSystem.Path.CurrentDataFile);
        SessionDateTime = name(end-14:end);
    end
elseif ischar(DataFile) || isstring(DataFile)
    load(DataFile);
    SessionDateTime = DataFile(end-18:end-4);
elseif isstruct(DataFile)
    SessionData = DataFile;

    % mismatch in time saved in .mat and the time used as file name
    SessionDateTime = strcat(datestr(SessionData.Info.SessionDate, 'yyyymmdd'), '_000000');
else
    disp('Error: Unknown input format. No further analysis can be performed.')
    return
end

if ~isfield(SessionData, 'SettingsFile')
    disp('Error: The selected file does not have the field "SettingsFile". No further Matching Analysis is performed.')
    AnalysisFigure = [];
    return
elseif ~isfield(SessionData.SettingsFile.GUIMeta, 'RiskType')
    disp('Error: The selected SessionFile may not be a TwoArmBanditVariant session. No further Matching Analysis is performed.')
    AnalysisFigure = [];
    return
elseif ~strcmpi(SessionData.SettingsFile.GUIMeta.RiskType.String{SessionData.SettingsFile.GUI.RiskType}, 'BlockFixHolding')
    disp('Error: The selected SessionData is not a Matching session. No further Matching Analysis is performed.')
    AnalysisFigure = [];
    return
end

AnalysisName = 'Matching_SS_B_ChoiceSymmetricQLearning_Diagnosis';

%% Load related data to local variabels
RatID = str2double(SessionData.Info.Subject);
if isnan(RatID)
    RatID = -1;
end
RatName = num2str(RatID);
% %%The following three lines doesn't not work, as the timestamp documented
% in the SessionData may not be the same as the one being used for saving
% Date = datestr(SessionData.Info.SessionDate, 'yyyymmdd');

nTrials = SessionData.nTrials;
if nTrials < 50
    disp('nTrial < 50. Impossible for analysis.')
    AnalysisFigure = [];
    return
end
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Baited = SessionData.Custom.TrialData.Baited(:, 1:nTrials);
IncorrectChoice = SessionData.Custom.TrialData.IncorrectChoice(1:nTrials);
NoDecision = SessionData.Custom.TrialData.NoDecision(1:nTrials);
NoTrialStart = SessionData.Custom.TrialData.NoTrialStart(1:nTrials);
BrokeFixation = SessionData.Custom.TrialData.BrokeFixation(1:nTrials);
EarlyWithdrawal = SessionData.Custom.TrialData.EarlyWithdrawal(1:nTrials);
StartNewTrial = SessionData.Custom.TrialData.StartNewTrial(1:nTrials);
SkippedFeedback = SessionData.Custom.TrialData.SkippedFeedback(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

SampleTime = SessionData.Custom.TrialData.SampleTime(1:nTrials);
MoveTime = SessionData.Custom.TrialData.MoveTime(1:nTrials);
FeedbackWaitingTime = SessionData.Custom.TrialData.FeedbackWaitingTime(1:nTrials);
% FeedbackDelay = SessionData.Custom.TrialData.FeedbackDelay(1:nTrials);
% FeedbackWaitingTime = rand(nTrials,1)*10; %delete this
% FeedbackWaitingTime = FeedbackWaitingTime';  %delete this
% FeedbackDelay = rand(nTrials,1)*10; %delete this
% FeedbackDelay= FeedbackDelay'; 

RewardProb = SessionData.Custom.TrialData.RewardProb(:, 1:nTrials);
LightLeft = SessionData.Custom.TrialData.LightLeft(1:nTrials);
LightLeftRight = [LightLeft; 1-LightLeft]; 
ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft]; 

BlockNumber = SessionData.Custom.TrialData.BlockNumber(:, 1:nTrials);
BlockTrialNumber = SessionData.Custom.TrialData.BlockTrialNumber(:, 1:nTrials);

% for files before April 2023, no DrinkingTime is available
try
    DrinkingTime = SessionData.Custom.TrialData.DrinkingTime(1:nTrials);
catch
    DrinkingTime = nan(1, nTrials);
end

LeftFeedbackDelayGraceTime = [];
RightFeedbackDelayGraceTime = [];
FirstDrinkingTime = [];
LatestRewardTimestamp = [];
for iTrial = 1:nTrials
    if ChoiceLeft(iTrial) == 1
        LeftFeedbackDelayGraceTime = [LeftFeedbackDelayGraceTime;...
                                      SessionData.RawEvents.Trial{iTrial}.States.LInGrace(:,2) -...
                                      SessionData.RawEvents.Trial{iTrial}.States.LInGrace(:,1)];
    elseif ChoiceLeft(iTrial) == 0
        RightFeedbackDelayGraceTime = [RightFeedbackDelayGraceTime;...
                                       SessionData.RawEvents.Trial{iTrial}.States.RInGrace(:,2) -...
                                       SessionData.RawEvents.Trial{iTrial}.States.RInGrace(:,1)];
    end
    
    FirstDrinkingTime = [FirstDrinkingTime SessionData.RawEvents.Trial{iTrial}.States.Drinking(1,1)];
    if iTrial == 1
        LatestRewardTimestamp(iTrial) = 0;
    elseif isnan(SessionData.RawEvents.Trial{iTrial-1}.States.Drinking(1,1))
        LatestRewardTimestamp(iTrial) = LatestRewardTimestamp(iTrial-1);
    else
        LatestRewardTimestamp(iTrial) = SessionData.RawEvents.Trial{iTrial-1}.States.Drinking(1,1) + SessionData.TrialStartTimestamp(iTrial-1);
    end
end
LatestRewardTime = SessionData.TrialStartTimestamp - LatestRewardTimestamp;

LeftFeedbackDelayGraceTime = LeftFeedbackDelayGraceTime(~isnan(LeftFeedbackDelayGraceTime))';
LeftFeedbackDelayGraceTime = LeftFeedbackDelayGraceTime(LeftFeedbackDelayGraceTime < SessionData.SettingsFile.GUI.FeedbackDelayGrace - 0.0001);
RightFeedbackDelayGraceTime = RightFeedbackDelayGraceTime(~isnan(RightFeedbackDelayGraceTime))';
RightFeedbackDelayGraceTime = RightFeedbackDelayGraceTime(RightFeedbackDelayGraceTime < SessionData.SettingsFile.GUI.FeedbackDelayGrace - 0.0001);

%% Common plots regardless of task design/ risk type
% create figure
AnalysisFigure = figure('Position', [   0    0 1191  842],... % DIN A3, 72 ppi (window will crop it to _ x 1024, same as disp resolution)
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateTime, '_Matching'),...
                        'MenuBar', 'none',...
                        'Resize', 'off');

% spacer for correct saving dimension
FrameAxes = axes(AnalysisFigure, 'Position', [0 0 1 1]); % spacer for correct saving dimension
set(FrameAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'w',...
    'YColor', 'w')

% Figure Info
FigureInfoAxes = axes(AnalysisFigure, 'Position', [0.01    0.96    0.48    0.01]);
set(FigureInfoAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'w',...
    'YColor', 'w')

FigureTitle = strcat(RatName, '_', SessionDateTime, '_', AnalysisName);

FigureTitleText = text(FigureInfoAxes, 0, 0,...
                       FigureTitle,...
                       'FontSize', 14,...
                       'FontWeight','bold',...
                       'Interpreter', 'none');

% colour palette
ColourPalette = CommonColourPalette();

%% Symmetric Q-Learning with Forgetting and Stickiness model
% model
if nargin < 2
    try
        Model = Matching_SS_B_ChoiceSymmetricQLearning_Model(SessionData);
    catch
        disp('Error: problem in modelling. N further analysis is possible')
        return
    end
end


%% extract info of the chains
ChainInitialParameters = [Model.ChainInitialParameters{:}];

BurnIn = Model.Prior.BurnIn;
nSample = Model.Prior.nSample;
ChainIdx = (BurnIn + 1):(BurnIn + nSample);

nChain = Model.Prior.nChain;

Chains = Model.Chains;
ConcatenatedChain = vertcat(Model.Chains{:});

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 2));
InverseTemperature = Values(Posterior == max(Posterior));

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 3));
ForgettingRate = Values(Posterior == max(Posterior));

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 4));
ChoiceStickiness = Values(Posterior == max(Posterior));

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 5));
ChoiceForgettingRate = Values(Posterior == max(Posterior));

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 6));
Bias = Values(Posterior == max(Posterior));

%% Chain progression
% learning rate chain
LearningRateChainAxes = axes(AnalysisFigure, 'Position', [0.01    0.61    0.28    0.08]);
hold(LearningRateChainAxes, 'on');

for iChain = 1:nChain
    LearningRateChainPlot(iChain) = plot(LearningRateChainAxes, ChainIdx, Chains{iChain}(:, 1),...
                                         'Color', ColourPalette.Session);

    LearningRateBurnInPlot(iChain) = plot(LearningRateChainAxes, [0, BurnIn + 1], [ChainInitialParameters(1, iChain), Chains{iChain}(1, 1)],...
                                          'Color', ColourPalette.Session,...
                                          'Marker', 'o',...
                                          'MarkerSize', 4);
    
end

ChainLegend = legend(LearningRateChainAxes, {'iChain', 'iBurnIn'},...
                     'Location', 'northeast',...
                     'NumColumns', 3);

set(LearningRateChainAxes,...
    'TickDir', 'out',...
    'XLim', [0, BurnIn + nSample],...
    'XTickLabel', {},...
    'YAxisLocation', 'right',...
    'FontSize', 10);
title(LearningRateChainAxes, 'Chain')
ylabel(LearningRateChainAxes, '\alpha')

% inverse temperature chain
InverseTemperatureChainAxes = axes(AnalysisFigure, 'Position', [0.01    0.50    0.28    0.08]);
hold(InverseTemperatureChainAxes, 'on');

for iChain = 1:nChain
    InverseTemperatureChainPlot(iChain) = plot(InverseTemperatureChainAxes, ChainIdx, Chains{iChain}(:, 2),...
                                               'Color', ColourPalette.Session);

    InverseTemperatureBurnInPlot(iChain) = plot(InverseTemperatureChainAxes, [0, BurnIn + 1], [ChainInitialParameters(2, iChain), Chains{iChain}(1, 2)],...
                                                'Color', ColourPalette.Session,...
                                                'Marker', 'o',...
                                                'MarkerSize', 4);
    
end

set(InverseTemperatureChainAxes,...
    'TickDir', 'out',...
    'XLim', [0, BurnIn + nSample],...
    'XTickLabel', {},...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(InverseTemperatureChainAxes, '\beta')

% forgetting rate chain
ForgettingRateChainAxes = axes(AnalysisFigure, 'Position', [0.01    0.39    0.28    0.08]);
hold(ForgettingRateChainAxes, 'on');

for iChain = 1:nChain
    ForgettingRateChainPlot(iChain) = plot(ForgettingRateChainAxes, ChainIdx, Chains{iChain}(:, 3),...
                                           'Color', ColourPalette.Session);

    ForgettingRateBurnInPlot(iChain) = plot(ForgettingRateChainAxes, [0, BurnIn + 1], [ChainInitialParameters(3, iChain), Chains{iChain}(1, 3)],...
                                            'Color', ColourPalette.Session,...
                                            'Marker', 'o',...
                                            'MarkerSize', 4);
    
end

set(ForgettingRateChainAxes,...
    'TickDir', 'out',...
    'XLim', [0, BurnIn + nSample],...
    'XTickLabel', {},...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(ForgettingRateChainAxes, '\gamma')

% choice stickiness chain
ChoiceStickinessChainAxes = axes(AnalysisFigure, 'Position', [0.01    0.28    0.28    0.08]);
hold(ChoiceStickinessChainAxes, 'on');

for iChain = 1:nChain
    ChoiceStickinessChainPlot(iChain) = plot(ChoiceStickinessChainAxes, ChainIdx, Chains{iChain}(:, 4),...
                                             'Color', ColourPalette.Session);

    ChoiceStickinessBurnInPlot(iChain) = plot(ChoiceStickinessChainAxes, [0, BurnIn + 1], [ChainInitialParameters(4, iChain), Chains{iChain}(1, 4)],...
                                              'Color', ColourPalette.Session,...
                                              'Marker', 'o',...
                                              'MarkerSize', 4);
    
end

set(ChoiceStickinessChainAxes,...
    'TickDir', 'out',...
    'XLim', [0, BurnIn + nSample],...
    'XTickLabel', {},...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(ChoiceStickinessChainAxes, '\phi')

% choice forgetting rate chain
ChoiceForgettingRateChainAxes = axes(AnalysisFigure, 'Position', [0.01    0.17    0.28    0.08]);
hold(ChoiceForgettingRateChainAxes, 'on');

for iChain = 1:nChain
    ChoiceForgettingRateChainPlot(iChain) = plot(ChoiceForgettingRateChainAxes, ChainIdx, Chains{iChain}(:, 5),...
                                                 'Color', ColourPalette.Session);

    ChoiceForgettingRateBurnInPlot(iChain) = plot(ChoiceForgettingRateChainAxes, [0, BurnIn + 1], [ChainInitialParameters(5, iChain), Chains{iChain}(1, 5)],...
                                                  'Color', ColourPalette.Session,...
                                                  'Marker', 'o',...
                                                  'MarkerSize', 4);
    
end

set(ChoiceForgettingRateChainAxes,...
    'TickDir', 'out',...
    'XLim', [0, BurnIn + nSample],...
    'XTickLabel', {},...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(ChoiceForgettingRateChainAxes, '\gamma_c')

% bias chain
BiasChainAxes = axes(AnalysisFigure, 'Position', [0.01    0.06    0.28    0.08]);
hold(BiasChainAxes, 'on');

for iChain = 1:nChain
    BiasChainPlot(iChain) = plot(BiasChainAxes, ChainIdx, Chains{iChain}(:, 6),...
                                 'Color', ColourPalette.Session);

    BiasBurnInPlot(iChain) = plot(BiasChainAxes, [0, BurnIn + 1], [ChainInitialParameters(6, iChain), Chains{iChain}(1, 6)],...
                                  'Color', ColourPalette.Session,...
                                  'Marker', 'o',...
                                  'MarkerSize', 4);
    
end

set(BiasChainAxes,...
    'TickDir', 'out',...
    'XLim', [0, BurnIn + nSample],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(BiasChainAxes, 'iSample')
ylabel(BiasChainAxes, 'Bias')

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')

%% kernel density of prior, likelihood, posterior
% learning rate
LearningRateDensityAxes = axes(AnalysisFigure, 'Position', [0.36    0.61    0.05    0.08]);
hold(LearningRateDensityAxes, 'on')

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 1));
LearningRateMAP = Values(Posterior == max(Posterior));

LearningRatePosteriorDensityPlot = plot(LearningRateDensityAxes, Posterior, Values,...
                                        'LineStyle', '-',...
                                        'Color', 'k');

Prior = betapdf(Values, Model.Prior.LearningRateAlpha, Model.Prior.LearningRateBeta);
LearningRatePriorDensityPlot = plot(LearningRateDensityAxes, Prior, Values,...
                                    'LineStyle', '--',...
                                    'Color', 'k');

Likelihood = Posterior ./ Prior;
LearningRateLikelihoodDensityPlot = plot(LearningRateDensityAxes, Likelihood, Values,...
                                         'LineStyle', '-',...
                                         'Color', ColourPalette.Session);

set(LearningRateDensityAxes,...
    'FontSize', 10,...
    'TickDir', 'out')
title(LearningRateDensityAxes, 'Probability density')

% inverse temperature
InverseTemperatureDensityAxes = axes(AnalysisFigure, 'Position', [0.36    0.50    0.05    0.08]);
hold(InverseTemperatureDensityAxes, 'on')

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 2));
InverseTemperatureMAP = Values(Posterior == max(Posterior));

InverseTemperaturePosteriorDensityPlot = plot(InverseTemperatureDensityAxes, Posterior, Values,...
                                              'LineStyle', '-',...
                                              'Color', 'k');

Prior = normpdf(Values, Model.Prior.InverseTemperatureMean, Model.Prior.InverseTemperatureSigma);
InverseTemperaturePriorDensityPlot = plot(InverseTemperatureDensityAxes, Prior, Values,...
                                          'LineStyle', '--',...
                                          'Color', 'k');

Likelihood = Posterior ./ Prior;
InverseTemperatureLikelihoodDensityPlot = plot(InverseTemperatureDensityAxes, Likelihood, Values,...
                                               'LineStyle', '-',...
                                               'Color', ColourPalette.Session);

set(InverseTemperatureDensityAxes,...
    'FontSize', 10,...
    'TickDir', 'out')

% forgetting rate
ForgettingRateDensityAxes = axes(AnalysisFigure, 'Position', [0.36    0.39    0.05    0.08]);
hold(ForgettingRateDensityAxes, 'on')

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 3));
ForgettingRateMAP = Values(Posterior == max(Posterior));

ForgettingRatePosteriorDensityPlot = plot(ForgettingRateDensityAxes, Posterior, Values,...
                                          'LineStyle', '-',...
                                          'Color', 'k');

Prior = betapdf(Values, Model.Prior.ForgettingRateAlpha, Model.Prior.ForgettingRateBeta);
ForgettingRatePriorDensityPlot = plot(ForgettingRateDensityAxes, Prior, Values,...
                                      'LineStyle', '--',...
                                      'Color', 'k');

Likelihood = Posterior ./ Prior;
ForgettingRateLikelihoodDensityPlot = plot(ForgettingRateDensityAxes, Likelihood, Values,...
                                           'LineStyle', '-',...
                                           'Color', ColourPalette.Session);

set(ForgettingRateDensityAxes,...
    'FontSize', 10,...
    'TickDir', 'out')

% choice stickiness
ChoiceStickinessDensityAxes = axes(AnalysisFigure, 'Position', [0.36    0.28    0.05    0.08]);
hold(ChoiceStickinessDensityAxes, 'on')

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 4));
ChoiceStickinessMAP = Values(Posterior == max(Posterior));

ChoiceStickinessPosteriorDensityPlot = plot(ChoiceStickinessDensityAxes, Posterior, Values,...
                                            'LineStyle', '-',...
                                            'Color', 'k');

Prior = normpdf(Values, Model.Prior.ChoiceStickinessMean, Model.Prior.ChoiceStickinessSigma);
ChoiceStickinessPriorDensityPlot = plot(ChoiceStickinessDensityAxes, Prior, Values,...
                                        'LineStyle', '--',...
                                        'Color', 'k');

Likelihood = Posterior ./ Prior;
ChoiceStickinessLikelihoodDensityPlot = plot(ChoiceStickinessDensityAxes, Likelihood, Values,...
                                             'LineStyle', '-',...
                                             'Color', ColourPalette.Session);

set(ChoiceStickinessDensityAxes,...
    'FontSize', 10,...
    'TickDir', 'out')

% choice forgetting rate
ChoiceForgettingRateDensityAxes = axes(AnalysisFigure, 'Position', [0.36    0.17    0.05    0.08]);
hold(ChoiceForgettingRateDensityAxes, 'on')

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 5));
ChoiceForgettingRateMAP = Values(Posterior == max(Posterior));

ChoiceForgettingRatePosteriorDensityPlot = plot(ChoiceForgettingRateDensityAxes, Posterior, Values,...
                                                'LineStyle', '-',...
                                                'Color', 'k');

Prior = betapdf(Values, Model.Prior.ChoiceForgettingRateAlpha, Model.Prior.ChoiceForgettingRateBeta);
ChoiceForgettingRatePriorDensityPlot = plot(ChoiceForgettingRateDensityAxes, Prior, Values,...
                                            'LineStyle', '--',...
                                            'Color', 'k');

Likelihood = Posterior ./ Prior;
ChoiceForgettingRateLikelihoodDensityPlot = plot(ChoiceForgettingRateDensityAxes, Likelihood, Values,...
                                                 'LineStyle', '-',...
                                                 'Color', ColourPalette.Session);

set(ChoiceForgettingRateDensityAxes,...
    'FontSize', 10,...
    'TickDir', 'out')

% bias
BiasDensityAxes = axes(AnalysisFigure, 'Position', [0.36    0.06    0.05    0.08]);
hold(BiasDensityAxes, 'on')

[Posterior, Values] = ksdensity(ConcatenatedChain(:, 6));
BiasMAP = Values(Posterior == max(Posterior));

BiasPosteriorDensityPlot = plot(BiasDensityAxes, Posterior, Values,...
                                'LineStyle', '-',...
                                'Color', 'k');

Prior = normpdf(Values, Model.Prior.BiasMean, Model.Prior.BiasSigma);
BiasPriorDensityPlot = plot(BiasDensityAxes, Prior, Values,...
                            'LineStyle', '--',...
                            'Color', 'k');

Likelihood = Posterior ./ Prior;
BiasLikelihoodDensityPlot = plot(BiasDensityAxes, Likelihood, Values,...
                                 'LineStyle', '-',...
                                 'Color', ColourPalette.Session);

set(BiasDensityAxes,...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(BiasDensityAxes, 'p.d.f.')

%% sample distribution of posterior
% learning rate
LearningRateSampleAxes = axes(AnalysisFigure, 'Position', [0.44    0.61    0.07    0.08]);
hold(LearningRateSampleAxes, 'on')

LearningRateSampleHistogram = histogram(LearningRateSampleAxes, ConcatenatedChain(:, 1), 30,...
                                        'FaceColor', 'w');

set(LearningRateSampleAxes,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Box', 'off',...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(LearningRateSampleAxes, '\alpha')
ylabel(LearningRateSampleAxes, 'Count')
title(LearningRateSampleAxes, 'Posterior sample')

% inverse temperature
InverseTemperatureSampleAxes = axes(AnalysisFigure, 'Position', [0.53    0.50    0.07    0.08]);
hold(InverseTemperatureSampleAxes, 'on')

InverseTemperatureSampleHistogram = histogram(InverseTemperatureSampleAxes, ConcatenatedChain(:, 2), 30,...
                                              'FaceColor', 'w');

set(InverseTemperatureSampleAxes,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Box', 'off',...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(InverseTemperatureSampleAxes, '\beta')

% forgetting rate
ForgettingRateSampleAxes = axes(AnalysisFigure, 'Position', [0.62    0.39    0.07    0.08]);
hold(ForgettingRateSampleAxes, 'on')

ForgettingRateSampleHistogram = histogram(ForgettingRateSampleAxes, ConcatenatedChain(:, 3), 30,...
                                          'FaceColor', 'w');

set(ForgettingRateSampleAxes,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Box', 'off',...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(ForgettingRateSampleAxes, '\gamma')

% choice stickiness
ChoiceStickinessSampleAxes = axes(AnalysisFigure, 'Position', [0.71    0.28    0.07    0.08]);
hold(ChoiceStickinessSampleAxes, 'on')

ChoiceStickinessSampleHistogram = histogram(ChoiceStickinessSampleAxes, ConcatenatedChain(:, 4), 30,...
                                            'FaceColor', 'w');

set(ChoiceStickinessSampleAxes,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Box', 'off',...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(ChoiceStickinessSampleAxes, '\phi')

% choice forgetting rate
ChoiceForgettingRateSampleAxes = axes(AnalysisFigure, 'Position', [0.80    0.17    0.07    0.08]);
hold(ChoiceForgettingRateSampleAxes, 'on')

ChoiceForgettingRateSampleHistogram = histogram(ChoiceForgettingRateSampleAxes, ConcatenatedChain(:, 5), 30,...
                                                'FaceColor', 'w');

set(ChoiceForgettingRateSampleAxes,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Box', 'off',...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(ChoiceForgettingRateSampleAxes, '\gamma_c')

% bias
BiasSampleAxes = axes(AnalysisFigure, 'Position', [0.89    0.06    0.07    0.08]);
hold(BiasSampleAxes, 'on')

BiasSampleHistogram = histogram(BiasSampleAxes, ConcatenatedChain(:, 6), 30,...
                                'FaceColor', 'w');

set(BiasSampleAxes,...
    'XAxisLocation', 'top',...
    'YAxisLocation', 'right',...
    'Box', 'off',...
    'FontSize', 10,...
    'TickDir', 'out')
xlabel(BiasSampleAxes, 'Bias')

%% sample joint posterior
BinValue = [];
ColourMap = gray;
ColourMap = flipud(ColourMap);
for i = 1:5
    for j = i+1:6
        JointSampleAxes(i, j) = axes(AnalysisFigure, 'Position', [0.35 + 0.09*i, 0.72 - 0.11*j, 0.07, 0.08]);
        
        JointSampleBinscatter(i, j) = binscatter(JointSampleAxes(i, j), ConcatenatedChain(:, i), ConcatenatedChain(:, j), [30, 30]);
        
        BinValue = [BinValue; JointSampleBinscatter(i, j).Values];

        set(JointSampleAxes(i, j),...
            'XAxisLocation', 'top',...
            'XTickLabel', {},...
            'YTickLabel', {},...
            'Box', 'off',...
            'FontSize', 10,...
            'TickDir', 'out')
        colormap(JointSampleAxes(i, j), ColourMap);
        colorbar(JointSampleAxes(i, j), 'off')
    end
end

CLim = [min(BinValue, [], 'all'), max(BinValue, [], 'all')];
for i = 1:5
    for j = i+1:6
        clim(JointSampleAxes(i, j), CLim);
    end
end

ColourBarAxes = axes(AnalysisFigure, 'Position', [0.44, 0.06, 0.43, 0]);
set(ColourBarAxes,...
    'XColor', 'w',...
    'YColor', 'w',...
    'FontSize', 10,...
    'XTick', []);
clim(ColourBarAxes, CLim);
colormap(ColourBarAxes, ColourMap);
ColourBar = colorbar(ColourBarAxes, 'location', 'southoutside');
ColourBar.Label.String = 'Count';

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