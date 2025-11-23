function AnalysisFigure = Matching_SS_BlockSwitch_LogOdds_Traces(DataFile, Model, FigureSize, AxeSize)
% for making an example of the task structure
% figure 'unit' set as 'inch' so that we know exactly the meters to pixels
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ April 2025

if nargin < 1
    global BpodSystem
    if isempty(BpodSystem) || isempty(BpodSystem.Data)
        [datafile, datapath] = uigetfile(OttLabDataServerFolderPath());
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
if nTrials < 200
    disp('nTrial < 200. Impossible for analysis.')
    AnalysisFigure = [];
    return
end

ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
if isempty(ChoiceLeft) || all(isnan(ChoiceLeft))
    disp('No choice made. Impossible for analysis.')
    AnalysisFigure = [];
    return
end

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

%% Initiatize figure
% create figure
if nargin < 3
    FigureSize = [0.5, 0.5, 7.1, 3.9];
end
AnalysisFigure = figure('Position', FigureSize,...
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateTime, '_Matching'),...
                        'MenuBar', 'none',...
                        'Resize', 'off',...
                        'Unit', 'inch',...
                        'Color', 'none');

set(AnalysisFigure,...
    'Position', FigureSize,...
    'unit', 'inch');

FrameAxes = axes(AnalysisFigure, 'Position', [0 0 1 1]); % spacer for correct saving dimension
set(FrameAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'none',...
    'YColor', 'none',...
    'Color', 'w')

% colour palette
ColourPalette = CommonColourPalette();

%% Symmetric Q-Learning with Forgetting and Stickiness model
% model
if nargin < 2
    try
        Model = Matching_SS_MLE_ChoiceSymmetricQLearning_Model(SessionData);
    catch
        disp('Error: problem in modelling. N further analysis is possible')
        return
    end
end

EstimatedParameters = [];
if isfield(Model, 'EstimatedParameters')
    EstimatedParameters = Model.EstimatedParameters;
elseif isfield(Model, 'MAPParameters')
    EstimatedParameters = Model.MAPParameters;
end

if all(isnan(EstimatedParameters))
    disp('Error: fail to run model');
    return
end

% extract parameters
LearningRate = EstimatedParameters(1); % alpha
InverseTemperature = EstimatedParameters(2); % beta
ForgettingRate = EstimatedParameters(3); % gamma
ChoiceStickiness = EstimatedParameters(4); % phi
ChoiceForgettingRate = EstimatedParameters(5); % gamma_c
Bias = EstimatedParameters(6);

% make prediction from Model (i.e. estimated parameters)
[~, Values] = ChoiceSymmetricQLearning(EstimatedParameters, nTrials, ChoiceLeft, Rewarded);
LeftValue = Values.LeftValue;
RightValue = Values.RightValue;
ChoiceMemory = Values.ChoiceMemory;

LogOdds = InverseTemperature * (LeftValue - RightValue) +...
          + ChoiceStickiness * ChoiceMemory + Bias;

PredictedLeftChoiceProb = 1 ./ (1 + exp(-LogOdds));
ModelResiduals = ChoiceLeft - PredictedLeftChoiceProb;
AbsModelResiduals = abs(ModelResiduals);

PredictedChoice = double(PredictedLeftChoiceProb >= 0.5);
PredictedChoice(isnan(ChoiceLeft)) = nan;

%% Total value
if nargin < 4
    AxeSize = [ 1.3, 1.1, 4.8, 2.4];
end

TracesAxes = axes(AnalysisFigure,...
                  'Position', AxeSize,...
                  'Units', 'inches',...
                  'Color', 'none');
set(TracesAxes,...
    'Position', AxeSize,...
    'Units', 'inches');
hold(TracesAxes, 'on');

if ~isempty(ChoiceLeft) && ~all(isnan(ChoiceLeft))
    % yyaxis(TracesAxes, 'left')

    idxTrial = 1:nTrials;
    RewardProbLeft = RewardProb(1,:);
    RewardProbRight = RewardProb(2,:);
    RewardProbLeftPlot = plot(TracesAxes, idxTrial, RewardProbLeft * 100,...
                              'LineStyle', '-',...
                              'Marker', 'none',...
                              'Color', 0.5 * [1, 1, 1],...
                              'LineWidth', 1);
    %{
    RewardProbRightPlot = plot(TracesAxes, idxTrial, RewardProbRight * 100,...
                               'LineStyle', '-',...
                               'Marker', 'none',...
                               'Color', ColourPalette.Right,...
                               'LineWidth', 1);
    %}

    BinWidth = 15;
    BlockChoicePlot = plot(TracesAxes, idxTrial, movmean(ChoiceLeft, BinWidth, 'omitnan') * 100,...
                           'LineStyle', '-',...
                           'Color', 'k',...
                           'Marker', 'none',...
                           'LineWidth', 1);
    % TrialLeftChoicePlot = plot(BlockSwitchAxes, idxTrial(ChoiceLeft==1), ChoiceLeft(ChoiceLeft==1) * 110 - 5,...
    %                            'LineStyle', 'none',...
    %                            'Marker', '.',...
    %                            'MarkerEdgeColor', sand, ...
    %                            'MarkerSize', 8);
    % TrialRightChoicePlot = plot(BlockSwitchAxes, idxTrial(ChoiceLeft==0), ChoiceLeft(ChoiceLeft==0) * 110 - 5,...
    %                             'LineStyle', 'none',...
    %                             'Marker', '.',...
    %                             'MarkerEdgeColor', turquoise,...
    %                             'MarkerSize', 8);
    
    set(TracesAxes,...
        'TickDir', 'out',...
        'XLim', [100, 500],...
        'XTick', 100:200:500,...
        'YColor', 'k',...
        'YLim', [0 120],...
        'YTick', [0, 100],...
        'FontSize', 24);
    xlabel(TracesAxes, 'Trial number')
    ylabel(TracesAxes, 'P(Left) (%)')
    % title('Block switching behaviour')
    %{
    yyaxis(TracesAxes, 'right')
    
    LogOddsPlot = plot(TracesAxes, 1:nTrials, movmean(LogOdds, 3),...
                       'LineStyle', '-',...
                       'Color', ColourPalette.RL,...
                       'LineWidth', 1);
    %}
    PredictedLeftChoiceProbPlot = plot(TracesAxes, 1:nTrials, movmean(PredictedLeftChoiceProb, BinWidth) * 100,...
                                       'LineStyle', '-',...
                                       'Color', ColourPalette.RL,...
                                       'LineWidth', 1);

    %LegendString = {'$\textsf{P(r)}_\textsf{L}$', '$\textsf{P(r)}_\textsf{R}$', '$\textsf{P(c=L)}$'};
    LegendString = {'P(reward)_L', 'Rat', 'Model'};
    
    %warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
    
    BlockSwitchLegend = legend(TracesAxes, LegendString,...
                               'Position', [ 1.4  2.8,  4.4,  0.6],...
                               'NumColumns', 3,...
                               'Units', 'inches',...
                               'Color', 'none',...
                               'Box', 'off');
    set(BlockSwitchLegend,... % 'Interpreter', 'latex',...
        'Position', [ 2.2  2.8,  3.6,  0.6],...
        'Units', 'inches');
    %{
    set(TracesAxes,...
        'YColor', ColourPalette.RL,...
        'YLim', [-3, 6],...
        'YTick', [0, 5],...
        'FontSize', 12);
    ylabel(TracesAxes,'Q_{L-R}')
    %}
end

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateTime, '_Matching_BlockSwitch_LogOdds_Traces.pdf'), 'ContentType', 'vector', 'Resolution', 300)
end