function AnalysisFigure = Matching_SS_BackgroundDA_TotalValue(DataFile, DataObject, Model, FigureSize, AxeSize)
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
if nargin < 4
    FigureSize = [0.2, 0.5, 4.5, 3.0];
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
    'Color', 'none')

% colour palette
ColourPalette = CommonColourPalette();

%% Symmetric Q-Learning with Forgetting and Stickiness model
% model
if nargin < 2
    [datafile, datapath] = uigetfile(OttLabDataServerFolderPath());
    load(fullfile(datapath, datafile));
elseif ischar(DataObject) || isstring(DataObject)
    load(DataObject);
else
    DataObject = DataObject;
end

try
    TE = DataObject.LoadProcessedPhotometryFile('TE');
catch
    disp('2nd input is not a PhotometryData.mat')
    return
end

PreAlignedData = DataObject.LoadProcessedPhotometryFile('Aligned', 'Channel', 'Green', 'Alignment', 'TimeChoice');
PreAlignedRawPhotometryData = PreAlignedData.Aligned;
PreAlignedBaselineData = PreAlignedData.AlignedBaseline;
if all(isnan(PreAlignedRawPhotometryData), 2)
    return
end

% Normalise to dF/F
NormalisedPreAlignedPhotometryData = (PreAlignedRawPhotometryData ./ PreAlignedBaselineData - 1) * 100;

% Extract dF/F from the interval of interested
PreAlignedTime = PreAlignedData.Time;
FocuseddF = NormalisedPreAlignedPhotometryData(:, PreAlignedTime >= -1.0 & PreAlignedTime <= -0.5);

% Estimate dF
EstimatingFunction = str2func('median');
try
    dFEstimation = EstimatingFunction(FocuseddF, 2);
catch
    disp('Incompatible estimating function for photometry data of nTrial x length of recording')
    return
end
dFEstimation = dFEstimation';

%% Symmetric Q-Learning with Forgetting and Stickiness model
% model
if nargin < 3
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

ChosenValue = LeftValue .* ChoiceLeft + RightValue .* (1 - ChoiceLeft);
UnchosenValue = LeftValue .* (1 - ChoiceLeft) + RightValue .* ChoiceLeft;
TotalValue = ChosenValue + UnchosenValue;

%% Background DA vs Total value
if nargin < 5
    AxeSize = [ 1.2, 0.8, 2.5, 2.0];
end

BackgroundDATotalValueAxes = axes(AnalysisFigure,...
                                  'Position', AxeSize,...
                                  'Units', 'inches',...
                                  'Color', 'none');
set(BackgroundDATotalValueAxes,...
    'Position', AxeSize,...
    'Units', 'inches');
hold(BackgroundDATotalValueAxes, 'on');

BackgroundDATotalValueScatter = scatter(BackgroundDATotalValueAxes, TotalValue, dFEstimation,...
                                        'Marker', '.',...
                                        'MarkerEdgeColor', 'k',...
                                        'SizeData', 1);

ValidIdx = ~isnan(ChosenValue) & ~isnan(dFEstimation);

[RValue, pValue] = corrcoef(TotalValue(ValidIdx), dFEstimation(ValidIdx));

pValueSymbol = '';
if pValue(1, 2) < 0.001
    pValueSymbol = '***';
elseif pValue(1, 2) < 0.01
    pValueSymbol = '**';
elseif pValue(1, 2) < 0.05
    pValueSymbol = '*';
end

BackgroundDATotalValueStatText = text(BackgroundDATotalValueAxes, 0.1, 0.2,...
                                      strcat(sprintf('R = %4.2f',...
                                                     RValue(1, 2)),...
                                             pValueSymbol),...
                                      'Units','inches',...
                                      'FontSize', 18,...
                                      'Color', 'k');

X = [ones(1, nTrials); TotalValue]';
b = X(ValidIdx, :) \ dFEstimation(ValidIdx)';

XData = [0, 0.6];
YData = XData * b(2) + b(1);

BackgroundDATotalValueRegressionLine = line(BackgroundDATotalValueAxes, XData, YData,...
                                            'Color', 'k',...
                                            'LineWidth', 1);

set(BackgroundDATotalValueAxes,...
    'TickDir', 'out',...
    'XLim', [0, 0.6],...
    'XTick', [0, 0.2, 0.4, 0.6],...
    'XColor', ColourPalette.Incorrect,...
    'YLim', [-10, 10],...
    'YTick', [-10, 0, 10],...
    'YColor', ColourPalette.Rewarded,...
    'FontSize', 18);
xlabel(BackgroundDATotalValueAxes, 'Reward rate (a.u.)')
ylabel(BackgroundDATotalValueAxes, sprintf('Baseline\ndopamine (%%)'))

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateTime, '_Matching_BackgroundDA_TotalValue.pdf'), 'ContentType', 'vector', 'Resolution', 300)
end