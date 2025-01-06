function TwoArmBanditVariant_Cued_MS_LeeGLM(DataFolderPath)
%{
MS = MultiSession
by Antonio Lee for AG Ott @HU Berlin

V1.0 20241115 dedicated script for analyzing multiple sessions with similar
settings (no checking is done). A folder is selected instead of a .mat file
(so that a bit of back and forth looking at the concatenated data and
individual session data is allowed)

The analysis is positioned to perform GLM on data from multiple session,
such that the statistical power is greater (esp. for high rewarding cue)
%}

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

AnalysisName = 'Cued_MultiSession_LeeGLM';

%% Check if all sessions are of the same SettingsFile
%{
% not sure how to best do it yet, as some settings are drawn in each trial
and displayed
for iSession = 1:length(DataHolder)
    SessionData = DataHolder{iSession};

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
end
%}
    
%% Common plots regardless of task design/ risk type
% colour palette for events (suitable for most colourblind people)
scarlet = [254, 60, 60]/255; % for incorrect sign, contracting with azure
denim = [31, 54, 104]/255; % mainly for unsuccessful trials
azure = [0, 162, 254]/255; % for rewarded sign

neon_green = [26, 255, 26]/255; % for NotBaited
neon_purple = [168, 12, 180]/255; % for SkippedBaited

sand = [225, 190 106]/255; % for left-right
turquoise = [64, 176, 166]/255;
LRPalette = [sand; turquoise];

% colour palette for cues: (1- P(r)) * 128 + 127
% P(0) = white; P(1) = smoky gray
% RewardProbCategories = unique(RewardProb);
% CuedPalette = ((1 - RewardProbCategories) * [128 128 128] + 127)/255;

%{
if p.TimelineView
    EarliestSessionDate = datetime(DataHolder{1}.Info.SessionDate);
    LatestSessionDate = datetime(DataHolder{end}.Info.SessionDate);
    FullDateRange = between(EarliestSessionDate, LatestSessionDate, 'Days');
    FullDateRangeChar = char(FullDateRange);
    ColourAdjustmentDenominator = str2double(FullDateRangeChar(1:end-1));
else
    ColourAdjustmentDenominator = size(DataHolder);
end
%}

%% Initiatize figure
% create figure
AnalysisFigure = figure('Position', [   0       0     595     420],... % DIN A4, 72 ppi
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

%% Analysis across trials
% Cue-sorted Move Time (MT) <-> Estimated reward rate
MTRewardRateAxes = axes(AnalysisFigure, 'Position', [0.01    0.58    0.18    0.30]);
hold(MTRewardRateAxes, 'on');

MTRewardRatePlot = {};

set(MTRewardRateAxes,...
    'TickDir', 'in',...
    'YLim', [-0.05 1],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(MTRewardRateAxes, 'Local reward rate')
ylabel(MTRewardRateAxes, 'Move Time (s)')

% Cue-sorted Move Time
CueSortedMTAxes = axes(AnalysisFigure, 'Position', [0.32    0.58    0.14    0.30]);
hold(CueSortedMTAxes, 'on');

set(CueSortedMTAxes,...
    'TickDir', 'in',...
    'XLim', [0 1],...
    'YLim', [-0.05 1],...
    'FontSize', 10);
xlabel(CueSortedMTAxes, 'Reward prob')

% Predicted MT
PredictedMTAxes = axes(AnalysisFigure, 'Position', [0.50    0.58    0.14    0.30]);
hold(PredictedMTAxes, 'on');

set(PredictedMTAxes,...
    'TickDir', 'in',...
    'XLim', [0 1],...
    'YLim', [-0.05 1],...
    'FontSize', 10);
xlabel(PredictedMTAxes, 'Reward prob')

% Move Time GLM Coeff.
MTGLMCoeffAxes = axes(AnalysisFigure, 'Position', [0.72    0.58    0.26    0.30]);
hold(MTGLMCoeffAxes, 'on');

set(MTGLMCoeffAxes,...
    'TickDir', 'in',...
    'YLim', [-0.05 1],...
    'YTick', [0 1],...
    'FontSize', 10);
% xlabel(MTGLMCoeffAxes, 'Regressor')
ylabel(MTGLMCoeffAxes, 'Coeff. \beta')

% Cue-sorted Time Investment (NotBaited Waiting Time) (TI) <-> Estimated reward rate
TIRewardRateAxes = axes(AnalysisFigure, 'Position', [0.01    0.12    0.18    0.30]);
hold(TIRewardRateAxes, 'on');

TIRewardRatePlot = {};

set(TIRewardRateAxes,...
    'TickDir', 'in',...
    'YLim', [-2 12],...
    'YTick', [0 5 10],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(TIRewardRateAxes, 'Local reward rate')
ylabel(TIRewardRateAxes, 'Invested Time (s)')

% Cue-sorted Time Investment (NotBaited Waiting Time)
CueSortedTIAxes = axes(AnalysisFigure, 'Position', [0.32    0.12    0.14    0.30]);
hold(CueSortedTIAxes, 'on');

set(CueSortedTIAxes,...
    'TickDir', 'in',...
    'XLim', [0 1],...
    'YLim', [-2 12],...
    'YTick', [0 5 10],...
    'FontSize', 10);
xlabel(CueSortedTIAxes, 'Reward prob')

% Cue-sorted Time Investment (NotBaited Waiting Time)
PredictedTIAxes = axes(AnalysisFigure, 'Position', [0.50    0.12    0.14    0.30]);
hold(PredictedTIAxes, 'on');

set(PredictedTIAxes,...
    'TickDir', 'in',...
    'XLim', [0 1],...
    'YLim', [-2 12],...
    'YTick', [0 5 10],...
    'FontSize', 10);
xlabel(PredictedTIAxes, 'Reward prob')

% Time Investment GLM Coeff.
TIGLMCoeffAxes = axes(AnalysisFigure, 'Position', [0.72    0.12    0.26    0.30]);
hold(TIGLMCoeffAxes, 'on');

set(TIGLMCoeffAxes,...
    'TickDir', 'in',...
    'YLim', [-2 12],...
    'YTick', [0 5 10],...
    'FontSize', 10);
xlabel(TIGLMCoeffAxes, 'Regressor')
ylabel(TIGLMCoeffAxes, 'Coeff. \beta')

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')

%% Aggregate data
SessionDateLabel = [];

AllMTChoiceLeft = [];
AllMT = [];
% AllMTSqrtZScore = []; % sqrt z-score is not a good estimation based on distribution
AllMTTrialRewardProb = [];
AllMTRewardedHistory = [];
AllMTConsumedWaterPercentage = [];
AllMTSessionNumber = [];

AllTIChoiceLeft = [];
AllTI = [];
AllTISqrtZScore = [];
AllTITrialRewardProb = [];
AllTIRewardedHistory = [];
AllTIConsumedWaterPercentage = [];
AllTISessionNumber = [];

for iSession = 1:length(DataHolder)
    % Import SessionData
    SessionData = DataHolder{iSession};
    SessionDateLabel = [SessionDateLabel, string(datestr(datetime(SessionData.Info.SessionDate), 'YYYYmmDD(ddd)'))];
    
    nTrials = SessionData.nTrials;
    if nTrials < 200
        disp(['Session ', num2str(iSession), ' has nTrial < 200. Impossible for analysis.'])
        continue
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
    
    RewardMagnitude = SessionData.Custom.TrialData.RewardMagnitude(:, 1:nTrials);
    TrialStartTimeStamp = SessionData.TrialStartTimestamp;
    TrialEndTimeStamp = SessionData.TrialEndTimestamp;

    idxTrial = 1:nTrials;
    SessionColor = [1, 1, 1] * 0.85; % ([1, 1, 1] - iSession / length(DataHolder)) * 0.9;
    
    %% insert data
    Chosen = ~isnan(ChoiceLeft);
    AllMTChoiceLeft = [AllMTChoiceLeft, ChoiceLeft(Chosen)];
    AllMT = [AllMT, MoveTime(Chosen)];
    
    % NotBaited Choice
    NotBaited = any(~Baited .* ChoiceLeftRight, 1) & (IncorrectChoice ~= 1);
    AllTIChoiceLeft = [AllTIChoiceLeft, ChoiceLeft(NotBaited)];
    AllTI = [AllTI, FeedbackWaitingTime(NotBaited)];
    
    % zscore(sqrt(TI))
    InvestedTimeSqrtZScore = zeros(size(ChoiceLeft));
    LeftNotBaitedIdx = find(ChoiceLeft == 1 & NotBaited == 1);
    InvestedTimeSqrtZScore(LeftNotBaitedIdx) = zscore(sqrt(FeedbackWaitingTime(LeftNotBaitedIdx)));
    RightNotBaitedIdx = find(ChoiceLeft == 0 & NotBaited == 1);
    InvestedTimeSqrtZScore(RightNotBaitedIdx) = zscore(sqrt(FeedbackWaitingTime(RightNotBaitedIdx)));

    AllTISqrtZScore = [AllTISqrtZScore, InvestedTimeSqrtZScore(NotBaited)];
    
    % trial reward prob
    if any(isnan(LightLeft)) % usually 2-arm task
        TrialRewardProb = max(RewardProb .* ChoiceLeftRight, [], 1);
    else % usually 1-arm task
        TrialRewardProb = max(RewardProb .* LightLeftRight, [], 1);
    end
    AllMTTrialRewardProb = [AllMTTrialRewardProb, TrialRewardProb(Chosen)];
    AllTITrialRewardProb = [AllTITrialRewardProb, TrialRewardProb(NotBaited)];

    % rewarded history (estimated with tau = 45s)
    RewardMagnitude = SessionData.Custom.TrialData.RewardMagnitude(:, 1:nTrials);
    RewardedMagnitude = sum(RewardMagnitude .* ChoiceLeftRight) .* Rewarded;
    RewardedMagnitude(isnan(RewardedMagnitude)) = 0;
    
    TrialStartTimestamp = SessionData.TrialStartTimestamp(:, 1:nTrials) - SessionData.TrialStartTimestamp(1);
    TrialTimeDuration = [0 diff(TrialStartTimestamp)];
    Tau = 5:5:100;
    
    RewardedHistory = 0;
    for iTrial = 1:nTrials-1
        RewardedHistory(iTrial+1) = RewardedHistory(iTrial) * exp(-TrialTimeDuration(iTrial+1)/Tau(9)) + RewardedMagnitude(iTrial);
    end
    
    RewardedHistory = RewardedHistory ./ SessionData.SettingsFile.GUI.RewardAmount; 
    
    AllMTRewardedHistory = [AllMTRewardedHistory, RewardedHistory(Chosen)];
    AllTIRewardedHistory = [AllTIRewardedHistory, RewardedHistory(NotBaited)];

    % water consumed
    ConsumedWaterPercentage = cumsum(RewardedMagnitude)./sum(RewardedMagnitude);
    AllMTConsumedWaterPercentage = [AllMTConsumedWaterPercentage, ConsumedWaterPercentage(Chosen)];
    AllTIConsumedWaterPercentage = [AllTIConsumedWaterPercentage, ConsumedWaterPercentage(NotBaited)];

    % session number
    AllMTSessionNumber = [AllMTSessionNumber, iSession * ones(size(find(Chosen)))];
    AllTISessionNumber = [AllTISessionNumber, iSession * ones(size(find(NotBaited)))];
    
    %% single-session data
    RewardProbCategories = unique(RewardProb);
    CuedPalette = ((1 - RewardProbCategories) * [128 128 128] + [96 96 96])/255;
    for iRewardProb = 1:length(RewardProbCategories)
        TargetRewardProb = RewardProbCategories(iRewardProb);
        Colour = CuedPalette(iRewardProb, :);
        
        % MT-Reward Rate
        ValidIdx = find(TrialRewardProb == TargetRewardProb & Chosen);
        XData = RewardedHistory(ValidIdx);
        YData = MoveTime(ValidIdx);

        [Coef, ~] = polyfit(XData, YData, 1);
        XVal = linspace(min(XData), max(XData), 3);
        YVal = XVal * Coef(1) + Coef(2);

        MTRewardRatePlot{end + 1} = plot(MTRewardRateAxes,...
                                         XVal,...
                                         YVal,...
                                         'LineStyle', '-',...
                                         'LineWidth', 0.1,...
                                         'Color', Colour);
        
        % TI-Reward Rate
        ValidIdx = find(TrialRewardProb == TargetRewardProb & NotBaited);
        XData = RewardedHistory(ValidIdx);
        YData = FeedbackWaitingTime(ValidIdx);

        [Coef, ~] = polyfit(XData, YData, 1);
        XVal = linspace(min(XData), max(XData), 3);
        YVal = XVal * Coef(1) + Coef(2);

        TIRewardRatePlot{end + 1} = plot(TIRewardRateAxes,...
                                         XVal,...
                                         YVal,...
                                         'LineStyle', '-',...
                                         'LineWidth', 0.1,...
                                         'Color', Colour);
    end
end

%% Pooled data analysis
RewardProbCategories = unique(AllMTTrialRewardProb');
CuedPalette = ((1 - RewardProbCategories) * [128 128 128])/255;
for iRewardProb = 1:length(RewardProbCategories)
    TargetRewardProb = RewardProbCategories(iRewardProb);
    Colour = CuedPalette(iRewardProb, :);
    
    % MT-Reward Rate
    ValidIdx = find(AllMTTrialRewardProb == TargetRewardProb);
    XData = AllMTRewardedHistory(ValidIdx);
    YData = AllMT(ValidIdx);

    [Coef, ~] = polyfit(XData, YData, 1);
    XVal = linspace(min(XData), max(XData), 3);
    YVal = XVal * Coef(1) + Coef(2);

    MTRewardRatePlot{end + 1} = plot(MTRewardRateAxes,...
                                     XVal,...
                                     YVal,...
                                     'LineStyle', '-',...
                                     'LineWidth', 2,...
                                     'Color', Colour);
    
    % TI-Reward Rate
    ValidIdx = find(AllTITrialRewardProb == TargetRewardProb);
    XData = AllTIRewardedHistory(ValidIdx);
    YData = AllTI(ValidIdx);

    [Coef, ~] = polyfit(XData, YData, 1);
    XVal = linspace(min(XData), max(XData), 3);
    YVal = XVal * Coef(1) + Coef(2);

    TIRewardRatePlot{end + 1} = plot(TIRewardRateAxes,...
                                     XVal,...
                                     YVal,...
                                     'LineStyle', '-',...
                                     'LineWidth', 2,...
                                     'Color', Colour);
end

% Left-Right Cue-Sorted MT
LeftMTSwarmchart = swarmchart(CueSortedMTAxes,...
                              AllMTTrialRewardProb(AllMTChoiceLeft == 1) - 0.05,...
                              AllMT(AllMTChoiceLeft == 1),...
                              'SizeData', 5,...
                              'Marker', '.',...
                              'MarkerEdgeColor', LRPalette(1,:),...
                              'XJitter', 'density',...
                              'XJitterWidth', 0.1);

RightMTSwarmchart = swarmchart(CueSortedMTAxes,...
                               AllMTTrialRewardProb(AllMTChoiceLeft == 0) + 0.05,...
                               AllMT(AllMTChoiceLeft == 0),...
                               'SizeData', 5,...
                               'Marker', '.',...
                               'MarkerEdgeColor', LRPalette(2,:),...
                               'XJitter', 'density',...
                               'XJitterWidth', 0.1);

% Left-Right Cue-Sorted TI
LeftTISwarmchart = swarmchart(CueSortedTIAxes,...
                              AllTITrialRewardProb(AllTIChoiceLeft == 1) - 0.05,...
                              AllTI(AllTIChoiceLeft == 1),...
                              'SizeData', 5,...
                              'Marker', '.',...
                              'MarkerEdgeColor', LRPalette(1,:),...
                              'XJitter', 'density',...
                              'XJitterWidth', 0.1);

RightTISwarmchart = swarmchart(CueSortedTIAxes,...
                               AllTITrialRewardProb(AllTIChoiceLeft == 0) + 0.05,...
                               AllTI(AllTIChoiceLeft == 0),...
                               'SizeData', 5,...
                               'Marker', '.',...
                               'MarkerEdgeColor', LRPalette(2,:),...
                               'XJitter', 'density',...
                               'XJitterWidth', 0.1);

% MT-GLM
X = [AllMTTrialRewardProb' == RewardProbCategories(1),...
     AllMTTrialRewardProb' == RewardProbCategories(end),...
     AllMTChoiceLeft' == 1,...
     AllMTChoiceLeft' == 0,...
     AllMTRewardedHistory',...
     AllMTConsumedWaterPercentage'];
LeeMTGLM = fitglm(X, AllMT,...
                  'y ~ x1:x3 + x2:x3 + x1:x4 + x2:x4 + x1:x5 + x2:x5 + x1:x6 + x2:x6',...
                  'Intercept', 'false'); 

MTGLMCoeffBar = bar(MTGLMCoeffAxes, LeeMTGLM.Coefficients.Estimate, 'w');

MTGLMCoeffErrorbar = errorbar(MTGLMCoeffAxes,...
                              LeeMTGLM.Coefficients.Estimate,...
                              LeeMTGLM.Coefficients.SE,...
                              'LineStyle', 'none',...
                              'Color', 'k');

XLabel = {'x_1*x_3', 'x_2*x_3',...
          'x_1*x_4', 'x_2*x_4',...
          'x_1*x_5', 'x_2*x_5',...
          'x_1*x_6', 'x_2*x_6'};

SignificantLevel = 0.05;
SignificanceIdx = find(LeeMTGLM.Coefficients.pValue < SignificantLevel);
CoeffSignificanceLine = line(MTGLMCoeffAxes,...
                             SignificanceIdx,...
                             zeros(size(SignificanceIdx)),...
                             'LineStyle', 'none',...
                             'Color', 'k',...
                             'Marker', '*');

MTGLMCoeffText = text(MTGLMCoeffAxes, 1, 0.8,...
                      sprintf("x_1 = P(r)_{low}, x_2 = P(r)_{high}\n" +...
                              "x_3 = c_L, x_4 = c_R\n" +...
                              "x_5 = RewardRate\n" +...
                              "x_6 = %%WaterIntake"),...
                      'FontSize', 8);

set(MTGLMCoeffAxes,...
    'XTick', 1:height(LeeMTGLM.Coefficients),...
    'XTickLabel', XLabel,...
    'XTickLabelRotation', 90)

% TI-GLM
X = [AllTITrialRewardProb' == RewardProbCategories(1),...
     AllTITrialRewardProb' == RewardProbCategories(end),...
     AllTIChoiceLeft' == 1,...
     AllTIChoiceLeft' == 0,...
     AllTIRewardedHistory',...
     AllTIConsumedWaterPercentage'];
LeeTIGLM = fitglm(X, AllTI,...
                  'y ~ x1:x3 + x2:x3 + x1:x4 + x2:x4 + x1:x5 + x2:x5 + x1:x6 + x2:x6',...
                  'Intercept', 'false'); 

TIGLMCoeffBar = bar(TIGLMCoeffAxes, LeeTIGLM.Coefficients.Estimate, 'w');

TIGLMCoeffErrorbar = errorbar(TIGLMCoeffAxes,...
                              LeeTIGLM.Coefficients.Estimate,...
                              LeeTIGLM.Coefficients.SE,...
                              'LineStyle', 'none',...
                              'Color', 'k');

SignificantLevel = 0.05;
SignificanceIdx = find(LeeTIGLM.Coefficients.pValue < SignificantLevel);
CoeffSignificanceLine = line(TIGLMCoeffAxes,...
                             SignificanceIdx,...
                             zeros(size(SignificanceIdx)),...
                             'LineStyle', 'none',...
                             'Color', 'k',...
                             'Marker', '*');

set(TIGLMCoeffAxes,...
    'XTick', 1:height(LeeTIGLM.Coefficients))

% Predicted MT
PredictedMT = LeeMTGLM.Fitted.Response;

LeftPredictedMTSwarmchart = swarmchart(PredictedMTAxes,...
                                       AllMTTrialRewardProb(AllMTChoiceLeft == 1) - 0.05,...
                                       PredictedMT(AllMTChoiceLeft == 1),...
                                       'SizeData', 5,...
                                       'Marker', '.',...
                                       'MarkerEdgeColor', LRPalette(1,:),...
                                       'XJitter', 'density',...
                                       'XJitterWidth', 0.1);

RightPredictedMTSwarmchart = swarmchart(PredictedMTAxes,...
                                        AllMTTrialRewardProb(AllMTChoiceLeft == 0) + 0.05,...
                                        PredictedMT(AllMTChoiceLeft == 0),...
                                        'SizeData', 5,...
                                        'Marker', '.',...
                                        'MarkerEdgeColor', LRPalette(2,:),...
                                        'XJitter', 'density',...
                                        'XJitterWidth', 0.1);

% Left-Right Cue-Sorted TI
PredictedTI = LeeTIGLM.Fitted.Response;

LeftPredictedTISwarmchart = swarmchart(PredictedTIAxes,...
                                       AllTITrialRewardProb(AllTIChoiceLeft == 1) - 0.05,...
                                       PredictedTI(AllTIChoiceLeft == 1),...
                                       'SizeData', 5,...
                                       'Marker', '.',...
                                       'MarkerEdgeColor', LRPalette(1,:),...
                                       'XJitter', 'density',...
                                       'XJitterWidth', 0.1);

RightPredictedTISwarmchart = swarmchart(PredictedTIAxes,...
                                        AllTITrialRewardProb(AllTIChoiceLeft == 0) + 0.05,...
                                        PredictedTI(AllTIChoiceLeft == 0),...
                                        'SizeData', 5,...
                                        'Marker', '.',...
                                        'MarkerEdgeColor', LRPalette(2,:),...
                                        'XJitter', 'density',...
                                        'XJitterWidth', 0.1);

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')

DataPath = strcat(DataFolderPath, '\', FigureTitle, '.png');
exportgraphics(AnalysisFigure, DataPath);

end % function