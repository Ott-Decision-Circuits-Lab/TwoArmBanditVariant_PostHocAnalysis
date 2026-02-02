function Matching_MS_LogOddsAnalysis(DataFolderPath)
%{
MS = MultiSession
First create on 20260127 by Antonio Lee for AG Ott @HU Berlin

to check if there are additional variables that affect on top of the
standard Lau-Glimcher log(odds) 
%}

%% load files
if nargin < 1
    DataFolderPath = uigetdir(OttLabDataServerFolderPath());
elseif ~ischar(DataFolderPath) && ~isstring(DataFolderPath)
    disp('Error: Unknown input format. No further analysis can be performed.')
    return
end

try
    load(fullfile(DataFolderPath, '\Selected_Data.mat'));
    % load(fullfile(DataFolderPath, '\Concatenated_Data.mat'));
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

AnalysisName = 'Matching_MS_LogOddsAnalysis';

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

%% Initiatize figure
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

%% Analysis across trials
SessionDateLabel = [];
nSessions = length(DataHolder);

%% Invested time -> LogOdds
% does Invested time influence StaySwitch
TIStaySwitchAxes = axes(AnalysisFigure, 'Position', [0.01, 0.76, 0.14, 0.18]);
hold(TIStaySwitchAxes, 'on');

AllValidTIStaySwitch = [];
AllValidInvestedTime = [];

set(TIStaySwitchAxes,...
    'TickDir', 'out',...
    'YLim', [0 100],...
    'YTick', [0 50 100],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(TIStaySwitchAxes, 'Time investment (s)')
% ylabel(TIStaySwitchAxes, 'Stayed (%)')

% does TI capture variance on top of log(odds)
TIGLMCoeffAxes = axes(AnalysisFigure, 'Position', [0.01, 0.52, 0.14, 0.18]);
hold(TIGLMCoeffAxes, 'on');

AllTIXData = [];
AllTIYData = [];

set(TIGLMCoeffAxes,...
    'TickDir', 'out',...
    'XLim', [0, 4],...
    'XTick', [1, 2, 3],...
    'XTickLabel', {'\beta_0', 'log(odds)', 'Time investment_{i-1}'},...
    'XTickLabelRotation', 90,...
    'YAxisLocation', 'right',...
    'FontSize', 10);
% ylabel(TIGLMCoeffAxes, '\beta')

%% RewardWaitingTime -> StaySwitch
% does RewardWaitingTime influence stay/switch
RWTStaySwitchAxes = axes(AnalysisFigure, 'Position', [0.18, 0.76, 0.14, 0.18]);
hold(RWTStaySwitchAxes, 'on');

AllValidRWTStaySwitch = [];
AllValidRewardWaitingTime = [];

set(RWTStaySwitchAxes,...
    'TickDir', 'out',...
    'YLim', [0 100],...
    'YTick', [0 50 100],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(RWTStaySwitchAxes, 'Reward waiting time (s)')
% ylabel(RWTStaySwitchAxes, 'Stayed (%)')

% does RewardWaitingTime capture variance on top of log(odds)
RWTGLMCoeffAxes = axes(AnalysisFigure, 'Position', [0.18, 0.52, 0.14, 0.18]);
hold(RWTGLMCoeffAxes, 'on');

AllRWTXData = [];
AllRWTYData = [];

set(RWTGLMCoeffAxes,...
    'TickDir', 'out',...
    'XLim', [0, 4],...
    'XTick', [1, 2, 3],...
    'XTickLabel', {'\beta_0', 'log(odds)', 'Reward waiting time_{i-1}'},...
    'XTickLabelRotation', 90,...
    'YAxisLocation', 'right',...
    'FontSize', 10);
% ylabel(RWTGLMCoeffAxes, '\beta')

%% RewardRate -> StaySwitch
% does RewardRate influence stay/switch
RRStaySwitchAxes = axes(AnalysisFigure, 'Position', [0.35, 0.76, 0.14, 0.18]);
hold(RRStaySwitchAxes, 'on');

AllValidRRStaySwitch = [];
AllValidRewardRate = [];

set(RRStaySwitchAxes,...
    'TickDir', 'out',...
    'YLim', [0 100],...
    'YTick', [0 50 100],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
xlabel(RRStaySwitchAxes, 'Reward rate')
ylabel(RRStaySwitchAxes, 'Stayed (%)')

% does RewardRate capture variance on top of log(odds)
RRGLMCoeffAxes = axes(AnalysisFigure, 'Position', [0.35, 0.52, 0.14, 0.18]);
hold(RRGLMCoeffAxes, 'on');

AllRRXData = [];
AllRRYData = [];

set(RRGLMCoeffAxes,...
    'TickDir', 'out',...
    'XLim', [0, 4],...
    'XTick', [1, 2, 3],...
    'XTickLabel', {'\beta_0', 'log(odds)', 'Reward rate'},...
    'XTickLabelRotation', 90,...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(RRGLMCoeffAxes, '\beta')

%% Plotting
for iSession = 1:length(DataHolder)
    % Import SessionData
    SessionData = DataHolder{iSession};
    SessionDateLabel = [SessionDateLabel, string(datestr(datetime(SessionData.Info.SessionDate), 'YYYYmmDD(ddd)'))];
    
    nTrials = SessionData.nTrials;
    if nTrials < 250
        disp(['Session ', num2str(iSession), ' has nTrial < 250. Impossible for analysis.'])
        continue
    end
    idxTrial = 1:nTrials;

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
    RewardProbCategories = unique(RewardProb);
    
    LightLeft = SessionData.Custom.TrialData.LightLeft(1:nTrials);
    LightLeftRight = [LightLeft; 1-LightLeft];
    TrialRewardProb = sum(RewardProb .* LightLeftRight, 1);
    
    ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft]; 
    ChoiceRewardProb = sum(RewardProb .* ChoiceLeftRight, 1);
    
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

    %% TI -> Stay/Switch
    % does TI influence stay/switch
    NotBaited = any(~Baited .* ChoiceLeftRight, 1) & (IncorrectChoice ~= 1);
    ValidNotBaited = [NotBaited(1:end-1), 0];
    
    Switched = [abs(ChoiceLeft(1:end-1) - ChoiceLeft(2:end)) == 1, false];
    Stayed = [abs(ChoiceLeft(1:end-1) - ChoiceLeft(2:end)) == 0, false];
    StaySwitch = nan(size(ChoiceLeft));
    StaySwitch(Switched == 1) = 0;
    StaySwitch(Stayed == 1) = 1;

    Valid = ValidNotBaited == 1 & ~isnan(StaySwitch);
    ValidInvestedTime = FeedbackWaitingTime(Valid);
    ValidTIStaySwitch = StaySwitch(Valid);

    Bin = linspace(0, 12, 6);
    [XData, YData, Error] = BinData(ValidInvestedTime, ValidTIStaySwitch, Bin);
    ValidData = ~isnan(XData) & ~isnan(YData) & ~isnan(Error);
    
    TIStaySwitchLine{iSession} = line(TIStaySwitchAxes,...
                                      'xdata', XData(ValidData),...
                                      'ydata', YData(ValidData) * 100,...
                                      'LineStyle', '-',...
                                      'LineWidth', 0.5,...
                                      'Color', ColourPalette.Session);
    
    AllValidTIStaySwitch = [AllValidTIStaySwitch, ValidTIStaySwitch];
    AllValidInvestedTime = [AllValidInvestedTime, ValidInvestedTime];

    % does TI capture variance on top of log(odds)
    Choices = 2 * ChoiceLeft - 1;
    
    Rewards = Rewarded == 1;
    Rewards = Rewards .* Choices;

    for iHistoryKernel = 1:5
        Choices(iHistoryKernel + 1, :) = [zeros(1, iHistoryKernel), Choices(1, 1:end - iHistoryKernel)];
        Rewards(iHistoryKernel + 1, :) = [zeros(1, iHistoryKernel), Rewards(1, 1:end - iHistoryKernel)];
    end

    Choices(isnan(Choices)) = 0;
    Rewards(isnan(Rewards)) = 0;

    XData = [Choices(2:end, :)', Rewards(2:end, :)'];
    YData = ChoiceLeft';
    LauGlimcherGLM = fitglm(XData, YData, 'distribution', 'binomial');
    LogOdds = LauGlimcherGLM.Fitted.LinearPredictor';
    
    % step-wise GLM
    LastFeedbackWaitingTime = [0, FeedbackWaitingTime(1:end-1)];
    LastFeedbackWaitingTime = LastFeedbackWaitingTime .* Choices(2, :);
    LastNotBaited = [false, NotBaited(1:end-1)];

    XData = [LogOdds(LastNotBaited)', LastFeedbackWaitingTime(LastNotBaited)'];
    YData = ChoiceLeft(LastNotBaited)';
    TIGLM = fitglm(XData, YData, 'distribution', 'binomial');
    TIGLMCoeff(:, iSession) = TIGLM.Coefficients.Estimate;

    AllTIXData = [AllTIXData; XData];
    AllTIYData = [AllTIYData; YData];
    
    %% RewardWaitingTime -> Stay/Switch
    % does RewardWaitingTime influence stay/switch
    ValidRewarded = [Rewarded(1:end-1), 0];
    
    Valid = ValidRewarded == 1 & ~isnan(StaySwitch);
    ValidRewardWaitingTime = FeedbackWaitingTime(Valid);
    ValidRWTStaySwitch = StaySwitch(Valid);

    Bin = linspace(0.5, 8, 5);
    [XData, YData, Error] = BinData(ValidRewardWaitingTime, ValidRWTStaySwitch, Bin);
    ValidData = ~isnan(XData) & ~isnan(YData) & ~isnan(Error);
    
    RWTStaySwitchLine{iSession} = line(RWTStaySwitchAxes,...
                                       'xdata', XData(ValidData),...
                                       'ydata', YData(ValidData) * 100,...
                                       'LineStyle', '-',...
                                       'LineWidth', 0.5,...
                                       'Color', ColourPalette.Session);
    
    AllValidRWTStaySwitch = [AllValidRWTStaySwitch, ValidRWTStaySwitch];
    AllValidRewardWaitingTime = [AllValidRewardWaitingTime, ValidRewardWaitingTime];

    % does RewardWaitingTime capture variance on top of log(odds)
    LastFeedbackWaitingTime = [0, FeedbackWaitingTime(1:end-1)];
    LastFeedbackWaitingTime = LastFeedbackWaitingTime .* Choices(2, :);
    LastRewarded = [false, Rewarded(1:end-1) == 1];

    XData = [LogOdds(LastRewarded)', LastFeedbackWaitingTime(LastRewarded)'];
    YData = ChoiceLeft(LastRewarded)';
    RWTGLM = fitglm(XData, YData, 'distribution', 'binomial');
    
    RWTGLMCoeff(:, iSession) = RWTGLM.Coefficients.Estimate;
    AllRWTXData = [AllRWTXData; XData];
    AllRWTYData = [AllRWTYData; YData];
    
    %% RewardRate -> Stay/Switch
    % does RewardRate influence stay/switch
    TrialDuration = diff(TrialStartTimeStamp);
    RewardRate = 0;
    for iTrial = 2:nTrials
        RewardRate(iTrial) = RewardRate(iTrial - 1)...
                                 * exp(-TrialDuration(iTrial - 1) /30)...
                                 + (Rewarded(iTrial - 1) == 1);
    end
    RewardRate = [RewardRate(2:end), 0];

    Valid = ~isnan(StaySwitch);
    ValidRewardRate = RewardRate(Valid);
    ValidRRStaySwitch = StaySwitch(Valid);

    Bin = linspace(0, 3, 10);
    [XData, YData, Error] = BinData(ValidRewardRate, ValidRRStaySwitch, Bin);
    ValidData = ~isnan(XData) & ~isnan(YData) & ~isnan(Error);
    
    RRStaySwitchLine{iSession} = line(RRStaySwitchAxes,...
                                      'xdata', XData(ValidData),...
                                      'ydata', YData(ValidData) * 100,...
                                      'LineStyle', '-',...
                                      'LineWidth', 0.5,...
                                      'Color', ColourPalette.Session);
    
    AllValidRRStaySwitch = [AllValidRRStaySwitch, ValidRRStaySwitch];
    AllValidRewardRate = [AllValidRewardRate, ValidRewardRate];

    % does RewardRate capture variance on top of log(odds)
    % RewardRate here is assigned with the sign of the current log(odds),
    % such that a +ve weight means it further supports model-conforming,
    % while -ve means it opposes the model
    RewardRate = [0, RewardRate(1:end-1)];
    RewardRate = RewardRate .* sign(LogOdds);

    XData = [LogOdds', RewardRate'];
    YData = ChoiceLeft';
    RRGLM = fitglm(XData, YData, 'distribution', 'binomial');
    
    RRGLMCoeff(:, iSession) = RRGLM.Coefficients.Estimate;
    AllRRXData = [AllRRXData; XData];
    AllRRYData = [AllRRYData; YData];
    
end

%% TI -> Stay/Switch
% does TI influence stay/switch
Bin = linspace(0, 12, 6);
[XData, YData, Error] = BinData(AllValidInvestedTime, AllValidTIStaySwitch, Bin);
ValidData = ~isnan(XData) & ~isnan(YData) & ~isnan(Error);

TIStaySwitchLine{iSession + 1} = line(TIStaySwitchAxes,...
                                      'xdata', XData(ValidData),...
                                      'ydata', YData(ValidData) * 100,...
                                      'LineStyle', '-',...
                                      'LineWidth', 0.5,...
                                      'Color', ColourPalette.Pooled);

% does TI capture variance on top of log(odds)
TIGLM = fitglm(AllTIXData, AllTIYData, 'distribution', 'binomial');

X = repmat(1:3, nSessions, 1);
X = reshape(X, 1, []);
TIGLMCoeff = reshape(TIGLMCoeff', 1, []);
SessionTIGLMCoeffSwarmchart = swarmchart(TIGLMCoeffAxes,...
                                         X,...
                                         TIGLMCoeff,...
                                         'Marker', '.',...
                                         'MarkerEdgeColor', ColourPalette.Session,...
                                         'XJitter', 'density',...
                                         'XJitterWidth', 0.4);

TICoeffBar = bar(TIGLMCoeffAxes, TIGLM.Coefficients.Estimate,...
                 'FaceColor', 'none');
TICoeffErrorbar = errorbar(TIGLMCoeffAxes,...
                           TIGLM.Coefficients.Estimate,...
                           TIGLM.Coefficients.SE,...
                           'LineStyle', 'none',...
                           'Color', 'k');

SignificantLevel = 0.05;
SignificanceIdx = find(TIGLM.Coefficients.pValue < SignificantLevel);
TISignificantLine = line(TIGLMCoeffAxes,...
                         SignificanceIdx,...
                         zeros(size(SignificanceIdx)),...
                         'LineStyle', 'none',...
                         'Color', 'k',...
                         'Marker', '*');

%% RWT -> Stay/Switch
% does RWT influence stay/switch
Bin = linspace(0.5, 8, 5);
[XData, YData, Error] = BinData(AllValidRewardWaitingTime, AllValidRWTStaySwitch, Bin);
ValidData = ~isnan(XData) & ~isnan(YData) & ~isnan(Error);

RWTStaySwitchLine{iSession + 1} = line(RWTStaySwitchAxes,...
                                       'xdata', XData(ValidData),...
                                       'ydata', YData(ValidData) * 100,...
                                       'LineStyle', '-',...
                                       'LineWidth', 0.5,...
                                       'Color', ColourPalette.Pooled);

% does TI capture variance on top of log(odds)
RWTGLM = fitglm(AllRWTXData, AllRWTYData, 'distribution', 'binomial');

X = repmat(1:3, nSessions, 1);
X = reshape(X, 1, []);
RWTGLMCoeff = reshape(RWTGLMCoeff', 1, []);
SessionRWTGLMCoeffSwarmchart = swarmchart(RWTGLMCoeffAxes,...
                                          X,...
                                          RWTGLMCoeff,...
                                          'Marker', '.',...
                                          'MarkerEdgeColor', ColourPalette.Session,...
                                          'XJitter', 'density',...
                                          'XJitterWidth', 0.4);

RWTCoeffBar = bar(RWTGLMCoeffAxes, RWTGLM.Coefficients.Estimate,...
                 'FaceColor', 'none');
RWTCoeffErrorbar = errorbar(RWTGLMCoeffAxes,...
                            RWTGLM.Coefficients.Estimate,...
                            RWTGLM.Coefficients.SE,...
                            'LineStyle', 'none',...
                            'Color', 'k');

SignificantLevel = 0.05;
SignificanceIdx = find(RWTGLM.Coefficients.pValue < SignificantLevel);
RWTSignificantLine = line(RWTGLMCoeffAxes,...
                          SignificanceIdx,...
                          zeros(size(SignificanceIdx)),...
                          'LineStyle', 'none',...
                          'Color', 'k',...
                          'Marker', '*');

%% RewardRate -> Stay/Switch
% does RewardRate influence stay/switch
Bin = linspace(0, 3, 10);
[XData, YData, Error] = BinData(AllValidRewardRate, AllValidRRStaySwitch, Bin);
ValidData = ~isnan(XData) & ~isnan(YData) & ~isnan(Error);

RRStaySwitchLine{iSession + 1} = line(RRStaySwitchAxes,...
                                      'xdata', XData(ValidData),...
                                      'ydata', YData(ValidData) * 100,...
                                      'LineStyle', '-',...
                                      'LineWidth', 0.5,...
                                      'Color', ColourPalette.Pooled);

% does TI influence stay/switch without mediation by log(odds)
RRGLM = fitglm(AllRRXData, AllRRYData', 'distribution', 'binomial');

X = repmat(1:3, nSessions, 1);
X = reshape(X, 1, []);
RRGLMCoeff = reshape(RRGLMCoeff', 1, []);
SessionRRGLMCoeffSwarmchart = swarmchart(RRGLMCoeffAxes,...
                                         X,...
                                         RRGLMCoeff,...
                                         'Marker', '.',...
                                         'MarkerEdgeColor', ColourPalette.Session,...
                                         'XJitter', 'density',...
                                         'XJitterWidth', 0.4);

RRCoeffBar = bar(RRGLMCoeffAxes, RRGLM.Coefficients.Estimate,...
                 'FaceColor', 'none');
RRCoeffErrorbar = errorbar(RRGLMCoeffAxes,...
                           RRGLM.Coefficients.Estimate,...
                           RRGLM.Coefficients.SE,...
                           'LineStyle', 'none',...
                           'Color', 'k');

SignificantLevel = 0.05;
SignificanceIdx = find(RRGLM.Coefficients.pValue < SignificantLevel);
RRSignificantLine = line(RRGLMCoeffAxes,...
                         SignificanceIdx,...
                         zeros(size(SignificanceIdx)),...
                         'LineStyle', 'none',...
                         'Color', 'k',...
                         'Marker', '*');

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')

DataPath = strcat(DataFolderPath, '\', FigureTitle, '.png');
exportgraphics(AnalysisFigure, DataPath);

DataPath = strcat(DataFolderPath, '\', FigureTitle, '.fig');
savefig(AnalysisFigure, DataPath);

end