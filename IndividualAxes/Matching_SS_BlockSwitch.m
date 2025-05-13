function AnalysisFigure = Matching_SS_BlockSwitch(DataFile, FigureSize, AxeSize)
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

carrot = [230, 97, 0]/255; % explore
violet = [93, 58, 155]/255; % exploit

% colour palette for cues: (1- P(r)) * 128 + 127
% P(0) = white; P(1) = smoky gray
RewardProbCategories = unique(RewardProb);
CuedPalette = ((1 - RewardProbCategories) * [128 128 128] + 127)/255;

% create figure
if nargin < 2
    FigureSize = [   2,    5,    5,  1.5];
end
AnalysisFigure = figure('Position', FigureSize,...
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateTime, '_Matching'),...
                        'MenuBar', 'none',...
                        'Resize', 'off',...
                        'Unit', 'inch',...
                        'Color', 'none');

set(AnalysisFigure, 'Position', FigureSize, 'unit', 'inch');

FrameAxes = axes(AnalysisFigure, 'Position', [0 0 1 1]); % spacer for correct saving dimension
set(FrameAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'none',...
    'YColor', 'none',...
    'Color', 'none')

%% Block switching behaviour across session
if nargin < 3
    AxeSize = [ 0.8, 0.55, 4.05,  0.90];
end

BlockSwitchAxes = axes(AnalysisFigure,...
                       'Position', AxeSize,...
                       'Units', 'inches',...
                       'Color', 'none');
set(BlockSwitchAxes, 'Position', AxeSize, 'Units', 'inches');

hold(BlockSwitchAxes, 'on');
if ~isempty(ChoiceLeft) && ~all(isnan(ChoiceLeft))
    idxTrial = 1:nTrials;
    RewardProbLeft = RewardProb(1,:);
    RewardProbRight = RewardProb(2,:);
    RewardProbLeftPlot = plot(BlockSwitchAxes, idxTrial, RewardProbLeft * 100,...
                              'LineStyle', '-',...
                              'Marker', 'none',...
                              'Color', sand,...
                              'LineWidth', 1);
    RewardProbRightPlot = plot(BlockSwitchAxes, idxTrial, RewardProbRight * 100,...
                               'LineStyle', '-',...
                               'Marker', 'none',...
                               'Color', turquoise,...
                               'LineWidth', 1);

    BinWidth = 10;
    ChoiceLeftSmoothed = smooth(ChoiceLeft, BinWidth, 'moving','omitnan'); %current bin width: 10 trials
    BlockChoicePlot = plot(BlockSwitchAxes, idxTrial, ChoiceLeftSmoothed * 100,...
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
    
    set(BlockSwitchAxes,...
        'TickDir', 'out',...
        'YLim', [-20 120],...
        'YTick', [0 50 100],...
        'FontSize', 10);
    xlabel(BlockSwitchAxes, 'i^{th} trial')
    ylabel(BlockSwitchAxes, 'Probability (%)')
    % title('Block switching behaviour')
    
    LegendString = {'$\textsf{P(r)}_\textsf{L}$', '$\textsf{P(r)}_\textsf{R}$', '$\bar\textsf{c}\textsf{=L}$'};
    
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
    
    BlockSwitchLegend = legend(BlockSwitchAxes, LegendString,...
                               'Position', [   1,  1.3,  2.7,  0.2],...
                               'NumColumns', 3,...
                               'Units', 'inches',...
                               'Color', 'none',...
                               'Box', 'off');
    set(BlockSwitchLegend,...
        'Interpreter', 'latex',...
        'Position', [   1,  1.3,  2.7,  0.2],...
        'Units', 'inches');
end
saveas(AnalysisFigure, strcat(RatName, '_', SessionDateTime, '_Matching_BlockSwitch'), 'svg')
end