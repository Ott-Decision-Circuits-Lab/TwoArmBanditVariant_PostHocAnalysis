function AnalysisFigure = TwoArmBanditVariant_Matching_MovingLauGlimcherGLM(DataFile)
% Matching Analysis Function
% Developed by Antonio Lee @ BCCN Berlin
% Version 2.0 ~ April 2024

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
AnalysisFigure = figure('Position', [   0    0 1191  842],... % DIN A3, 72 ppi (window will crop it to _ x 1024, same as disp resolution)
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateTime, '_Matching'),...
                        'MenuBar', 'none',...
                        'Resize', 'off');

FrameAxes = axes(AnalysisFigure, 'Position', [0 0 1 1]); % spacer for correct saving dimension
set(FrameAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'w',...
    'YColor', 'w')

%% Figure Info
FigureInfoAxes = axes(AnalysisFigure, 'Position', [0.01    0.96    0.48    0.01]);
set(FigureInfoAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'w',...
    'YColor', 'w')

FigureTitle = strcat(RatName, '_', SessionDateTime, '_Matching_MovingLauGlimcherGLM');

FigureTitleText = text(FigureInfoAxes, 0, 0,...
                       FigureTitle,...
                       'FontSize', 14,...
                       'FontWeight','bold',...
                       'Interpreter', 'none');

%% Block switching behaviour across session
BlockSwitchAxes = axes(AnalysisFigure, 'Position', [0.01    0.82    0.37    0.11]);
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
        'YAxisLocation', 'right',...
        'FontSize', 10);
    xlabel(BlockSwitchAxes, 'iTrial')
    ylabel(BlockSwitchAxes, 'Left Choices (%)')
    % title('Block switching behaviour')
end

%% Moving Lau-Glimcher model & Predicted Choice
try 
    % preparing the data for design matrix
    HistoryKernelSize = 5;

    Choices = [];
    Rewards = [];
    
    ChoiceLeftCoding = ChoiceLeft * 2 - 1;
    RewardedLeftCoding = Rewarded .* ChoiceLeftCoding;
    
    for i = 1:HistoryKernelSize
        Choices(:, i) = [zeros(i, 1); ChoiceLeftCoding(1:end-i)'];
        Rewards(:, i) = [zeros(i, 1); RewardedLeftCoding(1:end-i)'];
    end
    
    % concatenate to build design matrix X
    X = [Choices, Rewards];
    X(isnan(X)) = 0;
    
    LookingBack = 80;
    LookingFront = 80;
    Coefficients = [];
    for iWindow = 1:nTrials
        Window = max(1, iWindow - LookingBack):min(nTrials, iWindow + LookingFront);
        LauGlimcherGLM = fitglm(X(Window, :), ChoiceLeft(Window)',...
                                'distribution', 'binomial');
        Coefficients(:, iWindow) = LauGlimcherGLM.Coefficients.Estimate;
        
        % predict choices
        Index = min(iWindow, LookingBack+1);
        PredictedLeftChoiceProb(iWindow) = LauGlimcherGLM.Fitted.Response(Index);
        LogOdds(iWindow) = LauGlimcherGLM.Fitted.LinearPredictor(Index);   %logodds for both: left and right
        ModelResiduals(iWindow) = LauGlimcherGLM.Residuals.Raw(Index);
    end
    
    % for conservation of later calculations
    PredictedLeftChoiceProb = PredictedLeftChoiceProb';
    LogOdds = LogOdds';
    ModelResiduals = ModelResiduals';

    PredictedChoice = double(PredictedLeftChoiceProb>=0.5);
    PredictedChoice(isnan(ChoiceLeft)) = nan;
    SmoothedPredictedChoiceLeft = smooth(PredictedChoice, BinWidth, 'moving','omitnan');
    SmoothedPredictedChoicePlot = plot(BlockSwitchAxes, idxTrial, SmoothedPredictedChoiceLeft * 100, '-r', 'LineWidth', 0.2);
    
    ExploringTrial = find(abs(ChoiceLeft - PredictedLeftChoiceProb') >= 0.5);
    ExploitingTrial = find(abs(ChoiceLeft - PredictedLeftChoiceProb') < 0.5);
    PredictedExplorationChoicePlot = plot(BlockSwitchAxes, ExploringTrial, 110 - PredictedChoice(ExploringTrial) * 120,...
                                          'Color', 'none',...
                                          'Marker', '.',...
                                          'MarkerEdgeColor', carrot,...
                                          'MarkerSize', 5);

    PredictedExploitationChoicePlot = plot(BlockSwitchAxes, ExploitingTrial, PredictedChoice(ExploitingTrial) * 130 - 15,...
                                           'Color', 'none',...
                                           'Marker', '.',...
                                           'MarkerEdgeColor', violet,...
                                           'MarkerSize', 5);
    
    LegendString = {'P(r)_L', 'smoothed c', 'c_{explore}',...
                    % 'c_L', 'c_R',...
                    'P(r)_R', 'smoothed c hat', 'c_{exploit}'};
    BlockSwitchLegend = legend(BlockSwitchAxes, LegendString,...
                               'Position', [0.01    0.70    0.37    0.06],...
                               'NumColumns', 3);
    
catch
    disp('Error: fail to run model');
    return
end

%% psychometric
PsychometricAxes = axes(AnalysisFigure, 'Position', [0.23    0.56    0.15    0.11]);
hold(PsychometricAxes, 'on')

set(PsychometricAxes,...
    'FontSize', 10,...
    'XLim', [-5 5],...
    'YLim', [0, 100],...
    'YAxisLocation', 'right')
title(PsychometricAxes, 'Psychometric')
xlabel(PsychometricAxes, 'log(odds)')
ylabel(PsychometricAxes, 'Left Choices (%)')

% Choice Psychometric
ValidTrial = ~isnan(ChoiceLeft); % and EarlyWithdrawal is always 0
NotBaitedLogOdds = LogOdds(ValidTrial);
NotBaitedChoice = ChoiceLeft(ValidTrial);
dvbin = linspace(-5, 5, 10);
[xdata, ydata, error] = BinData(NotBaitedLogOdds, NotBaitedChoice, dvbin);
vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);

ChoicePsychometricErrorBar = errorbar(PsychometricAxes, xdata(vv), ydata(vv)*100, error(vv)*100,...
                                      'LineStyle', 'none',...
                                      'LineWidth', 1.5,...
                                      'Color', 'k',...
                                      'Marker', 'o',...
                                      'MarkerEdgeColor', 'k');

PsychometricGLM = fitglm(NotBaitedLogOdds, NotBaitedChoice(:), 'Distribution', 'binomial');
PsychometricGLMPlot = plot(PsychometricAxes, xdata, predict(PsychometricGLM, xdata)*100, '-', 'Color', [.5,.5,.5], 'LineWidth', 0.5);

%% Choice Coefficient of Moving Lau-Glimcher GLM
HistoryColor = [0; 0.2; 0.4; 0.6; 0.8] * [1, 1, 1];

ChoiceCoefficientAxes = axes(AnalysisFigure, 'Position', [0.01    0.42    0.37    0.11]);
hold(ChoiceCoefficientAxes, 'on');

ChoiceHistoryCoefficientPlot = plot(ChoiceCoefficientAxes, Coefficients(2:6, :)');

for iPlot = 1:length(ChoiceHistoryCoefficientPlot)
    ChoiceHistoryCoefficientPlot(iPlot).Color = HistoryColor(iPlot, :);
end

set(ChoiceCoefficientAxes,...
    'FontSize', 10,...
    'XLim', BlockSwitchAxes.XLim,...
    'XTickLabel', [],...
    'YLim', [-2.5 2.5],...
    'YAxisLocation', 'right')
ylabel(ChoiceCoefficientAxes, 'Choice Coeff.');

%% Reward Coefficient of Moving Lau-Glimcher GLM
RewardCoefficientAxes = axes(AnalysisFigure, 'Position', [0.01    0.29    0.37    0.11]);
hold(RewardCoefficientAxes, 'on');

RewardHistoryCoefficientPlot = plot(RewardCoefficientAxes, Coefficients(7:11, :)');

for iPlot = 1:length(RewardHistoryCoefficientPlot)
    RewardHistoryCoefficientPlot(iPlot).Color = HistoryColor(iPlot, :);
end

set(RewardCoefficientAxes,...
    'FontSize', 10,...
    'XLim', BlockSwitchAxes.XLim,...
    'XTickLabel', [],...
    'YLim', [-3 5],...
    'YAxisLocation', 'right')
ylabel(RewardCoefficientAxes, 'Reward Coeff.');

%% Intercept of Moving Lau-Glimcher GLM
InterceptAxes = axes(AnalysisFigure, 'Position', [0.01    0.21    0.37    0.06]);
hold(InterceptAxes, 'on');

IntecerptPlot = plot(InterceptAxes, Coefficients(1, :)');
IntecerptPlot.Color = HistoryColor(1, :);

ReferencePlot = plot(InterceptAxes, zeros(size(Coefficients(1, :)')), '--');
ReferencePlot.Color = HistoryColor(5, :);

set(InterceptAxes,...
    'FontSize', 10,...
    'XLim', BlockSwitchAxes.XLim,...
    'XTickLabel', [],...
    'YLim', [-1 1],...
    'YAxisLocation', 'right')
ylabel(InterceptAxes, 'Intercept');

%% Changes in angle between coefficients
AngleAxes = axes(AnalysisFigure, 'Position', [0.01    0.13    0.37    0.06]);
hold(AngleAxes, 'on');

for iAngle = 1:length(Coefficients)-1
    v1 = Coefficients(:, iAngle);
    v2 = Coefficients(:, iAngle + 1);
    CosAngles(iAngle) = dot(v1 / norm(v1), v2 / norm(v2));
end

CosAnglePlot = plot(AngleAxes, smooth(CosAngles));

set(AngleAxes,...
    'FontSize', 10,...
    'XLim', BlockSwitchAxes.XLim,...
    'XTickLabel', [],...
    'YLim', [0.99 1],...
    'YAxisLocation', 'right')
ylabel(AngleAxes, 'cos\theta');

%% Changes in angle between coefficients
MagnitudeAxes = axes(AnalysisFigure, 'Position', [0.01    0.05    0.37    0.06]);
hold(MagnitudeAxes, 'on');

for iMagnitude = 1:length(Coefficients)
    v1 = Coefficients(:, iMagnitude);
    Magnitudes(iMagnitude) = norm(v1);
end

MagnitudePlot = plot(MagnitudeAxes, Magnitudes);

set(MagnitudeAxes,...
    'FontSize', 10,...
    'XLim', BlockSwitchAxes.XLim,...
    'YLim', [0 10],...
    'YAxisLocation', 'right')
xlabel(MagnitudeAxes, 'iTrial');
ylabel(MagnitudeAxes, 'Magnitude');

%% Residual Histogram
ResidualHistogramAxes = axes(AnalysisFigure, 'Position', [0.04    0.56    0.15    0.11]);
ResidualHistogram = histogram(ResidualHistogramAxes, ModelResiduals, linspace(-1, 1, 11),...
                              'Normalization', 'pdf');
xlabel(ResidualHistogramAxes, 'Residuals')
ylabel(ResidualHistogramAxes, 'Probability density')
title(ResidualHistogramAxes, 'Histogram of residuals')


%{
%% Residual Histogram
ResidualLaggedAxes = axes(AnalysisFigure, 'Position', [0.23    0.36    0.15    0.11]);
ResidualLagged = plotResiduals(LauGlimcherGLM, 'lagged', 'Marker', '.', 'MarkerSize', 1);

%% Residual Histogram
ResidualFittedAxes = axes(AnalysisFigure, 'Position', [0.04    0.19    0.15    0.11]);
ResidualFitted = plotResiduals(LauGlimcherGLM, 'fitted', 'Marker', '.', 'MarkerSize', 1);

%% Residual Histogram
ResidualProbabilityAxes = axes(AnalysisFigure, 'Position', [0.23    0.19    0.15    0.11]);
ResidualProbability = plotResiduals(LauGlimcherGLM, 'Probability', 'Marker', '.', 'MarkerSize', 1);
%}

if ~all(isnan(FeedbackWaitingTime))
    %% Time Investment (TI) (only NotBaited Waiting Time) across session 
    TrialTIAxes = axes(AnalysisFigure, 'Position', [0.45    0.82    0.37    0.11]);
    hold(TrialTIAxes, 'on');
    
    NotBaited = any(~Baited .* ChoiceLeftRight, 1) & (IncorrectChoice ~= 1);
    Exploit = ChoiceLeft == (LogOdds'>0);
    
    ExploringTITrial = NotBaited & ~Exploit;
    ExploitingTITrial = NotBaited & Exploit;
    ExploringTI = FeedbackWaitingTime(ExploringTITrial);
    ExploitingTI = FeedbackWaitingTime(ExploitingTITrial);
    
    % NotBaited invested time per explore/exploit across session
    TrialExploringTIPlot = plot(TrialTIAxes, idxTrial(ExploringTITrial), ExploringTI,...
                                'Marker', '.',...
                                'MarkerSize', 4,...
                                'MarkerEdgeColor', carrot,...
                                'Color', 'none');
    
    TrialExploitingTIPlot = plot(TrialTIAxes, idxTrial(ExploitingTITrial), ExploitingTI,...
                                 'Marker', '.',...
                                 'MarkerSize', 4,...
                                 'MarkerEdgeColor', violet,...
                                 'Color', 'none');
    
    VevaiometryLegend = legend(TrialTIAxes, {'Explore', 'Exploit'},...
                               'Box', 'off',...
                               'Position', [0.75    0.90    0.05    0.03]);

    set(TrialTIAxes,...
        'TickDir', 'out',...
        'XLim', BlockSwitchAxes.XLim,...
        'YLim', [0, max(1, SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5)],...
        'YAxisLocation', 'right',...
        'FontSize', 10);
    ylabel(TrialTIAxes, 'Invested Time (s)')
    % title('Block switching behaviour')
    
    %% plot vevaiometric      
    VevaiometricAxes = axes(AnalysisFigure, 'Position', [0.89    0.82    0.10    0.11]);
    hold(VevaiometricAxes, 'on')
    
    ExploringLogOdds = LogOdds(ExploringTITrial);
    ExploitingLogOdds = LogOdds(ExploitingTITrial);
    
    ExploringTrialTIScatter = scatter(VevaiometricAxes, ExploringLogOdds, ExploringTI,...
                                      'Marker', '.',...
                                      'MarkerEdgeColor', carrot,...
                                      'SizeData', 18);
    
    ExploitingTrialTIScatter = scatter(VevaiometricAxes, ExploitingLogOdds, ExploitingTI,...
                                       'Marker', '.',...
                                       'MarkerEdgeColor', violet,...
                                       'SizeData', 18);

    [ExploreLineXData, ExploreLineYData] = Binvevaio(ExploringLogOdds, ExploringTI, 10);
    [ExploitLineXData, ExploitLineYData] = Binvevaio(ExploitingLogOdds, ExploitingTI, 10);
    
    ExplorePlot = plot(VevaiometricAxes, ExploreLineXData, ExploreLineYData,...
                       'Color', carrot,...
                       'LineWidth', 1);       
    
    ExploitPlot = plot(VevaiometricAxes, ExploitLineXData, ExploitLineYData,...
                       'Color', violet,...
                       'LineWidth', 1);
    
    set(VevaiometricAxes,...
        'FontSize', 10,...
        'XLim', [-5 5],...
        'YLim', [0 SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5])
    title(VevaiometricAxes, 'Vevaiometric');
    xlabel(VevaiometricAxes, 'log(odds)');
    % ylabel(VevaiometricAxes, 'Invested Time (s)');
    
    %% Psychometrics of NotBaited Choice with High- and Low-time investment (TI)
    % Time investment is limited to NotBaited trials
    TISortedPsychometricAxes = axes(AnalysisFigure, 'Position', [0.45    0.56    0.15    0.11]);
    hold(TISortedPsychometricAxes, 'on')
    
    TI = FeedbackWaitingTime(NotBaited);
    TImed = median(TI, "omitnan");
    HighTITrial = FeedbackWaitingTime>TImed & NotBaited;
    LowTITrial = FeedbackWaitingTime<=TImed & NotBaited;
    
    [xdata, ydata, error] = BinData(LogOdds(HighTITrial), ChoiceLeft(HighTITrial), dvbin);
    vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);

    HighTIErrorBar = errorbar(TISortedPsychometricAxes, xdata(vv), ydata(vv)*100, error(vv)*100,...
                              'LineStyle', 'none',...
                              'LineWidth', 1,...
                              'Marker', 'o',...
                              'MarkerFaceColor', 'none',...
                              'MarkerEdgeColor', 'k',...
                              'Color', 'k');

    [xdata, ydata, error] = BinData(LogOdds(LowTITrial), ChoiceLeft(LowTITrial), dvbin);
    vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);

    LowTIErrorBar = errorbar(TISortedPsychometricAxes, xdata(vv), ydata(vv)*100, error(vv)*100,...
                             'LineStyle', 'none',...
                             'LineWidth', 1,...
                             'Marker', 'o',...
                             'MarkerFaceColor', 'none',...
                             'MarkerEdgeColor', [0.5 0.5 0.5],...
                             'Color', [0.5 0.5 0.5]);
    
    HighTIGLM = fitglm(LogOdds(HighTITrial), ChoiceLeft(HighTITrial), 'Distribution', 'binomial');
    HighTIGLMPlot = plot(TISortedPsychometricAxes, xdata, predict(HighTIGLM, xdata)*100,...
                         'Marker', 'none',...
                         'Color', 'k',...
                         'LineWidth', 0.5);

    LowTIGLM = fitglm(LogOdds(LowTITrial), ChoiceLeft(LowTITrial), 'Distribution', 'binomial');
    LowTIGLMPlot = plot(TISortedPsychometricAxes, xdata, predict(LowTIGLM, xdata)*100,...
                        'Marker', 'none',...
                        'Color', [0.5, 0.5, 0.5],...
                        'LineWidth', 0.5);
    
    TIPsychometricLegend = legend(TISortedPsychometricAxes, {'High TI','Low TI'},...
                                  'Box', 'off',...
                                  'Position', [0.55    0.57    0.05    0.03]);
    
    set(TISortedPsychometricAxes,...
        'FontSize', 10,...
        'XLim', [-5 5],...
        'YLim', [0, 100])
    title(TISortedPsychometricAxes, 'TI Sorted Psychometric')
    xlabel(TISortedPsychometricAxes, 'log(odds)')
    % ylabel('Left Choices (%)')
    
    %% callibration plot
    CalibrationAxes = axes(AnalysisFigure, 'Position', [0.66    0.56    0.15    0.11]);
    hold(CalibrationAxes, 'on')
    
    Correct = Exploit(NotBaited); %'correct'
    edges = linspace(min(TI), max(TI), 8);
    
    [xdata, ydata, error] = BinData(TI, Correct, edges);
    vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);
    CalibrationErrorBar = errorbar(CalibrationAxes,...
                                   xdata(vv), ydata(vv)*100, error(vv),...
                                   'LineWidth', 2, ...
                                   'Color', 'k');
    
    set(CalibrationAxes,...
        'FontSize', 10,...
        'XLim', [0 SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5])
    % title(CalibrationHandle, 'Calibration');
    xlabel(CalibrationAxes, 'Invested Time (s)');
    ylabel(CalibrationAxes, 'Exploit Ratio (%)');
    
    %% Time Investment (TI) (only NotBaited Waiting Time) across session 
    LRTrialTIAxes = axes(AnalysisFigure, 'Position', [0.45    0.36    0.37    0.11]);
    hold(LRTrialTIAxes, 'on');
    
    % Smoothed NotBaited invested time per left/right across session
    LeftTITrial = NotBaited & ChoiceLeft==1;
    RightTITrial = NotBaited & ChoiceLeft==0;
    LeftTI = FeedbackWaitingTime(LeftTITrial);
    RightTI = FeedbackWaitingTime(RightTITrial);
    
    TrialLeftTIPlot = plot(LRTrialTIAxes, idxTrial(LeftTITrial), LeftTI,...
                           'LineStyle', 'none',...
                           'Marker', '.',...
                           'MarkerSize', 4,...
                           'MarkerEdgeColor', sand);
    
    TrialRightTIPlot = plot(LRTrialTIAxes, idxTrial(RightTITrial), RightTI,...
                            'LineStyle', 'none',...
                            'Marker', '.',...
                            'MarkerSize', 4,...
                            'MarkerEdgeColor', turquoise);
    
    LRVevaiometryLegend = legend(LRTrialTIAxes, {'Left', 'Right'},...
                                 'Box', 'off',...
                                 'Position', [0.75    0.44    0.05    0.03]);

    set(LRTrialTIAxes,...
        'TickDir', 'out',...
        'XLim', BlockSwitchAxes.XLim,...
        'YLim', [0, max(1, SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5)],...
        'YAxisLocation', 'right',...
        'FontSize', 10);
    ylabel(LRTrialTIAxes, 'Invested Time (s)')
    % title('Block switching behaviour')
    
    %% plot vevaiometric (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
    LRTIVevaiometricAxes = axes(AnalysisFigure, 'Position', [0.89    0.36    0.10    0.11]);
    hold(LRTIVevaiometricAxes, 'on')
    
    AbsModelResiduals = abs(ModelResiduals);
    
    LeftAbsResidual = AbsModelResiduals(LeftTITrial);
    RightAbsResidual = AbsModelResiduals(RightTITrial);
    
    LeftTIScatter = scatter(LRTIVevaiometricAxes, LeftAbsResidual, LeftTI,...
                            'Marker', '.',...
                            'MarkerEdgeColor', sand,...
                            'SizeData', 18);
    
    RightTIScatter = scatter(LRTIVevaiometricAxes, RightAbsResidual, RightTI,...
                             'Marker', '.',...
                             'MarkerEdgeColor', turquoise,...
                             'SizeData', 18);

    [LeftLineXData, LeftLineYData] = Binvevaio(LeftAbsResidual, LeftTI, 10);
    [RightLineXData, RightLineYData] = Binvevaio(RightAbsResidual, RightTI, 10);
    
    LeftTIPlot = plot(LRTIVevaiometricAxes, LeftLineXData, LeftLineYData,...
                      'Color', sand,...
                      'LineWidth', 1);       
    
    RightTIPlot = plot(LRTIVevaiometricAxes, RightLineXData, RightLineYData,...
                       'Color', turquoise,...
                       'LineWidth', 1);
    
    set(LRTIVevaiometricAxes,...
        'FontSize', 10,...
        'XLim', [0 1],...
        'YLim', [0 SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5])
    title(LRTIVevaiometricAxes, 'LRVevaiometric');
    xlabel(LRTIVevaiometricAxes, 'abs(Residuals)');
end

%% Move Time (MT) across session 
LRTrialMTAxes = axes(AnalysisFigure, 'Position', [0.45    0.19    0.37    0.11]);
hold(LRTrialMTAxes, 'on');

% Smoothed NotBaited invested time per left/right across session
LeftMT = MoveTime(ChoiceLeft==1);
RightMT = MoveTime(ChoiceLeft==0);

TrialLeftMTPlot = plot(LRTrialMTAxes, idxTrial(ChoiceLeft==1), LeftMT,...
                       'LineStyle', 'none',...
                       'Marker', '.',...
                       'MarkerSize', 4,...
                       'MarkerEdgeColor', sand);

TrialRightMTPlot = plot(LRTrialMTAxes, idxTrial(ChoiceLeft==0), RightMT,...
                        'LineStyle', 'none',...
                        'Marker', '.',...
                        'MarkerSize', 4,...
                        'MarkerEdgeColor', turquoise);

% LRMTVevaiometryLegend = legend(LRTrialMTAxes, {'Left', 'Right'},...
%                                'Box', 'off',...
%                                'Position', [0.75    0.35    0.05    0.03]);

set(LRTrialMTAxes,...
    'TickDir', 'out',...
    'XLim', BlockSwitchAxes.XLim,...
    'YLim', [0, 0.5],...
    'YAxisLocation', 'right',...
    'FontSize', 10);
ylabel(LRTrialMTAxes, 'Move Time (s)')
% title('Block switching behaviour')

%% plot vevaiometric (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
LRMTVevaiometricAxes = axes(AnalysisFigure, 'Position', [0.89    0.19    0.10    0.11]);
hold(LRMTVevaiometricAxes, 'on')

AbsModelResiduals = abs(ModelResiduals);

LeftAbsResidual = AbsModelResiduals(ChoiceLeft==1);
RightAbsResidual = AbsModelResiduals(ChoiceLeft==0);

LeftMTScatter = scatter(LRMTVevaiometricAxes, LeftAbsResidual, LeftMT,...
                        'Marker', '.',...
                        'MarkerEdgeColor', sand,...
                        'SizeData', 18);

RightMTScatter = scatter(LRMTVevaiometricAxes, RightAbsResidual, RightMT,...
                         'Marker', '.',...
                         'MarkerEdgeColor', turquoise,...
                         'SizeData', 18);

[LeftLineXData, LeftLineYData] = Binvevaio(LeftAbsResidual, LeftMT, 10);
[RightLineXData, RightLineYData] = Binvevaio(RightAbsResidual, RightMT, 10);

LeftMTPlot = plot(LRMTVevaiometricAxes, LeftLineXData, LeftLineYData,...
                  'Color', sand,...
                  'LineWidth', 1);       

RightMTPlot = plot(LRMTVevaiometricAxes, RightLineXData, RightLineYData,...
                   'Color', turquoise,...
                   'LineWidth', 1);

set(LRMTVevaiometricAxes,...
    'FontSize', 10,...
    'XLim', [0 1],...
    'YLim', [0 0.5])
% title(LRMTVevaiometricAxes, 'LRMTVevaiometric');
xlabel(LRMTVevaiometricAxes, 'abs(Residuals)');

DataFolder = OttLabDataServerFolderPath;
RatName = SessionData.Info.Subject;
% %%The following lines doesn't not work, as the timestamp documented
% in the SessionData may not be the same as the one being used for saving
% SessionDate = string(datetime(SessionData.Info.SessionDate), 'yyyyMMdd')';
% SessionTime = string(datetime(SessionData.Info.SessionStartTime_UTC), 'HHmmSS')';
% SessionDateTime = strcat(SessionDate, '_', SessionTime);
DataPath = strcat(DataFolder, RatName, '\bpod_session\', SessionDateTime, '\',...
                  RatName, '_TwoArmBanditVariant_', SessionDateTime, '_Matching_MovingLauGlimcherGLM.png');
exportgraphics(AnalysisFigure, DataPath);

DataPath = strcat(DataFolder, RatName, '\bpod_graph\',...
                  RatName, '_TwoArmBanditVariant_', SessionDateTime, '_Matching_MovingLauGlimcherGLM.png');
exportgraphics(AnalysisFigure, DataPath);

close(AnalysisFigure)

end % function