function AnalysisFigure = Matching_SS_SwitchStayGLM(DataFile)
% Matching Analysis Function
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ May 2025

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

AnalysisName = 'Matching_SS_SwitchStayGLM';

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

%% pre-filtered trials
PreFilteredTrial = ~isnan(ChoiceLeft);
ChoiceLeft = ChoiceLeft(PreFilteredTrial);
Baited = Baited(:, PreFilteredTrial);
IncorrectChoice = IncorrectChoice(PreFilteredTrial);
Rewarded = Rewarded(PreFilteredTrial);
FeedbackWaitingTime = FeedbackWaitingTime(PreFilteredTrial);

RewardProb = RewardProb(:, PreFilteredTrial);
ChoiceLeftRight = ChoiceLeftRight(:, PreFilteredTrial); 

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

%% Block switching behaviour across session
BlockSwitchAxes = axes(AnalysisFigure, 'Position', [0.01    0.82    0.37    0.11]);
hold(BlockSwitchAxes, 'on');
if ~isempty(ChoiceLeft) && ~all(isnan(ChoiceLeft))
    idxTrial = find(PreFilteredTrial == 1); % 1:nTrials;
    RewardProbLeft = RewardProb(1,:);
    RewardProbRight = RewardProb(2,:);
    RewardProbLeftPlot = plot(BlockSwitchAxes, idxTrial, RewardProbLeft * 100,...
                              'LineStyle', '-',...
                              'Marker', 'none',...
                              'Color', ColourPalette.Left,...
                              'LineWidth', 1);
    RewardProbRightPlot = plot(BlockSwitchAxes, idxTrial, RewardProbRight * 100,...
                               'LineStyle', '-',...
                               'Marker', 'none',...
                               'Color', ColourPalette.Right,...
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
    xlabel('iTrial')
    ylabel('Left Choices (%)')
    % title('Block switching behaviour')
end

%% Switch-stay-model & Predicted Choice
% in theory should just be the same as Lau Glimcher
try 
    Switched = abs(ChoiceLeft(2:end) - ChoiceLeft(1:end-1)) == 1; % need to exclude nan ~= 0 or nan ~= 1
    Stayed = abs(ChoiceLeft(2:end) - ChoiceLeft(1:end-1)) == 0;

    YData = nan(size(Switched));
    YData(Stayed) = 1;
    YData(Switched) = 0;
    YData = [nan, YData]; % first trial is unknown (philosophically explored)

    % build trial history kernels (n=5) last 5 trials
    Rewards = Rewarded == 1;
    Unrewards = Rewarded == 0;

    HistoryKernelSize = 5;
    SameChoice = zeros(length(YData), HistoryKernelSize);
    SameChoiceRewarded = zeros(length(YData), HistoryKernelSize);
    SameChoiceUnrewarded = zeros(length(YData), HistoryKernelSize);
    for j = 1:HistoryKernelSize
        Switched = [false(1, j), abs(ChoiceLeft(j:end-1) - ChoiceLeft(1:end-j)) == 1];
        Stayed = [false(1, j), abs(ChoiceLeft(j:end-1) - ChoiceLeft(1:end-j)) == 0];
        SameChoice(Stayed, j) = 1;
        SameChoice(Switched, j) = -1;
        
        SameChoiceRewarded(j+1:end, j) = Rewards(1:end-j) .* SameChoice(j+1:end, j)';
        SameChoiceUnrewarded(j+1:end, j) = Unrewards(1:end-j) .* SameChoice(j+1:end, j)';
        
    end
    
    % concatenate to build design matrix X
    XData = SameChoiceRewarded; % SameChoiceUnrewarded does not separately modulate, i.e. no counterfactual
    
    SwitchStayGLM = fitglm(XData, YData', 'distribution', 'binomial');
    
    % predict choices
    PredictedStayedProb = SwitchStayGLM.Fitted.Response;
    LogStayedOdds = SwitchStayGLM.Fitted.LinearPredictor;   %logodds for both: left and right
    
    PredictedStayed = double(PredictedStayedProb>=0.5);
    PredictedChoiceLeft = nan(size(ChoiceLeft));
    PredictedLeft = [false, abs(ChoiceLeft(1:end-1) - PredictedStayed(2:end)') == 0];
    PredictedRight = [false, abs(ChoiceLeft(1:end-1) - PredictedStayed(2:end)') == 1];
    PredictedChoiceLeft(PredictedLeft) = 1;
    PredictedChoiceLeft(PredictedRight) = 0;

    SmoothedPredictedChoiceLeft = smooth(PredictedChoiceLeft , BinWidth, 'moving','omitnan');
    SmoothedPredictedChoicePlot = plot(BlockSwitchAxes, idxTrial, SmoothedPredictedChoiceLeft * 100, '-r', 'LineWidth', 0.2);
    %{
    ExploringTrial = find(abs(YData - PredictedStayedProb') >= 0.5);
    ExploitingTrial = find(abs(YData - PredictedStayedProb') < 0.5);
    PredictedExplorationChoicePlot = plot(BlockSwitchAxes, ExploringTrial, 110 - PredictedStayed(ExploringTrial) * 120,...
                                          'Color', 'none',...
                                          'Marker', '.',...
                                          'MarkerEdgeColor', carrot,...
                                          'MarkerSize', 5);

    PredictedExploitationChoicePlot = plot(BlockSwitchAxes, ExploitingTrial, PredictedStayed(ExploitingTrial) * 130 - 15,...
                                           'Color', 'none',...
                                           'Marker', '.',...
                                           'MarkerEdgeColor', violet,...
                                           'MarkerSize', 5);
    
    LegendString = {'$\textsf{P(r)}_\textsf{L}$', '$\textsf{P(r)}_\textsf{R}$', '$\textsf{c}_\textsf{smooth}$',...
                    '$\hat\textsf{c}_\textsf{smooth}$', '$\hat\textsf{c}_\textsf{explore}$', '$\hat\textsf{c}_\textsf{exploit}$'};
    
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode')
    
    BlockSwitchLegend = legend(BlockSwitchAxes, LegendString,...
                               'Position', [0.01    0.70    0.37    0.06],...
                               'NumColumns', 3);
    set(BlockSwitchLegend,...
        'Interpreter', 'latex');
    %}
    model = true;
        
catch
    disp('error in running model');
    model = false;
    
end

if model
    %% psychometric
    PsychometricAxes = axes(AnalysisFigure, 'Position', [0.23    0.56    0.15    0.11]);
    hold(PsychometricAxes, 'on')
    
    set(PsychometricAxes,...
        'FontSize', 10,...
        'XLim', [-5 5],...
        'YLim', [0, 100],...
        'YAxisLocation', 'right')
    title(PsychometricAxes, 'Psychometric')
    xlabel(PsychometricAxes, 'log(odds_{stayed})')
    ylabel(PsychometricAxes, 'Stayed (%)')
    
    % Stayed Psychometric
    ValidTrial = ~isnan(YData); % and EarlyWithdrawal is always 0
    ValidLogOdds = LogStayedOdds(ValidTrial);
    ValidStayed = YData(ValidTrial);
    dvbin = linspace(-max(abs(ValidLogOdds)), max(abs(ValidLogOdds)), 10);
    [xdata, ydata, error] = BinData(ValidLogOdds, ValidStayed, dvbin);
    vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);
    
    StayedPsychometricErrorBar = errorbar(PsychometricAxes, xdata(vv), ydata(vv)*100, error(vv)*100,...
                                          'LineStyle', 'none',...
                                          'LineWidth', 1.5,...
                                          'Color', 'k',...
                                          'Marker', 'o',...
                                          'MarkerEdgeColor', 'k');
    
    PsychometricGLM = fitglm(ValidLogOdds, ValidStayed(:), 'Distribution', 'binomial');
    PsychometricGLMPlot = plot(PsychometricAxes, xdata, predict(PsychometricGLM, xdata)*100, '-', 'Color', [.5,.5,.5], 'LineWidth', 0.5);
    
    %% Coefficient of Stay-switch GLM
    ModelCoefficientAxes = axes(AnalysisFigure, 'Position', [0.04    0.56    0.15    0.11]);
    hold(ModelCoefficientAxes, 'on');
    
    set(ModelCoefficientAxes, 'FontSize', 10)
    xlabel(ModelCoefficientAxes, 'iTrial back');
    ylabel(ModelCoefficientAxes, 'Coeff.');
    title(ModelCoefficientAxes, 'GLM Fitted Coefficients')
    
    xdata = 1:HistoryKernelSize;
    ydataRewarded = SwitchStayGLM.Coefficients.Estimate(2:1+HistoryKernelSize);
    intercept = SwitchStayGLM.Coefficients.Estimate(1);

    SameChoiceRewardedCoeffPlot = plot(ModelCoefficientAxes, xdata, ydataRewarded', '-k');
    InterceptPlot = plot(ModelCoefficientAxes, xdata, intercept.*ones(size(xdata)), '--k');
    
    ModelCoefficientLegend = legend(ModelCoefficientAxes, {'Rewarded (same/diff from last =±1)', 'Intercept'},...
                                    'Position', [0.15    0.62    0.12    0.05],...
                                    'NumColumns', 1,...
                                    'Box', 'off');
    
    %% Residual Histogram
    ResidualHistogramAxes = axes(AnalysisFigure, 'Position', [0.04    0.36    0.15    0.11]);
    ResidualHistogram = plotResiduals(SwitchStayGLM, 'Histogram');
    
    %% Residual Histogram
    ResidualLaggedAxes = axes(AnalysisFigure, 'Position', [0.23    0.36    0.15    0.11]);
    ResidualLagged = plotResiduals(SwitchStayGLM, 'lagged', 'Marker', '.', 'MarkerSize', 1);
    
    %% Residual Histogram
    ResidualFittedAxes = axes(AnalysisFigure, 'Position', [0.04    0.19    0.15    0.11]);
    ResidualFitted = plotResiduals(SwitchStayGLM, 'fitted', 'Marker', '.', 'MarkerSize', 1);
    
    %% Residual Histogram
    ResidualProbabilityAxes = axes(AnalysisFigure, 'Position', [0.23    0.19    0.15    0.11]);
    ResidualProbability = plotResiduals(SwitchStayGLM, 'Probability', 'Marker', '.', 'MarkerSize', 1);
    
    %% Residual trend
    ResidualTrendAxes = axes(AnalysisFigure, 'Position', [0.04    0.06    0.34    0.09]);
    hold(ResidualTrendAxes, 'on')

    AbsModelResiduals = abs(SwitchStayGLM.Residuals.Raw);

    ResidualTrendPlot = plot(ResidualTrendAxes, idxTrial, AbsModelResiduals,...
                             'Color', 'k',...
                             'LineStyle', 'none',...
                             'Marker', '.',...
                             'MarkerSize', 1);
    
    MovMeanResidual = movmean(AbsModelResiduals, 10, 'omitnan');
    ResidualMovMeanPlot = plot(ResidualTrendAxes, idxTrial, MovMeanResidual,...
                               'Color', 'k',...
                               'LineStyle', '-',...
                               'Marker', 'none');
    
    MovStdResidual = movstd(AbsModelResiduals, 10, 'omitnan');
    ResidualMovStdPlot = plot(ResidualTrendAxes, idxTrial,  MovMeanResidual + MovStdResidual .* [1, -1],...
                              'Color', [0.5, 0.5, 0.5],...
                              'LineStyle', '-',...
                              'Marker', 'none');
    
    set(ResidualTrendAxes,...
        'Box', 'off',...
        'TickDir', 'out',...
        'FontSize', 10,...
        'YLim', [0, 1],...
        'YAxisLocation', 'right')
    xlabel(ResidualTrendAxes, 'iTrial')
    ylabel(ResidualTrendAxes, 'Abs(Residuals)')

    if ~all(isnan(FeedbackWaitingTime))
        %% Time Investment (TI) (only NotBaited Waiting Time) across session 
        TrialTIAxes = axes(AnalysisFigure, 'Position', [0.45    0.82    0.37    0.11]);
        hold(TrialTIAxes, 'on');
        
        NotBaited = any(~Baited .* ChoiceLeftRight, 1) & (IncorrectChoice ~= 1);
        
        StayedTITrial = NotBaited & Stayed;
        SwitchedTITrial = NotBaited & ~Stayed;
        StayedTI = FeedbackWaitingTime(StayedTITrial);
        SwitchedTI = FeedbackWaitingTime(SwitchedTITrial);
        
        % NotBaited invested time per explore/exploit across session
        TrialStayedTIPlot = plot(TrialTIAxes, idxTrial(StayedTITrial), StayedTI,...
                                 'Marker', '.',...
                                 'MarkerSize', 4,...
                                 'MarkerEdgeColor', ColourPalette.Explore,... % incorrect colour for now
                                 'Color', 'none');
        
        TrialSwitchedTIPlot = plot(TrialTIAxes, idxTrial(SwitchedTITrial), SwitchedTI,...
                                   'Marker', '.',...
                                   'MarkerSize', 4,...
                                   'MarkerEdgeColor', ColourPalette.Exploit,... % incorrect colour for now
                                   'Color', 'none');
        
        VevaiometryLegend = legend(TrialTIAxes, {'Stayed', 'Switched'},...
                                   'Box', 'off',...
                                   'Position', [0.75    0.90    0.05    0.03]);

        set(TrialTIAxes,...
            'TickDir', 'out',...
            'XLim', BlockSwitchAxes.XLim,...
            'YLim', [0, max(1, SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5)],...
            'YAxisLocation', 'right',...
            'FontSize', 10);
        ylabel('Invested Time (s)')
        % title('Block switching behaviour')
        
        %% plot vevaiometric      
        VevaiometricAxes = axes(AnalysisFigure, 'Position', [0.89    0.82    0.10    0.11]);
        hold(VevaiometricAxes, 'on')
        
        StayedLogStayedOdds = LogStayedOdds(StayedTITrial);
        SwitchedLogStayedOdds = LogStayedOdds(SwitchedTITrial);
        
        StayedTrialTIScatter = scatter(VevaiometricAxes, StayedLogStayedOdds, StayedTI,...
                                       'Marker', '.',...
                                       'MarkerEdgeColor', ColourPalette.Explore,...
                                       'SizeData', 18);
        %{
        SwitchedTrialTIScatter = scatter(VevaiometricAxes, SwitchedLogStayedOdds, SwitchedTI,...
                                         'Marker', '.',...
                                         'MarkerEdgeColor', ColourPalette.Exploit,...
                                         'SizeData', 18);
        %}
        [ExploreLineXData, ExploreLineYData] = Binvevaio(StayedLogStayedOdds, StayedTI, 10);
        % [ExploitLineXData, ExploitLineYData] = Binvevaio(SwitchedLogStayedOdds, SwitchedTI, 10);
        
        ExplorePlot = plot(VevaiometricAxes, ExploreLineXData, ExploreLineYData,...
                           'Color', ColourPalette.Explore,...
                           'LineWidth', 1);       
        %{
        ExploitPlot = plot(VevaiometricAxes, ExploitLineXData, ExploitLineYData,...
                           'Color', ColourPalette.Exploit,...
                           'LineWidth', 1);
        %}
        set(VevaiometricAxes,...
            'FontSize', 10,...
            'XLim', [-5 5],...
            'YLim', [0 SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5])
        title(VevaiometricAxes, 'Vevaiometric');
        xlabel(VevaiometricAxes, 'log(odds_{stayed})');
        % ylabel(VevaiometricAxes, 'Invested Time (s)');
        %{
        %% Psychometrics of NotBaited Choice with High- and Low-time investment (TI)
        % Time investment is limited to NotBaited trials
        TISortedPsychometricAxes = axes(AnalysisFigure, 'Position', [0.45    0.56    0.15    0.11]);
        hold(TISortedPsychometricAxes, 'on')
        
        TI = FeedbackWaitingTime(NotBaited);
        TImed = median(TI, "omitnan");
        HighTIStayedTrial = FeedbackWaitingTime>TImed & NotBaited & Stayed == 1;
        LowTIStayedTrial = FeedbackWaitingTime<=TImed & NotBaited & Stayed == 1;
        
        [xdata, ydata, error] = BinData(LogStayedOdds(HighTIStayedTrial), ChoiceLeft(HighTIStayedTrial), dvbin);
        vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);
    
        HighTIErrorBar = errorbar(TISortedPsychometricAxes, xdata(vv), ydata(vv)*100, error(vv)*100,...
                                  'LineStyle', 'none',...
                                  'LineWidth', 1,...
                                  'Marker', 'o',...
                                  'MarkerFaceColor', 'none',...
                                  'MarkerEdgeColor', 'k',...
                                  'Color', 'k');
    
        [xdata, ydata, error] = BinData(LogStayedOdds(LowTIStayedTrial), ChoiceLeft(LowTIStayedTrial), dvbin);
        vv = ~isnan(xdata) & ~isnan(ydata) & ~isnan(error);
    
        LowTIErrorBar = errorbar(TISortedPsychometricAxes, xdata(vv), ydata(vv)*100, error(vv)*100,...
                                 'LineStyle', 'none',...
                                 'LineWidth', 1,...
                                 'Marker', 'o',...
                                 'MarkerFaceColor', 'none',...
                                 'MarkerEdgeColor', [0.5 0.5 0.5],...
                                 'Color', [0.5 0.5 0.5]);
        
        HighTIGLM = fitglm(LogStayedOdds(HighTIStayedTrial), ChoiceLeft(HighTIStayedTrial), 'Distribution', 'binomial');
        HighTIGLMPlot = plot(TISortedPsychometricAxes, xdata, predict(HighTIGLM, xdata)*100,...
                             'Marker', 'none',...
                             'Color', 'k',...
                             'LineWidth', 0.5);
    
        LowTIGLM = fitglm(LogStayedOdds(LowTIStayedTrial), ChoiceLeft(LowTIStayedTrial), 'Distribution', 'binomial');
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
        xlabel('log(odds_{stayed})')
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
        %}
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
                               'MarkerEdgeColor', ColourPalette.Left);
        
        TrialRightTIPlot = plot(LRTrialTIAxes, idxTrial(RightTITrial), RightTI,...
                                'LineStyle', 'none',...
                                'Marker', '.',...
                                'MarkerSize', 4,...
                                'MarkerEdgeColor', ColourPalette.Right);
        
        LRVevaiometryLegend = legend(LRTrialTIAxes, {'Left', 'Right'},...
                                     'Box', 'off',...
                                     'Position', [0.75    0.44    0.05    0.03]);

        set(LRTrialTIAxes,...
            'TickDir', 'out',...
            'XLim', BlockSwitchAxes.XLim,...
            'YLim', [0, max(1, SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5)],...
            'YAxisLocation', 'right',...
            'FontSize', 10);
        ylabel('Invested Time (s)')
        % title('Block switching behaviour')
        
        %% plot vevaiometric (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
        LRTIVevaiometricAxes = axes(AnalysisFigure, 'Position', [0.89    0.36    0.10    0.11]);
        hold(LRTIVevaiometricAxes, 'on')
        %{
        AbsModelResiduals = abs(SwitchStayGLM.Residuals.Raw);
        
        LeftAbsResidual = AbsModelResiduals(LeftTITrial);
        RightAbsResidual = AbsModelResiduals(RightTITrial);
        %}
        LeftStayedTITrial = NotBaited & Stayed & ChoiceLeft == 1;
        LeftStayedTI = FeedbackWaitingTime(LeftStayedTITrial);
        LeftLogStayedOdds = LogStayedOdds(LeftStayedTITrial);

        LeftStayedTIScatter = scatter(LRTIVevaiometricAxes, LeftLogStayedOdds, LeftStayedTI,...
                                      'Marker', '.',...
                                      'MarkerEdgeColor', ColourPalette.Left,...
                                      'SizeData', 18);
        
        RightStayedTITrial = NotBaited & Stayed & ChoiceLeft == 0;
        RightStayedTI = FeedbackWaitingTime(RightStayedTITrial);
        RightLogStayedOdds = LogStayedOdds(RightStayedTITrial);
        
        RightStayedTIScatter = scatter(LRTIVevaiometricAxes, RightLogStayedOdds, RightStayedTI,...
                                       'Marker', '.',...
                                       'MarkerEdgeColor', ColourPalette.Right,...
                                       'SizeData', 18);
    
        [LeftLineXData, LeftLineYData] = Binvevaio(LeftLogStayedOdds, LeftStayedTI, 10);
        [RightLineXData, RightLineYData] = Binvevaio(RightLogStayedOdds, RightStayedTI, 10);
        
        LeftTIPlot = plot(LRTIVevaiometricAxes, LeftLineXData, LeftLineYData,...
                          'Color', ColourPalette.Left,...
                          'LineWidth', 1);       
        
        RightTIPlot = plot(LRTIVevaiometricAxes, RightLineXData, RightLineYData,...
                           'Color', ColourPalette.Right,...
                           'LineWidth', 1);
        
        set(LRTIVevaiometricAxes,...
            'FontSize', 10,...
            'XLim', [-5 5],...
            'YLim', [0 SessionData.SettingsFile.GUI.FeedbackDelayMax * 1.5])
        title(LRTIVevaiometricAxes, 'LRVevaiometric');
        xlabel(LRTIVevaiometricAxes, 'log(odds_{stayed})');
    end
    %{
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
                           'MarkerEdgeColor', ColourPalette.Left);
    
    TrialRightMTPlot = plot(LRTrialMTAxes, idxTrial(ChoiceLeft==0), RightMT,...
                            'LineStyle', 'none',...
                            'Marker', '.',...
                            'MarkerSize', 4,...
                            'MarkerEdgeColor', ColourPalette.Right);
    
    % LRMTVevaiometryLegend = legend(LRTrialMTAxes, {'Left', 'Right'},...
    %                                'Box', 'off',...
    %                                'Position', [0.75    0.35    0.05    0.03]);

    set(LRTrialMTAxes,...
        'TickDir', 'out',...
        'XLim', BlockSwitchAxes.XLim,...
        'YLim', [0, 0.5],...
        'YAxisLocation', 'right',...
        'FontSize', 10);
    ylabel('Move Time (s)')
    % title('Block switching behaviour')
    
    %% plot vevaiometric (L/R sorted residuals = abs(ChoiceLeft - P({ChoiceLeft}^))
    LRMTVevaiometricAxes = axes(AnalysisFigure, 'Position', [0.89    0.19    0.10    0.11]);
    hold(LRMTVevaiometricAxes, 'on')
    
    AbsModelResiduals = abs(SwitchStayGLM.Residuals.Raw);

    LeftAbsResidual = AbsModelResiduals(ChoiceLeft==1);
    RightAbsResidual = AbsModelResiduals(ChoiceLeft==0);
    
    LeftMTScatter = scatter(LRMTVevaiometricAxes, LeftAbsResidual, LeftMT,...
                            'Marker', '.',...
                            'MarkerEdgeColor', ColourPalette.Left,...
                            'SizeData', 18);
    
    RightMTScatter = scatter(LRMTVevaiometricAxes, RightAbsResidual, RightMT,...
                             'Marker', '.',...
                             'MarkerEdgeColor', ColourPalette.Right,...
                             'SizeData', 18);

    [LeftLineXData, LeftLineYData] = Binvevaio(LeftAbsResidual, LeftMT, 10);
    [RightLineXData, RightLineYData] = Binvevaio(RightAbsResidual, RightMT, 10);
    
    LeftMTPlot = plot(LRMTVevaiometricAxes, LeftLineXData, LeftLineYData,...
                    'Color', ColourPalette.Left,...
                    'LineWidth', 1);       
    
    RightMTPlot = plot(LRMTVevaiometricAxes, RightLineXData, RightLineYData,...
                     'Color', ColourPalette.Right,...
                     'LineWidth', 1);
    
    set(LRMTVevaiometricAxes,...
        'FontSize', 10,...
        'XLim', [0 1],...
        'YLim', [0 0.5])
    % title(LRMTVevaiometricAxes, 'LRMTVevaiometric');
    xlabel(LRMTVevaiometricAxes, 'abs(Residuals)');
    %}
end

DataFolder = OttLabDataServerFolderPath;
RatName = SessionData.Info.Subject;
% %%The following lines doesn't not work, as the timestamp documented
% in the SessionData may not be the same as the one being used for saving
% SessionDate = string(datetime(SessionData.Info.SessionDate), 'yyyyMMdd')';
% SessionTime = string(datetime(SessionData.Info.SessionStartTime_UTC), 'HHmmSS')';
% SessionDateTime = strcat(SessionDate, '_', SessionTime);
DataPath = strcat(DataFolder, RatName, '\bpod_session\', SessionDateTime, '\',...
                  RatName, '_TwoArmBanditVariant_', SessionDateTime, '_Matching_LauGlimcherGLM.png');
exportgraphics(AnalysisFigure, DataPath);

DataPath = strcat(DataFolder, RatName, '\bpod_graph\',...
                  RatName, '_TwoArmBanditVariant_', SessionDateTime, '_Matching_LauGlimcherGLM.png');
exportgraphics(AnalysisFigure, DataPath);

close(AnalysisFigure)
%}
end % function