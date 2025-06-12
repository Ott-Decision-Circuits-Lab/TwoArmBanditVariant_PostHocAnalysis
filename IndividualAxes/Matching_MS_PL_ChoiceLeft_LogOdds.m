function AnalysisFigure = Matching_MS_PL_ChoiceLeft_LogOdds(DataFolderPath, ModelsFilePath, FigureSize, AxeSize)
% for making an example of the task structure
% figure 'unit' set as 'inch' so that we know exactly the meters to pixels
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ April 2025

if nargin < 1
    DataFolderPath = uigetdir(OttLabDataServerFolderPath()); % end without '\'
    DataFolderPath = [DataFolderPath, '\'];
elseif ~ischar(DataFolderPath) && ~isstring(DataFolderPath)
    disp('Error: Unknown input format. No further analysis can be performed.')
    return
end

try
    load(fullfile(DataFolderPath, '\Selected_Data.mat')); % should be 'DataHolder' as variable name
    % load(fullfile(DataFolderPath, '\Concatenated_Data.mat'));
catch
    disp('Error: Selected DataFolderPath does not contain the required .mat for further steps.')
    return
end

SessionDateRange = DataFolderPath(end-17:end-1);
[~, RatName] = fileparts(fileparts(fileparts(fileparts(DataFolderPath))));

RatID = str2double(RatName);
if isnan(RatID)
    RatID = -1;
end
RatName = num2str(RatID);

AnalysisName = 'Matching_MS_PL_ChoiceLeft_LogOdds';

%% Bayesian Symmetric Q-Learning with Forgetting and Stickiness model
if nargin < 2
    [ModelsFileName, ModelsFilePath] = uigetfile(DataFolderPath);
    
    if ~strcmpi(ModelsFilePath, DataFolderPath)
        disp('Error: Models unlikely for the SessionDataset. Make sure they are in the same directory')
        return
    end

    ModelsFilePath = fullfile(ModelsFilePath, ModelsFileName);
end

try
    load(ModelsFilePath) % should be 'Models' as variable name
catch
    disp('Error: Problem in loading ModelsFilePath.')
    return
end

if exist('Models', 'var') == 0
    disp('Error: Loaded file is not a Models.')
    return
end

%% Load related data to local variabels
TimeInvestment = [];
TimeInvestmentZScore = [];
TimeInvestmentSqrtZScore = [];
Choices = [];

LogOdds = [];
ChosenValue = [];
UnchosenValue = [];
ChosenMemory = [];
LeftValue = [];
RightValue = [];
ChosenLogOdds = [];

for iSession = 1:length(DataHolder)
    SessionData = DataHolder{iSession};
    Model = Models{iSession};

    nTrials = SessionData.nTrials;
    TrialData = SessionData.Custom.TrialData;
    ChoiceLeft = TrialData.ChoiceLeft(1:nTrials);
    ChoiceLeftRight = [ChoiceLeft; 1 - ChoiceLeft];
    IncorrectChoice = TrialData.IncorrectChoice(1:nTrials);
    Baited = TrialData.Baited(:, 1:nTrials);
    NotBaited = any(~Baited .* ChoiceLeftRight, 1) & (IncorrectChoice ~= 1);
    FeedbackWaitingTime = TrialData.FeedbackWaitingTime(1:nTrials);
    Rewarded = TrialData.Rewarded(1:nTrials);
    
    if isfield(Model, 'Chains') % i.e. Bayesian
        Chain = vertcat(Model.Chains{:});
        [ProbDensity, Values] = ksdensity(Chain(:, 1));
        LearningRateMAPs = Values(ProbDensity == max(ProbDensity));
        [ProbDensity, Values] = ksdensity(Chain(:, 2));
        InverseTemperatureMAPs = Values(ProbDensity == max(ProbDensity));
        [ProbDensity, Values] = ksdensity(Chain(:, 3));
        ForgettingRateMAPs = Values(ProbDensity == max(ProbDensity));
        [ProbDensity, Values] = ksdensity(Chain(:, 4));
        ChoiceStickinessMAPs = Values(ProbDensity == max(ProbDensity));
        [ProbDensity, Values] = ksdensity(Chain(:, 5));
        ChoiceForgettingRateMAPs = Values(ProbDensity == max(ProbDensity));
        [ProbDensity, Values] = ksdensity(Chain(:, 6));
        BiasMAPs = Values(ProbDensity == max(ProbDensity));

        ParameterEstimates = [LearningRateMAPs, InverseTemperatureMAPs, ForgettingRateMAPs,...
                              ChoiceStickinessMAPs, ChoiceForgettingRateMAPs, BiasMAPs];
    else
        ParameterEstimates = Model.ParameterEstimates; % i.e. MLE
    end

    [NegLogDataLikelihood, Values] = ChoiceSymmetricQLearning(ParameterEstimates, nTrials, ChoiceLeft, Rewarded);

    SessionLogOdds = InverseTemperatureMAPs * (Values.LeftValue - Values.RightValue) +...
                     + ChoiceStickinessMAPs * Values.ChoiceMemory + BiasMAPs;
    
    SessionLeftValue = Values.LeftValue;
    SessionRightValue = Values.RightValue;
    SessionChoiceMemory = Values.ChoiceMemory;
    SessionChosenValue = SessionLeftValue .* ChoiceLeft + SessionRightValue .* (1 - ChoiceLeft);
    SessionUnchosenValue = SessionLeftValue .* (1 - ChoiceLeft) + SessionRightValue .* ChoiceLeft;
    SessionChosenMemory = SessionChoiceMemory .* ChoiceLeft - SessionChoiceMemory .* (1 - ChoiceLeft);
    
    SessionChosenLogOdds = InverseTemperatureMAPs * (SessionChosenValue - SessionUnchosenValue) +...
                           + ChoiceStickinessMAPs * SessionChosenMemory + BiasMAPs .* (2 * ChoiceLeft - 1);
    
    LeftValue = [LeftValue, SessionLeftValue];
    RightValue = [RightValue, SessionRightValue];
    LogOdds = [LogOdds, SessionLogOdds];
    ChosenValue = [ChosenValue, SessionChosenValue];
    UnchosenValue = [UnchosenValue, SessionUnchosenValue];
    ChosenMemory = [ChosenMemory, SessionChosenMemory];
    ChosenLogOdds = [ChosenLogOdds, SessionChosenLogOdds];


    SessionTimeInvestment = FeedbackWaitingTime;
    SessionTimeInvestment(~NotBaited) = nan;
    TimeInvestment = [TimeInvestment, SessionTimeInvestment];
    
    LeftTITrial = NotBaited & ChoiceLeft == 1;
    RightTITrial = NotBaited & ChoiceLeft == 0;
    LeftTI = FeedbackWaitingTime(LeftTITrial);
    RightTI = FeedbackWaitingTime(RightTITrial);
    
    [~, LeftTIMu, LeftTISigma] = zscore(LeftTI);
    [~, RightTIMu, RightTISigma] = zscore(RightTI);
    
    SessionTimeInvestmentZScore = sum((([FeedbackWaitingTime; FeedbackWaitingTime] - [LeftTIMu; RightTIMu]) ./ [LeftTISigma; RightTISigma]) .* ChoiceLeftRight, 1);
    TimeInvestmentZScore = [TimeInvestmentZScore, SessionTimeInvestmentZScore];

    [~, LeftTIMu, LeftTISigma] = zscore(sqrt(LeftTI));
    [~, RightTIMu, RightTISigma] = zscore(sqrt(RightTI));
    
    SessionTimeInvestmentSqrtZScore = sum(((sqrt([FeedbackWaitingTime; FeedbackWaitingTime]) - [LeftTIMu; RightTIMu]) ./ [LeftTISigma; RightTISigma]) .* ChoiceLeftRight, 1);
    TimeInvestmentSqrtZScore = [TimeInvestmentSqrtZScore, SessionTimeInvestmentSqrtZScore];

    Choices = [Choices, ChoiceLeft];
end

TotalValue = ChosenValue + UnchosenValue;
DiffValue = ChosenValue - UnchosenValue;

%% Initiatize figure
% create figure
if nargin < 4
    FigureSize = [0.2, 2.0, 4.1, 4.8];
end
AnalysisFigure = figure('Position', FigureSize,...
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateRange, '_Matching'),...
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

%% TI (z) vs Q
if nargin < 5
    AxeSize = [ 1.3, 1.3, 2.3, 2.8];
end

ChoiceLeftLogOddsAxes = axes(AnalysisFigure,...
                             'Position', AxeSize,...
                             'Units', 'inches',...
                             'Color', 'none');
set(ChoiceLeftLogOddsAxes,...
    'Position', AxeSize,...
    'Units', 'inches');
hold(ChoiceLeftLogOddsAxes, 'on');

Edges = linspace(-5, 5, 21);
BinCentres = movmean(Edges, 2, 'Endpoints', 'discard');
[LeftPDF, ~] = histcounts(LogOdds(Choices == 1), Edges, 'Normalization', 'pdf');
[RightPDF, ~] = histcounts(LogOdds(Choices == 0), Edges, 'Normalization', 'pdf');

ChoiceLeftProb = LeftPDF ./ (LeftPDF + RightPDF);
ChoiceLeftLogOddsLine = line(ChoiceLeftLogOddsAxes, BinCentres, ChoiceLeftProb * 100,...
                             'Color', 'k',...
                             'LineStyle', 'none',...
                             'Marker', 'o',...
                             'MarkerSize', 12);
%{
ValidIdx = ~isnan(ChoiceLeftProb);
Model = fit(BinCentres(ValidIdx)', ChoiceLeftProb(ValidIdx)', 'logistic');
PredictedChoiceLeftProb = feval(Model, BinCentres');

PredictedChoiceLeftLogOddsLine = line(ChoiceLeftLogOddsAxes, BinCentres, PredictedChoiceLeftProb * 100,...
                                      'Color', 'k',...
                                      'LineStyle', '-',...
                                      'Marker', 'none');

LegendText = {"data", "fitted"};
ChoiceLeftLogOddsLegend = legend(ChoiceLeftLogOddsAxes, LegendText,...
                                 'Position', [   2.4,  1.3,  1.2,  0.9],...
                                 'NumColumns', 1,...
                                 'Units', 'inches',...
                                 'Color', 'none',...
                                 'Box', 'off');

set(ChoiceLeftLogOddsLegend,...
    'Position', [   2.4,  1.3,  1.2,  0.9],...
    'Units', 'inches');
%}

nText = text(ChoiceLeftLogOddsAxes,...
             -4.5, 95,...
             sprintf(strcat('nTrial =', num2str(sum(~isnan(Choices))), '\nnRat = 1')),...
             'FontSize', 18);

set(ChoiceLeftLogOddsAxes,...
    'FontSize', 24,...
    'XLim', [-5, 5],...
    'YLim', [0, 100]);
xlabel(ChoiceLeftLogOddsAxes, 'Q_{L-R} (a.u.)')
ylabel(ChoiceLeftLogOddsAxes, 'P(c=L) (%)')


exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateRange, '_', AnalysisName, '.pdf'), 'ContentType', 'vector', 'Resolution', 300);
end