function AnalysisFigure = Matching_MS_PL_TimeInvestment_RewardRate_GLM(DataFolderPath, PhotometryDatasetFilePath, ModelsFilePath, FigureSize, AxeSize)
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

AnalysisName = 'Matching_MS_PL_TimeInvestment_RewardRate_GLM';

%{
%% PhotometryDataset
if nargin < 2
    AnimalDataFolder = fileparts(fileparts(fileparts(DataFolderPath)));
    PhotometryDatasetFilePath = fullfile(AnimalDataFolder, 'photometry', 'dataset');

    if ~isdir(PhotometryDatasetFilePath)
        disp('Error: No Photometry DatasetObject is found.')
        return
    end

    [PhotometryDatasetFileName, PhotometryDatasetFilePath] = uigetfile(PhotometryDatasetFilePath);
    
    if ~strcmpi(AnimalDataFolder, fileparts(fileparts(fileparts(PhotometryDatasetFilePath))))
        disp('Error: PhotometryDataset unlikely for the SessionDataset. Make sure they are in the adjacent directory')
        return
    end
    
    PhotometryDatasetFilePath = fullfile(PhotometryDatasetFilePath, PhotometryDatasetFileName);
end

try
    load(PhotometryDatasetFilePath); % should be 'DatasetObject' as variable name
catch
    disp('Error: Problem in loading PhotometryDatasetFilePath.')
    return
end

if exist('DatasetObject', 'var') == 0
    disp('Error: Loaded file is not a PhotometryDataset.')
    return
end
%}
%% Bayesian Symmetric Q-Learning with Forgetting and Stickiness model
if nargin < 3
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

%{
%% Filter out SessionData and Model for PhotometrySession
% not the most elegant way to filter session...ideally based on Photometry
% Dataset to select the SessionData from DataHolder
PhotometryDataHolder = {};
PhotometryModels = {};
for iSession = 1:length(DataHolder)
    if isfield(DataHolder{iSession}.Custom.SessionMeta, 'PhotometryValidation') && strcmpi(DataHolder{iSession}.Custom.SessionMeta.PhotometryBrainArea(end-1:end), 'PL') && DataHolder{iSession}.nTrials > 250
        PhotometryDataHolder = [PhotometryDataHolder, DataHolder(iSession)];
        PhotometryModels = [PhotometryModels, Models(iSession)];
    end
end

if length(PhotometryDataHolder) ~= height(DatasetObject.ProcessedPhotometryData)
    disp('Error: Mismatched nSessions between PhotometryDatasetObject and FilteredPhotometrySession from SessionData')
    return
end
%}

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
ChosenBias = [];

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
    
    SessionChosenLogOdds = SessionLogOdds .* ChoiceLeft - SessionLogOdds .* (1 - ChoiceLeft);
    SessionChosenBias = BiasMAPs .* ChoiceLeft - BiasMAPs .* (1 - ChoiceLeft);

    LeftValue = [LeftValue, SessionLeftValue];
    RightValue = [RightValue, SessionRightValue];
    LogOdds = [LogOdds, SessionLogOdds];
    ChosenValue = [ChosenValue, SessionChosenValue];
    UnchosenValue = [UnchosenValue, SessionUnchosenValue];
    ChosenMemory = [ChosenMemory, SessionChosenMemory];
    ChosenLogOdds = [ChosenLogOdds, SessionChosenLogOdds];
    ChosenBias = [ChosenBias, SessionChosenBias];

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
    SessionTimeInvestmentZScore(~NotBaited) = nan;
    TimeInvestmentZScore = [TimeInvestmentZScore, SessionTimeInvestmentZScore];

    [~, LeftTIMu, LeftTISigma] = zscore(sqrt(LeftTI));
    [~, RightTIMu, RightTISigma] = zscore(sqrt(RightTI));
    
    SessionTimeInvestmentSqrtZScore = sum(((sqrt([FeedbackWaitingTime; FeedbackWaitingTime]) - [LeftTIMu; RightTIMu]) ./ [LeftTISigma; RightTISigma]) .* ChoiceLeftRight, 1);
    SessionTimeInvestmentSqrtZScore(~NotBaited) = nan;
    TimeInvestmentSqrtZScore = [TimeInvestmentSqrtZScore, SessionTimeInvestmentSqrtZScore];

    Choices = [Choices, ChoiceLeft];
end

TotalValue = ChosenValue + UnchosenValue;
DiffValue = ChosenValue - UnchosenValue;

%% Initiatize figure
% create figure
if nargin < 4
    FigureSize = [0.2, 2.0, 5.8, 5.8];
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
    'Color', 'w')

% colour palette
ColourPalette = CommonColourPalette();

%% dF/F from DatasetObject
%{
ChannelName = 'GreenC';
AlignmentName = 'TimeChoice';
AlignmentWindow = [-2, 0];

[AlignedDemodData,...
 AlignedBaselineData,...
 AlignedTime] = DatasetObject.AlignData('Channel', ChannelName,...
                                        'Alignment', AlignmentName, ...
                                        'AlignmentWindow', AlignmentWindow);

% Prepare AlignData for plotting
NormalisedPreAlignedPhotometryData = (AlignedDemodData ./ AlignedBaselineData - 1) * 100;

% Estimate dF
EstimatingFunction = str2func('median');

try
    dFEstimation = EstimatingFunction(NormalisedPreAlignedPhotometryData , 2);
catch
    disp('Error: Incompatible estimating function for photometry data of nTrial x length of recording')
    return
end
dFEstimation = dFEstimation';
%}
%% TI  ~ Q + R
if nargin < 5
    AxeSize = [ 1.2, 2.2, 3.8, 2.8];
end

GLMCoeffAxes = axes(AnalysisFigure,...
                    'Position', AxeSize,...
                    'Units', 'inches',...
                    'Color', 'none');
set(GLMCoeffAxes,...
    'Position', AxeSize,...
    'Units', 'inches');
hold(GLMCoeffAxes, 'on');

XData = [ChosenLogOdds', TotalValue', ChosenValue', UnchosenValue'];
YData = TimeInvestmentZScore;
ValidIdx = ~isnan(ChosenValue) & ~isnan(YData);

GLM = fitglm(XData(ValidIdx, :), YData(ValidIdx));

GLMCoeffBar = bar(GLMCoeffAxes, GLM.Coefficients.Estimate,...
                  'FaceColor', 'none');
GLMCoeffErrorbar = errorbar(GLMCoeffAxes,...
                            GLM.Coefficients.Estimate,...
                            GLM.Coefficients.SE,...
                            'LineStyle', 'none',...
                            'Color', 'k');

pValues = GLM.Coefficients.pValue;
for iCoeff = 1:length(pValues)
    pValueSymbol = "";
    if pValues(iCoeff) < 0.001
        pValueSymbol = "***";
    elseif pValues(iCoeff) < 0.01
        pValueSymbol = "**";
    elseif pValues(iCoeff) < 0.05
        pValueSymbol = "*";
    end

    text(GLMCoeffAxes, iCoeff, 1.2, pValueSymbol,...
         'FontSize', 24,...
         'HorizontalAlignment', 'center')
end

Labels = {"Intercept", "Choice value", "Reward rate", "V_{chosen}", "V_{unchosen}"};

set(GLMCoeffAxes,...
    'FontSize', 24,...
    'XLim', [0, width(XData)+2],...
    'XTick', 1:width(XData)+1,...
    'XTickLabel', Labels,...
    'XTickLabelRotation', 90,...
    'YLim', [-0.5, 1.5],...
    'TickDir', 'out');
ylabel(GLMCoeffAxes, '\beta')
title(GLMCoeffAxes,...
      strcat('Time investment (z) ~ \betaX'))

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateRange, '_', AnalysisName, '.pdf'), 'ContentType', 'vector', 'Resolution', 300);
end