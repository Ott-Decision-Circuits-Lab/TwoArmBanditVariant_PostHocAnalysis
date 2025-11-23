function AnalysisFigure = Matching_MS_PL_BackgroundDA_GLM_Traces(DataFolderPath, PhotometryDatasetFilePath, ModelsFilePath, FigureSize, AxeSize)
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

AnalysisName = 'Matching_MS_PL_BackgroundDA_GLM_Traces';

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

%% Load related data to local variabels
LogOdds = [];
ChosenValue = [];
UnchosenValue = [];
ChosenMemory = [];
LeftValue = [];
RightValue = [];

for iPhotoSession = 1:length(PhotometryDataHolder)
    SessionData = PhotometryDataHolder{iPhotoSession};
    Model = PhotometryModels{iPhotoSession};

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
    
    LeftValue = [LeftValue, SessionLeftValue];
    RightValue = [RightValue, SessionRightValue];
    LogOdds = [LogOdds, SessionLogOdds];
    ChosenValue = [ChosenValue, SessionChosenValue];
    UnchosenValue = [UnchosenValue, SessionUnchosenValue];
    ChosenMemory = [ChosenMemory, SessionChosenMemory];
end

TotalValue = ChosenValue + UnchosenValue;

%% Initiatize figure
% create figure
if nargin < 4
    FigureSize = [0.2, 0.2, 5.6, 4.5];
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

%% dF/F from DatasetObject
ChannelName = 'GreenC';
AlignmentName = 'TimeChoice';
AlignmentWindow = [-2, 1];

[AlignedDemodData,...
 AlignedBaselineData,...
 AlignedTime] = DatasetObject.AlignData('Channel', ChannelName,...
                                        'Alignment', AlignmentName, ...
                                        'AlignmentWindow', AlignmentWindow);

% prepare AlignData for plotting
NormalisedPreAlignedPhotometryData = (AlignedDemodData ./ AlignedBaselineData - 1) * 100;

PhotometryDataValidity = sum(isnan(NormalisedPreAlignedPhotometryData), 2)./width(NormalisedPreAlignedPhotometryData) < 0.9;
PhotometryOutlier = sum(isoutlier(NormalisedPreAlignedPhotometryData), 2)./width(NormalisedPreAlignedPhotometryData) > 0.5;
PhotometryDataValidity = PhotometryDataValidity & ~PhotometryOutlier;
AlignmentDataValidity = ~isnan(DatasetObject.TrialEventTimeData.(AlignmentName)) & PhotometryDataValidity';

% sort
SortingBase = TotalValue;
[~, SortingIdx] = sort(SortingBase);
SortedValidIdx = SortingIdx(AlignmentDataValidity(SortingIdx) == 1);

SortedAlignedData = NormalisedPreAlignedPhotometryData(SortedValidIdx, :);

%% Background DA vs Total value
if nargin < 5
    AxeSize = [1.2, 1.2, 0.6 * diff(AlignmentWindow), 2.5];
end

BackgroundDATracesAxes = axes(AnalysisFigure,...
                        'Position', AxeSize,...
                        'Units', 'inches',...
                        'Color', 'none');
set(BackgroundDATracesAxes,...
    'Position', AxeSize,...
    'Units', 'inches');
hold(BackgroundDATracesAxes, 'on');

% GLM
Regressors = [TotalValue', ChosenValue', DatasetObject.TrialEventTimeData.Rewarded' == 1];
RegressorNames = {'Reward rate', 'V_{chosen}', 'Rewarded'};

nRegressor = width(Regressors);
YData = TimeChoicedFTrace;

CoefficientDynamics = [];
Upper95CI = []; % confidence interval
Lower95CI = [];
pValue = [];
for iTimePoint = 1:width(YData)
    LinearModel = fitlm(Regressors(AlignmentDataValidity, :), YData(AlignmentDataValidity, iTimePoint), 'RobustOpts', 'on');
    CoefficientDynamics(:, iTimePoint) = LinearModel.Coefficients.Estimate;
    pValue(:, iTimePoint) = LinearModel.Coefficients.pValue;

    CI = coefCI(LinearModel);
    Lower95CI(:, iTimePoint) = CI(:, 1);
    Upper95CI(:, iTimePoint) = CI(:, 2);
end
CoefficientDynamics = [CoefficientDynamics(2:end, :); CoefficientDynamics(1, :)];
pValue = [pValue(2:end, :); pValue(1, :)];
Lower95CI = [Lower95CI(2:end, :); Lower95CI(1, :)];
Upper95CI = [Upper95CI(2:end, :); Upper95CI(1, :)];

% plotting
ColourMap = [ColourPalette.Incorrect; ColourPalette.General; ColourPalette.Pooled; ColourPalette.Session];

for iRegressor = 1:(nRegressor + 1)
    CoefficientDynamicsLine(iRegressor) = line(BackgroundDATracesAxes,...
                                               'XData', AlignedTime,...
                                               'YData', smooth(CoefficientDynamics(iRegressor, :), 0.1),...
                                               'Color', ColourMap(iRegressor, :),...
                                               'LineStyle', '-',...
                                               'LineWidth', 2);
end

for iRegressor = 1:(nRegressor + 1) % separated to allow easy legend labelling
    Upper95CILine(iRegressor) = line(BackgroundDATracesAxes,...
                                     'XData', AlignedTime,...
                                     'YData', smooth(Upper95CI(iRegressor, :), 0.1),...
                                     'Color', ColourMap(iRegressor, :),...
                                     'LineStyle', ':',...
                                     'LineWidth', 1.5);
    
    Lower95CILine(iRegressor) = line(BackgroundDATracesAxes,...
                                     'XData', AlignedTime,...
                                     'YData', smooth(Lower95CI(iRegressor, :), 0.1),...
                                     'Color', ColourMap(iRegressor, :),...
                                     'LineStyle', ':',...
                                     'LineWidth', 1.5);
    
end

LegendText = [RegressorNames(:)', {'Intercept'}];

RegressordFLegend = legend(BackgroundDATracesAxes, LegendText,...
                           'Position', [4.5, 2.9, 1.3, 0.6],...
                           'NumColumns', 1,...
                           'Units', 'inches',...
                           'Color', 'w',...
                           'Box', 'on');

set(RegressordFLegend,...
    'Position', [3.1, 1.2, 2.2, 2.5],...
    'Units', 'inches');

% set colour
set(BackgroundDATracesAxes,...
    'Position', AxeSize,...
    'Units', 'inches',...
    'TickDir', 'out',...
    'XLim', AlignmentWindow,...
    'XTick', sort([AlignmentWindow, 0]),...
    'YLim', [-5, 15],...
    'YTick', [-5, 5, 15],...
    'FontSize', 24);
xlabel(BackgroundDATracesAxes, 'Time (s)')
ylabel(BackgroundDATracesAxes, '\beta')
title(BackgroundDATracesAxes, 'At choice')

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateRange, '_', AnalysisName, '_', AlignmentName, '.pdf'), 'ContentType', 'vector', 'Resolution', 300);
end