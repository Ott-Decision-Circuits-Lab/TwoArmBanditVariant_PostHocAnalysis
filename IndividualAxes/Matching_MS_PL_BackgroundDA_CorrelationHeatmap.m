function AnalysisFigure = Matching_MS_PL_BackgroundDA_CorrelationHeatmap(DataFolderPath, PhotometryDatasetFilePath, ModelsFilePath, FigureSize, AxeSize)
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

AnalysisName = 'Matching_MS_PL_BackgroundDA_CorrelationHeatmap';

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
ChosenLogOdds = [];

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
    
    SessionChosenLogOdds = InverseTemperatureMAPs * (SessionChosenValue - SessionUnchosenValue) +...
                           + ChoiceStickinessMAPs * SessionChosenMemory + BiasMAPs .* (2 * ChoiceLeft - 1);

    LeftValue = [LeftValue, SessionLeftValue];
    RightValue = [RightValue, SessionRightValue];
    LogOdds = [LogOdds, SessionLogOdds];
    ChosenValue = [ChosenValue, SessionChosenValue];
    UnchosenValue = [UnchosenValue, SessionUnchosenValue];
    ChosenMemory = [ChosenMemory, SessionChosenMemory];
    ChosenLogOdds = [ChosenLogOdds, SessionChosenLogOdds];
end

TotalValue = ChosenValue + UnchosenValue;
DiffValue = ChosenValue - UnchosenValue;

%% dF/F from DatasetObject
ChannelName = 'GreenC';
AlignmentName = 'TimeTrialStart';
AlignmentWindow = [0, 2];

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

dFEstimationZScore = [];
FirstTrialIdx = 1;
for iSession = 1:length(PhotometryDataHolder)
    SessionData = PhotometryDataHolder{iSession};
    nTrials = SessionData.nTrials;
    
    SessiondFEstimation = dFEstimation(FirstTrialIdx:FirstTrialIdx + nTrials - 1);
    SessiondFEstimationZScore = zscore(SessiondFEstimation);
    dFEstimationZScore = [dFEstimationZScore, SessiondFEstimationZScore];

    FirstTrialIdx = FirstTrialIdx + nTrials;
end

%% Initiatize figure
% create figure
if nargin < 4
    FigureSize = [0.2, 0.5, 8.0, 6.8];
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

%% Background DA vs Chosen value
if nargin < 5
    AxeSize = [ 2.1, 2.0, 4.0, 4.0];
end

BackgroundDACorrelationAxes = axes(AnalysisFigure,...
                                   'Position', AxeSize,...
                                   'Units', 'inches',...
                                   'Color', 'none');
set(BackgroundDACorrelationAxes,...
    'Position', AxeSize,...
    'Units', 'inches',...
    'YDir', 'reverse');
hold(BackgroundDACorrelationAxes, 'on');

Data = [dFEstimationZScore', TotalValue', ChosenValue', ChosenLogOdds'];
ValidIdx = ~isnan(ChosenLogOdds) & ~isoutlier(dFEstimationZScore);

[RValues, pValues] = corr(Data(ValidIdx, :));
IsUpper = logical(triu(ones(size(RValues))));
RValues(IsUpper) = NaN;
pValues(IsUpper) = NaN;

BackgroundDACorrelationImagesc = imagesc(BackgroundDACorrelationAxes, RValues(2:end, 1:end-1));
clim(BackgroundDACorrelationAxes, [0, 1]);

ColourMap = othercolor('Greys9', 128); % ColourMap, function from helper folder
colormap(BackgroundDACorrelationAxes, ColourMap);

ColourBar = colorbar(BackgroundDACorrelationAxes);
set(ColourBar,...
    'FontSize', 24,...
    'TickDirection', 'out')
ylabel(ColourBar, 'Pearson R')

Labels = {"dF/F (z)", "Reward rate", "Choice value"};
nVar = width(Data);

set(BackgroundDACorrelationAxes,...
    'FontSize', 24,...
    'XLim', [0.5, 0.5 + nVar-1],...
    'XTick', 1:nVar-1,...
    'XTickLabel', Labels(1:nVar-1),...
    'XTickLabelRotation', 90,...
    'YLim', [0.5, 0.5 + nVar-1],...
    'YTick', 1:nVar-1,...
    'YTickLabel', Labels(2:nVar));

set(BackgroundDACorrelationAxes,...
    'Position', AxeSize)

nVar = height(RValues);
for iRow = 1:nVar-1
    for jColumn = 1:nVar-1
        if isnan(RValues(iRow+1, jColumn))
            continue
        end

        pValueSymbol = '';
        if pValues(iRow+1, jColumn) < 0.001
            pValueSymbol = '***';
        elseif pValues(iRow+1, jColumn) < 0.01
            pValueSymbol = '**';
        elseif pValues(iRow+1, jColumn) < 0.05
            pValueSymbol = '*';
        end
        
        TextColour = [0, 0, 0];
        if RValues(iRow+1, jColumn) > 0.5
            TextColour = [1, 1, 1];
        end

        text(BackgroundDACorrelationAxes, jColumn, iRow,... %index for array is flipped than x-y coordinates
             strcat(sprintf('%4.2f',...
                            RValues(iRow+1, jColumn)),...
                    pValueSymbol),...
             'FontSize', 18,...
             'Color', TextColour,...
             'HorizontalAlignment', 'center');
    end
end

nText = text(BackgroundDACorrelationAxes, 2, 2,...
             sprintf(strcat('nTrial=', num2str(sum(ValidIdx)), '\nnRat=2')),...
             'FontSize', 18,...
             'Units', 'inches');

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateRange, '_', AnalysisName, '.pdf'), 'ContentType', 'vector', 'Resolution', 300);
end