function AnalysisFigure = Matching_MS_PL_TimeInvestment_RewardRate_BackgroundDA_Mediation(DataFolderPath, PhotometryDatasetFilePath, ModelsFilePath, FigureSize, AxeSize)
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

AnalysisName = 'Matching_MS_PL_TimeInvestment_RewardRate_BackgroundDA_Mediation';

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
    FigureSize = [0.2, 2.0, 8, 4];
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

dFEstimationZScore = [];
FirstTrialIdx = 1;
for iSession = 1:length(PhotometryDataHolder)
    SessionData = PhotometryDataHolder{iSession};
    nTrials = SessionData.nTrials;
    
    SessiondFEstimation = dFEstimation(FirstTrialIdx:FirstTrialIdx + nTrials - 1);
    ValidIdx = ~isnan(SessiondFEstimation);

    SessiondFEstimationZScore = nan(size(SessiondFEstimation));
    SessiondFEstimationZScore(ValidIdx) = zscore(SessiondFEstimation(ValidIdx));

    dFEstimationZScore = [dFEstimationZScore, SessiondFEstimationZScore];

    FirstTrialIdx = FirstTrialIdx + nTrials;
end

%% axes structure
if nargin < 5
    AxeSize = [ 0, 0, 8, 4];
end

MediationAxes = axes(AnalysisFigure,...
                     'Position', AxeSize,...
                     'Units', 'inches',...
                     'Color', 'none');
set(MediationAxes,...
    'Position', AxeSize,...
    'Units', 'inches',...
    'XLim', [0, 8],...
    'XColor', 'none',...
    'YLim', [0, 4],...
    'YColor', 'none');
hold(MediationAxes, 'on');

% circles
CircleRadius = 1;
XCenter = [1, 1];
MCenter = XCenter + [2.6, 2.0];
YCenter = XCenter + [5.2, 0];

XDataAnnotation = annotation('ellipse',...
                             'Units', 'inches');
MDataAnnotation = annotation('ellipse',...
                             'Units', 'inches');
YDataAnnotation = annotation('ellipse',...
                             'Units', 'inches');
set(XDataAnnotation,...
    'Position', [XCenter(1)-CircleRadius, XCenter(2)-CircleRadius, 2*CircleRadius, 2*CircleRadius],...
    'LineWidth', 2);
set(MDataAnnotation,...
    'Position', [MCenter(1)-CircleRadius, MCenter(2)-CircleRadius, 2*CircleRadius, 2*CircleRadius],...
    'LineWidth', 2);
set(YDataAnnotation,...
    'Position', [YCenter(1)-CircleRadius, YCenter(2)-CircleRadius, 2*CircleRadius, 2*CircleRadius],...
    'LineWidth', 2);

% arrows
XYArrowAnnotation = annotation('arrow',...
                               'Units', 'inches');
XMArrowAnnotation = annotation('arrow',...
                               'Units', 'inches');
MYArrowAnnotation = annotation('arrow',...
                               'Units', 'inches');
MMArrowAnnotation = annotation('arrow',...
                               'Units', 'inches');
XYCenterCoordinateDiff = YCenter - XCenter;
XYCenterDistanceDiff = sqrt(XYCenterCoordinateDiff(1)^2 + XYCenterCoordinateDiff(2)^2);
XYArrowCoordinateBegin = XCenter + XYCenterCoordinateDiff * CircleRadius / XYCenterDistanceDiff;
XYArrowCoordinateEnd = XCenter + XYCenterCoordinateDiff * (1 - CircleRadius / XYCenterDistanceDiff);
set(XYArrowAnnotation,...
    'Position', [XYArrowCoordinateBegin(1), XYArrowCoordinateBegin(2), XYArrowCoordinateEnd(1)-XYArrowCoordinateBegin(1), XYArrowCoordinateEnd(2)-XYArrowCoordinateBegin(2)],...
    'LineWidth', 2,...
    'HeadStyle', 'cback2',...
    'HeadLength', 20,...
    'HeadWidth', 20);

XMCenterCoordinateDiff = MCenter - XCenter;
XMCenterDistanceDiff = sqrt(XMCenterCoordinateDiff(1)^2 + XMCenterCoordinateDiff(2)^2);
XMArrowCoordinateBegin = XCenter + XMCenterCoordinateDiff * CircleRadius / XMCenterDistanceDiff;
XMArrowCoordinateEnd = XCenter + XMCenterCoordinateDiff * (1 - CircleRadius / XMCenterDistanceDiff);
set(XMArrowAnnotation,...
    'Position', [XMArrowCoordinateBegin(1), XMArrowCoordinateBegin(2), XMArrowCoordinateEnd(1)-XMArrowCoordinateBegin(1), XMArrowCoordinateEnd(2)-XMArrowCoordinateBegin(2)],...
    'LineWidth', 2,...
    'HeadStyle', 'cback2',...
    'HeadLength', 20,...
    'HeadWidth', 20);

MYCenterCoordinateDiff = MCenter - YCenter;
MYCenterDistanceDiff = sqrt(MYCenterCoordinateDiff(1)^2 + MYCenterCoordinateDiff(2)^2);
MYArrowCoordinateBegin = MCenter - MYCenterCoordinateDiff * CircleRadius / MYCenterDistanceDiff;
MYArrowCoordinateEnd = MCenter - MYCenterCoordinateDiff * (1 - CircleRadius / MYCenterDistanceDiff);
set(MYArrowAnnotation,...
    'Position', [MYArrowCoordinateBegin(1), MYArrowCoordinateBegin(2), MYArrowCoordinateEnd(1)-MYArrowCoordinateBegin(1), MYArrowCoordinateEnd(2)-MYArrowCoordinateBegin(2)],...
    'LineWidth', 2,...
    'HeadStyle', 'cback2',...
    'HeadLength', 20,...
    'HeadWidth', 20);

set(MMArrowAnnotation,...
    'Position', [XMArrowCoordinateEnd(1), XMArrowCoordinateEnd(2), MYArrowCoordinateBegin(1)-XMArrowCoordinateEnd(1), MYArrowCoordinateBegin(2)-XMArrowCoordinateEnd(2)],...
    'LineWidth', 2,...
    'HeadStyle', 'cback2',...
    'HeadLength', 20,...
    'HeadWidth', 20);

%% variables
% text
XDataName = 'Reward rate';
MDataName = 'Dopamine';
YDataName = sprintf('Time\ninvestment');
XData = PhotometryTotalValue;
MData = TimeChoicedFEstimationZScore;
YData = PhotometryTimeInvestmentZScore;

XDataText = text(MediationAxes, XCenter(1), XCenter(2), XDataName,...
                 'FontSize', 24,...
                 'HorizontalAlignment', 'center');
MDataText = text(MediationAxes, MCenter(1), MCenter(2), MDataName,...
                 'FontSize', 24,...
                 'HorizontalAlignment', 'center');
YDataText = text(MediationAxes, YCenter(1), YCenter(2), YDataName,...
                 'FontSize', 24,...
                 'HorizontalAlignment', 'center');

% Y ~ M + X
bLM = fitlm([MData', XData'], YData');
b = bLM.Coefficients.Estimate(2);
bError = bLM.Coefficients.SE(2);
bTextPos = (MYArrowCoordinateBegin + MYArrowCoordinateEnd)/2;
bString = string(strcat('\beta_b', sprintf('=%4.2f±%4.2f', b, bError)));
bText = text(MediationAxes, bTextPos(1), bTextPos(2), bString,...
             'FontSize', 18,...
             'HorizontalAlignment', 'left',...
             'VerticalAlignment', 'bottom',...
             'Interpreter', 'tex');

cPrime = bLM.Coefficients.Estimate(3);
cPrimeError = bLM.Coefficients.SE(3);
cPrimeTextPos = (XYArrowCoordinateBegin + XYArrowCoordinateEnd)/2;
cPrimeString = string(strcat('\beta_{c\prime}', sprintf('=%4.2f±%4.2f', cPrime, cPrimeError)));
cPrimeText = text(MediationAxes, cPrimeTextPos(1), cPrimeTextPos(2), cPrimeString,...
                  'FontSize', 18,...
                  'HorizontalAlignment', 'center',...
                  'VerticalAlignment', 'top',...
                  'Interpreter', 'tex');

% M ~ X
aLM = fitlm(XData', MData');
a = aLM.Coefficients.Estimate(2);
aError = aLM.Coefficients.SE(2);
aTextPos = (XMArrowCoordinateBegin + XMArrowCoordinateEnd)/2;
aString = string(strcat('\beta_a', sprintf('=%4.2f±%4.2f', a, aError)));
aText = text(MediationAxes, aTextPos(1), aTextPos(2), aString,...
             'FontSize', 18,...
             'HorizontalAlignment', 'right',...
             'VerticalAlignment', 'bottom',...
             'Interpreter', 'tex');

% a*b
ab = a*b;
abError = abs(ab * sqrt((aError/a)^2 + (bError/b)^2));
abTextPos = (XMArrowCoordinateEnd + MYArrowCoordinateBegin)/2;
abString = string(strcat('\beta_a\beta_b', sprintf('=%4.2f±%4.2f', ab, abError)));
abText = text(MediationAxes, abTextPos(1), abTextPos(2), abString,...
              'FontSize', 18,...
              'HorizontalAlignment', 'center',...
              'VerticalAlignment', 'bottom',...
              'Interpreter', 'tex');

% c
cLM = fitlm(XData', YData');
c = cLM.Coefficients.Estimate(2);
cError = cLM.Coefficients.SE(2);
cTextPos = (XYArrowCoordinateBegin + XYArrowCoordinateEnd)/2;
cString = string(strcat('\beta_c', sprintf('=%4.2f±%4.2f', c, cError)));
cText = text(MediationAxes, cTextPos(1), cTextPos(2), cString,...
             'FontSize', 18,...
             'HorizontalAlignment', 'center',...
             'VerticalAlignment', 'bottom',...
             'Interpreter', 'tex');

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateRange, '_', AnalysisName, '.pdf'), 'ContentType', 'vector', 'Resolution', 300);
end