function AnalysisFigure = Matching_MS_RunLengthAnalysis()
%{
MS = MultiSession
First create on 20250725 by Antonio Lee for AG Ott @HU Berlin

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

AnalysisName = 'Matching_MS_RunLengthAnalysis';

%% Initiatize figure
% create figure
AnalysisFigure = figure('Position', [   0,    0,  842,  595],... % DIN A4, 72 ppi
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

%% Analysis across sessions
nStayedTrials = 10; % e.g. R71 has 53 session; 10 continuous run of same choice is only 29 occurrence

Counts = zeros(1, nStayedTrials);
Samples = [];
for iStayedTrial = 1:nStayedTrials

    String1 = replace(num2str([0, ones(1, iStayedTrial), 0]), ' ', '');
    String2 = replace(num2str([1, zeros(1, iStayedTrial), 1]), ' ', '');
    
    for iSession = 1:length(DataHolder)
        SessionData = DataHolder{iSession};

        nTrials = SessionData.nTrials;
        TrialData = SessionData.Custom.TrialData;
        
        ChoiceLeft = TrialData.ChoiceLeft(1:nTrials);
        Rewarded = TrialData.Rewarded(1:nTrials);
        NoTrialStart = TrialData.NoTrialStart(1:nTrials);
        BrokeFixation = TrialData.BrokeFixation(1:nTrials);
        
        ChoiceLeftChar = replace(num2str(ChoiceLeft == 1), ' ', '');
        ChoiceRightChar = replace(num2str(ChoiceLeft == 0), ' ', '');
        
        Index1 = strfind(ChoiceLeftChar, String1); % may include cases [NaN/Right Lefts NaN/Right]
        Index2 = strfind(ChoiceRightChar, String2); % may include cases [Right NaNs/Lefts Right]
        Index = intersect(Index1, Index2);

        Counts(iStayedTrial) = Counts(iStayedTrial) + length(Index);

        Index1 = strfind(ChoiceLeftChar, String2); % may include cases [NaN/Left Rights NaN/Left]
        Index2 = strfind(ChoiceRightChar, String1); % may include cases [Left NaNs/Rights Left]
        Index = intersect(Index1, Index2);

        Counts(iStayedTrial) = Counts(iStayedTrial) + length(Index);
        
    end

    Samples = [Samples, ones(1, Counts(iStayedTrial)) * iStayedTrial];
end

Probabilities = Counts ./ sum(Counts);

%% plot counts distribution
CountAxes = axes(AnalysisFigure, 'Position', [0.10, 0.10, 0.27, 0.31]);
hold(CountAxes, 'on');

CountBar = bar(CountAxes, Probabilities,...
               'EdgeColor', 'k',...
               'FaceColor', 'none');

set(CountAxes,...
    'TickDir', 'out',...
    'XLim', [0, nStayedTrials + 1],...
    'XTick', [1 5 10],...
    'FontSize', 12);
xlabel(CountAxes, 'Run length')
ylabel(CountAxes, 'Probability')

%% exponential fits
nGeometrics = 4;

NegLogLikelihood = zeros(1, nGeometrics);
for iGeometric = 1:nGeometrics
    Model = MultiGeometrics_Model(iGeometric, Samples);
    Models{iGeometric} = Model;
    
    EstimatedParameters = Model.EstimatedParameters;
    NegLogLikelihood(iGeometric) = Model.MinNegLogDataLikelihood;
    
    WeightIdx = 2 * (1:(iGeometric-1));
    if iGeometric > 1
        Weights = EstimatedParameters(WeightIdx);
        Weights = [Weights, 1 - sum(Weights)];
    else
        Weights = 1;
    end

    ProbIdx = 2 * (1:iGeometric) - 1;
    TransitionProbs = EstimatedParameters(ProbIdx);
    
    PredictedRunLengthProb = [];
    for iStayedTrial = 1:nStayedTrials
        PredictedRunLengthProb(iStayedTrial, 1:iGeometric)...
            = Weights...
              .* TransitionProbs...
              .* (1 - TransitionProbs) .^ (iStayedTrial - 1);

    end

    PredictedRunLengthProb = sum(PredictedRunLengthProb, 2);
    
    ColourPalette = CommonColourPalette(1./iGeometric);
    PredictedRunLengthProbLine{iGeometric}...
        = line(CountAxes, 1:nStayedTrials, PredictedRunLengthProb,...
               'Color', ColourPalette.RewardProbDark);
end

LegendString = [{'Data'}, strcat(string(1:nGeometrics), 'geometric(s)')];
CountLegend = legend(CountAxes, LegendString,...
                     'Position', [0.19, 0.24, 0.17, 0.17],...
                     'NumColumns', 1);

%% likelihood
LikelihoodAxes = axes(AnalysisFigure, 'Position', [0.10, 0.55, 0.27, 0.31]);
hold(LikelihoodAxes, 'on');

MaxLikelihoodLine = line(LikelihoodAxes, 1:nGeometrics, -NegLogLikelihood,...
                         'Color', 'k',...
                         'Marker', 'none');

set(LikelihoodAxes,...
    'TickDir', 'out',...
    'XLim', [0, nGeometrics + 1],...
    'FontSize', 12);
xlabel(LikelihoodAxes, 'Number of geometrics')
ylabel(LikelihoodAxes, 'log(likelihood)')

%% 2 geometrics
Model2Axes = axes(AnalysisFigure, 'Position', [0.40, 0.55, 0.27, 0.31]);
hold(Model2Axes, 'on');

LegendString = {'', ''};
Model = Models{2};
for iGeometric = 1:2
    Prob = Model.EstimatedParameters(iGeometric * 2 - 1);
    if iGeometric == 1
        Weight = Model.EstimatedParameters(2);
    else
        Weight = 1 - Model.EstimatedParameters(2);
    end

    PredictedRunLengthProb = Weight...
                             * Prob...
                             * (1 - Prob) .^ ((1:nStayedTrials) - 1);

    ColourPalette = CommonColourPalette(1./iGeometric);
    GeometricsLine{iGeometric}...
        = line(Model2Axes, 1:nStayedTrials, PredictedRunLengthProb,...
               'Color', ColourPalette.RewardProbDark);

    LegendString{iGeometric} = sprintf('p=%4.2f; y=%4.2f*p*(1-p)^x', Prob, Weight);
end

Model2Legend = legend(Model2Axes, LegendString,...
                     'Position', [0.42, 0.30, 0.25, 0.09],...
                     'NumColumns', 1);

set(Model2Axes,...
    'TickDir', 'out',...
    'XLim', [0, nStayedTrials + 1],...
    'XTick', [1 5 10],...
    'YLim', CountAxes.YLim,...
    'YTick', CountAxes.YTick,...
    'YTickLabel', {},...
    'FontSize', 12);
xlabel(Model2Axes, 'Run length')
title(Model2Axes, 'Sum of 2 geometrics')

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')
end % function