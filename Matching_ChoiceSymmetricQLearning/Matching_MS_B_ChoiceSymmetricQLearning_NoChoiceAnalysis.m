function AnalysisFigure = Matching_MS_B_ChoiceSymmetricQLearning_NoChoiceAnalysis(DataFolderPath)
% MS = MultiSession
% B = Bayesian <- Prior using simulation & MCMC (Hamiltonian MC) sampling
% from prior to get marginal posterior
% Matching Analysis Function
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ Jan 2025
% Model iteration see the end of script

%% load files
if nargin < 1
    DataFolderPath = uigetdir(OttLabDataServerFolderPath());
elseif ~ischar(DataFolderPath) && ~isstring(DataFolderPath)
    disp('Error: Unknown input format. No further analysis can be performed.')
    return
end

try
    load(fullfile(DataFolderPath, '\Selected_Data.mat'));
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

AnalysisName = 'Matching_MS_B_ChoiceSymmetricQLearning';

%% Hierarchaical Symmetric Q-Learning with Forgetting and Stickiness model
try
    load(fullfile(DataFolderPath, strcat('\', AnalysisName, '.mat')));
catch
    disp('Error: no models are found')
    return
end

if ~exist('Models', 'var')
    disp('Error: Loaded data is not a Models')
    return
end

%%
nSessions = length(DataHolder);

NoTrialStartCount = zeros(1, nSessions);
DiscountedNoTrialStartNLL = zeros(1, nSessions);
UndiscountedNoTrialStartNLL = zeros(1, nSessions);

NoDecisionCount = zeros(1, nSessions);
DiscountedNoDecisionNLL = zeros(1, nSessions);
UndiscountedNoDecisionNLL = zeros(1, nSessions);

for iSession = 1:nSessions

    % Import SessionData
    SessionData = DataHolder{iSession};

    nTrials = SessionData.nTrials;
    if nTrials < 200
        disp(['Session ', num2str(iSession), ' has nTrial < 200. Impossible for analysis.'])
        continue
    end

    %%
    ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
    Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);
    
    %% Get Values from Model
    Model = Models{iSession};
    Chain = vertcat(Model.Chains{:});
    
    [ProbDensity, Values] = ksdensity(Chain(:, 1));
    LearningRateMAP = Values(ProbDensity == max(ProbDensity));
    
    [ProbDensity, Values] = ksdensity(Chain(:, 2));
    InverseTemperatureMAP = Values(ProbDensity == max(ProbDensity));
    
    [ProbDensity, Values] = ksdensity(Chain(:, 3));
    ForgettingRateMAP = Values(ProbDensity == max(ProbDensity));
    
    [ProbDensity, Values] = ksdensity(Chain(:, 4));
    ChoiceStickinessMAP = Values(ProbDensity == max(ProbDensity));
    
    [ProbDensity, Values] = ksdensity(Chain(:, 5));
    ChoiceForgettingRateMAP = Values(ProbDensity == max(ProbDensity));
    
    [ProbDensity, Values] = ksdensity(Chain(:, 6));
    BiasMAP = Values(ProbDensity == max(ProbDensity));
    
    MAPEstimates = [LearningRateMAP, InverseTemperatureMAP, ForgettingRateMAP,...
                    ChoiceStickinessMAP, ChoiceForgettingRateMAP, BiasMAP];
    
    [~, Values] = ChoiceSymmetricQLearning(MAPEstimates, nTrials, ChoiceLeft, Rewarded);
    
    LeftValue = Values.LeftValue;
    RightValue = Values.RightValue;
    ChoiceMemory = Values.ChoiceMemory;
    ChoiceLeftLogOdds = Values.ChoiceLeftLogOdds;
    
    ChoiceLeftProb = 1 ./ (1 + exp(-ChoiceLeftLogOdds));

    %% check if Values should be discounted when NoTrialStart
    NoTrialStart = SessionData.Custom.TrialData.NoTrialStart(1:nTrials);
    Ndx = strfind(NoTrialStart, [0, 1, 0]);
    
    for iNdx = 1:length(Ndx)
        Idx = Ndx(iNdx);
        if ChoiceLeft(Idx + 2) == 0
            DiscountedNoTrialStartNLL(iSession) = DiscountedNoTrialStartNLL(iSession) - log(ChoiceLeftProb(Idx + 2));
            UndiscountedNoTrialStartNLL(iSession) = UndiscountedNoTrialStartNLL(iSession) - log(ChoiceLeftProb(Idx + 1));
        elseif ChoiceLeft(Idx + 2) == 1
            DiscountedNoTrialStartNLL(iSession) = DiscountedNoTrialStartNLL(iSession) - log((1 - ChoiceLeftProb(Idx + 2)));
            UndiscountedNoTrialStartNLL(iSession) = UndiscountedNoTrialStartNLL(iSession) - log((1 - ChoiceLeftProb(Idx + 1)));
        end
    end
    
    NoTrialStartCount(iSession) =  length(Ndx);
    
    %% check if Values should be discounted when NoDecision
    NoDecision = SessionData.Custom.TrialData.NoDecision(1:nTrials);
    Ndx = strfind(NoDecision, [0, 1, 0]);
    
    for iNdx = 1:length(Ndx)
        Idx = Ndx(iNdx);
        if ChoiceLeft(Idx + 2) == 0
            DiscountedNoDecisionNLL(iSession) = DiscountedNoDecisionNLL(iSession) - log(ChoiceLeftProb(Idx + 2));
            UndiscountedNoDecisionNLL(iSession) = UndiscountedNoDecisionNLL(iSession) - log(ChoiceLeftProb(Idx + 1));
        elseif ChoiceLeft(Idx + 2) == 1
            DiscountedNoDecisionNLL(iSession) = DiscountedNoDecisionNLL(iSession) - log((1 - ChoiceLeftProb(Idx + 2)));
            UndiscountedNoDecisionNLL(iSession) = UndiscountedNoDecisionNLL(iSession) - log((1 - ChoiceLeftProb(Idx + 1)));
        end
    end
    
    NoDecisionCount(iSession) =  length(Ndx);
    
end

%% create figure
% create figure
AnalysisFigure = figure('Position', [    0,    0,  842,  ],... % DIN A4, 72 ppi
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

%%
NoTrialStartAxes = axes(AnalysisFigure, 'Position', [0.05, 0.06, 0.85, 0.85]);
hold(NoTrialStartAxes, 'on')

yyaxis(NoTrialStartAxes, 'left')

SessionNoTrialStartNLLPlot = plot(NoTrialStartAxes,...
                                [0, 1],...
                                [DiscountedNoTrialStartNLL', UndiscountedNoTrialStartNLL'],...
                                'Color', ColourPalette.Session,...
                                'Marker', 'o');

[h,p,ci,stats] = ttest(DiscountedNoTrialStartNLL, UndiscountedNoTrialStartNLL);
SignificanceText = text(0.5, 2, sprintf('p=%4.2g', p));
set(SignificanceText,...
    'HorizontalAlignment', 'center')

set(NoTrialStartAxes,...
    'FontSize', 18,...
    'XLim', [-0.5, 1.5],...
    'XTick', [0, 1],...
    'XTickLabel', {'Extra discounting', 'No discounting'})

title(NoTrialStartAxes, 'NoTrialStart Trial')
ylabel(NoTrialStartAxes, 'log(likelihood)')

yyaxis(NoTrialStart, 'right')

TrialDiscountedNoTrialStartNLL = sum(DiscountedNoTrialStartNLL) / sum(NoTrialStartCount);
TrialUndiscountedNoTrialStartNLL = sum(UndiscountedNoTrialStartNLL) / sum(NoTrialStartCount);

NoTrialStartNLLPlot = plot(NoTrialStartAxes,...
                           [0, 1],...
                           [TrialDiscountedNoTrialStartNLL', TrialUndiscountedNoTrialStartNLL'],...
                           'Color', ColourPalette.Pooled,...
                           'Marker', 'o');

ylabel(NoTrialStartAxes, 'log(likelihood)_{trial average}')

%%
NoDecisionAxes = axes(AnalysisFigure, 'Position', [0.01    0.75    0.15    0.19]);
hold(NoDecisionAxes, 'on')

yyaxis(NoDecisionAxes, 'left')

SessionNoDecisionNLLPlot = plot(NoDecisionAxes,...
                                [0, 1],...
                                [DiscountedNoDecisionNLL', UndiscountedNoDecisionNLL'],...
                                'Color', ColourPalette.Session,...
                                'Marker', 'o');

[h,p,ci,stats] = ttest(DiscountedNoDecisionNLL, UndiscountedNoDecisionNLL);
SignificanceText = text(0.5, 2, sprintf('p=%4.2g', p));
set(SignificanceText,...
    'HorizontalAlignment', 'center')

set(NoDecisionAxes,...
    'FontSize', 18,...
    'XLim', [-0.5, 1.5],...
    'XTick', [0, 1],...
    'XTickLabel', {'Extra discounting', 'No discounting'})
title(NoDecisionAxes, 'NoDecision Trial')
ylabel(NoDecisionAxes, 'log(likelihood)')

yyaxis(NoDecision, 'right')

TrialDiscountedNoDecisionNLL = sum(DiscountedNoDecisionNLL) / sum(NoDecisionCount);
TrialUndiscountedNoDecisionNLL = sum(UndiscountedNoDecisionNLL) / sum(NoDecisionCount);

NoDecisionNLLPlot = plot(NoDecisionAxes,...
                         [0, 1],...
                         [TrialDiscountedNoDecisionNLL', TrialUndiscountedNoDecisionNLL'],...
                         'Color', ColourPalette.Pooled,...
                         'Marker', 'o');

ylabel(NoDecisionAxes, 'log(likelihood)_{trial average}')

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')
end