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

%% Check if all sessions are of the same SettingsFile
%{
% not sure how to best do it yet, as some settings are drawn in each trial
and displayed
for iSession = 1:length(DataHolder)
    SessionData = DataHolder{iSession};

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
end
%}

%% Initiatize figure
% create figure
AnalysisFigure = figure('Position', [   0       0    1191     842],... % DIN A3, 72 ppi
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
end

Probabilities = Counts ./ sum(Counts);

%% exponential fits
nExponentials = 2;

NegLogLikelihood = zeros(1, nExponentials);
for iExponential = 1:nExponentials
    [Model, Stats] = fit([1:nStayedTrials]', Probabilities', strcat('exp', num2str(iExponential)));
end

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')
end % function