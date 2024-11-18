function TwoArmBanditVariant_Matching_MultiSession_KatahiraGLME(DataFolderPath)
%{
First create on 20240519 by Antonio Lee for AG Ott @HU Berlin

V1.0 20240519 dedicated script for analyzing multiple sessions with similar
settings (no checking is done). A folder is selected instead of a .mat file
(so that a bit of back and forth looking at the concatenated data and
individual session data is allowed)

The analysis is positioned to look at the data using Katahira models with
mixed effect, i.e. sessions are considered together. Use to confirm no
interaction between the reward history of 2 same consecutive choices on the
next one, and thus proving assymetric learning rate is unnecessary
%}

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

AnalysisName = 'Matching_MultiSession_KatahiraGLME';

%% Initiatize figure
% create figure
AnalysisFigure = figure('Position', [   0       0     595     420],... % DIN A5, 72 ppi
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
FigureInfoAxes = axes(AnalysisFigure, 'Position', [0.01    0.98    0.98    0.01]);
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

%% Analysis across sessions
% Coefficient of Katahira GLME
ModelCoefficientAxes = axes(AnalysisFigure, 'Position', [0.12    0.12    0.80    0.76]);
hold(ModelCoefficientAxes, 'on');

Last2RewardedHistory = [];
Last1RewardedHistory = [];
SessionNumber = [];
Stayed = [];

set(ModelCoefficientAxes,...
    'FontSize', 10)
xlabel(ModelCoefficientAxes, 'iTrial back');
ylabel(ModelCoefficientAxes, 'Coeff.');
title(ModelCoefficientAxes, 'GLME Fitted Coefficients')

%% Plotting
for iSession = 1:length(DataHolder)
    % Import SessionData
    SessionData = DataHolder{iSession};
    SessionDateLabel = [SessionDateLabel, string(datestr(datetime(SessionData.Info.SessionDate), 'YYYYmmDD(ddd)'))];
    
    nTrials = SessionData.nTrials;
    if nTrials < 200
        disp(['Session ', num2str(iSession), ' has nTrial < 200. Impossible for analysis.'])
        continue
    end
    
    ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
    Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);
    
    %% Concatenate data
    ValidIdx = find(~isnan(ChoiceLeft(1:end-2)) &...
                    ChoiceLeft(1:end-2) == ChoiceLeft(2:end-1) &...
                    ~isnan(ChoiceLeft(3:end)));
    
    Last2RewardedHistory = [Last2RewardedHistory, Rewarded(ValidIdx)];
    Last1RewardedHistory = [Last1RewardedHistory, Rewarded(ValidIdx + 1)];
    SessionNumber = [SessionNumber, iSession * ones(size(ValidIdx))];
    Stayed = [Stayed, ChoiceLeft(ValidIdx) == ChoiceLeft(ValidIdx + 2)];
    
end

disp('YOu aRE a bEAutIFul HUmaN BeiNG, saID anTOniO.')

%% Katahira model
Data = [Last2RewardedHistory', Last1RewardedHistory', SessionNumber', Stayed'];
Data = array2table(Data);

KatahiraModel = fitglme(Data, 'X4 ~ 1 + X1 * X2 + (1 + X1 * X2 | X3)',...
                        'Distribution', 'Binomial');

%% Average across sessions
SessionColor = [1, 1, 1] * 0.85; % ([1, 1, 1] - iSession / length(DataHolder)) * 0.9;

% Coefficient of Katahira GLME
Estimate = KatahiraModel.Coefficients.Estimate;
SE = KatahiraModel.Coefficients.SE;
PValue = KatahiraModel.Coefficients.pValue;
Upper95CI = KatahiraModel.Coefficients.Upper;
Lower95CI = KatahiraModel.Coefficients.Lower;

%{
DataPath = strcat(DataFolderPath, '\', FigureTitle, '.png');
exportgraphics(AnalysisFigure, DataPath);
%}
end % function