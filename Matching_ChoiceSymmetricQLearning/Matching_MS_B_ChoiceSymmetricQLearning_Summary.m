function Matching_MS_B_ChoiceSymmetricQLearning_Summary(DataFolderPath)
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

SessionDateRange = DataFolderPath(end-16:end);
[~, RatName] = fileparts(fileparts(fileparts(DataFolderPath)));

RatID = str2double(RatName);
if isnan(RatID)
    RatID = -1;
end
RatName = num2str(RatID);

AnalysisName = 'Matching_MS_B_ChoiceSymmetricQLearning';

FigureTitle = strcat(RatName, '_', SessionDateRange, '_', AnalysisName);
DataPath = strcat(DataFolderPath, '\', FigureTitle, '.pdf');

%% 
AnalysisFigure = Matching_MS_B_ChoiceSymmetricQLearning_Visualisation(DataFolderPath);
exportgraphics(AnalysisFigure, DataPath, 'Append', true);

close(AnalysisFigure)

%%
AnalysisFigure = Matching_MS_B_ChoiceSymmetricQLearning_Visualisation(DataFolderPath);
exportgraphics(AnalysisFigure, DataPath, 'Append', true);

close(AnalysisFigure)

%%
load(fullfile(DataFolderPath, '\Selected_Data.mat'));
load(fullfile(DataFolderPath, strcat('\', AnalysisName, '.mat')));

for iSession = 1:length(DataHolder)
    SessionData = DataHolder{iSession};
    Model = Models{iSession};

    AnalysisFigure = Matching_SS_B_ChoiceSymmetricQLearning_Visualisation(SessionData, Model);
    exportgraphics(AnalysisFigure, DataPath, 'Append', true);
    
    close(AnalysisFigure)

    AnalysisFigure = Matching_SS_B_ChoiceSymmetricQLearning_Diagnosis(SessionData, Model);
    exportgraphics(AnalysisFigure, DataPath, 'Append', true);
    
    close(AnalysisFigure)
end

end