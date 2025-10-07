function Model = Matching_SS_MLE_ChoiceForaging_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [0.00, 2, 0.00, 0.0];
UpperBound = [1.00, 10, 1, 1];

% Free parameters
% 20250708 tested with simulation that works well as initial parameters
LearningRate = 0.30; % alpha
InverseTemperature = 3; % beta
Threshold = 0.6; % theta
ForgettingRate = 0.15; % gamma

InitialParameters = [LearningRate, InverseTemperature, Threshold, ForgettingRate];

CalculateMLE = @(Parameters) ChoiceForaging(Parameters, nTrials, ChoiceLeft, Rewarded);

try
    [EstimatedParameters, MinNegLogDataLikelihood] =...
        fmincon(CalculateMLE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
catch
    disp('Error: fail to run model');
    EstimatedParameters = [];
    MinNegLogDataLikelihood = nan;
end

Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.InitialParameters = InitialParameters;

Model.EstimatedParameters = EstimatedParameters;
Model.MinNegLogDataLikelihood = MinNegLogDataLikelihood;

end % end function