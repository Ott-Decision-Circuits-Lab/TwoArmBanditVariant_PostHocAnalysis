function Model = Matching_SS_MLE_ChoiceSymmetricRLearning_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [0.10, 0.00, -5, 0.05, -1.5, 0.5, -0.5];
UpperBound = [0.45, 0.30, 5, 0.20, -0.5, 1, 0.5];

% Free parameters
% 20241220 tested with simulation that works well as initial parameters
LearningRate = 0.25; % alpha
RLearningRate = 0.25; % alpha_r
InverseTemperature = 0; % beta
ForgettingRate = 0.15; % gamma
ChoiceStickiness = -1; % phi
ChoiceForgettingRate = 1; % gamma_c
Bias = 0;

InitialParameters = [LearningRate, RLearningRate, InverseTemperature, ForgettingRate, ChoiceStickiness, ChoiceForgettingRate, Bias];

CalculateMLE = @(Parameters) ChoiceSymmetricRLearning(Parameters, nTrials, ChoiceLeft, Rewarded);

try
    [EstimatedParameters, MinNegLogDataLikelihood] =...
        fmincon(CalculateMLE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
catch
    disp('Error: fail to run model');
    EstimatedParameters = [];
    MinNegLogDataLikelihood = nan;
end

Model.EstimatedParameters = EstimatedParameters;
Model.MinNegLogDataLikelihood = MinNegLogDataLikelihood;

end % end function