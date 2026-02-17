function Model = Matching_SS_MLE_ChoiceSymmetricQLearningLastC_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [0.10, 2, 0.05, -2, -0.5];
UpperBound = [0.45, 10, 0.25, 2, 0.5];

% Free parameters
% 20241220 tested with simulation that works well as initial parameters
LearningRate = 0.25; % alpha
InverseTemperature = 8; % beta
ForgettingRate = 0.15; % gamma
ChoiceStickiness = 0; % phi
Bias = 0;

InitialParameters = [LearningRate, InverseTemperature, ForgettingRate, ChoiceStickiness, Bias];

CalculateMLE = @(Parameters) ChoiceSymmetricQLearningLastC(Parameters, nTrials, ChoiceLeft, Rewarded);

try
    [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
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
Model.Grad = Grad;
Model.Hessian = Hessian;
try
    Model.ParameterStandardError = sqrt(diag(inv(Hessian)))';
catch
    Model.ParameterStandardError = nan(size(EstimatedParameters));
end

end % end function