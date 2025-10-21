function Model = Matching_SS_MLE_LauGlimcherGLM_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
% LowerBound = [0.10, 6, 0.05, -1.5, 0.5, -0.5];
% UpperBound = [0.45, 10, 0.25, 0.5, 1, 0.5];

% Free parameters
% 20251001 tested with simulation that works well as initial parameters
Bias = 0; % beta_0
RewardCoeff = randn(1, 5) - 0.5; % beta_r_i
ChoiceCoeff = randn(1, 5) - 0.5; % beta_c_i

InitialParameters = [Bias, RewardCoeff, ChoiceCoeff];

CalculateMLE = @(Parameters) LauGlimcherGLM(Parameters, nTrials, ChoiceLeft, Rewarded);

try
    [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
        fmincon(CalculateMLE, InitialParameters, [], [], [], [], [], []); %LowerBound, UpperBound);
catch
    disp('Error: fail to run model');
    EstimatedParameters = [];
    MinNegLogDataLikelihood = nan;
end

% Model.LowerBound = LowerBound;
% Model.UpperBound = UpperBound;
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