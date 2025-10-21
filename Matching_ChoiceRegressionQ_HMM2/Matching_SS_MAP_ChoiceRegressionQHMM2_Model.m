function Model = Matching_SS_MAP_ChoiceRegressionQHMM2_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
% LowerBound = [0.10, 6, 0.05, -1.5, 0.5, -0.5];
% UpperBound = [0.45, 10, 0.25, 0.5, 1, 0.5];

% Free parameters
% 20251001 from Ashwood or Venditto initial parameters
try
    load(fullfile(DataFolderPath, strcat('\Matching_ConcatMS_MAP_ChoiceRegressionQHMM2.mat')));
catch
    disp('Error: A subject-level parameter estimated from a concatenated session is needed for initializing parameters. No further model can be estimated. ')
    return
end

SubjectConcatSessionEstimatedParameter = Model.EstimatedParameters;

for iInitialCond = 1:10
    
GLMWeights = randn(2, 4) * 1; % 2 states x 4 weights (QCoeff, reward, choice, bias)
TransitionMatrix = 

InitialParameters = [Bias, QCoeff, Minus1ChoiceCoeff, Minus1RewardCoeff, Minus2ChoiceCoeff, Minus2RewardCoeff];

CalculateMLE = @(Parameters) ChoiceRegressionQ(Parameters, nTrials, ChoiceLeft, Rewarded);

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