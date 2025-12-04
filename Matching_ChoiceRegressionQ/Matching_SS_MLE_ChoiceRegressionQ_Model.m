function Model = Matching_SS_MLE_ChoiceRegressionQ_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [-5, 0, -5, -5, -5, -5];
UpperBound = [5, 1, 5, 5, 5, 5];

% Free parameters
% 20251001 tested with simulation that works well as initial parameters
Bias = 0; % beta_0
QCoeff = 0; % gamma
Minus1ChoiceCoeff = 0; % beta_c1
Minus1RewardCoeff = 0; % beta_r1
Minus2ChoiceCoeff = 0; % beta_c2
Minus2RewardCoeff = 0; % beta_r2

CalculateMLE = @(Parameters) ChoiceRegressionQ(Parameters, nTrials, ChoiceLeft, Rewarded);

Model = struct();
Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.MinNegLogDataLikelihood = 0;

for iInitialCond = 1:10
    % 20250708 tested with simulation that works well as initial parameters
    Bias = randn(); % beta_0
    QCoeff = rand(); % gamma
    Minus1ChoiceCoeff = randn(); % beta_c1
    Minus1RewardCoeff = randn(); % beta_r1
    Minus2ChoiceCoeff = randn(); % beta_c2
    Minus2RewardCoeff = randn(); % beta_r2
    
    InitialParameters = [Bias, QCoeff, Minus1ChoiceCoeff, Minus1RewardCoeff, Minus2ChoiceCoeff, Minus2RewardCoeff];
    
    try
        [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
            fmincon(CalculateMLE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
    catch
        disp('Error: fail to run model');
        EstimatedParameters = [];
        MinNegLogDataLikelihood = nan;
    end
    
    if Model.MinNegLogDataLikelihood < MinNegLogDataLikelihood
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
    end
end % end function