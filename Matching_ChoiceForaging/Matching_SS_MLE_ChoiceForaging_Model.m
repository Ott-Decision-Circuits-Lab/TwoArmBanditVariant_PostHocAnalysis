function Model = Matching_SS_MLE_ChoiceForaging_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [0.00, 2, 0.00, 0.0];
UpperBound = [1.00, 10, 1, 1];

% Free parameters
CalculateMLE = @(Parameters) ChoiceForaging(Parameters, nTrials, ChoiceLeft, Rewarded);

Model = struct();
Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.MinNegLogDataLikelihood = 0;

for iInitialCond = 1:10
    % 20250708 tested with simulation that works well as initial parameters
    LearningRate = rand() / 2; % alpha
    InverseTemperature = rand() * 5; % beta
    Threshold = rand() * 2; % theta
    ForgettingRate = rand() / 2; % gamma
    
    InitialParameters = [LearningRate, InverseTemperature, Threshold, ForgettingRate];

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
end
end % end function