function Model = Matching_SS_MLE_ChoiceSymmetricQLearning_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [0.05, -10, 0.05, -2, 0.5, -2];
UpperBound = [0.65,  20, 0.45,  1,   1,  2];

CalculateMLE = @(Parameters) ChoiceSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded);

Model = struct();
Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.MinNegLogDataLikelihood = Inf;

for iInitialCond = 1:20
    InitialParameters...
        = LowerBound + rand * (UpperBound - LowerBound);
    
    try
        [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
            fmincon(CalculateMLE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
    catch
        disp('Error: fail to run model');
        EstimatedParameters = [];
        MinNegLogDataLikelihood = nan;
    end
    
    if Model.MinNegLogDataLikelihood > MinNegLogDataLikelihood
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