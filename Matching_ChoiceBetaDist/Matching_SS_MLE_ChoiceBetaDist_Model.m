function Model = Matching_SS_MLE_ChoiceBetaDist_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
LowerBound = [0, 0, 0, 0, 0];
UpperBound = [10, 1, 1, 1, 1];

% Free parameters
CalculateMLE = @(Parameters) ChoiceBetaDist(Parameters, nTrials, ChoiceLeft, Rewarded);

Model = struct();
Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.MinNegLogDataLikelihood = 0;

for iInitialCond = 1:10
    % 20250708 tested with simulation that works well as initial parameters
    AlphaLearningStepSize = rand * 5; % \delta_\alpha
    AlphaForgettingRate = rand / 5; % \gamma_\alpha
    RLearningRate = rand / 5; % \alpha_R
    RForgettingRate = rand / 5; % \gamma_R
    Bias = rand; %
    
    InitialParameters = [AlphaLearningStepSize, AlphaForgettingRate, RLearningRate, RForgettingRate, Bias];

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