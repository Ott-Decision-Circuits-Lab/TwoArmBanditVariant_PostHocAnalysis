function Model = Matching_ConcatMS_MLE_ChoiceRegressionQHMM2_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

% Parametric estimation
% LowerBound = [0.10, 6, 0.05, -1.5, 0.5, -0.5];
% UpperBound = [0.45, 10, 0.25, 0.5, 1, 0.5];

% Free parameters
Model = struct();
Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.MinNegLogDataLikelihood = 0;

for iInitialCond = 1:10    
    StateInitialProbability = randn(1, 2); % pi_1, pi_2
    StateTransitionProbability = randn(2, 2);
    GLMCoeff = randn(2, 4) * 5; % 2 states x 4 weights (QCoeff, reward, choice, bias)
    
    InitialParameters = [StateInitialProbability, reshape(StateTransitionProbability, 1, []), reshape(GLMCoeff, 1, [])];
    
    MinNegLogDataLikelihood = 0;
    ImprovedNegLogDataLikelihood = 1;
    nIterations = 1;
    if ImprovedNegLogDataLikelihood > 10^-5 && nIterations < 300
        %% calculate new expected log(odds)/Q
        
        %% Baum-Welch algorithm (E-step)
        % forward
        % backward
        % gamma
        % xi
        %% Baum-Welch algorihm (M-step)
        % state initial probability
        % state transition probability
        % GLM coeff fit
        % 
        NewMinNegLogDataLikelihood = abc;
        if NewMinNegLogDataLikelihood - MinNegLogDataLikelihood > 0
            ImprovedNegLogDataLikelihood = NewMinNegLogDataLikelihood - MinNegLogDataLikelihood;
            MinNegLogDataLikelihood = NewMinNegLogDataLikelihood;
        end
        
        nIterations = nIterations + 1;
    end
    CalculateMLE = @(Parameters) ConcatChoiceRegressionQHMM2(Parameters, nTrials, ChoiceLeft, Rewarded);
    
    try
        [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
            fmincon(CalculateMLE, InitialParameters, [], [], [], [], [], []); %LowerBound, UpperBound);
    catch
        disp('Error: fail to run model');
        EstimatedParameters = [];
        MinNegLogDataLikelihood = nan;
    end
    
    if Model.MinNegLogDataLikelihood < MinNegLogDataLikelihood
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