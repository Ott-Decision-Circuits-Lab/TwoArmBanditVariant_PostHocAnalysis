function Prior = ChoiceBeliefStateHMM_InitializePrior(nStates, InitialGLMMedian, InitialGLMMAD)
Prior.nStates = nStates;

Prior.InitialStateProbabilityAlpha...
    = ones(nStates, 1) + normrnd(0, 0.1, nStates, 1);

Prior.StateTransitionMatrixAlpha...
    = eye(nStates) * 0.9 + normrnd(0.2, 0.1, nStates, nStates);

Prior.StrategyWeightMu...
    = repmat(InitialGLMMedian(2:2:8)', nStates, 1) + normrnd(0, 0.1, nStates, 4);

Prior.StrategyWeightSigma...
    = repmat(InitialGLMMAD(2:2:8)', nStates, 1) + abs(normrnd(0, 0.1, nStates, 4));

% using method-of-moment to generate alpha and beta
TempMedian = InitialGLMMedian(1:2:7);
TempMAD = InitialGLMMAD(1:2:7);
Prior.StrategyParameterAlpha...
    = TempMedian...
        .* (TempMedian .* (1 - TempMedian)...
        ./ TempMAD - 1);
Prior.StrategyParameterBeta...
    = (1 - TempMedian)...
        .* (TempMedian .* (1 - TempMedian)...
        ./ TempMAD - 1); 

for iPara = 1:4
    if Prior.StrategyParameterAlpha(iPara) < 1 || Prior.StrategyParameterBeta(iPara) < 1
        Prior.StrategyParameterAlpha(iPara) = 1 + normrnd(0, 0.1);
        Prior.StrategyParameterBeta(iPara) = 1 + normrnd(0, 0.1);
    end
end
Prior.StrategyParameterAlpha = Prior.StrategyParameterAlpha';
Prior.StrategyParameterBeta = Prior.StrategyParameterBeta';

Prior.StrategyWeightL1Lambda = ones(1, nStates);
end