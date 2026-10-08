function Parameters = ChoiceBeliefStateHMM_InitializeParameters(Prior, nSessions)
% initialize the first estimate of session-level parameters based on Prior
%% extract prior
nStates = Prior.nStates;
InitialStateProbabilityAlpha = Prior.InitialStateProbabilityAlpha; % nStates
StateTransitionMatrixAlpha = Prior.StateTransitionMatrixAlpha; % nStates x nStates
StrategyWeightMu = Prior.StrategyWeightMu; % nStates * 4
StrategyWeightSigma = Prior.StrategyWeightSigma; % nStates * 4
StrategyParameterAlpha = Prior.StrategyParameterAlpha; % 4
StrategyParameterBeta = Prior.StrategyParameterBeta; % 4

%% initialize
Parameters = cell(1, nSessions);
for iSession = 1:nSessions
    TempStruct = struct();

    TempStruct.InitialStateProbability...
        = dirrnd(InitialStateProbabilityAlpha, 1)'; % nState x 1 x nSessions
    
    TempStruct.StateTransitionMatrix = nan(nStates, nStates);
    TempStruct.StrategyWeight = nan(nStates, 4);
    for iState = 1:nStates
        TempStruct.StateTransitionMatrix(iState, :)...
            = dirrnd(StateTransitionMatrixAlpha(iState, :), 1); % nState_i x nState_j x nSessions
        
        TempStruct.StrategyWeight(iState, :)...
            = normrnd(StrategyWeightMu(iState, :), StrategyWeightSigma(iState, :)); % nState x 4
    end
    
    TempStruct.StrategyParameter...
	    = betarnd(StrategyParameterAlpha, StrategyParameterBeta); % 4

    Parameters{iSession} = TempStruct;
end
end

