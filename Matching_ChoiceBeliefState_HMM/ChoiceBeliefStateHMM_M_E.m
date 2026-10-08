function Prior = ChoiceBeliefStateHMM_M_E(Parameters, DataHolder, Prior)
%% extract parameters
nSessions = length(DataHolder);
nStates = Prior.nStates;

% concatenate strategy-related parameters
AllStrategyWeight = zeros(nStates, 4, nSessions);
AllStrategyParameter = zeros(4, nSessions);
for iSession = 1:nSessions
    AllStrategyWeight(:, :, iSession) = Parameters{iSession}.StrategyWeight;
    AllStrategyParameter(:, iSession) = Parameters{iSession}.StrategyParameter;
end

%% calculate expected value of the prior parameters
Prior.StrategyWeightMu = mean(AllStrategyWeight, 3);
Prior.StrategyWeightSigma = std(AllStrategyWeight, 0, 3);

% estimate using method-of-moment
TempMean = mean(AllStrategyParameter, 2);
TempStd = var(AllStrategyParameter, 0, 2);
Prior.StrategyParameterAlpha...
    = TempMean...
        .* (TempMean .* (1 - TempMean)...
        ./ TempStd - 1);
Prior.StrategyParameterBeta...
    = (1 - TempMean)...
        .* (TempMean .* (1 - TempMean)...
        ./ TempStd - 1);

end

