function [LogPosterior, Values]...
    = ChoiceBeliefStateHMM_MAP(Estimands, SessionParameters, SessionData, Prior, SessionExpectedStateInfo)
%{
calculate full ECLL even though this only requires optimizing
strategy-related parameters
%}

%% extract info
nStates = Prior.nStates;

% reshape and convert Parameters
StrategyWeight = Estimands(1:end-4);
StrategyWeight = reshape(StrategyWeight, 4, nStates)';
StrategyParameter = Estimands(end-3:end);
StrategyParameter = 1 ./ (1 + exp(-StrategyParameter));

% extract SessionData
TrialData = SessionData.Custom.TrialData;
ChoiceLeft = TrialData.ChoiceLeft;
Rewarded = TrialData.Rewarded;

% filter data (NoChoice -> No Emission Prob)
Valid = ~isnan(ChoiceLeft) & ~isnan(Rewarded);
nTrials = sum(Valid);
ChoiceLeft = ChoiceLeft(Valid);
Rewarded = Rewarded(Valid);

%% calculate data log(likelihood)
Gamma = SessionExpectedStateInfo.Gamma;
DataEmissionProb = nan(nStates, nTrials);
for iState = 1:nStates
    Parameters = [StrategyWeight(iState, :), StrategyParameter];
    Parameters = Parameters([5, 1, 6, 2, 7, 3, 8, 4]);
    [~, Values] = ChoiceBeliefState(Parameters, nTrials, ChoiceLeft, Rewarded);

    ChoiceLeftProb = Values.ChoiceLeftProb;
    DataEmissionProb(iState, :)...
        = ChoiceLeftProb .* ChoiceLeft...
            + (1 - ChoiceLeftProb) .* (1 - ChoiceLeft);
end

LogDataLikelihood = sum(Gamma .* log(DataEmissionProb), 'all');

%% calculate state log(likelihood)
InitialStateProbability = SessionParameters.InitialStateProbability;
StateTransitionMatrix = SessionParameters.StateTransitionMatrix;
Zeta = SessionExpectedStateInfo.Zeta;

LogStateLikelihood...
    = sum(Gamma(:, 1) .* log(InitialStateProbability), 'all')...
        + sum(Zeta .* permute(repmat(log(StateTransitionMatrix), 1, 1, nTrials - 1), [1, 3, 2]), 'all');

%% calculate log prior
LogParameterLikelihood = 0;

% InitialStateProbability
InitialStateProbabilityAlpha = Prior.InitialStateProbabilityAlpha; % nStates
InitialStateProbabilityPDF = dirpdf(InitialStateProbability, InitialStateProbabilityAlpha);
LogParameterLikelihood...
    = LogParameterLikelihood...
        + log(InitialStateProbabilityPDF);

% StateTransitionMatrix
StateTransitionMatrixAlpha = Prior.StateTransitionMatrixAlpha; % nStates x nStates
for iState = 1:nStates
    StateTransitionMatrixPDF = dirpdf(StateTransitionMatrix(iState, :), StateTransitionMatrixAlpha(iState, :));
    LogParameterLikelihood...
        = LogParameterLikelihood...
            + log(StateTransitionMatrixPDF);

end

% StrategyWeight
StrategyWeightMu = Prior.StrategyWeightMu; % nStates * 4
StrategyWeightSigma = Prior.StrategyWeightSigma; % nStates * 4
for iState = 1:nStates
    StrategyWeightPDF = pdf('Normal', StrategyWeight(iState, :), StrategyWeightMu, StrategyWeightSigma);
    for jState = 1:nStates
        LogParameterLikelihood...
            = LogParameterLikelihood...
                + log(StrategyWeightPDF(jState));
    end
end

% StrategyParameter
StrategyParameterAlpha = Prior.StrategyParameterAlpha; % 4
StrategyParameterBeta = Prior.StrategyParameterBeta; % 4
StrategyParameterPDF = pdf('Beta', StrategyParameter, StrategyParameterAlpha', StrategyParameterBeta');
for iPara = 1:4
    LogParameterLikelihood...
        = LogParameterLikelihood...
            + log(StrategyParameterPDF(iPara));
end

%% calculate log posterior & gradients
LogPosterior = LogDataLikelihood +...
                 LogStateLikelihood +...
                 LogParameterLikelihood;

end % end function