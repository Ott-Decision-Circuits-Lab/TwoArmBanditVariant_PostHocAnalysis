function [LogPosterior, SessionParameters]...
    = ChoiceBeliefStateHMM_CalculateSessionECLL(Estimand, SessionParameters, SessionData, Prior, SessionExpectedStateInfo)
%{
calculate ECLL
%}

%% un-parameterize Estimand (strategy-related Prior)
nStates = Prior.nStates;
SessionParameters.StrategyWeight = Estimand(1:4*nStates);
SessionParameters.StrategyWeight = reshape(SessionParameters.StrategyWeight, nStates, []);
SessionParameters.StrategyParameter = Estimand(4*nStates + 1:4*nStates + 4);
SessionParameters.StrategyParameter = log(exp(SessionParameters.StrategyParameter) + 1);

%% extract info
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
    Parameters = [SessionParameters.StrategyWeight(iState, :), SessionParameters.StrategyParameter];
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

%% calculate log(parameter | prior)
LogParameterLikelihood = 0;

% InitialStateProbability
InitialStateProbabilityAlpha = Prior.InitialStateProbabilityAlpha; % nStates
InitialStateProbabilityPDF = dirpdf(InitialStateProbability, InitialStateProbabilityAlpha);
LogParameterLikelihood...
    = LogParameterLikelihood...
        + log(InitialStateProbabilityPDF);

% StateTransitionMatrix
StateTransitionMatrixAlpha = Prior.StateTransitionMatrixAlpha; % nStates x nStates
StateTransitionMatrixPDF = dirpdf(StateTransitionMatrix, StateTransitionMatrixAlpha);
LogParameterLikelihood...
    = LogParameterLikelihood...
        + sum(log(StateTransitionMatrixPDF));

% StrategyWeight
StrategyWeightMu = Prior.StrategyWeightMu; % nStates * 4
StrategyWeightSigma = Prior.StrategyWeightSigma; % nStates * 4
StrategyWeightPDF = pdf('Normal', SessionParameters.StrategyWeight, StrategyWeightMu, StrategyWeightSigma);
LogParameterLikelihood...
    = LogParameterLikelihood...
        + sum(log(StrategyWeightPDF), 'all');

% StrategyParameter
StrategyParameterAlpha = Prior.StrategyParameterAlpha; % 4
StrategyParameterBeta = Prior.StrategyParameterBeta; % 4
StrategyParameterPDF = pdf('Beta', SessionParameters.StrategyParameter, StrategyParameterAlpha, StrategyParameterBeta);
LogParameterLikelihood...
    = LogParameterLikelihood...
        + sum(log(StrategyParameterPDF), 'all');

%% calculate penalty for StrategyWeight using L1 regularization (Laplace)
LogPenalty = 0;
for iState = 1:nStates
    StrategyWeightL1Lambda = Prior.StrategyWeightL1Lambda(iState);
    PenaltyLogPDF...
        = sum(abs(SessionParameters.StrategyWeight(iState, :))) / StrategyWeightL1Lambda...
            - 4 * log(2 * StrategyWeightL1Lambda);

    LogPenalty...
        = LogPenalty + PenaltyLogPDF;
end

%% calculate log posterior & gradients
LogPosterior = LogDataLikelihood +...
                 LogStateLikelihood +...
                 LogParameterLikelihood +...
                 LogPenalty;

end % end function