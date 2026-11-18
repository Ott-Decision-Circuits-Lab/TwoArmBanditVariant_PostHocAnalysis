function [ECLL, Parameters, Prior]...
    = ChoiceBeliefStateHMM_M(Parameters, DataHolder, Prior, ExpectedStateInfo)
%{
maximize ECLL by optimizing Session strategy-related parameters
%}

%% update state-related prior using method of moments
% these priors cannot be optimized as the session parameters are from the
% E-step, thus here they are just summarized statistically using
% method-of-moment
nStates = Prior.nStates;
nSessions = length(DataHolder);
AllInitialStateProbability = zeros(nStates, nSessions);
AllStateTransitionMatrix = zeros(nStates, nStates, nSessions);
for iSession = 1:nSessions
    InitialStateProbability = Parameters{iSession}.InitialStateProbability;
    StateTransitionMatrix = Parameters{iSession}.StateTransitionMatrix;
    
    AllInitialStateProbability(:, iSession) = InitialStateProbability;
    AllStateTransitionMatrix(:, :, iSession) = StateTransitionMatrix;
end

TempMean = mean(AllInitialStateProbability, 2);
TempVar = var(AllInitialStateProbability, 0, 2);
Prior.InitialStateProbabilityAlpha = TempMean .* (TempMean .* (1 - TempMean) ./ TempVar - 1);
for iState = 1:nStates
    TempMean = mean(AllStateTransitionMatrix(iState, :, :), 3);
    TempVar = var(AllStateTransitionMatrix(iState, :, :), 0, 3);
    Prior.StateTransitionMatrixAlpha(iState, :) = TempMean .* (TempMean .* (1 - TempMean) ./ TempVar - 1);
end

%% maximize MAP by optimizing parameter estimates
% get initial values of estimate (which are the strategy-related prior)
StrategyWeightMu = Prior.StrategyWeightMu; % nStates * 4
StrategyWeightSigma = Prior.StrategyWeightSigma; % nStates * 4
StrategyParameterAlpha = Prior.StrategyParameterAlpha; % 4
StrategyParameterBeta = Prior.StrategyParameterBeta; % 4
StrategyWeightL1Lambda = Prior.StrategyWeightL1Lambda; % nStates

% reparameterization and linearization
StrategyWeightSigma = log(exp(StrategyWeightSigma) - 1);
StrategyParameterAlpha = log(exp(StrategyParameterAlpha) - 1);
StrategyParameterBeta = log(exp(StrategyParameterBeta) - 1);
StrategyWeightL1Lambda = log(exp(StrategyWeightL1Lambda) - 1);

InitialEstimate...
    = [reshape(StrategyWeightMu, 1, []),...
       reshape(StrategyWeightSigma, 1, []),...
       StrategyParameterAlpha,...
       StrategyParameterBeta,...
       StrategyWeightL1Lambda];

% set up optimization
UpperBound = ones(1, 9*nStates + 8) * 40;
LowerBound = ones(1, 9*nStates + 8) * -40;
Options.Display = 'off';

LogPosterior = @(Estimand)...
    ChoiceBeliefStateHMM_CalculateECLL(Estimand, Parameters, DataHolder, Prior, ExpectedStateInfo);

FinalEstimate = fmincon(LogPosterior, InitialEstimate, [], [], [], [], LowerBound, UpperBound); %, [], Options);

% get best estimate of Prior and Parameter
[ECLL, Parameters, Prior]...
    = ChoiceBeliefStateHMM_CalculateECLL(FinalEstimate, Parameters, DataHolder, Prior, ExpectedStateInfo);

end

%{
nSessions = length(DataHolder);
for iSession = 1:nSessions
    %% extract session info
    SessionParameters = Parameters{iSession};
    SessionData = DataHolder{iSession};
    SessionExpectedStateInfo = ExpectedStateInfo{iSession};
    
    % linearize and reparameterize strategy-related SessionParameters
    StrategyWeight = SessionParameters.StrategyWeight;
    StrategyParameter = SessionParameters.StrategyParameter;

    StrategyParameter = log(StrategyParameter ./ (1 - StrategyParameter));
    InitialParameters = [reshape(StrategyWeight', 1, []), StrategyParameter]; % nStates * 4 + 4

    %% formulate optimization problem
    LogPosteriorPDF = @(Estimands)...
        ChoiceBeliefStateHMM_MAP(Estimands, SessionParameters, SessionData,...
                                 Prior, SessionExpectedStateInfo);
    %{
    Sampler = hmcSampler(LogPosteriorPDF, InitialParameters,...
                         'CheckGradient', false,...
                         'UseNumericalGradient', true);

    % MAP estimated by L-BFGS
    [MAPParameters, FitInfo] = estimateMAP(Sampler,...
                                           'VerbosityLevel', 0);

    NewECLL = FitInfo.Objective(end);
    %}
    % if any(isnan(MAPParameters)) || isnan(NewECLL) || isinf(NewECLL)
    UpperBound = ones(size(InitialParameters)) * 20;
    LowerBound = ones(size(InitialParameters)) * -20;
    Options.Display = 'off';
    LogPosteriorPDF = @(Estimands)...
        -ChoiceBeliefStateHMM_MAP(Estimands, SessionParameters, SessionData,...
                                 Prior, SessionExpectedStateInfo);

    [MAPParameters, NewECLL] = fmincon(LogPosteriorPDF, InitialParameters, [], [], [], [], LowerBound, UpperBound, [], Options);
    % end

    %% update SessionParameters
    StrategyWeight = MAPParameters(1:end-4);
    StrategyWeight = reshape(StrategyWeight, 4, [])';
    StrategyParameter = MAPParameters(end-3:end);
    StrategyParameter = 1 ./ (1 + exp(-StrategyParameter));

    Parameters{iSession}.StrategyWeight = StrategyWeight;
    Parameters{iSession}.StrategyParameter = StrategyParameter;

    %% eject ECLL
    ECLL = ECLL + NewECLL;

end
%}