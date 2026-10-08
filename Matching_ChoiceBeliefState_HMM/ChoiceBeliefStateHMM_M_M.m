function [ECLL, Parameters]...
    = ChoiceBeliefStateHMM_M_M(Parameters, DataHolder, Prior, ExpectedStateInfo)
%{
maximize ECLL by optimizing Session strategy-related parameters
%}

%% maximize MAP by optimizing parameter estimates
ECLL = 0;
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
        UpperBound = ones(size(InitialParameters)) * 40;
        LowerBound = ones(size(InitialParameters)) * -40;
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