function Model = ChoiceBeliefStateHMM_EM(DataHolder, Prior)
%{
Basically Baum-Welch algorithm
Data cannot have NoChoice
%}

%% initialise Session-level parameters based on prior
nSessions = length(DataHolder);
Parameters = ChoiceBeliefStateHMM_InitializeParameters(Prior, nSessions);

%%
Model.InitialParameters = Parameters;
Model.Prior = Prior;
ImprovedLoss = Inf;
ECLL = Inf; % expected complete data log likelihood
while ImprovedLoss > 1e-5 && length(ECLL) < 301
    %% E-step
    %{
    using forward-backward algorithm to calculate the expected value of
    state and state transition, and update state-related parameters and
    prior
    %}
    [ExpectedStateInfo, Prior] = ChoiceBeliefStateHMM_E(Parameters, DataHolder, Prior);

    %% M-step
    %{
    maximizing ECLL by optimizing strategy-related parameters and
    prior based on the expected state(s) they are in using another EM
    %}
    [NewECLL, Parameters, Prior]...
        = ChoiceBeliefStateHMM_M(Parameters, DataHolder, Prior, ExpectedStateInfo);
    
    %% update loop
    ImprovedLoss = ECLL(end) - NewECLL;
    ECLL(end + 1) = NewECLL;

end

%% report
ECLL = ECLL(2:end);
if ImprovedLoss < 1e-5
    Flag = 1;
elseif length(ECLL) == 300
    Flag = 0;
end

Posterior = Model.EstimatedParameters;
ExpectedStateInfo
Parameters

end