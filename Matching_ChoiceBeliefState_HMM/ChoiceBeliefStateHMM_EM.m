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
    try
        [ExpectedStateInfo, Prior] = ChoiceBeliefStateHMM_E(Parameters, DataHolder, Prior);
    catch
        disp('Error: problem at E step')
        error('')
    end

    %% M-step
    %{
    maximizing ECLL by optimizing strategy-related parameters and
    prior based on the expected state(s) they are in using another EM
    %}
    try
        [NewECLL, Parameters, Prior]...
            = ChoiceBeliefStateHMM_M(Parameters, DataHolder, Prior, ExpectedStateInfo);
    catch
        disp('Error: problem at M step')
        error('')
    end

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

Model.ECLL = ECLL;
Model.EstimatedParameters = Parameters;
Model.Posterior = Prior;
Model.ExpectedStateInfo = ExpectedStateInfo;

end