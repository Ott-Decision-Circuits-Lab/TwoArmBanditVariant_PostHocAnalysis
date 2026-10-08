function [NewECLL, Parameters, Prior]...
        = ChoiceBeliefStateHMM_M(Parameters, DataHolder, Prior, ExpectedStateInfo)
%{
EM within M-Step for optimizing both prior and session estimates for
strategy-related parameters
%}

ImprovedLoss = Inf;
ECLL = Inf;
while ImprovedLoss > 1e-5 && length(ECLL) < 301
    %% E-Step
    % calculate expected value of strategy-related prior
    Prior = ChoiceBeliefStateHMM_M_E(Parameters, DataHolder, Prior);
    
    %% M-Step
    % maximize ECLL by optimizing strategy-related parameters
    [NewECLL, Parameters]...
        = ChoiceBeliefStateHMM_M_M(Parameters, DataHolder, Prior, ExpectedStateInfo);
    
    ImprovedLoss = ECLL(end) - NewECLL;
    ECLL(end+1) = NewECLL;
    
end

end