function [NewECLL, Parameters, Prior]...
        = ChoiceBeliefStateHMM_M(Parameters, DataHolder, Prior, ExpectedStateInfo)
%{
EM within M-Step for optimizing both prior and session estimates for
strategy-related parameters
%}

ImprovedLoss = Inf;
ECLL = Inf;
while ImprovedLoss > 1e-2 && length(ECLL) < 301
    %% E-Step
    % calculate expected value of strategy-related prior
    try
        Prior = ChoiceBeliefStateHMM_M_E(Parameters, DataHolder, Prior);
    catch
        disp('Error: problem at M_E step')
        error('')
    end

    %% M-Step
    % maximize ECLL by optimizing strategy-related parameters
    try
        [NewECLL, Parameters]...
            = ChoiceBeliefStateHMM_M_M(Parameters, DataHolder, Prior, ExpectedStateInfo);
    catch
        disp('Error: problem at M_M step')
        error('')
    end

    ImprovedLoss = ECLL(end) - NewECLL;
    ECLL(end+1) = NewECLL;
    
end

end