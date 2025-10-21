function [NegLogDataLikelihood, Values] = ChoiceRegressionQHMM2(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

InitialStateDistribution = Parameters(1:2); % pi_1, pi_2
StateTransitionMatrix = Parameters(3:6); % A_{i,j}
GLMCoeff = Parameters(7:14); % w_{k, 1:4}

% Q-value initialisation
ExpectedChoiceLeftLogOdds = 0; %E(P(c=L))
NegLogDataLikelihood = 0;

% Likelihood iteration
for iTrial = 1:nTrials % <- -1 as raw ChoiceLeft has one pre-allocated nan
    if iTrial > 1
        LastChoice = CurrentChoice;
        LastRewarded = CurrentRewarded;
    else
        LastChoice = 0;
        LastRewarded = 0;
    end

    if ChoiceLeft(iTrial) == 1
        CurrentChoice = 1;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ExpectedChoiceLeftLogOdds(iTrial))));

    elseif ChoiceLeft(iTrial) == 0
        CurrentChoice = -1;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(ExpectedChoiceLeftLogOdds(iTrial))));

    else
        CurrentChoice = 0;
        NegLogDataLikelihood = NegLogDataLikelihood;

    end

    if Rewarded(iTrial) == 1
        CurrentRewarded = 1 * CurrentChoice;
    elseif Rewarded(iTrial) == 0
        CurrentRewarded = -1 * CurrentChoice;
    else
        CurrentRewarded = 0;
    end

    ExpectedChoiceLeftLogOdds(iTrial + 1) = Bias...
                                    + QCoeff * ExpectedChoiceLeftLogOdds(iTrial)...
                                    + Minus1ChoiceCoeff * CurrentChoice...
                                    + Minus1RewardCoeff * CurrentRewarded...
                                    + Minus2ChoiceCoeff * LastChoice...
                                    + Minus2RewardCoeff * LastRewarded;
    
end % end for-loop

Values.ChoiceLeftLogOdds = ExpectedChoiceLeftLogOdds(1:nTrials);

end % end function