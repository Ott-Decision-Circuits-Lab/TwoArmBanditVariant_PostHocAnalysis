function [NegLogDataLikelihood, Values] = ChoiceRegressionQ(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

Bias = Parameters(1); % beta_0
QCoeff = Parameters(2); % gamma
Minus1ChoiceCoeff = Parameters(3); % beta_c1
Minus1RewardCoeff = Parameters(4); % beta_r1
Minus2ChoiceCoeff = Parameters(5); % beta_c2
Minus2RewardCoeff = Parameters(6); % beta_r2

% Q-value initialisation
ChoiceLeftLogOdds = 0;
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
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLeftLogOdds(iTrial))));

    elseif ChoiceLeft(iTrial) == 0
        CurrentChoice = -1;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(ChoiceLeftLogOdds(iTrial))));

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

    ChoiceLeftLogOdds(iTrial + 1) = Bias...
                                    + QCoeff * ChoiceLeftLogOdds(iTrial)...
                                    + Minus1ChoiceCoeff * CurrentChoice...
                                    + Minus1RewardCoeff * CurrentRewarded...
                                    + Minus2ChoiceCoeff * LastChoice...
                                    + Minus2RewardCoeff * LastRewarded;
    
end % end for-loop

Values.ChoiceLeftLogOdds = ChoiceLeftLogOdds(1:nTrials);

end % end function