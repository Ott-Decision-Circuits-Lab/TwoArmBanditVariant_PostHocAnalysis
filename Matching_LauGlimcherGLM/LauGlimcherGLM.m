function [NegLogDataLikelihood, Values] = LauGlimcherGLM(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

Bias = Parameters(1); % beta_0
RewardCoeff = Parameters(2:6); % beta_r_i
ChoiceCoeff = Parameters(7:11); % beta_c_i

% value initialisation
NegLogDataLikelihood = 0;

Choice = (ChoiceLeft == 1) - (ChoiceLeft == 0);
Reward = (Rewarded == 1) .* Choice;

Choices = zeros(5, nTrials);
Rewards = zeros(5, nTrials);
for i = 1:5
    Choices(i, :) = [zeros(1, i), Choice(1:end-i)];
    Rewards(i, :) = [zeros(1, i), Reward(1:end-i)];
end

ChoiceLeftLogOdds = Bias + RewardCoeff * Rewards + ChoiceCoeff * Choices;

% Likelihood iteration
for iTrial = 1:nTrials % <- -1 as raw ChoiceLeft has one pre-allocated nan
    if ChoiceLeft(iTrial) == 1
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLeftLogOdds(iTrial))));

    elseif ChoiceLeft(iTrial) == 0
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(ChoiceLeftLogOdds(iTrial))));

    else
        NegLogDataLikelihood = NegLogDataLikelihood;

    end
end % end for-loop

Values.ChoiceLeftLogOdds = ChoiceLeftLogOdds(1:nTrials);

end % end function