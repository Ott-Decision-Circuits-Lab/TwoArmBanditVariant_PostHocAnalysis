function [NegLogDataLikelihood, Values] = ChoiceForagingRewardRate(Parameters, nTrials, ChoiceLeft, Rewarded, TrialTimeDuration)
% Parameters import
% fixed variable

TrialTimeDuration = [TrialTimeDuration, 20]; % add a dummy timeDuration

LearningRate = Parameters(1); % alpha
InverseTemperature = Parameters(2); % beta
BackgroundRewardRateWeight = Parameters(3); % phi
ForgettingRate = Parameters(4); % gamma_V
BackgroundRewardRateDecayConstant = Parameters(5); % tau
BasalBackgroundRewardRate = Parameters(6); % theta

% Foraging-value initialisation
ExploitingValue = BackgroundRewardRateWeight;
NegLogDataLikelihood = 0;
LastValidChoice = NaN;
BackgroundRewardRate = 0;

% Likelihood iteration
for iTrial = 1:nTrials
    ExploitingLogOdds = InverseTemperature * (ExploitingValue(iTrial) - BackgroundRewardRateWeight * (BackgroundRewardRate(iTrial) + BasalBackgroundRewardRate));

    ExploitingValue(iTrial + 1) = (1 - ForgettingRate) * ExploitingValue(iTrial);
    if Rewarded(iTrial) == 1
        ExploitingValue(iTrial + 1) = ExploitingValue(iTrial + 1) + LearningRate * (1 - ExploitingValue(iTrial));
    elseif Rewarded(iTrial) == 0
        ExploitingValue(iTrial + 1) = ExploitingValue(iTrial + 1) + LearningRate * (0 - ExploitingValue(iTrial));
    end

    BackgroundRewardRate(iTrial + 1) = (1 - BackgroundRewardRateDecayConstant) * BackgroundRewardRate(iTrial) +...
                                       Rewarded(iTrial) == 1;

    if isnan(LastValidChoice) % initial trials without choice
        Exploited(iTrial) = nan;
        ChoiceLogOdds(iTrial) = nan;
        ChoiceLeftLogOdds(iTrial) = nan;

        if ~isnan(ChoiceLeft(iTrial)) % first trial with choice
            Exploited(iTrial) = 0;
            ChoiceLogOdds(iTrial) = 0;
            ChoiceLeftLogOdds(iTrial) = 0;
            NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds(iTrial))));
            LastValidChoice = ChoiceLeft(iTrial);

        end

    elseif abs(ChoiceLeft(iTrial) - LastValidChoice) == 0 % if exploited
        Exploited(iTrial) = 1;
        ChoiceLogOdds(iTrial) = ExploitingLogOdds;

        if LastValidChoice == 1
            ChoiceLeftLogOdds(iTrial) = ExploitingLogOdds;
        elseif LastValidChoice == 0
            ChoiceLeftLogOdds(iTrial) = -ExploitingLogOdds;
        end
        
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds(iTrial))));
        LastValidChoice = ChoiceLeft(iTrial);
        
    elseif abs(ChoiceLeft(iTrial) - LastValidChoice) == 1
        Exploited(iTrial) = 0;
        ChoiceLogOdds(iTrial) = ExploitingLogOdds;
        if LastValidChoice == 1
            ChoiceLeftLogOdds(iTrial) = ExploitingLogOdds;
        elseif LastValidChoice == 0
            ChoiceLeftLogOdds(iTrial) = -ExploitingLogOdds;
        end
        
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(ChoiceLogOdds(iTrial))));
        LastValidChoice = ChoiceLeft(iTrial);

    else % subsequent trials without choice
        Exploited(iTrial) = nan;
        ChoiceLogOdds(iTrial) = nan;
        ChoiceLeftLogOdds(iTrial) = nan;
        NegLogDataLikelihood = NegLogDataLikelihood;
        
    end

end % end for-loop

Values.ExploitingValue = ExploitingValue(1:nTrials);
Values.Exploited = Exploited(1:nTrials);
Values.ChoiceExploitLogOdds = ChoiceLogOdds(1:nTrials);
Values.ChoiceLeftLogOdds = ChoiceLeftLogOdds(1:nTrials);
Values.BackgroundRewardRate = BackgroundRewardRate(1:nTrials);

end % end function