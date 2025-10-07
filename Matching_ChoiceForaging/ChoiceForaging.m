function [NegLogDataLikelihood, Values] = ChoiceForaging(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

LearningRate = Parameters(1); % alpha
InverseTemperature = Parameters(2); % beta
Threshold = Parameters(3); % theta
ForgettingRate = Parameters(4); % gamma

% Foraging-value initialisation
ExploitingValue = Threshold;
NegLogDataLikelihood = 0;
LastValidChoice = NaN;

% Likelihood iteration
for iTrial = 1:nTrials
    ExploitingLogOdds = InverseTemperature * (ExploitingValue(iTrial) - Threshold);

    ExploitingValue(iTrial + 1) = (1 - ForgettingRate) * ExploitingValue(iTrial);

    if Rewarded(iTrial) == 1
        ExploitingValue(iTrial + 1) = ExploitingValue(iTrial + 1) + LearningRate * (1 - ExploitingValue(iTrial));
    elseif Rewarded(iTrial) == 0
        ExploitingValue(iTrial + 1) = ExploitingValue(iTrial + 1) + LearningRate * (0 - ExploitingValue(iTrial));
    end

    if isnan(LastValidChoice) % initial trials without choice
        Exploited(iTrial) = nan;
        ChoiceLogOdds(iTrial) = nan;
        ChoiceLeftLogOdds(iTrial) = nan;

        if ~isnan(ChoiceLeft(iTrial)) % first trial with choice
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
        ChoiceLogOdds(iTrial) = -ExploitingLogOdds;
        if LastValidChoice == 1
            ChoiceLeftLogOdds(iTrial) = ExploitingLogOdds;
        elseif LastValidChoice == 0
            ChoiceLeftLogOdds(iTrial) = -ExploitingLogOdds;
        end
        
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds(iTrial))));
        LastValidChoice = ChoiceLeft(iTrial);

    else % subsequent trials without choice
        Exploited(iTrial) = nan;
        ChoiceLogOdds(iTrial) = nan;
        ChoiceLeftLogOdds(iTrial) = nan;
        NegLogDataLikelihood = NegLogDataLikelihood;
        
    end

end % end for-loop

Values.ExploitingValue = ExploitingValue;
Values.Exploited = Exploited;
Values.ChoiceExploitLogOdds = ChoiceLogOdds;
Values.ChoiceLeftLogOdds = ChoiceLeftLogOdds;

end % end function