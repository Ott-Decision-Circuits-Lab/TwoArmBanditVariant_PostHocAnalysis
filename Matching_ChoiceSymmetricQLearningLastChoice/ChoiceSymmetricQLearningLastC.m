function [NegLogDataLikelihood, Values] = ChoiceSymmetricQLearningLastC(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

LearningRate = Parameters(1); % alpha
InverseTemperature = Parameters(2); % beta
ForgettingRate = Parameters(3); % gamma
ChoiceStickiness = Parameters(4); % phi
Bias = Parameters(5);

% Q-value initialisation
LeftValue = 0;
RightValue = 0;
ChoiceMemory = 0;
NegLogDataLikelihood = 0;

% Likelihood iteration
for iTrial = 1:nTrials % <- -1 as raw ChoiceLeft has one pre-allocated nan
    LogOdds(iTrial) = InverseTemperature * (LeftValue(iTrial) - RightValue(iTrial)) +...
                      + ChoiceStickiness * ChoiceMemory(iTrial) + Bias;
    
    LeftValue(iTrial + 1) = (1 - ForgettingRate) * LeftValue(iTrial);
    RightValue(iTrial + 1) = (1 - ForgettingRate) * RightValue(iTrial);
    ChoiceMemory(iTrial + 1) =  0;

    if ChoiceLeft(iTrial) == 1
        ChoiceLogOdds = LogOdds(iTrial);
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        if Rewarded(iTrial) == 1
            LeftValue(iTrial + 1) = LeftValue(iTrial + 1) + LearningRate * (1 - LeftValue(iTrial));
        else
            LeftValue(iTrial + 1) = LeftValue(iTrial + 1) + LearningRate * (0 - LeftValue(iTrial));
        end
        
        ChoiceMemory(iTrial + 1) =  1;

    elseif ChoiceLeft(iTrial) == 0
        ChoiceLogOdds = -LogOdds(iTrial);
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        if Rewarded(iTrial) == 1
            RightValue(iTrial + 1) = RightValue(iTrial + 1) + LearningRate * (1 - RightValue(iTrial));
        else
            RightValue(iTrial + 1) = RightValue(iTrial + 1) + LearningRate * (0 - RightValue(iTrial));
        end
        
        ChoiceMemory(iTrial + 1) =  -1;

    else
        NegLogDataLikelihood = NegLogDataLikelihood;
        
        % Un-do the forgetting (?)
        % LeftValue(iTrial + 1) = LeftValue(iTrial);
        % RightValue(iTrial + 1) = RightValue(iTrial);
        % ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial); 
    end
end % end for-loop

Values.LeftValue = LeftValue(1:nTrials);
Values.RightValue = RightValue(1:nTrials);
Values.ChoiceMemory = ChoiceMemory(1:nTrials);
Values.ChoiceLeftLogOdds = LogOdds(1:nTrials);

end % end function