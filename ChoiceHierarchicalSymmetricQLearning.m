function [NegLogDataLikelihood, Values] = ChoiceHierarchicalSymmetricQLearning(Parameters)

% Parameters import
% fixed variable
global nTrials ChoiceLeft Rewarded

LearningRate = Parameters(1); % alpha
InverseTemperature = Parameters(2); % beta
ForgettingRate = Parameters(3); % gamma
ChoiceStickiness = Parameters(4); % phi
ChoiceForgettingRate = Parameters(5); % c_gamma
Bias = 0; %Parameters(6);

% Q-value initialisation
LeftValue = 0.25;
RightValue = 0.25;
ChoiceMemory = 0;
NegLogDataLikelihood = 0;

% Likelihood iteration
for iTrial = 1:nTrials-1 % <- -1 as raw ChoiceLeft has one pre-allocated nan
    LogOdds = InverseTemperature * (LeftValue(iTrial) - RightValue(iTrial)) +...
              + ChoiceStickiness * ChoiceMemory(iTrial) + Bias;
    
    LeftValue(iTrial + 1) = (1 - ForgettingRate) * LeftValue(iTrial);
    RightValue(iTrial + 1) = (1 - ForgettingRate) * RightValue(iTrial);
    ChoiceMemory(iTrial + 1) = (1 - ChoiceForgettingRate) * ChoiceMemory(iTrial);

    if ChoiceLeft(iTrial) == 1
        ChoiceLogOdds = LogOdds;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        if Rewarded(iTrial) == 1
            LeftValue(iTrial + 1) = LeftValue(iTrial + 1) + LearningRate * (1 - LeftValue(iTrial));
        else
            LeftValue(iTrial + 1) = LeftValue(iTrial + 1) + LearningRate * (0 - LeftValue(iTrial));
        end

        ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial + 1) + ChoiceForgettingRate;
        
    elseif ChoiceLeft(iTrial) == 0
        ChoiceLogOdds = -LogOdds;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        if Rewarded(iTrial) == 1
            RightValue(iTrial + 1) = RightValue(iTrial + 1) + LearningRate * (1 - RightValue(iTrial));
        else
            RightValue(iTrial + 1) = RightValue(iTrial + 1) + LearningRate * (0 - RightValue(iTrial));
        end

        ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial + 1) - ChoiceForgettingRate;
        
    else
        NegLogDataLikelihood = NegLogDataLikelihood;
        
        % Un-do the forgetting (?)
        % LeftValue(iTrial + 1) = LeftValue(iTrial);
        % RightValue(iTrial + 1) = RightValue(iTrial);
        % ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial); 
    end
end % end for-loop

Values.LeftValue = LeftValue;
Values.RightValue = RightValue;
Values.ChoiceMemory = ChoiceMemory;

end % end function