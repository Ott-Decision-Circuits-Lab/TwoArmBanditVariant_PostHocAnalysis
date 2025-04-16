function [NegLogDataLikelihood, Values] = ChoiceSymmetricRLearning(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

LearningRate = Parameters(1); % alpha
RLearningRate = Parameters(2); % alpha_r
InverseTemperature = Parameters(3); % beta
ForgettingRate = Parameters(4); % gamma
ChoiceStickiness = Parameters(5); % phi
ChoiceForgettingRate = Parameters(6); % gamma_c
Bias = 0; %Parameters(7);

% Q-value initialisation
RValue = 0.25;
LeftValue = 0.25;
RightValue = 0;
ChoiceMemory = 0;
NegLogDataLikelihood = 0;

% Likelihood iteration
for iTrial = 1:nTrials-1 % <- -1 as raw ChoiceLeft has one pre-allocated nan
    LogOdds = (8 + InverseTemperature * RValue(iTrial)) * (LeftValue(iTrial) - RightValue(iTrial)) +...
              + ChoiceStickiness * ChoiceMemory(iTrial) + Bias;
    
    RValue(iTrial + 1) = (1 - ForgettingRate) * RValue(iTrial);
    LeftValue(iTrial + 1) = (1 - ForgettingRate) * LeftValue(iTrial);
    RightValue(iTrial + 1) = (1 - ForgettingRate) * RightValue(iTrial);
    ChoiceMemory(iTrial + 1) = (1 - ChoiceForgettingRate) * ChoiceMemory(iTrial);
    
    if Rewarded(iTrial) == 1
        RDelta = 1 - RValue(iTrial);
    else
        RDelta = 0 - RValue(iTrial);
    end

    RValue(iTrial + 1) = RValue(iTrial + 1) + RLearningRate * RDelta;

    if ChoiceLeft(iTrial) == 1
        ChoiceLogOdds = LogOdds;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        LeftValue(iTrial + 1) = LeftValue(iTrial + 1) + LearningRate * (RDelta - LeftValue(iTrial));

        ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial + 1) + ChoiceForgettingRate;
        
    elseif ChoiceLeft(iTrial) == 0
        ChoiceLogOdds = -LogOdds;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        RightValue(iTrial + 1) = RightValue(iTrial + 1) + LearningRate * (RDelta - RightValue(iTrial));

        ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial + 1) - ChoiceForgettingRate;
        
    else
        NegLogDataLikelihood = NegLogDataLikelihood;
        
        % Un-do the forgetting (?)
        % LeftValue(iTrial + 1) = LeftValue(iTrial);
        % RightValue(iTrial + 1) = RightValue(iTrial);
        % ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial); 
    end
end % end for-loop

Values.RValue = RValue;
Values.LeftValue = LeftValue;
Values.RightValue = RightValue;
Values.ChoiceMemory = ChoiceMemory;

end % end function