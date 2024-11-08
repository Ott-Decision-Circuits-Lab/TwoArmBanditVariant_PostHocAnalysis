function NegLogDataLikelihood = ChoiceAsymmetricQLearning(Parameters)
global nTrials
global ChoiceLeft
global Rewarded

% Parameters import
PosLearningRate = Parameters(1); % alpha_pos
NegLearningRate = Parameters(2); % alpha_neg
InverseTemperature = Parameters(3); % beta
Bias = 0; %Parameters(4);

% Q-value initialisation
LeftValue = 0.5;
RightValue = 0.5;
NegLogDataLikelihood = 0;

% Likelihood iteration
for iTrial = 1:nTrials
    if ChoiceLeft(iTrial) == 1
        ChoiceLogOdds = InverseTemperature * (LeftValue(iTrial) - RightValue(iTrial)) + Bias;
        NegLogDataLikelihood = NegLogDataLikelihood - 1 ./ (1 + exp(-ChoiceLogOdds));
        
        if Rewarded(iTrial) == 1
            LeftValue(iTrial+1) = LeftValue(iTrial) + PosLearningRate * (1 - LeftValue(iTrial));
            RightValue(iTrial+1) = RightValue(iTrial);
        else
            LeftValue(iTrial+1) = LeftValue(iTrial) + NegLearningRate * (0 - LeftValue(iTrial));
            RightValue(iTrial+1) = RightValue(iTrial);
        end
    elseif ChoiceLeft(iTrial) == 0
        ChoiceLogOdds = InverseTemperature * (RightValue(iTrial) - LeftValue(iTrial)) + Bias;
        NegLogDataLikelihood = NegLogDataLikelihood - 1 ./ (1 + exp(-ChoiceLogOdds));
        
        if Rewarded(iTrial) == 1
            LeftValue(iTrial+1) = LeftValue(iTrial);
            RightValue(iTrial+1) = RightValue(iTrial) + PosLearningRate * (1 - RightValue(iTrial));
        else
            LeftValue(iTrial+1) = LeftValue(iTrial);
            RightValue(iTrial+1) = RightValue(iTrial) + NegLearningRate * (0 - RightValue(iTrial));
        end
    else
        NegLogDataLikelihood = NegLogDataLikelihood;
        
        LeftValue(iTrial+1) = LeftValue(iTrial);
        RightValue(iTrial+1) = RightValue(iTrial);
    end
end % end for-loop
end % end function