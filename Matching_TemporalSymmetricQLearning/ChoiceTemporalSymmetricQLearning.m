function [NegLogDataLikelihood, Values] = ChoiceTemporalSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded, InterChoiceInterval, FeedbackWaitingTime)

% Parameters import
% fixed variable

LearningRate = Parameters(1); % alpha
InverseTemperature = Parameters(2); % beta
ForgettingTau = Parameters(3); % tau
ChoiceStickiness = Parameters(4); % phi
ChoiceForgettingTau = Parameters(5); % tau_c
Bias = 0; %Parameters(6);

% Q-value initialisation
LeftValue = 0; % at TimeChoice/TimeNoDecision/TimeNoTrialStart
RightValue = 0;
ChoiceMemory = 0;
NegLogDataLikelihood = 0;

% Likelihood iteration
for iTrial = 1:nTrials-1 % <- -1 as raw ChoiceLeft has one pre-allocated nan
    LogOdds = InverseTemperature * (LeftValue(iTrial) - RightValue(iTrial)) +...
              + ChoiceStickiness * ChoiceMemory(iTrial) + Bias;
    
    LeftValue(iTrial + 1) = LeftValue(iTrial) * exp(- InterChoiceInterval(iTrial) ./ ForgettingTau);
    RightValue(iTrial + 1) = RightValue(iTrial) * exp(- InterChoiceInterval(iTrial) ./ ForgettingTau);
    ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial) * exp(- InterChoiceInterval(iTrial) ./ ChoiceForgettingTau);

    if ChoiceLeft(iTrial) == 1
        ChoiceLogOdds = LogOdds;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        if Rewarded(iTrial) == 1
            LeftValue(iTrial + 1) = LeftValue(iTrial + 1) +...
                                    LearningRate *...
                                    (1 - LeftValue(iTrial) * exp(- FeedbackWaitingTime(iTrial) ./ ForgettingTau)) *...
                                    exp(- (InterChoiceInterval(iTrial) - FeedbackWaitingTime(iTrial)) ./ ForgettingTau);
        else
            LeftValue(iTrial + 1) = LeftValue(iTrial + 1) +...
                                    LearningRate *...
                                    (0 - LeftValue(iTrial) * exp(- FeedbackWaitingTime(iTrial) ./ ForgettingTau)) *...
                                    exp(- (InterChoiceInterval(iTrial) - FeedbackWaitingTime(iTrial)) ./ ForgettingTau);
        end

        ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial + 1) + exp(-(InterChoiceInterval(iTrial) - FeedbackWaitingTime(iTrial)) ./ ChoiceForgettingTau);
        
    elseif ChoiceLeft(iTrial) == 0
        ChoiceLogOdds = -LogOdds;
        NegLogDataLikelihood = NegLogDataLikelihood - log(1 ./ (1 + exp(-ChoiceLogOdds)));
        
        if Rewarded(iTrial) == 1
            RightValue(iTrial + 1) = RightValue(iTrial + 1) +...
                                     LearningRate *...
                                     (1 - RightValue(iTrial) * exp(-FeedbackWaitingTime(iTrial) ./ ForgettingTau)) *...
                                     exp(- (InterChoiceInterval(iTrial) -FeedbackWaitingTime(iTrial)) ./ ForgettingTau);
        else
            RightValue(iTrial + 1) = RightValue(iTrial + 1) +...
                                     LearningRate *...
                                     (0 - RightValue(iTrial) * exp(-FeedbackWaitingTime(iTrial) ./ ForgettingTau)) *...
                                     exp(-(InterChoiceInterval(iTrial) - FeedbackWaitingTime(iTrial)) ./ ForgettingTau);
        end
        
        ChoiceMemory(iTrial + 1) =  ChoiceMemory(iTrial + 1) - exp(-(InterChoiceInterval(iTrial) - FeedbackWaitingTime(iTrial)) ./ ChoiceForgettingTau);
        
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