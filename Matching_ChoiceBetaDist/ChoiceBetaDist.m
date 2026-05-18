function [NegLogDataLikelihood, Values] = ChoiceBetaDist(Parameters, nTrials, ChoiceLeft, Rewarded)

% Parameters import
% fixed variable

AlphaLearningStepSize = Parameters(1); % \delta_\alpha
AlphaForgettingRate = Parameters(2); % \gamma_\alpha
RLearningRate = Parameters(3); % \alpha_R
RForgettingRate = Parameters(4); % \gamma_R
Bias = Parameters(5); %

% Foraging-value initialisation
AlphaInf = 1;
LeftAlpha = 1;
RightAlpha = 1;
LeftBeta = 1;
RightBeta = 1;
SignalDetectionThreshold = 0;

NegLogDataLikelihood = 0;
ChoiceLeftProb = [];
ChoiceLeftLogOdds = [];

% Likelihood iteration
for iTrial = 1:nTrials
    % decision from values
    LeftValue(iTrial) = betacdf(SignalDetectionThreshold(iTrial),...
                                LeftAlpha(iTrial),...
                                LeftBeta(iTrial),...
                                'upper');

    RightValue(iTrial) = betacdf(SignalDetectionThreshold(iTrial),...
                                 RightAlpha(iTrial),...
                                 RightBeta(iTrial),...
                                 'upper');
    
    ChoiceLeftProb(iTrial)...
        = Bias * LeftValue(iTrial)...
            ./ (Bias * LeftValue(iTrial) + (1 - Bias) * RightValue(iTrial));

    ChoiceLeftLogOdds(iTrial) = log(ChoiceLeftProb(iTrial) ./ (1 - ChoiceLeftProb(iTrial)));
    
    % outcome-value conversion
    LeftAlpha(iTrial + 1)...
        = LeftAlpha(iTrial)...
            - AlphaForgettingRate...
            * (LeftAlpha(iTrial) - AlphaInf);
    LeftBeta(iTrial + 1)...
        = LeftBeta(iTrial)...
            - AlphaForgettingRate...
            * (LeftBeta(iTrial) - AlphaInf);
    RightAlpha(iTrial + 1)...
        = RightAlpha(iTrial)...
            - AlphaForgettingRate...
            * (RightAlpha(iTrial) - AlphaInf);
    RightBeta(iTrial + 1)...
        = RightBeta(iTrial)...
            - AlphaForgettingRate...
            * (RightBeta(iTrial) - AlphaInf);

    if ChoiceLeft(iTrial) == 1
        if Rewarded(iTrial) == 1
            LeftAlpha(iTrial + 1) = LeftAlpha(iTrial + 1)...
                                      + AlphaLearningStepSize;
        else
            LeftBeta(iTrial + 1) = LeftBeta(iTrial + 1)...
                                   + AlphaLearningStepSize;
        end

        NegLogDataLikelihood = NegLogDataLikelihood - log(ChoiceLeftProb(iTrial));

    elseif  ChoiceLeft(iTrial) == 0
        if Rewarded(iTrial) == 1
            RightAlpha(iTrial + 1) = RightAlpha(iTrial + 1)...
                                       + AlphaLearningStepSize;
        else
            RightBeta(iTrial + 1) = RightBeta(iTrial + 1)...
                                    + AlphaLearningStepSize;
        end

        NegLogDataLikelihood = NegLogDataLikelihood - log(1 - ChoiceLeftProb(iTrial));

    else
        NegLogDataLikelihood = NegLogDataLikelihood;
    end
    
    SignalDetectionThreshold(iTrial + 1)...
        = (1 - RForgettingRate)...
            * SignalDetectionThreshold(iTrial);

    if ~isnan(Rewarded(iTrial))
        SignalDetectionThreshold(iTrial + 1)...
            = SignalDetectionThreshold(iTrial + 1)...
                +  RLearningRate...
                * (Rewarded(iTrial) - SignalDetectionThreshold(iTrial));
    end
    
end % end for-loop

Values.LeftValue = LeftValue(1:nTrials);
Values.RightValue = RightValue(1:nTrials);
Values.ChoiceLeftProb = ChoiceLeftProb(1:nTrials);
Values.ChoiceLeftLogOdds = ChoiceLeftLogOdds(1:nTrials);

Values.LeftAlpha = LeftAlpha(1:nTrials);
Values.LeftBeta = LeftBeta(1:nTrials);
Values.RightAlpha = RightAlpha(1:nTrials);
Values.RightBeta = RightBeta(1:nTrials);
Values.SignalDetectionThreshold = SignalDetectionThreshold(1:nTrials);

end % end function