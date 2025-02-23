function [LogPosterior, GradLogPosterior] = CalculateLogPosteriorChoiceSymmetricQ(Parameters, SessionData, Prior)

%% unpack data
nTrials = SessionData.nTrials;
if nTrials < 50
    disp('nTrial < 50. Impossible for analysis.')
    LogPosterior = nan;
    GradLogPosterior = nan;
    return
end

ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

%% calculate log likelihood & gradients
[NegLogDataLikelihood, ~] = ChoiceSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded);
GradLogLikelihood = zeros(1, length(Parameters));

for iParameter = 1:length(Parameters)
    ThetaPlus = Parameters(iParameter) * 1.01;
    
    NewParameters = Parameters;
    NewParameters(iParameter) = ThetaPlus;

    [NegLogDataLikelihoodPlus, ~] = ChoiceSymmetricQLearning(NewParameters, nTrials, ChoiceLeft, Rewarded);
    
    ThetaMinus = Parameters(iParameter) * 0.99;
    NewParameters(iParameter) = ThetaMinus;
    
    [NegLogDataLikelihoodMinus, ~] = ChoiceSymmetricQLearning(NewParameters, nTrials, ChoiceLeft, Rewarded);

    GradLogLikelihood(iParameter) = - ((NegLogDataLikelihoodPlus - NegLogDataLikelihood) ./ (ThetaPlus - Parameters(iParameter)) +...
                                      (NegLogDataLikelihoodMinus - NegLogDataLikelihood) ./ (ThetaMinus - Parameters(iParameter))) ./ 2;
end

%% calculate log prior & gradients
LearningRatePrior = pdf('Beta', Parameters(1), Prior.LearningRateAlpha, Prior.LearningRateBeta);
LogLearningRatePrior = log(LearningRatePrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
GradLogLearningRatePrior = (Prior.LearningRateAlpha - 1) ./ Parameters(1) +...
                               - (Prior.LearningRateBeta - 1) ./ (1 - Parameters(1));

InverseTemperaturePrior = pdf('Normal', Parameters(2), Prior.InverseTemperatureMean, Prior.InverseTemperatureSigma);
LogInverseTemperaturePrior = log(InverseTemperaturePrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
GradLogInverseTemperaturePrior = -(Parameters(2) - Prior.InverseTemperatureMean) ./ Prior.InverseTemperatureSigma.^2;

ForgettingRatePrior = pdf('Beta', Parameters(3), Prior.ForgettingRateAlpha, Prior.ForgettingRateBeta);
LogForgettingRatePrior = log(ForgettingRatePrior);
GradLogForgettingRatePrior = (Prior.ForgettingRateAlpha - 1) ./ Parameters(3) +...
                               - (Prior.ForgettingRateBeta - 1) ./ (1 - Parameters(3));

ChoiceStickinessPrior = pdf('Normal', Parameters(4), Prior.ChoiceStickinessMean, Prior.ChoiceStickinessSigma);
LogChoiceStickinessPrior = log(ChoiceStickinessPrior);
GradLogChoiceStickinessPrior = -(Parameters(4) - Prior.ChoiceStickinessMean) ./ Prior.ChoiceStickinessSigma.^2;

ChoiceForgettingRatePrior = pdf('Beta', Parameters(5), Prior.ChoiceForgettingRateAlpha, Prior.ChoiceForgettingRateBeta);
LogChoiceForgettingRatePrior = log(ChoiceForgettingRatePrior);
GradLogChoiceForgettingRatePrior = (Prior.ChoiceForgettingRateAlpha - 1) ./ Parameters(5) +...
                                       - (Prior.ChoiceForgettingRateBeta - 1) ./ (1 - Parameters(5));

BiasPrior = pdf('Normal', Parameters(6), Prior.BiasMean, Prior.BiasSigma);
LogBiasPrior = log(BiasPrior);
GradLogBiasPrior = -(Parameters(6) - Prior.BiasMean) ./ Prior.BiasSigma.^2;

%% calculate log posterior & gradients
LogPosterior = - NegLogDataLikelihood +...
                 LogLearningRatePrior +...
                 LogInverseTemperaturePrior +...
                 LogForgettingRatePrior +...
                 LogChoiceStickinessPrior +...
                 LogChoiceForgettingRatePrior +...
                 LogBiasPrior;

GradLogPosterior = GradLogLikelihood +...
                   [GradLogLearningRatePrior,...
                    GradLogInverseTemperaturePrior,...
                    GradLogForgettingRatePrior,...
                    GradLogChoiceStickinessPrior,...
                    GradLogChoiceForgettingRatePrior,...
                    GradLogBiasPrior];

end