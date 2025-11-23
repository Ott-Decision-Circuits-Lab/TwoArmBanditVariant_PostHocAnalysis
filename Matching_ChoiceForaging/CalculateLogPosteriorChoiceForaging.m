function [LogPosterior, GradLogPosterior] = CalculateLogPosteriorChoiceForaging(Parameters, SessionData, Prior)

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
[NegLogDataLikelihood, ~] = ChoiceForaging(Parameters, nTrials, ChoiceLeft, Rewarded);
GradLogLikelihood = zeros(1, length(Parameters));

for iParameter = 1:length(Parameters)
    ThetaPlus = Parameters(iParameter) * 1.01;
    
    NewParameters = Parameters;
    NewParameters(iParameter) = ThetaPlus;

    [NegLogDataLikelihoodPlus, ~] = ChoiceForaging(NewParameters, nTrials, ChoiceLeft, Rewarded);
    
    ThetaMinus = Parameters(iParameter) * 0.99;
    NewParameters(iParameter) = ThetaMinus;
    
    [NegLogDataLikelihoodMinus, ~] = ChoiceForaging(NewParameters, nTrials, ChoiceLeft, Rewarded);

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

ThresholdPrior = pdf('Normal', Parameters(3), Prior.ThresholdMean, Prior.ThresholdSigma);
LogThresholdPrior = log(ThresholdPrior);
GradLogThresholdPrior = -(Parameters(3) - Prior.ThresholdMean) ./ Prior.ThresholdSigma.^2;

ForgettingRatePrior = pdf('Beta', Parameters(4), Prior.ForgettingRateAlpha, Prior.ForgettingRateBeta);
LogForgettingRatePrior = log(ForgettingRatePrior);
GradLogForgettingRatePrior = (Prior.ForgettingRateAlpha - 1) ./ Parameters(4) +...
                               - (Prior.ForgettingRateBeta - 1) ./ (1 - Parameters(4));

%% calculate log posterior & gradients
LogPosterior = - NegLogDataLikelihood +...
                 LogLearningRatePrior +...
                 LogInverseTemperaturePrior +...
                 LogForgettingRatePrior +...
                 LogThresholdPrior;

GradLogPosterior = GradLogLikelihood +...
                   [GradLogLearningRatePrior,...
                    GradLogInverseTemperaturePrior,...
                    GradLogForgettingRatePrior,...
                    GradLogThresholdPrior];

end