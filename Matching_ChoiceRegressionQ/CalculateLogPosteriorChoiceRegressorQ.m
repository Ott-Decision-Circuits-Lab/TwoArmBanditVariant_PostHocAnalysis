function [LogPosterior, GradLogPosterior] = CalculateLogPosteriorChoiceRegressorQ(Parameters, SessionData, Prior)

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
[NegLogDataLikelihood, ~] = ChoiceRegressorQ(Parameters, nTrials, ChoiceLeft, Rewarded);
GradLogLikelihood = zeros(1, length(Parameters));

for iParameter = 1:length(Parameters)
    ThetaPlus = Parameters(iParameter) * 1.01;
    
    NewParameters = Parameters;
    NewParameters(iParameter) = ThetaPlus;

    [NegLogDataLikelihoodPlus, ~] = ChoiceRegressorQ(NewParameters, nTrials, ChoiceLeft, Rewarded);
    
    ThetaMinus = Parameters(iParameter) * 0.99;
    NewParameters(iParameter) = ThetaMinus;
    
    [NegLogDataLikelihoodMinus, ~] = ChoiceRegressorQ(NewParameters, nTrials, ChoiceLeft, Rewarded);

    GradLogLikelihood(iParameter) = - ((NegLogDataLikelihoodPlus - NegLogDataLikelihood) ./ (ThetaPlus - Parameters(iParameter)) +...
                                      (NegLogDataLikelihoodMinus - NegLogDataLikelihood) ./ (ThetaMinus - Parameters(iParameter))) ./ 2;
end

%% calculate log prior & gradients
BiasPrior = pdf('Normal', Parameters(1), Prior.BiasMean, Prior.BiasSigma);
LogBiasPrior = log(BiasPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
GradLogBiasPrior = -(Parameters(1) - Prior.BiasMean) ./ Prior.BiasSigma.^2;

QCoeffPrior = pdf('Beta', Parameters(2), Prior.QCoeffAlpha, Prior.QCoeffBeta);
LogQCoeffPrior = log(QCoeffPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
GradLogQCoeffPrior = (Prior.QCoeffAlpha - 1) ./ Parameters(2) +...
                         - (Prior.QCoeffBeta - 1) ./ (1 - Parameters(2));

Minus1ChoicePrior = pdf('Normal', Parameters(3), Prior.Minus1ChoiceMean, Prior.Minus1ChoiceSigma);
LogMinus1ChoicePrior = log(Minus1ChoicePrior);
GradLogMinus1ChoicePrior = -(Parameters(3) - Prior.Minus1ChoiceMean) ./ Prior.Minus1ChoiceSigma.^2;

Minus1RewardPrior = pdf('Normal', Parameters(4), Prior.Minus1RewardMean, Prior.Minus1RewardSigma);
LogMinus1RewardPrior = log(Minus1RewardPrior);
GradLogMinus1RewardPrior = -(Parameters(4) - Prior.Minus1RewardMean) ./ Prior.Minus1RewardSigma.^2;

Minus2ChoicePrior = pdf('Normal', Parameters(5), Prior.Minus2ChoiceMean, Prior.Minus2ChoiceSigma);
LogMinus2ChoicePrior = log(Minus2ChoicePrior);
GradLogMinus2ChoicePrior = -(Parameters(5) - Prior.Minus2ChoiceMean) ./ Prior.Minus2ChoiceSigma.^2;

Minus2RewardPrior = pdf('Normal', Parameters(6), Prior.Minus2RewardMean, Prior.Minus2RewardSigma);
LogMinus2RewardPrior = log(Minus2RewardPrior);
GradLogMinus2RewardPrior = -(Parameters(6) - Prior.Minus2RewardMean) ./ Prior.Minus2RewardSigma.^2;

%% calculate log posterior & gradients
LogPosterior = - NegLogDataLikelihood +...
                 LogBiasPrior +...
                 LogQCoeffPrior +...
                 LogMinus1ChoicePrior +...
                 LogMinus1RewardPrior +...
                 LogMinus2ChoicePrior +...
                 LogMinus2RewardPrior;

GradLogPosterior = GradLogLikelihood +...
                   [GradLogBiasPrior,...
                    GradLogQCoeffPrior,
                    GradLogMinus1ChoicePrior,...
                    GradLogMinus1RewardPrior,...
                    GradLogMinus2ChoicePrior,...
                    GradLogMinus2RewardPrior];

end