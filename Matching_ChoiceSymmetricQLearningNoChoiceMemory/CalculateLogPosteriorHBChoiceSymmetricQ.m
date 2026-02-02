function LogPosterior = CalculateLogPosteriorHBChoiceSymmetricQ(Parameters, DataHolder, HyperPrior)
%% calculate log prior & gradients i.e. P(theta_i|THETA)
LearningRatePrior = pdf('Beta', Parameters(1), Parameters(2) * Parameters(3), (1 - Parameters(2)) * Parameters(3));
LogLearningRatePrior = log(LearningRatePrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
% GradLogLearningRatePrior = (Parameters(2) * Parameters(3) - 1) ./ Parameters(1) +...
%                                - ((1 - Parameters(2)) * Parameters(3) - 1) ./ (1 - Parameters(1));

InverseTemperaturePrior = pdf('Normal', Parameters(4), Parameters(5), sqrt(1 / Parameters(6)));
LogInverseTemperaturePrior = log(InverseTemperaturePrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
% GradLogInverseTemperaturePrior = -(Parameters(4) - Parameters(5)) * Parameters(6);

ForgettingRatePrior = pdf('Beta', Parameters(7), Parameters(8) * Parameters(9), (1 - Parameters(8)) * Parameters(9));
LogForgettingRatePrior = log(ForgettingRatePrior);
% GradLogForgettingRatePrior = (Parameters(8) * Parameters(9) - 1) ./ Parameters(7) +...
%                                - ((1 - Parameters(8)) * Parameters(9) - 1) ./ (1 - Parameters(7));

ChoiceStickinessPrior = pdf('Normal', Parameters(10), Parameters(11), sqrt(1 / Parameters(12)));
LogChoiceStickinessPrior = log(ChoiceStickinessPrior);
% GradLogChoiceStickinessPrior = -(Parameters(10) - Parameters(11)) * Parameters(12);

ChoiceForgettingRatePrior = pdf('Beta', Parameters(13), Parameters(14) * Parameters(15), (1 - Parameters(14)) * Parameters(15));
LogChoiceForgettingRatePrior = log(ChoiceForgettingRatePrior);
% GradLogChoiceForgettingRatePrior = (Parameters(14) * Parameters(15) - 1) ./ Parameters(13) +...
%                                        - ((1 - Parameters(14)) * Parameters(15) - 1) ./ (1 - Parameters(13));

BiasPrior = pdf('Normal', Parameters(16), Parameters(17), sqrt(1 / Parameters(18)));
LogBiasPrior = log(BiasPrior);
% GradLogBiasPrior = -(Parameters(16) - Parameters(17)) * Parameters(18);

LogPrior = LogLearningRatePrior...
               + LogInverseTemperaturePrior...
               + LogForgettingRatePrior...
               + LogChoiceStickinessPrior...
               + LogChoiceForgettingRatePrior...
               + LogBiasPrior;

%% calculate log likelihood of hyper-prior, i.e. P(THETA)
LearningRateMuPrior = pdf('Beta', Parameters(2), HyperPrior.LearningRateMuAlpha, HyperPrior.LearningRateMuBeta);
LogLearningRateMuPrior = log(LearningRateMuPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
% GradLogLearningRateMuPrior = (HyperPrior.LearningRateMuAlpha - 1) ./ Parameters(2)...
%                                 - (HyperPrior.LearningRateMuBeta - 1) ./ (1 - Parameters(2));

LearningRateKappaPrior = pdf('Gamma', Parameters(3), HyperPrior.LearningRateKappaAlpha, HyperPrior.LearningRateKappaTheta);
LogLearningRateKappaPrior = log(LearningRateKappaPrior);
% GradLogLearningRateKappaPrior = (HyperPrior.LearningRateKappaAlpha - 1) ./ Parameters(3)...
%                                     - HyperPrior.LearningRateKappaTheta;

InverseTemperatureMeanPrior = pdf('Normal', Parameters(5), HyperPrior.InverseTemperatureMeanMu, HyperPrior.InverseTemperatureMeanSigma);
LogInverseTemperatureMeanPrior = log(InverseTemperatureMeanPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
% GradLogInverseTemperatureMeanPrior = -(Parameters(5) - HyperPrior.InverseTemperatureMeanMu) ./ HyperPrior.InverseTemperatureMeanSigma.^2;

InverseTemperaturePrecisionPrior = pdf('Gamma', Parameters(6), HyperPrior.InverseTemperaturePrecisionAlpha, HyperPrior.InverseTemperaturePrecisionTheta);
LogInverseTemperaturePrecisionPrior = log(InverseTemperaturePrecisionPrior);
% GradLogInverseTemperaturePrecisionPrior = (HyperPrior.InverseTemperaturePrecisionAlpha - 1) ./ Parameters(6)...
%                                               - HyperPrior.InverseTemperaturePrecisionTheta;

ForgettingRateMuPrior = pdf('Beta', Parameters(8), HyperPrior.ForgettingRateMuAlpha, HyperPrior.ForgettingRateMuBeta);
LogForgettingRateMuPrior = log(ForgettingRateMuPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
% GradLogForgettingRateMuPrior = (HyperPrior.ForgettingRateMuAlpha - 1) ./ Parameters(8)...
%                                    - (HyperPrior.ForgettingRateMuBeta - 1) ./ (1 - Parameters(8));

ForgettingRateKappaPrior = pdf('Gamma', Parameters(9), HyperPrior.ForgettingRateKappaAlpha, HyperPrior.ForgettingRateKappaTheta);
LogForgettingRateKappaPrior = log(ForgettingRateKappaPrior);
% GradLogForgettingRateKappaPrior = (HyperPrior.ForgettingRateKappaAlpha - 1) ./ Parameters(9)...
%                                       - HyperPrior.ForgettingRateKappaTheta;

ChoiceStickinessMeanPrior = pdf('Normal', Parameters(11), HyperPrior.ChoiceStickinessMeanMu, HyperPrior.ChoiceStickinessMeanSigma);
LogChoiceStickinessMeanPrior = log(ChoiceStickinessMeanPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
% GradLogChoiceStickinessMeanPrior = -(Parameters(11) - HyperPrior.ChoiceStickinessMeanMu) ./ HyperPrior.ChoiceStickinessMeanSigma.^2;

ChoiceStickinessPrecisionPrior = pdf('Gamma', Parameters(12), HyperPrior.ChoiceStickinessPrecisionAlpha, HyperPrior.ChoiceStickinessPrecisionTheta);
LogChoiceStickinessPrecisionPrior = log(ChoiceStickinessPrecisionPrior);
% GradLogChoiceStickinessPrecisionPrior = (HyperPrior.ChoiceStickinessPrecisionAlpha - 1) ./ Parameters(12)...
%                                             - HyperPrior.ChoiceStickinessPrecisionTheta;

ChoiceForgettingRateMuPrior = pdf('Beta', Parameters(14), HyperPrior.ChoiceForgettingRateMuAlpha, HyperPrior.ChoiceForgettingRateMuBeta);
LogChoiceForgettingRateMuPrior = log(ChoiceForgettingRateMuPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
% GradLogChoiceForgettingRateMuPrior = (HyperPrior.ChoiceForgettingRateMuAlpha - 1) ./ Parameters(14)...
%                                          - (HyperPrior.ChoiceForgettingRateMuBeta - 1) ./ (1 - Parameters(14));

ChoiceForgettingRateKappaPrior = pdf('Gamma', Parameters(15), HyperPrior.ChoiceForgettingRateKappaAlpha, HyperPrior.ChoiceForgettingRateKappaTheta);
LogChoiceForgettingRateKappaPrior = log(ChoiceForgettingRateKappaPrior);
% GradLogChoiceForgettingRateKappaPrior = (HyperPrior.ChoiceForgettingRateKappaAlpha - 1) ./ Parameters(15)...
%                                             - HyperPrior.ChoiceForgettingRateKappaTheta;

BiasMeanPrior = pdf('Normal', Parameters(17), HyperPrior.BiasMeanMu, HyperPrior.BiasMeanSigma);
LogBiasMeanPrior = log(BiasMeanPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
% GradLogBiasMeanPrior = -(Parameters(17) - HyperPrior.BiasMeanMu) ./ HyperPrior.BiasMeanSigma.^2;

BiasPrecisionPrior = pdf('Gamma', Parameters(18), HyperPrior.BiasPrecisionAlpha, HyperPrior.BiasPrecisionTheta);
LogBiasPrecisionPrior = log(BiasPrecisionPrior);
% GradLogBiasPrecisionPrior = (HyperPrior.BiasPrecisionAlpha - 1) ./ Parameters(18)...
%                                 - HyperPrior.BiasPrecisionTheta;

LogHyperPrior = LogLearningRateMuPrior...
                    + LogLearningRateKappaPrior...
                    + LogInverseTemperatureMeanPrior...
                    + LogInverseTemperaturePrecisionPrior...
                    + LogForgettingRateMuPrior...
                    + LogForgettingRateKappaPrior...
                    + LogChoiceStickinessMeanPrior...
                    + LogChoiceStickinessPrecisionPrior...
                    + LogChoiceForgettingRateMuPrior...
                    + LogChoiceForgettingRateKappaPrior...
                    + LogBiasMeanPrior...
                    + LogBiasPrecisionPrior;

%% unpack data and calculate datalikelihood, i.e. P(y|theta_i)
Thetas = Parameters(3 * (1:6) -2);
LogPosterior = 0;
% GradLogPosterior = zeros(size(Parameters));
for iSession = 1:length(DataHolder)
    SessionData = DataHolder{iSession};
    nTrials = SessionData.nTrials;

    if nTrials < 50
        fprintf('Session #%2.0f has only %3.0f trials. Too little for analysis.\n', iSession, nTrials)
        continue
    end
    
    ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
    Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

    [NegLogDataLikelihood, ~] = ChoiceSymmetricQLearning(Thetas, nTrials, ChoiceLeft, Rewarded);
    LogPosterior = LogPosterior - NegLogDataLikelihood;

    %{
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
    %}
end

%% calculate log posterior & gradients
LogPosterior = LogPosterior +...
               LogPrior +...
               LogHyperPrior;

%{
GradLogPosterior = GradLogLikelihood +...
                   [GradLogLearningRateMuPrior,...
                    GradLogInverseTemperaturePrior,...
                    GradLogForgettingRatePrior,...
                    GradLogChoiceStickinessPrior,...
                    GradLogChoiceForgettingRatePrior,...
                    GradLogBiasPrior];
%}
end