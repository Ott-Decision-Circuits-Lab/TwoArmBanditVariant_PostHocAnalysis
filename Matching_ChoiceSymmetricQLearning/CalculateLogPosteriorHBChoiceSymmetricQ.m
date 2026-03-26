function LogPosterior = CalculateLogPosteriorHBChoiceSymmetricQ(Parameters, DataHolder, HyperPrior)
% convert Parameters back from real number space to designated space
HyperParameters = Parameters(1:12);
HyperParameters([1, 5, 9]) = 1 ./ (1 + exp(-HyperParameters([1, 5, 9])));
HyperParameters([2, 4, 6, 8, 10, 12]) = exp(HyperParameters([2, 4, 6, 8, 10, 12]));

SessionParameters = Parameters(13:end);
SessionParameters = reshape(SessionParameters, 6, []);
% SessionParameters([1, 3, 5], :) = 1 ./ (1 + exp(-SessionParameters([1, 3, 5], :)));

%% P(HyperPrior) ~ f(Hypermeters)
% calculate log likelihood of hyper-prior, i.e. P(THETA)
LearningRateMuPrior = betapdf(HyperParameters(1), HyperPrior.LearningRateMuAlpha, HyperPrior.LearningRateMuBeta);
LogLearningRateMuPrior = log(LearningRateMuPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)

LearningRateKappaPrior = gampdf(HyperParameters(2), HyperPrior.LearningRateKappaAlpha, HyperPrior.LearningRateKappaTheta);
LogLearningRateKappaPrior = log(LearningRateKappaPrior);

InverseTemperatureMeanPrior = normpdf(HyperParameters(3), HyperPrior.InverseTemperatureMeanMu, HyperPrior.InverseTemperatureMeanSigma);
LogInverseTemperatureMeanPrior = log(InverseTemperatureMeanPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)

InverseTemperaturePrecisionPrior = gampdf(HyperParameters(4), HyperPrior.InverseTemperaturePrecisionAlpha, HyperPrior.InverseTemperaturePrecisionTheta);
LogInverseTemperaturePrecisionPrior = log(InverseTemperaturePrecisionPrior);

ForgettingRateMuPrior = betapdf(HyperParameters(5), HyperPrior.ForgettingRateMuAlpha, HyperPrior.ForgettingRateMuBeta);
LogForgettingRateMuPrior = log(ForgettingRateMuPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)

ForgettingRateKappaPrior = gampdf(HyperParameters(6), HyperPrior.ForgettingRateKappaAlpha, HyperPrior.ForgettingRateKappaTheta);
LogForgettingRateKappaPrior = log(ForgettingRateKappaPrior);

ChoiceStickinessMeanPrior = normpdf(HyperParameters(7), HyperPrior.ChoiceStickinessMeanMu, HyperPrior.ChoiceStickinessMeanSigma);
LogChoiceStickinessMeanPrior = log(ChoiceStickinessMeanPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)

ChoiceStickinessPrecisionPrior = gampdf(HyperParameters(8), HyperPrior.ChoiceStickinessPrecisionAlpha, HyperPrior.ChoiceStickinessPrecisionTheta);
LogChoiceStickinessPrecisionPrior = log(ChoiceStickinessPrecisionPrior);

ChoiceForgettingRateMuPrior = betapdf(HyperParameters(9), HyperPrior.ChoiceForgettingRateMuAlpha, HyperPrior.ChoiceForgettingRateMuBeta);
LogChoiceForgettingRateMuPrior = log(ChoiceForgettingRateMuPrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)

ChoiceForgettingRateKappaPrior = gampdf(HyperParameters(10), HyperPrior.ChoiceForgettingRateKappaAlpha, HyperPrior.ChoiceForgettingRateKappaTheta);
LogChoiceForgettingRateKappaPrior = log(ChoiceForgettingRateKappaPrior);

BiasMeanPrior = normpdf(HyperParameters(11), HyperPrior.BiasMeanMu, HyperPrior.BiasMeanSigma);
LogBiasMeanPrior = log(BiasMeanPrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
% GradLogBiasMeanPrior = -(Parameters(17) - HyperPrior.BiasMeanMu) ./ HyperPrior.BiasMeanSigma.^2;

BiasPrecisionPrior = gampdf(HyperParameters(12), HyperPrior.BiasPrecisionAlpha, HyperPrior.BiasPrecisionTheta);
LogBiasPrecisionPrior = log(BiasPrecisionPrior);

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
LogPosterior =  LogHyperPrior;

for iSession = 1:length(DataHolder)
    SessionData = DataHolder{iSession};
    nTrials = SessionData.nTrials;

    if nTrials < 50
        fprintf('Session #%2.0f has only %3.0f trials. Too little for analysis.\n', iSession, nTrials)
        continue
    end
    
    %% log data posterior, i.e. P(Y_i | theta_i)
    ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
    Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);
    
    Parameters = SessionParameters(:, iSession);
    [NegLogDataLikelihood, ~] = ChoiceSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded);
    
    %% generate and calculate log prior & gradients i.e. P(theta_i|THETA)
    LearningRate = Parameters(1); % betarnd(Parameters(1) * Parameters(2), (1 - Parameters(1)) * Parameters(2));
    LearningRatePrior = betapdf(LearningRate, HyperParameters(1) * HyperParameters(2), (1 - HyperParameters(1)) * HyperParameters(2));
    LogLearningRatePrior = log(LearningRatePrior); % = -log beta(alpha, beta) + (alpha-1) * logx + (beta-1) * log(1-x)
    % GradLogLearningRatePrior = (Parameters(2) * Parameters(3) - 1) ./ Parameters(1) +...
    %                                - ((1 - Parameters(2)) * Parameters(3) - 1) ./ (1 - Parameters(1));
    
    InverseTemperature = Parameters(2); % normrnd(Parameters(3), sqrt(1 / Parameters(4)));
    InverseTemperaturePrior = normpdf(InverseTemperature, HyperParameters(3), sqrt(1 / HyperParameters(4)));
    LogInverseTemperaturePrior = log(InverseTemperaturePrior); % = -log sigma - 0.5 * log (2*pi) - (x-mu).^2 ./ (2 * sigma^2)
    % GradLogInverseTemperaturePrior = -(Parameters(4) - Parameters(5)) * Parameters(6);
    
    ForgettingRate = Parameters(3); % betarnd(Parameters(5) * Parameters(6), (1 - Parameters(5)) * Parameters(6));
    ForgettingRatePrior = betapdf(ForgettingRate, HyperParameters(5) * HyperParameters(6), (1 - HyperParameters(5)) * HyperParameters(6));
    LogForgettingRatePrior = log(ForgettingRatePrior);
    % GradLogForgettingRatePrior = (Parameters(8) * Parameters(9) - 1) ./ Parameters(7) +...
    %                                - ((1 - Parameters(8)) * Parameters(9) - 1) ./ (1 - Parameters(7));
    
    ChoiceStickiness = Parameters(4); % normrnd(Parameters(7), sqrt(1 / Parameters(8)));
    ChoiceStickinessPrior = normpdf(ChoiceStickiness, HyperParameters(7), sqrt(1 / HyperParameters(8)));
    LogChoiceStickinessPrior = log(ChoiceStickinessPrior);
    % GradLogChoiceStickinessPrior = -(Parameters(10) - Parameters(11)) * Parameters(12);
    
    ChoiceForgettingRate = Parameters(5); % betarnd(Parameters(9) * Parameters(10), (1 - Parameters(9)) * Parameters(10));
    ChoiceForgettingRatePrior = betapdf(ChoiceForgettingRate, HyperParameters(9) * HyperParameters(10), (1 - HyperParameters(9)) * HyperParameters(10));
    LogChoiceForgettingRatePrior = log(ChoiceForgettingRatePrior);
    % GradLogChoiceForgettingRatePrior = (Parameters(14) * Parameters(15) - 1) ./ Parameters(13) +...
    %                                        - ((1 - Parameters(14)) * Parameters(15) - 1) ./ (1 - Parameters(13));
    
    Bias = Parameters(6); % normpdf(Parameters(11), sqrt(1 / Parameters(12)));
    BiasPrior = normpdf(Bias, HyperParameters(11), sqrt(1 / HyperParameters(12)));
    LogBiasPrior = log(BiasPrior);
    % GradLogBiasPrior = -(Parameters(16) - Parameters(17)) * Parameters(18);
    
    LogPrior = LogLearningRatePrior...
                   + LogInverseTemperaturePrior...
                   + LogForgettingRatePrior...
                   + LogChoiceStickinessPrior...
                   + LogChoiceForgettingRatePrior...
                   + LogBiasPrior;

    LogPosterior = LogPosterior - NegLogDataLikelihood + LogPrior;
    
end
end