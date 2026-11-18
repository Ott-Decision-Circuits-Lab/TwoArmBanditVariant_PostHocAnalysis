function [ECLL, Parameters, Prior]...
    = ChoiceBeliefStateHMM_CalculateECLL(Estimand, Parameters, DataHolder, Prior, ExpectedStateInfo)

%% un-parameterize Estimand (strategy-related Prior)
nStates = Prior.nStates;
StrategyWeightMu = Estimand(1:4*nStates);
Prior.StrategyWeightMu = reshape(StrategyWeightMu, nStates, []);

StrategyWeightSigma = Estimand(4*nStates + 1:8*nStates);
StrategyWeightSigma = log(exp(StrategyWeightSigma) + 1);
Prior.StrategyWeightSigma = reshape(StrategyWeightSigma, nStates, []);

StrategyParameterAlpha = Estimand(8*nStates + 1:8*nStates + 4);
StrategyParameterAlpha = log(exp(StrategyParameterAlpha) + 1);
Prior.StrategyWeightSigma = StrategyParameterAlpha;

StrategyParameterBeta = Estimand(8*nStates + 5:8*nStates + 8);
StrategyParameterBeta = log(exp(StrategyParameterBeta) + 1);
Prior.StrategyParameterBeta = StrategyParameterBeta;

StrategyWeightL1Lambda = Estimand(8*nStates + 9:9*nStates + 8);
StrategyWeightL1Lambda = log(exp(StrategyWeightL1Lambda) + 1);
Prior.StrategyWeightL1Lambda = StrategyWeightL1Lambda;

%%
ECLL = 0;
nSessions = length(DataHolder);
for iSession = 1:nSessions
    SessionParameters = Parameters{iSession};
    SessionData = DataHolder{iSession};
    SessionExpectedStateInfo = ExpectedStateInfo{iSession};
    
    nTrials = SessionData.nTrials;
    if nTrials < 50
        continue
    end
    % get initial values of estimate (which are the strategy-related parameter)
    StrategyWeight = randn(nStates, 4) .* Prior.StrategyWeightSigma + Prior.StrategyWeightMu;
    StrategyParameter = betarnd(Prior.StrategyParameterAlpha, Prior.StrategyParameterBeta);
    
    % reparameterization and linearization
    StrategyParameter = log(exp(StrategyParameter) - 1);
    InitialEstimate...
        = [reshape(StrategyWeight, 1, []),...
           StrategyParameter];

    % set up optimization
    UpperBound = ones(1, 4*nStates + 4) * 40;
    LowerBound = ones(1, 4*nStates + 4) * -40;
    Options.Display = 'off';
    
    LogPosterior = @(Estimand)...
        ChoiceBeliefStateHMM_CalculateSessionECLL(Estimand, SessionData, SessionParameters, Prior, SessionExpectedStateInfo);
    
    FinalEstimate = fmincon(LogPosterior, InitialEstimate, [], [], [], [], LowerBound, UpperBound, [], Options);
    
    % get best estimate of Parameter
    [SessionECLL, SessionParameters]...
        = ChoiceBeliefStateHMM_CalculateSessionECLL(FinalEstimate, SessionData, SessionParameters, Prior, SessionExpectedStateInfo);
    
    % insert result
    Parameters{iSession} = SessionParameters;
    ECLL = ECLL + SessionECLL;
end
%{
%% unpack data
nTrials = SessionData.nTrials;
if nTrials < 50
    disp('nTrial < 50. Impossible for analysis.')
    LogPosterior = nan;
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
%}