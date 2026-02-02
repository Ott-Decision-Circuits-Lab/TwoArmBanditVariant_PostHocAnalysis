function Model = Matching_MS_HB_ChoiceSymmetricQLearning_Model(DataHolder, HyperPrior)
% MS = MultiSession
% HB = Hierarchical Bayesian <- find hyperparameters using simulation & MCMC (Hamiltonian MC) sampling
% from hyperparameters to get prior
% Matching Analysis Function
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ Dez 2025
% Model iteration see the end of script

%% hyper-prior
% set initial point of MCMC around simulation results
SamplerInitialHyperParameters = [0.25, 16,... % LearningRate: Mu, Kappa
                                 8, 1,...     % InverseTemperature: Mean, Sigma
                                 1/6, 12,...  % ForgettingRate: Mu, Kappa
                                 -1, 1,...    % ChoiceStickiness: Mean, Sigma
                                 5/8, 16,...  % ChoiceForgettingRate: Mu, Kappa
                                 0, 1]';      % Bias: Mean, Sigma

LogPosteriorPDF = @(Parameters) CalculateLogPosteriorHBChoiceSymmetricQ(Parameters, DataHolder, HyperPrior);
HyperPriorSampler = hmcSampler(LogPosteriorPDF, SamplerInitialHyperParameters,...
                               'CheckGradient', false,...
                               'UseNumericalGradient', true);

% MAP estimated by L-BFGS
disp(strcat('Estimating MAP from -', datestr(datetime('now'))))
[MAPParameters, FitInfo] = estimateMAP(HyperPriorSampler,...
                                       'VerbosityLevel', 0);

EstimationSuccess = true;
if all(isnan(MAPParameters))
    MAPParameters = SamplerInitialHyperParameters;
    EstimationSuccess = false;
end

% Tune sampler to improve acceptance ratio
%{
- set initial point of MCMC at MAP estimates
- Tutorial suggests NumStepSizeTuningIterations as 50, but takes too much time
- 'MassVectorTuningMethod',  <- 'hessian' can speed up, but chain got
stuck, also MassVector may get into Nan
- 'NumStepsLimit', 100 <- limit NumStep to a lower range so that it won't tune the step size too small
- ideal acceptance ratio = 0.65, but 0.3 is also okay
%}
[HyperPriorSampler, Info] = tuneSampler(HyperPriorSampler,...
                                        'Start', MAPParameters,...
                                        'NumStepSizeTuningIterations', 10,...
                                        'NumStepsLimit', 100,...
                                        'VerbosityLevel', 1,...
                                        'NumPrint', 1);

for iChain = 1:HyperPrior.nChain
    InitialParameters = MAPParameters;
    BetaIdx = [1, 2, 7, 8, 13, 14];
    GammaIdx = 3 * (1:6);

    InitialParameters(BetaIdx) = log(InitialParameters(BetaIdx) ./ (1 - InitialParameters(BetaIdx)));
    InitialParameters(GammaIdx) = sqrt(InitialParameters(GammaIdx));

    InitialParameters = InitialParameters + randn(size(InitialParameters));
    InitialParameters(BetaIdx) = 1 ./ (1 + exp(-InitialParameters(BetaIdx)));
    InitialParameters(GammaIdx) = InitialParameters(GammaIdx) .^ 2;
    
    ChainInitialParameters{iChain} = InitialParameters;
    
    Chains{iChain} = drawSamples(HyperPriorSampler,...
                                 'Start', ChainInitialParameters{iChain},...
                                 'Burnin', HyperPrior.BurnIn,...
                                 'NumSamples', HyperPrior.nSample,...
                                 'VerbosityLevel', 1,...
                                 'NumPrint', 500);
end

Diags = diagnostics(PriorSampler, Chains);

% Parameters = median(Chain);
% nTrials = SessionData.nTrials;
% ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
% Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);
% 
% [NegLogDataLikelihood, Values] = ChoiceSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded);

Model.HyperPrior = HyperPrior;
Model.SamplerInitialHyperParameters = SamplerInitialHyperParameters;
Model.Sampler = PriorSampler;
Model.SamplerTuningInfo = Info;
Model.MAPParameters = MAPParameters;
Model.FitInfo = FitInfo;
Model.EstimationFlag = EstimationSuccess;
Model.ChainInitialParameters = ChainInitialParameters;
Model.Chains = Chains;
Model.Diags = Diags;
% Model.ParametersMedian = Parameters;
% Model.NegLogDataLikelihood = NegLogDataLikelihood;
% Model.PredictedValues = Values;
end % function