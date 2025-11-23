function Model = Matching_SS_B_ChoiceForaging_Model(SessionData, Prior)
% set initial point of MCMC around simulation results for MAP estimation
SamplerInitialParameters = [0.45, 4, 0.12, 0.16]'; 

LogPosteriorPDF = @(Parameters) CalculateLogPosteriorChoiceForaging(Parameters, SessionData, Prior);
Sampler = hmcSampler(LogPosteriorPDF, SamplerInitialParameters,...
                     'CheckGradient', false,...
                     'UseNumericalGradient', true);

% MAP estimated by L-BFGS
[MAPParameters, FitInfo] = estimateMAP(Sampler,...
                                       'VerbosityLevel', 0);

EstimationSuccess = true;
if all(isnan(MAPParameters))
    MAPParameters = [0.45, 4, 0.12, 0.16]';
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
[Sampler, Info] = tuneSampler(Sampler,...
                              'Start', MAPParameters,...
                              'NumStepSizeTuningIterations', 10,...
                              'NumStepsLimit', 100,...
                              'VerbosityLevel', 1,...
                              'NumPrint', 1);

for iChain = 1:Prior.nChain
    InitialParameters = MAPParameters;
    InitialParameters([1, 4) = log(InitialParameters([1, 4]) ./ (1 - InitialParameters([1, 4])));
    InitialParameters = InitialParameters + randn(4, 1);
    InitialParameters([1, 4]) = 1 ./ (1 + exp(-InitialParameters([1, 4])));
    
    ChainInitialParameters{iChain} = InitialParameters;
    
    Chains{iChain} = drawSamples(Sampler,...
                                'Start', ChainInitialParameters{iChain},...
                                'Burnin', Prior.BurnIn,...
                                'NumSamples', Prior.nSample,...
                                'VerbosityLevel', 1,...
                                'NumPrint', 500);
end

Diags = diagnostics(Sampler, Chains);

% Parameters = median(Chain);
% nTrials = SessionData.nTrials;
% ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
% Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);
% 
% [NegLogDataLikelihood, Values] = ChoiceSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded);

Model.Prior = Prior;
Model.SamplerInitialParameters = SamplerInitialParameters;
Model.Sampler = Sampler;
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

end