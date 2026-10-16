function Models = Matching_MS_MAP_ChoiceBeliefStateHMM_Model(DataHolder, Models)
% concatenate SessionData in long-format
nSessions = length(DataHolder);

%% initial weights for GLM using 1-state/regular belief state models
EstimatedParameter = nan(8, nSessions);
ModelValidity = false(1, nSessions);
for iSession = 1:nSessions
    Model = Models{iSession};
    if isempty(Model) || ~isfield(Model, 'EstimatedParameters') || all(isnan(Model.EstimatedParameters))
        fprintf('Warning: Invalid model for session number %2.0f', iSession)
        continue
    end
    
    EstimatedParameter(:, iSession) = Model.EstimatedParameters;
    ModelValidity(iSession) = true;
end

InitialGLMMedian = median(EstimatedParameter(:, ModelValidity), 2);
InitialGLMMAD = mad(EstimatedParameter(:, ModelValidity), 1, 2);

%% run GLM-HMM with 2:5 states using GLM weight to inform prior
Models = cell(1, 4);
for nStates = 2:5
    MinECLL = Inf;
    for iInitialization = 1:20
        %% initialise prior using 1-state estimate with randomness
        Prior.nStates = nStates;
        
        Prior.InitialStateProbabilityAlpha...
            = ones(nStates, 1) + normrnd(0, 0.1, nStates, 1);

        Prior.StateTransitionMatrixAlpha...
            = eye(nStates) * 0.9 + normrnd(0.2, 0.1, nStates, nStates);

        Prior.StrategyWeightMu...
            = repmat(InitialGLMMedian(2:2:8)', nStates, 1) + normrnd(0, 0.1, nStates, 4);

        Prior.StrategyWeightSigma...
            = repmat(InitialGLMMAD(2:2:8)', nStates, 1) + abs(normrnd(0, 0.1, nStates, 4));
        
        % using method-of-moment to generate alpha and beta
        TempMedian = InitialGLMMedian(1:2:7);
        TempMAD = InitialGLMMAD(1:2:7);
        Prior.StrategyParameterAlpha...
            = TempMedian...
                .* (TempMedian .* (1 - TempMedian)...
                ./ TempMAD - 1);
        Prior.StrategyParameterBeta...
            = (1 - TempMedian)...
                .* (TempMedian .* (1 - TempMedian)...
                ./ TempMAD - 1); 
        
        for iPara = 1:4
            if Prior.StrategyParameterAlpha(iPara) < 1 || Prior.StrategyParameterBeta(iPara) < 1
                Prior.StrategyParameterAlpha(iPara) = 1 + normrnd(0, 0.1);
                Prior.StrategyParameterBeta(iPara) = 1 + normrnd(0, 0.1);
            end
        end
        Prior.StrategyParameterAlpha = Prior.StrategyParameterAlpha';
        Prior.StrategyParameterBeta = Prior.StrategyParameterBeta';
        
        %% run E-M algorithm to optimise prior
        try
            Model = ChoiceBeliefStateHMM_EM(DataHolder, Prior);
        catch
            disp('Error: fail to run model');
            continue
        end
        
        if MinECLL > Model.ECLL(end)
            MinECLL = Model.ECLL(end);
            Models{nStates - 1} = Model;
        end
    end
end