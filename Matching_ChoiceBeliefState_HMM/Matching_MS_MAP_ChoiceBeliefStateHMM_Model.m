function Models = Matching_MS_MAP_ChoiceBeliefStateHMM_Model(DataHolder, Models)
% concatenate SessionData in long-format
nSessions = length(DataHolder);

%% initial weights for GLM using belief state models fitted independently by MLE
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

%% run GLM-HMM with 2:5 states using GLM weight to inform prior (not 1 state as dirpdf doesn't work)
Models = cell(1, 4);
for nStates = 2:5
    MinECLL = Inf;
    for iInitialization = 1:20
        %% initialise prior using 1-state estimate with randomness
        Prior = ChoiceBeliefStateHMM_InitializePrior(nStates, InitialGLMMedian, InitialGLMMAD);
        
        %% run E-M algorithm to optimise prior
        try
            Model = ChoiceBeliefStateHMM_EM(DataHolder, Prior);
        catch
            disp('Error: fail to run model');
            continue
        end
        
        if MinECLL > Model.ECLL(end)
            MinECLL = Model.ECLL(end);
            Models{nStates} = Model;
        end
    end
end