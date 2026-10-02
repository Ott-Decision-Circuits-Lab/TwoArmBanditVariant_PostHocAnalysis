function [EstimatedParameter, MinNegLogDataLikelihood, Values, Grad, Hessian, NegLogDataLikelihood]...
    = ChoiceBeliefStateHMM_EM(InitialParameter, Data, Prior)
%{
Basically Baum-Welch algorithm
Data cannot have NoChoice
%}

%% extract parameters
nStates = InitialParameter.nStates;
InitialStateProbability = InitialParameter.InitialStateProbability; % nState x 1
StateTransitionMatrix = InitialParameter.StateTransitionMatrix; % nState_i x nState_j
StateStrategyWeight = InitialParameter.StateStrategyWeight; % nState x 4
StrategyParameter = InitialParameter.StrategyParameter; % 4

%% extract data
nTrials = height(Data);
ChoiceLeft = Data.ChoiceLeft;
Rewarded = Data.Rewarded;

SessionIdx = Data.SessionIdx;
nSessions = max(SessionIdx);
TrialIdx = Data.TrialIdx;

%%
ImprovedLoss = Inf;
NegLogDataLikelihood = Inf;
while ImprovedLoss > 1e-5 && length(NegLogDataLikelihood) < 301
    %% E-step
    %{
    evaluate marginal posterior of hidden state gamma(t, j)=P(z_t=j|D, theta) 
    and joint posterior of consecutive states zeta(t, i, j)=P(z_t=i, z_t+1=j| D, theta)
    using forward-backward algorithm
    alpha_i = joint posterior of choice_1:i and state_i given all X before
    beta_i = joint posterior of choice_i+1:end given all X after and state_i 
    after
    %}
    
    StateEmissionProbability = nan(nStates, nTrials);
    for iState = 1:nStates
        Parameters = [StrategyParameter, StateStrategyWeight(iState, :)];
        Parameters = Parameters([1, 5, 2, 6, 3, 7, 4, 8]);
        
        DataProb = nan(1, nTrials);
        for iSession = 1:nSessions
            IsSession = SessionIdx == iSession;
            SessionnTrials = max(TrialIdx(IsSession));
            SessionChoiceLeft = ChoiceLeft(IsSession);
            SessionRewarded = Rewarded(IsSession);
            [~, Values] = ChoiceBeliefState(Parameters, SessionnTrials, SessionChoiceLeft, SessionRewarded);

            SessionChoiceLeftProb = Values.ChoiceLeftProb;
            SessionDataProb...
                = SessionChoiceLeftProb .* SessionChoiceLeft...
                    + (1 - SessionChoiceLeftProb) .* (1 - SessionChoiceLeft);
            DataProb(IsSession) = SessionDataProb;
        end

        StateEmissionProbability(iState, :) = DataProb;
    end

    % forward pass
    Alpha = nan(nStates, nTrials);
    for iTrial = 1:nTrials
        if TrialIdx(iTrial) == 1
            Alpha(:, iTrial)...
                = InitialStateProbability...
                    .* StateEmissionProbability(:, iTrial);
        else
            Alpha(:, iTrial)...
                = (Alpha(:, iTrial - 1)' * StateTransitionMatrix)'...
                    .* StateEmissionProbability(:, iTrial);
        end
    end

    % backward pass
    Beta = nan(nStates, nTrials);
    for iTrial = nTrials:-1:1
        if iTrial == nTrials || TrialIdx(iTrial + 1) == 1
            Beta(:, iTrial) = 1;
        else
            Beta(:, iTrial)...
                = (Beta(:, iTrial + 1)' / StateTransitionMatrix)'...
                    .* StateEmissionProbability(:, iTrial + 1);
        end
    end
    
    % calculate gamma
    Gamma = Alpha .* Beta; % nStates x nTrials
    Gamma = Gamma ./ sum(Gamma, 1); % <- could be numerically unstable

    % calculate zeta
    nTransitions = sum(TrialIdx ~= 1);
    Zeta = nan(nStates, nTransitions, nStates);
    for iSession = 1:nSessions
        IsSession = SessionIdx == iSession;
        SessionnTransitions = sum(IsSession) - 1;
        SessionStateEmissionProbability = StateEmissionProbability(:, IsSession);

        SessionAlpha = Alpha(:, IsSession);
        SessionBeta = Beta(:, IsSession);
        SessionZeta...
            = repmat(SessionAlpha(:, 1:end-1), 1, 1, nStates)...
                .* permute(...
                    repmat(...
                    SessionBeta(:, 2:end) .* SessionStateEmissionProbability(:, 2:end),...
                    1, 1, nStates),...
                    [3, 2, 1])...
                .* permute(repmat(StateTransitionMatrix, 1, 1, SessionnTransitions), [1, 3, 2]);
        
        SessionZeta = SessionZeta ./ sum(SessionZeta, [1, 3]);  % <- could be numerically unstable

        FirstIdx = find(IsSession, 1, 'first');
        LastIdx = find(IsSession, 1, 'last');
        Idx = (FirstIdx - iSession + 1):(LastIdx - iSession);
        Zeta(:, Idx, :) = SessionZeta;
    end

    %% M-step
    %{
    maximizing log(likelihood) by optimizing StateStrategyWeight and
    StrategyParameter based on the expected state(s) they are in
    %}
    
    % update pi and transition matrix by new expectation
    InitialStateProbability...
        = sum(Gamma(TrialIdx==1, :), 2) ./ nSessions;

    StateTransitionMatrix...
        = sum(Zeta, 2) ./ nTransitions;
    
    % format optimation problem
    StateInfo.Gamma = Gamma;
    StateInfo.Zeta = Zeta;
    StateInfo.InitialStateProbability = InitialStateProbability;
    StateInfo.StateTransitionMatrix = StateTransitionMatrix;
    
    CalculateECLL = @(Parameters) ChoiceBeliefStateHMM(Parameters, Data, StateInfo, Prior);
    
    % update Parameters by new state information
    Model = struct();
    Model.LowerBound = LowerBound;
    Model.UpperBound = UpperBound;
    Model.MinNegLogDataLikelihood = Inf;
    for iInitialCond = 1:20
        InitialParameters...
            = LowerBound + rand * (UpperBound - LowerBound);
        
        try
            [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
                fmincon(CalculateECLL, InitialParameters, [], [], [], [], LowerBound, UpperBound);
        catch
            disp('Error: fail to run model');
            EstimatedParameters = [];
            MinNegLogDataLikelihood = nan;
        end
        
        if Model.MinNegLogDataLikelihood > MinNegLogDataLikelihood
            Model.LowerBound = LowerBound;
            Model.UpperBound = UpperBound;
            Model.InitialParameters = InitialParameters;
            
            Model.EstimatedParameters = EstimatedParameters;
            Model.MinNegLogDataLikelihood = MinNegLogDataLikelihood;
            Model.Grad = Grad;
            Model.Hessian = Hessian;
            try
                Model.ParameterStandardError = sqrt(diag(inv(Hessian)))';
            catch
                Model.ParameterStandardError = nan(size(EstimatedParameters));
            end
        end
    end
    
    %% update loop
    ImprovedLoss = NegLogDataLikelihood(end) - NewMinNegLogDataLikelihood;
    NegLogDataLikelihood(end + 1) = NewMinNegLogDataLikelihood;
    
end

%% report
NegLogDataLikelihood = NegLogDataLikelihood(2:end);
if ImprovedLoss < 1e-5
    Flag = 1;
elseif length(NegLogDataLikelihood) == 300
    Flag = 0;
end

end