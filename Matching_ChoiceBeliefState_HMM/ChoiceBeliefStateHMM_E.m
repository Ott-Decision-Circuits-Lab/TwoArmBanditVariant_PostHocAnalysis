function [ExpectedStateInfo, Prior] = ChoiceBeliefStateHMM_E(Parameters, DataHolder, Prior)
%{
evaluate marginal posterior of hidden state gamma(t, j)=P(z_t=j|D, theta) 
and joint posterior of consecutive states zeta(t, i, j)=P(z_t=i, z_t+1=j| D, theta)
using forward-backward algorithm
alpha_i = joint posterior of choice_1:i and state_i given all X before
beta_i = joint posterior of choice_i+1:end given all X after and state_i 
after
%}

%% extract parameters
nStates = Prior.nStates;
nSessions = length(DataHolder);

%% iterate through Session
ExpectedStateInfo = cell(1, nSessions);
AllInitialStateProbability = zeros(nStates, nSessions);
AllStateTransitionMatrix = zeros(nStates, nStates, nSessions);
for iSession = 1:nSessions
    %% preparing information
    % extract SessionData
    SessionData = DataHolder{iSession};
    TrialData = SessionData.Custom.TrialData;
    ChoiceLeft = TrialData.ChoiceLeft;
    Rewarded = TrialData.Rewarded;
    
    % filter data (NoChoice -> No Emission Prob)
    Valid = ~isnan(ChoiceLeft) & ~isnan(Rewarded);
    nTrials = sum(Valid);
    ChoiceLeft = ChoiceLeft(Valid);
    Rewarded = Rewarded(Valid);

    % extract SessionParameters
    InitialStateProbability = Parameters{iSession}.InitialStateProbability;
    StateTransitionMatrix = Parameters{iSession}.StateTransitionMatrix;
    StrategyParameter = Parameters{iSession}.StrategyParameter;
    StrategyWeight = Parameters{iSession}.StrategyWeight;
    
    % calculate DataEmissionProb
    DataEmissionProb = nan(nStates, nTrials);
    for iState = 1:nStates
        StrategyParameters = [StrategyParameter, StrategyWeight(iState, :)];
        StrategyParameters = StrategyParameters([1, 5, 2, 6, 3, 7, 4, 8]);
        
        [~, Values] = ChoiceBeliefState(StrategyParameters, nTrials, ChoiceLeft, Rewarded);
        ChoiceLeftProb = Values.ChoiceLeftProb;
        DataEmissionProb(iState, :)...
            = ChoiceLeftProb .* ChoiceLeft...
                + (1 - ChoiceLeftProb) .* (1 - ChoiceLeft);

    end
    
    %% forward-backward algorithm
    % forward pass
    Alpha = nan(nStates, nTrials);
    for iTrial = 1:nTrials
        if iTrial == 1
            Alpha(:, iTrial)...
                = InitialStateProbability...
                    .* DataEmissionProb(:, iTrial);
        else
            Alpha(:, iTrial)...
                = (Alpha(:, iTrial - 1)' * StateTransitionMatrix)'...
                    .* DataEmissionProb(:, iTrial);
        end
        Alpha(:, iTrial) = Alpha(:, iTrial) ./ sum(Alpha(:, iTrial)); % numerically stablility
    end
    
    % backward pass
    Beta = nan(nStates, nTrials);
    for iTrial = nTrials:-1:1
        if iTrial == nTrials
            Beta(:, iTrial) = 1;
        else
            Beta(:, iTrial)...
                = (StateTransitionMatrix * Beta(:, iTrial + 1))...
                    .* DataEmissionProb(:, iTrial + 1);
        end
        Beta(:, iTrial) = Beta(:, iTrial) ./ sum(Beta(:, iTrial)); % numerically stablility
    end

    % calculate gamma
    Gamma = Alpha .* Beta; % nStates x nTrials
    Gamma = Gamma ./ sum(Gamma, 1); % <- numerically unstable if only normalized here

    % calculate zeta
    Zeta...
        = repmat(Alpha(:, 1:end-1), 1, 1, nStates)...
            .* permute(...
                repmat(...
                Beta(:, 2:end) .* DataEmissionProb(:, 2:end),...
                1, 1, nStates),...
                [3, 2, 1])...
            .* permute(repmat(StateTransitionMatrix, 1, 1, nTrials-1), [1, 3, 2]);
    Zeta = Zeta ./ sum(Zeta, [1, 3]);  % <- numerically unstable if only normalized here
    
    %% format result and update state parameters
    ExpectedStateInfo{iSession}.Alpha = Alpha;
    ExpectedStateInfo{iSession}.Beta = Beta;
    ExpectedStateInfo{iSession}.Gamma = Gamma;
    ExpectedStateInfo{iSession}.Zeta = Zeta;
    
    InitialStateProbability = Gamma(:, 1);
    StateTransitionMatrix = mean(Zeta, 2);
    
    Parameters{iSession}.InitialStateProbability = InitialStateProbability;
    Parameters{iSession}.StateTransitionMatrix = StateTransitionMatrix;
    AllInitialStateProbability(:, iSession) = InitialStateProbability;
    AllStateTransitionMatrix(:, :, iSession) = StateTransitionMatrix;

end

% update state-related prior using method of moments
TempMean = mean(AllInitialStateProbability, 2);
TempVar = var(AllInitialStateProbability, 0, 2);
Prior.InitialStateProbabilityAlpha = TempMean .* (TempMean .* (1 - TempMean) ./ TempVar - 1);
for iState = 1:nStates
    TempMean = mean(AllStateTransitionMatrix(iState, :, :), 3);
    TempVar = var(AllStateTransitionMatrix(iState, :, :), 0, 3);
    Prior.StateTransitionMatrixAlpha(iState, :) = TempMean .* (TempMean .* (1 - TempMean) ./ TempVar - 1);
end

end