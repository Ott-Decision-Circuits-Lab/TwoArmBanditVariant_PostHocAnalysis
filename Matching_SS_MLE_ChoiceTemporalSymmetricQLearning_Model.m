function Model = Matching_SS_MLE_ChoiceTemporalSymmetricQLearning_Model(SessionData)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

ChoiceTime = [];
for iTrial = 1:nTrials
    if isnan(ChoiceLeft(iTrial))
        ChoiceTime(iTrial) = SessionData.RawEvents.Trial{iTrial}.States.ITI(1, 1);
    elseif ChoiceLeft(iTrial) == 1
        ChoiceTime(iTrial) = SessionData.RawEvents.Trial{iTrial}.States.StartLIn(1, 1);
    elseif ChoiceLeft(iTrial) == 0
        ChoiceTime(iTrial) = SessionData.RawEvents.Trial{iTrial}.States.StartRIn(1, 1);
    end
end

AbsChoiceTime = ChoiceTime + SessionData.TrialStartTimestamp;

InterChoiceInterval = diff(AbsChoiceTime);
FeedbackWaitingTime = SessionData.Custom.TrialData.FeedbackWaitingTime;

% Parametric estimation
LowerBound = [0.10, 6, 5, -1.5, 5, -0.5];
UpperBound = [0.45, 10, 600, 0.5, 600, 0.5];

% Free parameters
% 20241220 tested with simulation that works well as initial parameters
LearningRate = 0.25; % alpha
InverseTemperature = 8; % beta
ForgettingTau = 60; % tau
ChoiceStickiness = -1; % phi
ChoiceForgettingTau = 60; % tau_c
Bias = 0;

InitialParameters = [LearningRate, InverseTemperature, ForgettingTau, ChoiceStickiness, ChoiceForgettingTau, Bias];

CalculateMLE = @(Parameters) ChoiceTemporalSymmetricQLearning(Parameters, nTrials, ChoiceLeft, Rewarded, InterChoiceInterval, FeedbackWaitingTime);

try
    [EstimatedParameters, MinNegLogDataLikelihood] =...
        fmincon(CalculateMLE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
catch
    disp('Error: fail to run model');
    EstimatedParameters = [];
    MinNegLogDataLikelihood = nan;
end

Model.EstimatedParameters = EstimatedParameters;
Model.MinNegLogDataLikelihood = MinNegLogDataLikelihood;

end % end function