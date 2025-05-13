function QLearningTIModel = Matching_SS_MSE_TIQLearning_Model(SessionData, QLearningModel)
nTrials = SessionData.nTrials;
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft(1:nTrials);
Rewarded = SessionData.Custom.TrialData.Rewarded(1:nTrials);

if isfield(QLearningModel, 'EstimatedParameters')
    EstimatedParameters = QLearningModel.EstimatedParameters;
elseif isfield(QLearningModel, 'MAPParameters')
    EstimatedParameters = QLearningModel.MAPParameters;
end

[~, Values] = ChoiceSymmetricQLearning(EstimatedParameters, nTrials, ChoiceLeft, Rewarded);

LeftValue = Values.LeftValue;
RightValue = Values.RightValue;

ChosenValue = LeftValue .* ChoiceLeft + RightValue .* (1 - ChoiceLeft);
UnchosenValue = LeftValue .* (1 - ChoiceLeft) + RightValue .* ChoiceLeft;

Baited = SessionData.Custom.TrialData.Baited(:, 1:nTrials);
ChoiceLeftRight = [ChoiceLeft; 1 - ChoiceLeft];
IncorrectChoice = SessionData.Custom.TrialData.IncorrectChoice(1:nTrials);
NotBaited = any(~Baited .* ChoiceLeftRight, 1) & (IncorrectChoice ~= 1);

FeedbackWaitingTime = SessionData.Custom.TrialData.FeedbackWaitingTime(1:nTrials);
InvestedTime = FeedbackWaitingTime(NotBaited);

% Parametric estimation
LowerBound = [0, 0, -4, -4];
UpperBound = [0.5, 1, 8, 8];

% Free parameters
% 20250427 tested with simulation that works well as initial parameters
BackgroundScaling = 0.2; % beta_1
BackgroundIntercept = 0.5; % beta_0
TauScaling = 4; % beta_3
TauIntercept = 2; % beta_2

InitialParameters = [BackgroundScaling, BackgroundIntercept, TauScaling, TauIntercept];
CalculateMSE = @(Parameters) TIQLearning(Parameters, ChosenValue(NotBaited), UnchosenValue(NotBaited), InvestedTime);

try
    [EstimatedParameters, MinMeanSquaredError] =...
        fmincon(CalculateMSE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
catch
    disp('Error: fail to run model');
    EstimatedParameters = [];
    MinMeanSquaredError = nan;
end

QLearningTIModel.LowerBound = LowerBound;
QLearningTIModel.UpperBound = UpperBound;
QLearningTIModel.InitialParameters = InitialParameters;

QLearningTIModel.EstimatedParameters = EstimatedParameters;
QLearningTIModel.MinMeanSquaredError = MinMeanSquaredError;

end % end function