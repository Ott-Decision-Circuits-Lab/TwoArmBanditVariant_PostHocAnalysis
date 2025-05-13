function [MeanSquaredError, PredictedInvestedTime] = TIQLearning(Parameters, ChosenValue, UnchosenValue, InvestedTime)

% Parameters import
% fixed variable

BackgroundScaling = Parameters(1); % beta_1
BackgroundIntercept = Parameters(2); % beta_0
TauScaling = Parameters(3); % beta_3
TauIntercept = Parameters(4); % beta_2

% MSE algorithm
DiffValue = ChosenValue - UnchosenValue;
TotalValue = ChosenValue + UnchosenValue;

PredictedInvestedTime = - 1.5 * log((BackgroundScaling .* TotalValue + BackgroundIntercept) ./ (TauScaling .* DiffValue + TauIntercept));
% 
% 
% PredictedInvestedTime = -(TauScaling .* DiffValue + TauIntercept) ./ log(BackgroundScaling .* TotalValue + BackgroundIntercept);

MeanSquaredError = mean((PredictedInvestedTime - InvestedTime) .^ 2);

end % end function