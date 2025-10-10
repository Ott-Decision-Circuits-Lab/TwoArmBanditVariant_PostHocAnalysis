function NegLogDataLikelihood = MultiGeometrics(Parameters, Samples)
% Parameters import
% fixed variable

NegLogDataLikelihood = nan;
if Parameters > 10 % i.e 6 exponentials
    disp('Error: More than 5 exponentials to fit.')
    return
end

% initialisation
NegLogDataLikelihood = 0;

% Likelihood iteration
iGeometric = (length(Parameters) + 1) / 2;

WeightIdx = 2 * (1:(iGeometric-1));
if iGeometric > 1
    Weights = Parameters(WeightIdx);
    Weights = [Weights, 1 - sum(Weights)];
else
    Weights = 1;
end

ProbIdx = 2 * (1:iGeometric) - 1;
TransitionProbs = Parameters(ProbIdx);

for iSample = 1:length(Samples)
    NegLogDataLikelihood...
        = NegLogDataLikelihood...
          - log(sum(Weights...
                    .* TransitionProbs...
                    .* (1 - TransitionProbs) .^ (Samples(iSample) - 1)));
end % end for-loop

end % end function