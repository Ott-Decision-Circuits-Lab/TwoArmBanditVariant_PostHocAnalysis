function Model = MultiGeometrics_Model(nGeometrics, Samples)

% Parametric estimation
LowerBound = zeros(1, 2 * nGeometrics - 1);
UpperBound = ones(1, 2 * nGeometrics - 1);

% Free parameters
% 20250708 tested with simulation that works well as initial parameters
CalculateMLE = @(Parameters) MultiGeometrics(Parameters, Samples);

Model = struct();
Model.Samples = Samples;
Model.LowerBound = LowerBound;
Model.UpperBound = UpperBound;
Model.MinNegLogDataLikelihood = 0;

for iInitialCond = 1:10
    InitialParameters = rand(1, 2 * nGeometrics - 1);
    try
        [EstimatedParameters, MinNegLogDataLikelihood, ~, ~, ~, Grad, Hessian] =...
            fmincon(CalculateMLE, InitialParameters, [], [], [], [], LowerBound, UpperBound);
    catch
        disp('Error: fail to run model');
        EstimatedParameters = [];
        MinNegLogDataLikelihood = nan;
    end
    
    if Model.MinNegLogDataLikelihood < MinNegLogDataLikelihood
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
end % end function