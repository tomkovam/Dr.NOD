function [tablePredictors, nPredictors] = computeSignificantUnivariablePredictors(tableDataForBMM, univariableCutoff, maxPredictors)



tablePredictors = table();
tablePredictors.predictor = tableDataForBMM.Properties.VariableNames';
nPredictors = size(tablePredictors, 1);

tablePredictors.univariablePValue = NaN*ones(nPredictors, 1);
tablePredictors.explainedDeviance = NaN*ones(nPredictors, 1);
tablePredictors.estimate = NaN*ones(nPredictors, 1);

for iCol = 1:nPredictors
    try
        mdl =  fitglm(tableDataForBMM(:,[iCol,end]), 'linear', 'Distribution', 'poisson', 'DispersionFlag', true);
        tablePredictors.univariablePValue(iCol) = mdl.coefTest;
        tablePredictors.explainedDeviance(iCol) = mdl.Rsquared.Deviance;
        tablePredictors.estimate(iCol) = mdl.Coefficients.Estimate(end);
    catch
        warning('PROBLEM in column %s', tablePredictors.predictor{iCol});
    end
end

tablePredictors.isResponseVariable = false(nPredictors, 1); tablePredictors.isResponseVariable(end) = true;
tablePredictors.isUsed = tablePredictors.univariablePValue<univariableCutoff;

if (exist('maxPredictors', 'var'))
    lstPredictors = find(tablePredictors.isUsed); nPredictors = size(tablePredictors, 1);
    if (length(lstPredictors)>maxPredictors)
        isUsed = tablePredictors.isUsed;
        tmp = tablePredictors; tmp.iPredictor = (1:nPredictors)';
        tmp = sortrows(tmp,'explainedDeviance','descend');
        lstPredictors = sort(tmp.iPredictor(1:maxPredictors)); % mutation_frequency will be always selected and will be the last
        if (lstPredictors(end) ~= nPredictors), error('Mutation frequency not selected.'); end
        tablePredictors.isUsed = false(nPredictors, 1); tablePredictors.isUsed(lstPredictors) = true;
        fprintf('%d predictors originally predicted, %d top-explained-deviance selected: %s\n', sum(isUsed), sum(tablePredictors.isUsed), strjoin(tablePredictors.predictor(lstPredictors), '|'));
    end
end