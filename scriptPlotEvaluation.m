% scriptPlotEvaluation.m
clear; clc; close all; addpath(genpath('code/')); rng(1);
%% Evaluations of the background mutagenesis model
%%
univariableCutoff = 0.001;
maxPredictors = 20;
qtlMF = 0.99;
imagesPath = ['results/evaluation/', 'eval_summary_', num2str(maxPredictors), '_', num2str(1e4*qtlMF),'/']; createDir(imagesPath);
%% Input file with all the parameters
inputPropertiesFile = 'inputParameters.properties';
%%
[tableTissues, sProperties] = loadParameters(inputPropertiesFile);
runAgain = sProperties.runAgain; tailDirection = sProperties.tailDirection; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; mutTypeName = sProperties.mutTypeName; nGencodeGenes = sProperties.nGencodeGenes;
%%
lstMinCADD_PHRED = [0, 10]; 
nBinCADD = length(lstMinCADD_PHRED);
nTissues = size(tableTissues, 1);
cResults = cell(nTissues, nBinCADD);
for iTissue = 1:nTissues
    tissueName = tableTissues.tissue{iTissue};
    tissueNameSV = tableTissues.tissueSV{iTissue};
    biosampleABC = tableTissues.biosampleABC{iTissue};
    for iBinCADD = 1:nBinCADD
        minCADD_PHRED = lstMinCADD_PHRED(iBinCADD);
        fprintf('\n=================== %s %d ===================\n', tissueName, minCADD_PHRED);
        sProperties.minCADD_PHRED = minCADD_PHRED;
        exclusionType = 'excludePOLE_MSI';
        sProperties.exclusionType = exclusionType;

        results = evaluateBackgroundMutagenesisModel(runAgain, tissueName, tissueNameSV, biosampleABC, minCADD_PHRED, exclusionType, sProperties, qtlMF, univariableCutoff, maxPredictors);
        cResults{iTissue, iBinCADD} = results;
    end
end
toc
%%
iTissue = 4;
tissueName = tableTissues.tissue{iTissue};
tissueNameSV = tableTissues.tissueSV{iTissue};
biosampleABC = tableTissues.biosampleABC{iTissue};
exclusionType = 'includeAll';
sProperties.exclusionType = exclusionType;
%
minCADD_PHRED = lstMinCADD_PHRED(1);
sProperties.minCADD_PHRED = minCADD_PHRED;
resultsCRC_allSamples = evaluateBackgroundMutagenesisModel(runAgain, tissueName, tissueNameSV, biosampleABC, minCADD_PHRED, exclusionType, sProperties, qtlMF, univariableCutoff, maxPredictors);
%
minCADD_PHRED = lstMinCADD_PHRED(2);
sProperties.minCADD_PHRED = minCADD_PHRED;
resultsCRC_allSamples_highCADD = evaluateBackgroundMutagenesisModel(runAgain, tissueName, tissueNameSV, biosampleABC, minCADD_PHRED, exclusionType, sProperties, qtlMF, univariableCutoff, maxPredictors);
%%
iBinCADD = 1;
fig = createMaximisedFigure(1, [0 0 40 20]);
nR = 2; nC = 4; xS = 0.8; yS = 0.8; xB = 0.05; yB = 0.05; xM = -0.01; yM = -0.01;
for iTissue = 1:nTissues
    results = cResults{iTissue, iBinCADD};
    nSizes = cResults{iTissue, iBinCADD}.nSizes;
    nFolds = cResults{iTissue, iBinCADD}.nFolds;
    lstSizes = cResults{iTissue, iBinCADD}.lstSizes;
    s = cResults{iTissue, iBinCADD}.s;
    matWindowsFold_training = NaN*ones(nSizes, nFolds);
    matWindowsFold_unseen = NaN*ones(nSizes, nFolds);
    labels = cell(nSizes, 1);
    for iSize = 1:nSizes
        typeName = sprintf('enhancerWindows%d', iSize);
        matWindowsFold_training(iSize,:) = 100*s.(typeName).rPearson_perFold_training.^2;
        matWindowsFold_unseen(iSize,:) = 100*s.(typeName).rPearson_perFold_unseen.^2;
        labels{iSize} = sprintf('10^{%g}', log10(lstSizes(iSize)));
    end
    myGeneralSubplot(nR,nC,iTissue,xS,yS,xB,yB,xM,yM); hold on;
    plot(median(matWindowsFold_unseen, 2), '-', 'LineWidth', 2);
    boxplot(matWindowsFold_unseen', 'Color', 'k');
    set(gca, 'XTickLabel', log10(lstSizes), 'TickLabelInterpreter', 'tex', 'FontSize', 14); 
    xlabel('Size (log_{10} bp)'); ylabel('Explained variance (%)'); ylim([0, 100]); grid on;
    title(tableTissues.tissuePrint{iTissue});
end
mySaveAs(fig, imagesPath, 'Fig1.png');
%% Windows: explained variance increases with increasing size
iBinCADD = 1;
matGeneLevel = NaN*ones(nTissues, nFolds);
tableTissues.medianSize = NaN*ones(nTissues, 1);
tableTissues.nSamplesIncluded = NaN*ones(nTissues, 1);
tableTissues.labelSamples = cell(nTissues, 1);
tableTissues.labelMedianSize = cell(nTissues, 1);
for iTissue = 1:nTissues
    s = cResults{iTissue, iBinCADD}.s;
    matGeneLevel(iTissue,:) = 100*s.geneLevel.rPearson_perFold_unseen.^2;
    tableTissues.nSamplesIncluded(iTissue) = sum(~cResults{iTissue, iBinCADD}.tableSamples.isExcluded);
    tableTissues.labelSamples{iTissue} = sprintf('%d s.', tableTissues.nSamplesIncluded(iTissue));
    tableTissues.medianSize(iTissue) = median(cResults{iTissue, iBinCADD}.tableGenesNasserExpressed.nPositionsInEnhancers); % s.geneLevel.nPositionsInEnhancersCADD/3
    tableTissues.labelMedianSize{iTissue} = [num2sepNumStr(tableTissues.medianSize(iTissue)), ' bp'];
end
%
tableShermanT4 = readtable('data/other/41587_2022_1353_MOESM3_ESM.xlsx', 'sheet', 'T4-Tiled region R2'); tableShermanT4.Properties.VariableNames = {'cohort', 'Mbp', 'kbp10', 'kbp10_auto', 'NBR', 'nSamples', 'nSNVs'};
lstRows = [22, 13, 23, 32, 31, 30, 16, 24, 5]; lstCols = 3:5;
matShermanT4 = 100*table2array(tableShermanT4(lstRows,lstCols));
%% Gene-level: explained variance ranges between 20-60%, comparable with previously published values.
fig = createMaximisedFigure(2); 
nR = 2; nC = 1; xS = 0.8; yS = 0.8; xB = 0.12; yB = 0.17; xM = -0.05; yM = -0.02;
colour_basic = [0,.5,1];
colour_CRC = [.25,.25,.75]; % [1,0,.5];
fontSizeSamples = 12;

myGeneralSubplot(nR,nC,1,.772,yS,xB,yB,xM,yM); hold on;
boxplot(matGeneLevel', 'Color', 'k'); 
h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    h1 = patch(get(h(j),'XData'),get(h(j),'YData'),colour_basic,'FaceAlpha',.5);
end
set(gca, 'XTickLabel', tableTissues.tissuePrint, 'TickLabelInterpreter', 'tex', 'FontSize', 16); grid on; box off;
vCRC = 100*resultsCRC_allSamples.s.geneLevel.rPearson_perFold_unseen.^2;
h2 = boxchart(4*ones(nFolds, 1), vCRC, 'BoxFaceColor', colour_CRC);
xShift = .3;
text(4+xShift, median(vCRC), sprintf('%d s.', sum(~resultsCRC_allSamples.tableSamples.isExcluded)), 'Color', colour_CRC, 'FontSize', fontSizeSamples);
xValues = (1:nTissues);
text(xValues+xShift, median(matGeneLevel, 2), tableTissues.labelSamples, 'Color', colour_basic, 'FontSize', fontSizeSamples);
text(xValues, 110+0*xValues, tableTissues.labelMedianSize, 'Color', .3*[1,1,1], 'FontSize', fontSizeSamples, 'HorizontalAlignment', 'center');
ylabel('Explained variance (%)'); xlim(.5+[0,nTissues]);
set(gca, 'YTick', 0:20:100);
ylim([0, 100]);
hL = legend([h1, h2], {'Main analysis', 'CRC including hypermutated samples'}, 'Location', 'NorthEast');
hL.Position(1) = .95 - 1.2*hL.Position(3);
hL.Position(2) = hL.Position(2) + .01; 

myGeneralSubplot(nR,nC,2,.87,yS,xB,yB,xM,yM); hold on; cmapDig = lines(3);
% bar(matShermanT4, 'EdgeColor', 'none')
xValues = repmat(.22*[-1, 0, 1], length(lstRows), 1) + repmat((1:length(lstRows))', 1, length(lstCols));
for iCol = 1:3
    plot(xValues(:,iCol), matShermanT4(:,iCol), 'o', 'MarkerEdgeColor', (1+cmapDig(iCol,:))/2, 'MarkerFaceColor', cmapDig(iCol,:), 'MarkerSize', 10);
end
grid on; box off;
set(gca, 'XTick', 1:length(lstRows), 'XTickLabel', strrep(tableShermanT4.cohort(lstRows), '_', ' '), 'FontSize', 16); ylabel('Explained variance (%)'); grid on;
xValues = 1:length(lstRows);
text(xValues, 95+0*xValues, num2str(tableShermanT4.nSamples(lstRows), '%d s.'), 'HorizontalAlignment', 'center', 'FontSize', fontSizeSamples, 'Color', .3*[1,1,1]);
hL = legend({'Dig Epi', 'Dig Epi&auto', 'Dig NBR'}, 'Location', 'East');
% hL.Position(1) = hL.Position(1) + .03;
hL.Position(1) = .95 - hL.Position(3);
hL.Position(2) = hL.Position(2) + .02;

fontSizeLetters = 26;
dim = [.01 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.01 .55 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
mySaveAs(fig, imagesPath, 'Fig2_new.png');
%% Predictors: which predictors get selected, which are independent predictors, overlapping predictors across folds and between 0-CADD and 10-CADD.
lstPredictors = cResults{1,1}.tablePredictors.predictor;
nPredictors = length(lstPredictors);
matPredictorsTissues_percSelected = NaN*ones(nPredictors, nTissues);
matPredictorsTissues_percIndependentPredictor = NaN*ones(nPredictors, nTissues);
matPredictorsTissues_percSelected_highCADD = NaN*ones(nPredictors, nTissues);
matPredictorsTissues_percIndependentPredictor_highCADD = NaN*ones(nPredictors, nTissues);
lstCols = 6+(1:nFolds);

matOverlapTissuesFoldPairs = NaN*ones(nTissues, nFolds^2);
matOverlapTissuesFolds_all_vs_highCADD = NaN*ones(nTissues, nFolds);

for iTissue = 1:nTissues
    matPredictorIsUsed = table2array(cResults{iTissue, 1}.tablePredictors(:,lstCols));
    matPredictorIsUsed_highCADD = table2array(cResults{iTissue, 2}.tablePredictors(:,lstCols));
    matPredictorsTissues_percSelected(:,iTissue) = 100*mean(matPredictorIsUsed, 2);
    matPredictorsTissues_percSelected_highCADD(:,iTissue) = 100*mean(matPredictorIsUsed_highCADD, 2);

    tmp = cell(2,1);
    for iBinCADD = 1:2
        tmp{iBinCADD}.matPredictorsFoldsMultiP = NaN*ones(nPredictors, nFolds);
        for iFold = 1:nFolds
            rowNames = cResults{iTissue, iBinCADD}.tableFolds.mdl{iFold}.Coefficients.Row;
            pValues = cResults{iTissue, iBinCADD}.tableFolds.mdl{iFold}.Coefficients.pValue;
            [isOK, iPredictor] = ismember(rowNames, lstPredictors);
            tmp{iBinCADD}.matPredictorsFoldsMultiP(iPredictor(isOK),iFold) = pValues(isOK);
        end
    end
    multivariableCutoff = 0.05;
    matPredictorsTissues_percIndependentPredictor(:,iTissue) = 100*mean(tmp{1}.matPredictorsFoldsMultiP < multivariableCutoff, 2);
    matPredictorsTissues_percIndependentPredictor_highCADD(:,iTissue) = 100*mean(tmp{2}.matPredictorsFoldsMultiP < multivariableCutoff, 2);

    matOverlap = NaN*ones(nFolds); lstUsedPredictors = 1:nPredictors-1;
    for iFold = 1:nFolds
        for jFold = 1:iFold-1
            matOverlap(iFold, jFold) = 100*mean(matPredictorIsUsed(lstUsedPredictors,iFold) == matPredictorIsUsed(lstUsedPredictors,jFold));
        end
    end
    %tableTissues.medianOverlapFolds(iTissue) = median(matOverlap(:), 'omitnan');
    matOverlapTissuesFoldPairs(iTissue,:) = matOverlap(:);

    for iFold = 1:nFolds
        matOverlapTissuesFolds_all_vs_highCADD(iTissue, iFold) = 100*mean(matPredictorIsUsed(lstUsedPredictors, iFold) == matPredictorIsUsed_highCADD(lstUsedPredictors, iFold));
    end
end
tablePredictors = table(lstPredictors);
tablePredictors.meanSelected = mean(matPredictorsTissues_percSelected, 2);
tablePredictors.meanSelected(end) = NaN; % here we ignore the response variable
tmp = matPredictorsTissues_percIndependentPredictor; tmp(isnan(tmp)) = 0;
tablePredictors.meanIndependent = mean(tmp, 2, 'omitnan');

tablePredictors.meanSelected_highCADD = mean(matPredictorsTissues_percSelected_highCADD, 2);
tablePredictors.meanSelected_highCADD(end) = NaN; % here we ignore the response variable
tmp = matPredictorsTissues_percIndependentPredictor_highCADD; tmp(isnan(tmp)) = 0;
tablePredictors.meanIndependent_highCADD = mean(tmp, 2, 'omitnan');

tablePredictors.printName = strrep(strrep(strrep(strrep(strrep(tablePredictors.lstPredictors, 'nPositionsInEnhancers', 'positions'), 'mean_mfInFlanks', 'flanking mutation frequency'), 'mean_replicationTiming', 'replication timing'), 'mean_GC', 'GC content'), 'mean_baseActivity', 'base activity');
%%
lstSuffixSave = {'_allCADD', '_highCADD'};
for iBinCADD = 1:2
    if (iBinCADD == 1)
        current_matPredictorsTissues_percSelected = matPredictorsTissues_percSelected;
        current_matPredictorsTissues_percIndependentPredictor = matPredictorsTissues_percIndependentPredictor;
        suffix = '';
    else
        current_matPredictorsTissues_percSelected = matPredictorsTissues_percSelected_highCADD;
        current_matPredictorsTissues_percIndependentPredictor = matPredictorsTissues_percIndependentPredictor_highCADD;
        suffix = '_highCADD';
    end
    fig = createMaximisedFigure(iBinCADD, [0 0 40 20]);
    nR = 1; nC = 2; xS = 0.8; yS = 0.9; xB = 0.12; yB = 0.17; xM = 0; yM = 0;
    cmap = flipud(lbmap(100, 'RedBlue'));

    myGeneralSubplot(nR,nC,1,xS,yS,xB,yB,xM,yM); hold on;
    for iPredictor = 1:nPredictors
        for iTissue = 1:nTissues
            percSelected = current_matPredictorsTissues_percSelected(iPredictor, iTissue);
            percIndep = current_matPredictorsTissues_percIndependentPredictor(iPredictor, iTissue);
            if (percSelected>0)
                colour = .5*[1,1,1];
                if (percIndep>0)
                    colour = cmap(round(percIndep), :);
                end
                plot(iTissue, iPredictor, 'o', 'MarkerFaceColor', colour, 'MarkerEdgeColor', 'none', 'MarkerSize', round(percSelected/10)+1);
            end
        end
    end
    xlim(.5+[0,nTissues+1.5]); ylim(.5+[0,nPredictors-1]);
    set(gca, 'TickLength', [0 0], 'XTick', 1:nTissues, 'XTickLabel', tableTissues.tissuePrint, 'YTick', 1:nPredictors, 'YTickLabel', tablePredictors.printName, 'TickLabelInterpreter', 'none');

    xVal = nTissues+1;
    yVal = nPredictors - 1;
    text(xVal, yVal, 'Selected');
    for iFold = 1:nFolds
        sizeValue = 100*iFold/nFolds; % for sizeValue = 10:15:100
        yVal = yVal - 1;
        plot(xVal, yVal, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', 'MarkerSize', round(sizeValue/10)+1);
        text(xVal+.3, yVal, sprintf('%d/%d folds', iFold, nFolds)); % text(xVal+.3, yVal, sprintf('%d %%', sizeValue));
    end
    yVal = yVal - 3;
    text(xVal, yVal, 'Independent');
    for iFold = 0:nFolds
        colourValue = 100*iFold/nFolds;
        yVal = yVal - 1;
        if (colourValue == 0)
            colour = .5*[1,1,1];
        else
            colour = cmap(round(colourValue), :);
        end
        plot(xVal, yVal, 'o', 'MarkerFaceColor', colour, 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
        text(xVal+.3, yVal, sprintf('%d/%d folds', iFold, nFolds), 'Color', colour); 
    end
    %
    myGeneralSubplot(nR,nC,2,xS,yS,xB,yB,xM,yM); hold on;
    plot(tablePredictors.(['meanSelected', suffix]), tablePredictors.(['meanIndependent', suffix]), 'o', 'MarkerFaceColor', colour_basic)
    xValues = tablePredictors.(['meanSelected', suffix]); yValues = tablePredictors.(['meanIndependent', suffix]); labels = tablePredictors.printName; fontSize= 12;
    [xValuesText, yValuesText] = labelRepelSimple(xValues, yValues, labels, fontSize);
    text(xValuesText, yValuesText, labels, 'Interpreter','none');

    set(gca, 'FontSize', 16);
    xlabel('Average selected (%)'); % Percentage of folds that the predictor has been selected, averaged across all tissues.
    ylabel('Average independent (%)'); % Percentage of folds that the predictor has been independent (p<0.05 in multivariable model), averaged across all tissues.
    mySaveAs(fig, imagesPath, ['Fig3_', lstSuffixSave{iBinCADD}, '.png']);
end
%% What is the overlap between selected predictors across the folds?
fig = createMaximisedFigure(4, [0 0 25 15]);
boxplot(matOverlapTissuesFoldPairs', 'Color', 'k');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    h1 = patch(get(h(j),'XData'),get(h(j),'YData'),colour_basic,'FaceAlpha',.5);
end
set(gca, 'XTickLabel', tableTissues.tissuePrint, 'TickLabelInterpreter', 'tex', 'FontSize', 16); grid on; box off;
ylabel({'Overlap of selected predictors', 'between folds (%)'}); xlim(.5+[0,nTissues]); ylim([0, 100]);
mySaveAs(fig, imagesPath, 'Fig4.png');
%% What is the overlap between selected predictors in all-CADD vs. high-CADD models? (one value per fold)
fig = createMaximisedFigure(5, [0 0 25 15]);
boxplot(matOverlapTissuesFolds_all_vs_highCADD', 'Color', 'k');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    h1 = patch(get(h(j),'XData'),get(h(j),'YData'),colour_basic,'FaceAlpha',.5);
end
set(gca, 'XTickLabel', tableTissues.tissuePrint, 'TickLabelInterpreter', 'tex', 'FontSize', 16); grid on; box off;
ylabel({'Overlap of selected predictors', 'between all-CADD vs. high-CADD models (%)'}); xlim(.5+[0,nTissues]); ylim([0, 100]);
mySaveAs(fig, imagesPath, 'Fig5.png');
%% Gene-level high-CADD: explained variance ranges between 6-33%, somewhat lower than in all-CADD case, as expected given the lower number of mutations.
iBinCADD = 2;
matGeneLevel_highCADD = NaN*ones(nTissues, nFolds);
for iTissue = 1:nTissues
    s = cResults{iTissue, iBinCADD}.s;
    matGeneLevel_highCADD(iTissue,:) = 100*s.geneLevel.rPearson_perFold_unseen.^2;
end
fig = createMaximisedFigure(6, [0 0 25 15]); hold on;
boxplot(matGeneLevel_highCADD', 'Color', 'k'); 
h = findobj(gca,'Tag','Box'); 
for j=1:length(h)
    h1 = patch(get(h(j),'XData'),get(h(j),'YData'),colour_basic,'FaceAlpha',.5);
end
set(gca, 'XTickLabel', tableTissues.tissuePrint, 'TickLabelInterpreter', 'tex', 'FontSize', 16); grid on; box off;
vCRC = 100*resultsCRC_allSamples_highCADD.s.geneLevel.rPearson_perFold_unseen.^2;
h2 = boxchart(4*ones(nFolds, 1), vCRC, 'BoxFaceColor', colour_CRC);
xShift = .3;
text(4+xShift, median(vCRC), sprintf('%d s.', sum(~resultsCRC_allSamples_highCADD.tableSamples.isExcluded)), 'Color', colour_CRC, 'FontSize', fontSizeSamples);
xValues = (1:nTissues);
text(xValues+xShift, median(matGeneLevel_highCADD, 2), tableTissues.labelSamples, 'Color', colour_basic, 'FontSize', fontSizeSamples);
ylabel('Explained variance high-CADD (%)'); xlim(.5+[0,nTissues]);
hL = legend([h1, h2], {'Main analysis', 'CRC including hypermutated samples'}, 'Location', 'NorthOutside');
mySaveAs(fig, imagesPath, 'Fig6_new.png');
%% Using the actual final model (used throughout the rest of the paper) and its sizeEffectM (observed/expected) estimates, 
% we check that the candidate hits are not due to the region (surrounding genes) having high sizeEffectM due to non-overlapping samples.
% This would represent a situation when the model wrongly underestimates the mutation rate in that region. 
% In other words, a situation when the model thinks that these regions are under positive selection,
% while they have just increased background mutation rate that is not captured by the model.
[~, ~, sResults] = loadData1();

for iBinCADD = 2%1:2
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        results = cResults{iTissue, iBinCADD};
        if (~isequal(results.tableGenesNasserExpressed.geneName, sResults{iTissue}.geneName)), error('Gene names do not match.'); end
        if (max(results.tableGenes_pValues.pM_fullModel_SNVs_highCADD - sResults{iTissue}.pM)>0 && ~isequal(results.tableGenes_pValues.pM_fullModel_SNVs_highCADD, sResults{iTissue}.pM)), error('P-values pM do not match.'); end
        isCandidate = sResults{iTissue}.isCandidate;

        tmp_tableGenes = results.tableGenesNasserExpressed(:,{'chrNumeric', 'pos0_GENCODE', 'pos1_GENCODE', 'geneName', 'nPositionsInEnhancers', 'iGencode', 'geneType'});
        tmp_tableGenes.sizeEffectM = results.tableGenes_pValues.eM_fullModel_SNVs_highCADD; % fold-change observed/expected
        tmp_tableGenes.pM = results.tableGenes_pValues.pM_fullModel_SNVs_highCADD;
        tmp_tableGenes.nMutSamples = results.tableGenes_pValues.nMutSamplesInEnhancers_SNVs_highCADD;
        tmp_tableGenes.nTheoreticalMutations = results.tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF;
        tmp_tableGenes.isCandidate = isCandidate;
        tmp_tableGenes.mf = tmp_tableGenes.nMutSamples ./ tmp_tableGenes.nTheoreticalMutations;
        tmp_tableGenes.expected_mf = results.stats.SNVs_highCADD.expected_mf;
        tmp_tableGenes.observed_mf = results.stats.SNVs_highCADD.observed_mf;
        tmp_tableGenes.foldChange = results.stats.SNVs_highCADD.observed_mf ./ results.stats.SNVs_highCADD.expected_mf;
        tmp_tableGenes.log2foldChange = log(tmp_tableGenes.foldChange);
        tmp_tableGenes.excess_mutSamples = (results.stats.SNVs_highCADD.observed_mf - results.stats.SNVs_highCADD.expected_mf) .* tmp_tableGenes.nTheoreticalMutations;
        tmp_tableGenes.excess_mf = (results.stats.SNVs_highCADD.observed_mf - results.stats.SNVs_highCADD.expected_mf);
        [tmp_tableGenes, perm] = sortrows(tmp_tableGenes,{'chrNumeric', 'pos0_GENCODE', 'pos1_GENCODE'},'ascend');
        matGenesSamplesNMut_SNVs_highCADD = results.matGenesSamplesNMut_SNVs_highCADD(perm,~results.tableSamples.isExcluded);
        maxDistance = 1e5; 
        matToPlot = [];
        jGene = 1;
        lstGenes = find(tmp_tableGenes.isCandidate)';
        labels = tmp_tableGenes.geneName(lstGenes);
        isAboveThreshold = false(length(labels), 1);
        for iGene = lstGenes
            tmp_tableGenes.nMutSamples_unique = sum((matGenesSamplesNMut_SNVs_highCADD - matGenesSamplesNMut_SNVs_highCADD(iGene,:))>0, 2);
            tmp_tableGenes.mf_unique = tmp_tableGenes.nMutSamples_unique ./ tmp_tableGenes.nTheoreticalMutations;
            tmp_tableGenes.fc_unique = tmp_tableGenes.mf_unique ./ tmp_tableGenes.expected_mf;
            tmp_tableGenes.distance = min(abs([...
                tmp_tableGenes.pos0_GENCODE(iGene) - tmp_tableGenes.pos0_GENCODE, tmp_tableGenes.pos0_GENCODE(iGene) - tmp_tableGenes.pos1_GENCODE, ...
                tmp_tableGenes.pos0_GENCODE - tmp_tableGenes.pos1_GENCODE(iGene), tmp_tableGenes.pos1_GENCODE - tmp_tableGenes.pos1_GENCODE(iGene)]), [], 2);
            isNearbyGene = tmp_tableGenes.chrNumeric == tmp_tableGenes.chrNumeric(iGene) & tmp_tableGenes.distance < maxDistance;
            isOK = isNearbyGene & tmp_tableGenes.iGencode ~= tmp_tableGenes.iGencode(iGene);
            matToPlot = [matToPlot; tmp_tableGenes.fc_unique(isOK), jGene + 0*tmp_tableGenes.fc_unique(isOK)];
            if (median(tmp_tableGenes.fc_unique(isOK))>2)
                isAboveThreshold(jGene) = true;
                fprintf('%s %s: median fc_unique=%.2fx\n', tissueName, labels{jGene}, median(tmp_tableGenes.fc_unique(isOK)));
            end
            jGene = jGene + 1;
        end
        sResults{iTissue}.matToPlot = matToPlot;
        sResults{iTissue}.labels = labels;
        sResults{iTissue}.genesAboveThreshold = labels(isAboveThreshold);
        sResults{iTissue}.nMutSamples = tmp_tableGenes.nMutSamples(lstGenes);
        sResults{iTissue}.expected_mutSamples = tmp_tableGenes.expected_mf(lstGenes).*tmp_tableGenes.nTheoreticalMutations(lstGenes);
        sResults{iTissue}.excess_mutSamples = tmp_tableGenes.excess_mutSamples(lstGenes);
        sResults{iTissue}.expected_mf = tmp_tableGenes.expected_mf(lstGenes);
        sResults{iTissue}.observed_mf = tmp_tableGenes.observed_mf(lstGenes);
        sResults{iTissue}.excess_mf = tmp_tableGenes.excess_mf(lstGenes);
        sResults{iTissue}.foldChange = tmp_tableGenes.foldChange(lstGenes);
    end
    %%
    fig = createMaximisedFigure(2);
    nR = 4; nC = 2; iS = 1; xS = 0.9; yS = 0.45; xB = 0.08; yB = 0.08; xM = -0.02; yM = -0.05;
    for iTissue = 2:nTissues
        labels = sResults{iTissue}.labels;
        matToPlot = sResults{iTissue}.matToPlot;
        myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;
        xLimVal = 0.5+[0,length(labels)];
        plot(xLimVal, 2*[1,1], '-', 'LineWidth', 2, 'Color', .5*[1,1,1]);
        lstAboveThreshold = find(ismember(labels, sResults{iTissue}.genesAboveThreshold));
        isOK = ismember(matToPlot(:,2), lstAboveThreshold);
        colour = [0, .25, 1]; %[0.8500    0.3250    0.0980]; %[1,.25,0]
        hB = boxchart(matToPlot(~isOK,2), matToPlot(~isOK,1), 'MarkerStyle', 'none', 'BoxFaceColor', colour);
        nRows = size(matToPlot, 1);
        xValues = matToPlot(:,2) + randn(nRows, 1)/10-1/20;
        yValues = matToPlot(:,1);
        plot(xValues, yValues, 'o', 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
        colour = [1, 0, 0];
        plot(xValues(isOK), yValues(isOK), 'o', 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
        hB = boxchart(matToPlot(isOK,2), matToPlot(isOK,1), 'MarkerStyle', 'none', 'BoxFaceColor', colour);
        set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'FontSize', 12);
        ylabel({'Flanking genes', 'observed/expected'});
        title(tableTissues.tissuePrint{iTissue});
        drawnow;
    end
    mySaveAs(fig, imagesPath, ['Fig_posthoc_CADD', num2str(lstMinCADD_PHRED(iBinCADD)), '_distance', num2str(maxDistance), '.png']);
    %%
    fig = createMaximisedFigure(2, [0 0 35 15]);
    iTissue = 1; labels = sResults{iTissue}.labels;
    matToPlot = sResults{iTissue}.matToPlot;
    axes('Position', [.08, .25, .9, .65]); hold on;
    xLimVal = 0.5+[0,length(labels)];
    plot(xLimVal, 2*[1,1], '-', 'LineWidth', 2, 'Color', .5*[1,1,1]);
    lstAboveThreshold = find(ismember(labels, sResults{iTissue}.genesAboveThreshold));
    isOK = ismember(matToPlot(:,2), lstAboveThreshold);
    colour = [0, .25, 1]; 
    hB = boxchart(matToPlot(~isOK,2), matToPlot(~isOK,1), 'MarkerStyle', 'none', 'BoxFaceColor', colour);
    nRows = size(matToPlot, 1);
    xValues = matToPlot(:,2) + randn(nRows, 1)/10-1/20;
    yValues = matToPlot(:,1);
    plot(xValues, yValues, 'o', 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
    colour = [1, 0, 0];
    plot(xValues(isOK), yValues(isOK), 'o', 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
    hB = boxchart(matToPlot(isOK,2), matToPlot(isOK,1), 'MarkerStyle', 'none', 'BoxFaceColor', colour);
    set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'FontSize', 12);
    ylabel('Flanking genes observed/expected');
    title(tableTissues.tissuePrint{iTissue});
    drawnow;
    mySaveAs(fig, imagesPath, ['Fig_posthoc_CADD', num2str(lstMinCADD_PHRED(iBinCADD)), '_distance', num2str(maxDistance), '_blood.png']);
    %%
    for iType = 1:2
        fig = createMaximisedFigure(2);
        nR = 4; nC = 2; iS = 1; xS = 0.9; yS = 0.65; xB = 0.05; yB = 0.05; xM = -0.02; yM = -0.05;
        for iTissue = 2:nTissues
            labels = sResults{iTissue}.labels;
            if (iType == 1)
                vToPlot = [sResults{iTissue}.nMutSamples, sResults{iTissue}.expected_mutSamples, sResults{iTissue}.excess_mutSamples];
                yLabelText = 'Mutated samples';
                suffix = 'mutSamples';
            else
                vToPlot = [sResults{iTissue}.observed_mf, sResults{iTissue}.expected_mf, sResults{iTissue}.excess_mf]/tableTissues.nSamplesIncluded(iTissue);
                yLabelText = 'Mutation frequency';
                suffix = 'mf';
            end
            myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;
            bar(vToPlot, 'EdgeColor', 'none');
            xLimVal = 0.5+[0,length(labels)]; xlim(xLimVal);
            set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels);
            ylabel(yLabelText);
            title(tableTissues.tissuePrint{iTissue});
            drawnow;
            if (iTissue == nTissues)
                hL = legend({'observed', 'expected', 'excess'}, 'Location', 'East');
                hL.Position(1) = hL.Position(1) + 2*hL.Position(3);
            end
        end
        mySaveAs(fig, imagesPath, ['Fig_excess_CADD', num2str(lstMinCADD_PHRED(iBinCADD)), '_', suffix, '.png']);
    end
end
