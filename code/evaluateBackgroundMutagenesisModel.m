function results = evaluateBackgroundMutagenesisModel(runAgain, tissueName, tissueNameSV, biosampleABC, minCADD_PHRED, exclusionType, sProperties, qtlMF, univariableCutoff, maxPredictors)
%% Evaluates the gene-level background mutagenesis models on various datasets

if (~exist('univariableCutoff', 'var'))
    univariableCutoff = 0.001;
end
if (~exist('maxPredictors', 'var'))
    maxPredictors = 20;
end

plotFigures = false;
% runAgain = true;
typeName = 'SNVs_highCADD';
enhancerAnalysis = sProperties.enhancerAnalysis;

suffix = [tissueName, '_', biosampleABC, '_', enhancerAnalysis]; % evaluateBMM_test_
if (plotFigures)
    imagesPath = ['results/evaluation/', 'eval_', num2str(minCADD_PHRED), '_', num2str(maxPredictors), '_', num2str(1e4*qtlMF),'/']; createDir(imagesPath);
end

fileNameEvaluation = ['save/evaluateModel/evaluation_', suffix, '_', num2str(minCADD_PHRED), '_', exclusionType, '_max', num2str(maxPredictors), '_', num2str(1e4*qtlMF), '.mat']; createDir(fileparts(fileNameEvaluation));
if (runAgain || ~exist(fileNameEvaluation, 'file'))
    t1 = tic;
    fprintf('Computing pM %s...\n', suffix);
    %%
    levelOutputArguments = 3;
    [tableGenesNasserExpressed, tableGenes_pValues, stats, tableSamples, ~, ... % levelOutputArguments = 1
        ~, ~, ~, ~, matUniqueEnhancersGenes, ~, ... % levelOutputArguments = 2
        matGenesSamplesNMut_SNVs_highCADD, ~, ~, tableGenes_annotations, tableGenes_mean_trinucleotdies, tableUniqueEnhancers, ...
        ~, ~, ~, ~, tableUE_mean_trinucleotdies, ~, tableUE_annotations_hyperUE, ~] = ... % levelOutputArguments = 3
        computeMainAnalysis(false, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
    %% Definition of the 6 folds in the cross-validation is shared for all experiments and is split by chromosomes to have approximately equal number of enhancers/genes each.
    iFoldPerChromosome = [1,1, 2,2,2, 3,3,3, 4,4,4, 5,5,5,5,5, 6,6,6,6,6,6, NaN, NaN]';
    nFolds = max(iFoldPerChromosome);
    %% Preparation of gene-level data
    nMutSamplesPerRow = tableGenes_pValues.(['nMutSamplesInEnhancers_', typeName]);
    [tableDataForBMM, tableGenes_annotations] = prepare_tableDataForBMM(nMutSamplesPerRow, tableGenes_annotations, tableGenes_mean_trinucleotdies);
    s = struct();
    s.geneLevel.nRows = size(tableGenes_annotations, 1);
    s.geneLevel.chrNumeric = tableGenesNasserExpressed.chrNumeric;
    s.geneLevel.tableDataForBMM = tableDataForBMM;
    s.geneLevel.nPositionsInEnhancersCADD = tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF; % s.geneLevel.tableDataForBMM.nPositionsInEnhancers;
    %% Preparation of enhancer-level data
    nMutSamplesPerRow = tableUE_annotations_hyperUE.(['nMutSamplesInEnhancers_', typeName]);
    [tableDataForBMM, tableUE_annotations_hyperUE] = prepare_tableDataForBMM(nMutSamplesPerRow, tableUE_annotations_hyperUE, tableUE_mean_trinucleotdies);
    s.enhancerLevel.nRows = size(tableUE_annotations_hyperUE, 1);
    s.enhancerLevel.chrNumeric = tableUniqueEnhancers.chrNumeric;
    s.enhancerLevel.tableDataForBMM = tableDataForBMM;
    s.enhancerLevel.nPositionsInEnhancersCADD = tableUE_annotations_hyperUE.nTheoreticalMutations_PHRED_geqCUTOFF;
    %% Preparation of enhancer-windows data (for every enhancer row as above, we merge this and the subsequent enhancers together in a window - such that it's size is <= maxSize; so these represent sliding windows over enhancers)
    lstSizes = 10.^(3.5:.5:6.5); nSizes = length(lstSizes);
    for iSize = 1:nSizes
        maxSize = lstSizes(iSize);
        typeName = sprintf('enhancerWindows%d', iSize);
        [tableDataForBMM_averaged, tableDataForBMM_averaged_info] = slidingWindowEnhancers(tableUniqueEnhancers, tableUE_annotations_hyperUE.nTheoreticalMutations_PHRED_geqCUTOFF, tableDataForBMM, nMutSamplesPerRow, maxSize);
        s.(typeName).nRows = size(tableDataForBMM_averaged_info, 1);
        s.(typeName).chrNumeric = tableUniqueEnhancers.chrNumeric;
        s.(typeName).tableDataForBMM = tableDataForBMM_averaged;
        s.(typeName).nPositionsInEnhancersCADD = tableDataForBMM_averaged_info.nTheoreticalMutations_PHRED_geqCUTOFF;
    end
    %%
    clear nMutSamplesPerRow tableDataForBMM tableDataForBMM_averaged tableDataForBMM_averaged_info
    %%
    lstTypes = [{'geneLevel', 'enhancerLevel'}, strcat('enhancerWindows', strrep(cellstr(num2str((1:nSizes)')), ' ', ''))']; nTypes = length(lstTypes);
    for iType = 1:nTypes
        typeName = lstTypes{iType};
        if (~isinf(qtlMF))
            mutationFrequencyCutoffValue = quantile(s.(typeName).tableDataForBMM.mutationFrequency, qtlMF);
        else
            mutationFrequencyCutoffValue = Inf;
        end
        s.(typeName).isRelevantRow = ...
            s.(typeName).tableDataForBMM.mutationFrequency<mutationFrequencyCutoffValue & ...
            ~isnan(s.(typeName).tableDataForBMM.mutationFrequency) & ...
            s.(typeName).chrNumeric<23;
            %  REMOVED: & nPositionsInEnhancers>500 & nPositionsInEnhancers<=maxLength
            s.(typeName).truth = s.(typeName).tableDataForBMM.mutationFrequency.*s.(typeName).nPositionsInEnhancersCADD; 
            s.(typeName).truth(~s.(typeName).isRelevantRow) = NaN;
            s.(typeName).predictionsTraining = NaN*ones(s.(typeName).nRows, 1);
            s.(typeName).predictionsUnseen = NaN*ones(s.(typeName).nRows, 1);
            s.(typeName).rPearson_perFold_training = NaN*ones(nFolds, 1);
            s.(typeName).rPearson_perFold_unseen = NaN*ones(nFolds, 1);
            s.(typeName).rSpearman_perFold_training = NaN*ones(nFolds, 1);
            s.(typeName).rSpearman_perFold_unseen = NaN*ones(nFolds, 1);
    end
    clear mutationFrequencyCutoffValue
    %% First we select predictors based on the entire dataset, just for our information
    tablePredictors = computeSignificantUnivariablePredictors(s.geneLevel.tableDataForBMM(s.geneLevel.isRelevantRow,:), univariableCutoff, maxPredictors);
    enhancersPerChromosome = histcounts(s.geneLevel.chrNumeric, 1:25)';
    tableChromosomes = table(iFoldPerChromosome, enhancersPerChromosome);
    tableFolds = grpstats(tableChromosomes, 'iFoldPerChromosome', 'sum');
    %% In each fold, we build the gene-level model and use it to predict gene-level, enhancer-level, and enhancer-windows predictions on the training and unseen chromosomes
    for iFold = 1:nFolds
        %% The same gene-level model is used in all evaluations
        isThisFold = ismember(s.geneLevel.chrNumeric, find(iFoldPerChromosome == iFold));
        isTrainingFold = s.geneLevel.isRelevantRow & ~isThisFold;
        %isTestingFold = s.geneLevel.isRelevantRow & isThisFold;
        tablePredictorsThisFold = computeSignificantUnivariablePredictors(s.geneLevel.tableDataForBMM(isTrainingFold,:), univariableCutoff, maxPredictors);
        tablePredictors.(['isUsedFold', num2str(iFold)]) = tablePredictorsThisFold.isUsed;
        mdl = fitglm(s.geneLevel.tableDataForBMM(isTrainingFold, tablePredictorsThisFold.isUsed), 'linear', 'Distribution', 'poisson', 'DispersionFlag', true);
        tableFolds.mdl{iFold} = mdl;
        %% Evaluate mdl on training and unseen (test) data (gene-level, enhancer-level, enhancer-windows1, 2, ..., nSizes)
        for iType = 1:nTypes
            typeName = lstTypes{iType};
            isThisFold = ismember(s.(typeName).chrNumeric, find(iFoldPerChromosome == iFold));
            
            % unseen (test) chromosomes
            isTestingFold = s.(typeName).isRelevantRow & isThisFold;
            mult = s.(typeName).nPositionsInEnhancersCADD(isTestingFold);
            ytrueTest = s.(typeName).truth(isTestingFold);
            ypredTest = predict(mdl, s.(typeName).tableDataForBMM(isTestingFold,:)).*mult;
            s.(typeName).predictionsUnseen(isTestingFold) = ypredTest;
            s.(typeName).rPearson_perFold_unseen(iFold) = corr(ypredTest, ytrueTest, 'rows', 'pairwise', 'type', 'Pearson');
            s.(typeName).rSpearman_perFold_unseen(iFold) = corr(ypredTest, ytrueTest, 'rows', 'pairwise', 'type', 'Spearman');
            
            % training chromosomes
            isTrainingFold = s.(typeName).isRelevantRow & ~isThisFold;
            mult = s.(typeName).nPositionsInEnhancersCADD(isTrainingFold);
            ytrueTest = s.(typeName).truth(isTrainingFold);
            ypredTest = predict(mdl, s.(typeName).tableDataForBMM(isTrainingFold,:)).*mult;
            s.(typeName).predictionsTraining(isTrainingFold) = ypredTest;
            s.(typeName).rPearson_perFold_training(iFold) = corr(ypredTest, ytrueTest, 'rows', 'pairwise', 'type', 'Pearson');
            s.(typeName).rSpearman_perFold_training(iFold) = corr(ypredTest, ytrueTest, 'rows', 'pairwise', 'type', 'Spearman');
        end
    end
    clear isThisFold isTrainingFold isTestingFold mult ytrueTest ypredTest typeName iType iFold mdl tablePredictorsThisFold maxSize
    %%
    if (plotFigures)
        matWindowsFold_training = NaN*ones(nSizes, nFolds);
        matWindowsFold_unseen = NaN*ones(nSizes, nFolds);
        labels = cell(nSizes, 1);
        for iSize = 1:nSizes
            typeName = sprintf('enhancerWindows%d', iSize);
            matWindowsFold_training(iSize,:) = 100*s.(typeName).rPearson_perFold_training.^2;
            matWindowsFold_unseen(iSize,:) = 100*s.(typeName).rPearson_perFold_unseen.^2;
            %matWindowsFold_unseen(iSize,:) = 100*s.(typeName).rSpearman_perFold_unseen;
            labels{iSize} = sprintf('10^{%g}', log10(lstSizes(iSize)));
        end
        fig = createMaximisedFigure(1); hold on;
        plot(median(matWindowsFold_unseen, 2), '-', 'LineWidth', 2);
        boxplot(matWindowsFold_unseen', 'Color', 'k');
        set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex', 'FontSize', 16); grid on;
        xlabel('Size (bp)'); ylabel('Explained variance (%)'); ylim([0, 100]);
        sgtitle(tissueName);
        mySaveAs(fig, imagesPath, ['Fig1_',tissueName,'.png']);
        %%
        matWindowsFold_unseen = NaN*ones(2, nFolds);
        labels = cell(2, 1);
        matWindowsFold_unseen(1,:) = 100*s.enhancerLevel.rPearson_perFold_unseen.^2;
        matWindowsFold_unseen(2,:) = 100*s.geneLevel.rPearson_perFold_unseen.^2;
        labels{1} = sprintf('enhancer-levels (%.0f bp)', median(s.enhancerLevel.nPositionsInEnhancersCADD/3));
        labels{2} = sprintf('gene-levels (%.0f bp)', median(s.geneLevel.nPositionsInEnhancersCADD/3));
        fig = createMaximisedFigure(2); hold on;
        plot(median(matWindowsFold_unseen, 2), '-', 'LineWidth', 2);
        boxplot(matWindowsFold_unseen', 'Color', 'k');
        set(gca, 'XTickLabel', labels, 'TickLabelInterpreter', 'tex', 'FontSize', 16); grid on;
        ylabel('Explained variance (%)'); ylim([0, 100]);
        sgtitle(tissueName);
        mySaveAs(fig, imagesPath, ['Fig2_',tissueName,'.png']);
        %%
        typeName = 'geneLevel'; % geneLevel enhancerLevel enhancerWindows5
        yValuesPredAll = s.(typeName).predictionsUnseen;%./s.(typeName).nPositionsInEnhancersCADD;
        %yValuesPredAll = s.(typeName).predictionsTraining./s.(typeName).nPositionsInEnhancersCADD;
        yValuesTruthAll = s.(typeName).truth;%./s.(typeName).nPositionsInEnhancersCADD;
        xValuesAll = (tableUniqueEnhancers.min_pos0 + tableUniqueEnhancers.max_pos1)/2;
        chrNumber = s.(typeName).chrNumeric;

        %     corr(s.(typeName).predictionsUnseen./s.(typeName).nPositionsInEnhancersCADD, s.(typeName).truth./s.(typeName).nPositionsInEnhancersCADD, 'rows', 'pairwise', 'type', 'Pearson')
        %     corr(s.(typeName).predictionsUnseen, s.(typeName).truth, 'rows', 'pairwise', 'type', 'Pearson')

        tableLiterature = readtable(sProperties.TABLE_LITERATURE);
        lstHits = tableLiterature.geneSymbol(strcmpi(tableLiterature.tissues_candidate, tissueName));
        labelsAll = tableGenesNasserExpressed.geneName;
        %labelsAll = tableUniqueEnhancers.name;


        fig = createMaximisedFigure(3);

        for iChr = 1:22
            subplot(11,2,iChr); hold on;
            isOK = chrNumber == iChr;
            xValues = xValuesAll(isOK);
            yValuesPred = yValuesPredAll(isOK);
            yValuesTruth = yValuesTruthAll(isOK);
            labels = labelsAll(isOK);

            colourPredictions = [1,0,1];
            colourTruth = [0,0,1];
            plot(xValues, yValuesPred, '.', 'Color', (1+colourPredictions/2)/2);
            plot(xValues, yValuesTruth, '.', 'Color', (1+colourTruth/2)/2);

            mfFold = yValuesTruth./yValuesPred;

            minFold = 2.5;
            isProblematic = (mfFold > minFold) | isnan(yValuesPred);
            se = strel('disk',3); % Any stretches of < 5 consequtive truths will get removed
            isProblematicE = imopen(isProblematic, se);
            plot(xValues(isProblematicE), yValuesTruth(isProblematicE), 'o', 'Color', (1+colourTruth)/2);
            plot(xValues(isProblematicE), yValuesPred(isProblematicE), 'o', 'Color', (1+colourPredictions)/2);
            %         text(3e6+xValues(isProblematicE), yValuesTruth(isProblematicE), labels(isProblematicE), 'HorizontalAlignment', 'left');
            
            isHit = ismember(labels, lstHits);
            isProblematicHit = isProblematicE & isHit;
            text(1e6+xValues(isProblematicHit), yValuesTruth(isProblematicHit), labels(isProblematicHit), 'HorizontalAlignment', 'left', 'Color', 'r');
            isProblematicHit = ~isProblematicE & isHit;
            text(1e6+xValues(isProblematicHit), yValuesTruth(isProblematicHit), labels(isProblematicHit), 'HorizontalAlignment', 'left', 'Color', 'k');
            ylabel(sprintf('chr%d', iChr));
            %         nSteps1 = 10000; nSteps2 = 10000; degree = 35;
            %
            %         [x1, y1] = mySmoothing(xValues', yValuesPred', nSteps1, nSteps2, degree); colour = colourPredictions;
            %         plot(x1, y1, '-', 'MarkerFaceColor', (1+colour)/2, 'Color', colour, 'MarkerSize', 7, 'LineWidth', 1.5);
            %
            %         [x1, y1] = mySmoothing(xValues', yValuesTruth', nSteps1, nSteps2, degree); colour = colourTruth;
            %         plot(x1, y1, '-', 'MarkerFaceColor', (1+colour)/2, 'Color', colour, 'MarkerSize', 7, 'LineWidth', 1.5);
            drawnow;
        end
        sgtitle(tissueName);
        mySaveAs(fig, imagesPath, ['Fig3_',tissueName,'.png']);
        %%
        lstCols = 6:size(tablePredictors, 2);
        matPredictors = table2array(tablePredictors(:,lstCols));
        %imagesc((matPredictors)); colorbar;

        fig = createMaximisedFigure(4); colormap(flipud(lbmap(100, 'RedBlue')));
        matOverlap = NaN*ones(nFolds+1);
        for iPredictor = 1:nFolds+1
            for jPredictor = 1:nFolds+1
                matOverlap(iPredictor, jPredictor) = 100*mean(matPredictors(:,iPredictor) == matPredictors(:,jPredictor));
            end
        end
        imagesc((matOverlap)); colorbar; caxis([50, 100]);
        for iPredictor = 1:nFolds+1
            for jPredictor = 1:nFolds+1
                text(jPredictor, iPredictor, sprintf('%.0f%%', matOverlap(iPredictor, jPredictor)), 'HorizontalAlignment', 'center');
            end
        end
        labels = strrep(tablePredictors.Properties.VariableNames(lstCols), 'isUsed', '');
        set(gca, 'XTickLabel', labels, 'YTickLabels', labels, 'FontSize', 16);
        sgtitle(tissueName);
        mySaveAs(fig, imagesPath, ['Fig4_',tissueName,'.png']);
        clear fig
    end
    %% At the moment, we save all the data, to make sure we don't miss something that would need recomputing
    toc(t1)
    save(fileNameEvaluation);
end

fprintf('Loading %s...\n', fileNameEvaluation);
results = load(fileNameEvaluation);
% close fig;
if (exist('fig', 'var'))
    close(fig);
end

if (~isfield(results, 'matUniqueEnhancersGenes'))
    error('matUniqueEnhancersGenes does not exist');
end

if (~isfield(results, 'matGenesSamplesNMut_SNVs_highCADD'))
    error('matGenesSamplesNMut_SNVs_highCADD does not exist');
end