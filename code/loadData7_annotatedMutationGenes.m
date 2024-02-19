function [dataSupTables, tableGencodeGenes, tableMutations_candidate] = loadData7_annotatedMutationGenes(tableGencodeGenes, tableMutations_candidate, matGeneGencodeIsCandidateMut, dataDepMap, sProperties)
%% Annotates genes and candidate mutations with literature, closest gene per candidate mutation, promoter/close/distant candidate mutations etc.

saveFileData = [sProperties.DIRECTORY_SAVE, '/main/data7_annotatedMutationGenes_',sProperties.exclusionType,'.mat'];
if (~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
    %%
    tableGencodeGenes.isCandidateSolid = tableGencodeGenes.iTissue > 1;
    nGenes = size(tableGencodeGenes, 1);
    tableLiterature = readtable(sProperties.TABLE_LITERATURE); 
    lstGenesStrongSupport = tableLiterature.geneSymbol(tableLiterature.literatureEvidenceOncogene>2 | tableLiterature.literatureEvidenceTSG>1);
    [isOK, indexLiterature] = ismember(tableGencodeGenes.geneSymbol, tableLiterature.geneSymbol);
    tableGencodeGenes.literatureEvidenceOncogene = NaN*ones(nGenes, 1);
    tableGencodeGenes.literatureEvidenceTSG = NaN*ones(nGenes, 1);
    tableGencodeGenes.literatureEvidenceOncogene(isOK) = tableLiterature.literatureEvidenceOncogene(indexLiterature(isOK));
    tableGencodeGenes.literatureEvidenceTSG(isOK) = tableLiterature.literatureEvidenceTSG(indexLiterature(isOK));
    %%
    tableNasserAll = readtable(sProperties.GENES_ABC); % 'data/genes/Nasser2021_genesTSS.formatted.txt'
    tableNasserAll.Properties.VariableNames = {'nBiosamplesAll', 'chromosome', 'geneSymbol', 'TSS'};
    [isOK, index] = ismember(tableGencodeGenes.geneSymbol, tableNasserAll.geneSymbol);
    tableGencodeGenes.TSS = NaN*ones(nGenes, 1);
    tableGencodeGenes.TSS(isOK) = tableNasserAll.TSS(index(isOK));
    [isOK, index] = ismember(tableNasserAll.geneSymbol, tableGencodeGenes.geneSymbol);
    tableNasserAll.isInGencode = isOK;
    tableNasserAll.geneType = cell(size(tableNasserAll, 1), 1); tableNasserAll.geneType(:) = {''};
    tableNasserAll.geneType(isOK) = tableGencodeGenes.geneType2(index(isOK));

    nMuts = size(tableMutations_candidate, 1);
    tableMutations_candidate.geneSymbol_closestGene = cell(nMuts, 1); tableMutations_candidate.geneSymbol_closestGene(:) = {''};
    tableMutations_candidate.distance_closestGene = NaN*ones(nMuts, 1);
    tableMutations_candidate.isCloserToAnotherGene = false(nMuts, 1);
    tableMutations_candidate.geneSymbol_closestProteinCodingGene = cell(nMuts, 1); tableMutations_candidate.geneSymbol_closestProteinCodingGene(:) = {''};
    tableMutations_candidate.distance_closestProteinCodingGene = NaN*ones(nMuts, 1);
    tableMutations_candidate.isCloserToAnotherProteinCodingGene = false(nMuts, 1);
    tableMutations_candidate.distance_thisGene = NaN*ones(nMuts, 1); % If more candidates, then the minimal distance
    for iMut = 1:nMuts
        tmp = tableNasserAll.TSS;
        tmp(~strcmp(tableNasserAll.chromosome, tableMutations_candidate.chr{iMut})) = Inf; % This is to compare genes on this chromosome only
        [minDistance, iGene] = min(abs(tmp - tableMutations_candidate.pos1(iMut)));
        tableMutations_candidate.geneSymbol_closestGene{iMut} = tableNasserAll.geneSymbol{iGene};
        tableMutations_candidate.distance_closestGene(iMut) = minDistance;
        lstTargetCandidateGenes = strsplit(tableMutations_candidate.candidateGenes{iMut}, ' ');
        tableMutations_candidate.isCloserToAnotherGene(iMut) = ~ismember(tableNasserAll.geneSymbol{iGene}, lstTargetCandidateGenes);


        tmp(~strcmp(tableNasserAll.geneType, 'protein_coding')) = Inf; % This is to remove other than protein-coding genes (or genes outside GENCODE)
        [minDistance, iGene] = min(abs(tmp - tableMutations_candidate.pos1(iMut)));
        tableMutations_candidate.geneSymbol_closestProteinCodingGene{iMut} = tableNasserAll.geneSymbol{iGene};
        tableMutations_candidate.distance_closestProteinCodingGene(iMut) = minDistance;
        lstTargetCandidateGenes = strsplit(tableMutations_candidate.candidateGenes{iMut}, ' ');
        tableMutations_candidate.isCloserToAnotherProteinCodingGene(iMut) = ~ismember(tableNasserAll.geneSymbol{iGene}, lstTargetCandidateGenes);


        tmp = tableNasserAll.TSS(ismember(tableNasserAll.geneSymbol, lstTargetCandidateGenes));
        tableMutations_candidate.distance_thisGene(iMut) = min(abs(tmp - tableMutations_candidate.pos1(iMut)));
    end
    %%
    lstCandidateGenes_iGene = find(tableGencodeGenes.isIn_matGeneGencodeIsCandidateMut);

    nMutationGenePairs = sum(matGeneGencodeIsCandidateMut(:));
    tableMutationGenePairs = table();
    tableMutationGenePairs.iMutationCandidate = NaN*ones(nMutationGenePairs, 1);
    tableMutationGenePairs.geneSymbol_thisGene = cell(nMutationGenePairs, 1);
    tableMutationGenePairs.geneSymbol_closestGene = cell(nMutationGenePairs, 1);
    tableMutationGenePairs.geneSymbol_closestProteinCodingGene = cell(nMutationGenePairs, 1);
    tableMutationGenePairs.distance_thisGene = NaN*ones(nMutationGenePairs, 1);
    tableMutationGenePairs.distance_closestGene = NaN*ones(nMutationGenePairs, 1);
    tableMutationGenePairs.distance_closestProteinCodingGene = NaN*ones(nMutationGenePairs, 1);
    tableMutationGenePairs.isCloserToAnotherGene = false(nMutationGenePairs, 1);
    tableMutationGenePairs.isCloserToAnotherProteinCodingGene = false(nMutationGenePairs, 1);

    iPair = 0;
    for iCandidateGene = 1:size(matGeneGencodeIsCandidateMut, 1)
        iGene = lstCandidateGenes_iGene(iCandidateGene);
        geneSymbol = tableGencodeGenes.geneSymbol{iGene};

        isRelMut = matGeneGencodeIsCandidateMut(iCandidateGene,:)';
        for iMut = find(isRelMut)'
            iPair = iPair + 1;
            tableMutationGenePairs.iMutationCandidate(iPair) = iMut;
            tableMutationGenePairs.geneSymbol_thisGene{iPair} = geneSymbol;
            tableMutationGenePairs.geneSymbol_closestGene{iPair} = tableMutations_candidate.geneSymbol_closestGene{iMut};
            tableMutationGenePairs.geneSymbol_closestProteinCodingGene{iPair} = tableMutations_candidate.geneSymbol_closestProteinCodingGene{iMut};
            tableMutationGenePairs.distance_thisGene(iPair) = abs(tableMutations_candidate.pos1(iMut) - tableGencodeGenes.TSS(iGene));
            tableMutationGenePairs.distance_closestGene(iPair) = tableMutations_candidate.distance_closestGene(iMut);
            tableMutationGenePairs.distance_closestProteinCodingGene(iPair) = tableMutations_candidate.distance_closestProteinCodingGene(iMut);
            tableMutationGenePairs.isCloserToAnotherGene(iPair) = ~strcmp(tableMutationGenePairs.geneSymbol_thisGene{iPair}, tableMutationGenePairs.geneSymbol_closestGene{iPair});
            tableMutationGenePairs.isCloserToAnotherProteinCodingGene(iPair) = ~strcmp(tableMutationGenePairs.geneSymbol_thisGene{iPair}, tableMutationGenePairs.geneSymbol_closestProteinCodingGene{iPair});
        end
    end

    tableMutationGenePairs.pos0 = tableMutations_candidate.pos0(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.pos1 = tableMutations_candidate.pos1(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.isHighCADD = tableMutations_candidate.isHighCADD(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.isExcluded = tableMutations_candidate.isExcluded(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.iTissue = tableMutations_candidate.iTissue(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.chr = tableMutations_candidate.chr(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.context = tableMutations_candidate.context(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.ref = tableMutations_candidate.ref(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.alt = tableMutations_candidate.alt(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.CADD_PHRED = tableMutations_candidate.CADD_PHRED(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.isPlusStrand = tableMutations_candidate.isPlusStrand(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.isInFunSeq2 = tableMutations_candidate.isInFunSeq2(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.FunSeq2_motif_analysis = tableMutations_candidate.FunSeq2_motif_analysis(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.isMOTIFBR = tableMutations_candidate.isMOTIFBR(tableMutationGenePairs.iMutationCandidate);
    tableMutationGenePairs.isMOTIFG = tableMutations_candidate.isMOTIFG(tableMutationGenePairs.iMutationCandidate);    
    tableMutationGenePairs.tissue = dataDepMap.tableTissuesWithPancancer_DepMap.tissuePrint(tableMutationGenePairs.iTissue);
    %%
    tableGencodeGenes.nMutSamples_v2 = NaN*ones(nGenes, 1); 
    tableGencodeGenes.nMutSamplesHighCADD_v2 = NaN*ones(nGenes, 1); 

    tableGencodeGenes.nMutations = NaN*ones(nGenes, 1);
    tableGencodeGenes.nMutationsHighCADD = NaN*ones(nGenes, 1);
    tableGencodeGenes.pMutationsHighCADD_promoter = NaN*ones(nGenes, 1);
    tableGencodeGenes.pMutationsHighCADD_distant = NaN*ones(nGenes, 1);
    tableGencodeGenes.pMutationsHighCADD_isCloserToAnotherGene = NaN*ones(nGenes, 1);
    tableGencodeGenes.DepMap_meanDependencyTissueMatched = NaN*ones(nGenes, 1);
    tableGencodeGenes.DepMap_pDependentTissueMatched = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_nMut = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_pMutTFBS = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_pMutTFBS_break = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_pMutTFBS_gain = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_nMutSamples = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_pMutSamplesTFBS = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_pMutSamplesTFBS_break = NaN*ones(nGenes, 1);
    tableGencodeGenes.FunSeq2_pMutSamplesTFBS_gain = NaN*ones(nGenes, 1);



    for iCandidateGene = 1:size(matGeneGencodeIsCandidateMut, 1)
        iGene = lstCandidateGenes_iGene(iCandidateGene);
        geneSymbol = tableGencodeGenes.geneSymbol{iGene};

        isRelSNV = ~tableMutations_candidate.isIndel & ~tableMutations_candidate.isExcluded & matGeneGencodeIsCandidateMut(iCandidateGene,:)';
        lstSamples = unique(tableMutations_candidate.iSample(isRelSNV));
        tableGencodeGenes.nMutations(iGene) = sum(isRelSNV);
        tableGencodeGenes.nMutSamples_v2(iGene) = length(lstSamples);

        lstMuts = find(tableMutations_candidate.isHighCADD & isRelSNV);
        lstSamples = unique(tableMutations_candidate.iSample(lstMuts));
        tableGencodeGenes.nMutationsHighCADD(iGene) = length(lstMuts);
        tableGencodeGenes.nMutSamplesHighCADD_v2(iGene) = length(lstSamples);

        tableGencodeGenes.pMutationsHighCADD_promoter(iGene) = 100*mean(abs(tableMutations_candidate.pos1(lstMuts) - tableGencodeGenes.TSS(iGene)) <= 250);     % 250 bp
        tableGencodeGenes.pMutationsHighCADD_distant(iGene) = 100*mean(abs(tableMutations_candidate.pos1(lstMuts) - tableGencodeGenes.TSS(iGene)) > 20e3);      % 20 kbp
        tableGencodeGenes.pMutationsHighCADD_isCloserToAnotherGene(iGene) = 100*mean(~strcmp(tableMutations_candidate.geneSymbol_closestGene(lstMuts), geneSymbol));  % the closest gene (its TSS) is not this gene

        iTissue = tableGencodeGenes.iTissue(iGene);
        if (~isnan(iTissue))
            isDepMapThisGene = strcmp(dataDepMap.lstGenesDepMap, geneSymbol);
            if (sum(isDepMapThisGene)==1)
                tableGencodeGenes.DepMap_meanDependencyTissueMatched(iGene) = dataDepMap.matDepMatGenesTissues_average(isDepMapThisGene,iTissue);
                tableGencodeGenes.DepMap_pDependentTissueMatched(iGene) = dataDepMap.matDepMatGenesTissues_pAbove50(isDepMapThisGene,iTissue);
            end
        end

        lstMuts = find(tableMutations_candidate.isInFunSeq2 & isRelSNV);
        lstSamples = unique(tableMutations_candidate.iSample(lstMuts));
        isOK = tableMutations_candidate.isMOTIFBR(lstMuts)==1 | tableMutations_candidate.isMOTIFG(lstMuts)==1;
        tableGencodeGenes.FunSeq2_nMut(iGene) = length(lstMuts);
        tableGencodeGenes.FunSeq2_nMutSamples(iGene) = length(lstSamples);
        tableGencodeGenes.FunSeq2_pMutTFBS(iGene) = 100*mean(isOK);
        tableGencodeGenes.FunSeq2_pMutSamplesTFBS(iGene) = 100*length(unique(tableMutations_candidate.iSample(isOK)))/length(lstSamples);

        isOK = tableMutations_candidate.isMOTIFBR(lstMuts)==1;
        tableGencodeGenes.FunSeq2_pMutTFBS_break(iGene) = 100*mean(isOK);
        tableGencodeGenes.FunSeq2_pMutSamplesTFBS_break(iGene) = 100*length(unique(tableMutations_candidate.iSample(isOK)))/length(lstSamples);

        isOK = tableMutations_candidate.isMOTIFG(lstMuts)==1;
        tableGencodeGenes.FunSeq2_pMutTFBS_gain(iGene) = 100*mean(isOK);
        tableGencodeGenes.FunSeq2_pMutSamplesTFBS_gain(iGene) = 100*length(unique(tableMutations_candidate.iSample(isOK)))/length(lstSamples);
    end
    %% Check that mutatedSample counts match
    isOK = tableGencodeGenes.isCandidate;
    tableGencodeGenesCandidates = tableGencodeGenes(isOK,:);
    if (max(abs(tableGencodeGenesCandidates.nMutSamples - tableGencodeGenesCandidates.nMutSamples_v2))>0), error('Mutated sample counts do not match.'); end
    if (max(abs(tableGencodeGenesCandidates.nMutSamplesHighCADD - tableGencodeGenesCandidates.nMutSamplesHighCADD_v2))>0), error('Mutated sample counts do not match.'); end
    tableGencodeGenesCandidates.nMutSamples_v2 = [];
    tableGencodeGenesCandidates.nMutSamplesHighCADD_v2 = [];
    %%
    lstColGenes = {'chromosome', 'pos0', 'pos1', 'geneSymbol', 'geneNameGencode', 'strand', 'tissuePrint', 'isUP', 'pM', 'pE', 'qCombined', 'sizeEffectM', 'sizeEffectE', ...
        'nMutSamples', 'nMutSamplesHighCADD', 'nMutSamplesHighCADD_hasRNA', 'nMutations', 'nMutationsHighCADD', 'pMutationsHighCADD_promoter', 'pMutationsHighCADD_distant', 'pMutationsHighCADD_isCloserToAnotherGene', ...
        'isDriver', 'isONCOGENE', 'isTSG', 'typePrognostic_tissueMatched', 'pValuePrognostic_tissueMatched', 'isFavorable_tissueMatched', 'literatureEvidenceOncogene', 'literatureEvidenceTSG', ...
        'DepMap_meanDependencyTissueMatched', 'DepMap_pDependentTissueMatched', 'FunSeq2_nMut', 'FunSeq2_pMutTFBS', 'FunSeq2_pMutTFBS_break', 'FunSeq2_pMutTFBS_gain', 'FunSeq2_nMutSamples', 'FunSeq2_pMutSamplesTFBS', 'FunSeq2_pMutSamplesTFBS_break', 'FunSeq2_pMutSamplesTFBS_gain'};
    lstColGenesBlood = {'chromosome', 'pos0', 'pos1', 'geneSymbol', 'geneNameGencode', 'strand', 'tissuePrint', 'isUP', 'pM', 'pE', 'qCombined', 'sizeEffectM', 'sizeEffectE', ...
        'nMutSamples', 'nMutSamplesHighCADD', 'nMutSamplesHighCADD_hasRNA', 'nMutations', 'nMutationsHighCADD', 'pMutationsHighCADD_promoter', 'pMutationsHighCADD_distant', 'pMutationsHighCADD_isCloserToAnotherGene', ...
        'isDriver', 'isONCOGENE', 'isTSG', 'DepMap_meanDependencyTissueMatched', 'DepMap_pDependentTissueMatched'};
    lstColGenesDLBCL = {'chromosome', 'pos0', 'pos1', 'geneSymbol', 'geneNameGencode', 'strand', 'tissuePrint', 'isUP', 'pM', 'pE', 'qCombined', 'sizeEffectM', 'sizeEffectE', ...
        'nMutSamples', 'nMutSamplesHighCADD', 'nMutSamplesHighCADD_hasRNA', 'nMutations', 'nMutationsHighCADD', 'pMutationsHighCADD_promoter', 'pMutationsHighCADD_distant', 'pMutationsHighCADD_isCloserToAnotherGene', ...
        'isDriver', 'isONCOGENE', 'isTSG'}; % , 'DepMap_meanDependencyTissueMatched', 'DepMap_pDependentTissueMatched'
    lstColTissues = {'iTissue', 'tissuePrint', 'biosampleABC', 'tissuePrognostic', 'nSamplesWGS', 'nSamplesWGSandRNA', 'PCAWG_projects', 'nUniqueEnhancers', 'nEnhancerGenePairs', ...
        'nUsedGenes', 'nDriverUpregulatedGenes', 'nDriverDownregulatedGenes', 'nCandidates', 'nCandidates_CDG_observed', 'nCandidates_CDG_expected', 'enrichmentCDG', 'pFisherCDG'};
    lstColsMutationGenePairs = {'chr', 'pos0', 'pos1', 'ref', 'alt', 'tissue', 'geneSymbol_thisGene', 'distance_thisGene', 'distance_closestGene', 'distance_closestProteinCodingGene', 'isCloserToAnotherGene', 'isCloserToAnotherProteinCodingGene', ...
    'isHighCADD', 'CADD_PHRED', 'context', 'isPlusStrand', 'isInFunSeq2', 'isMOTIFBR', 'isMOTIFG'}; 
    %%
    dataSupTables.tableGencodeGenesCandidates = tableGencodeGenesCandidates;
    dataSupTables.tableMutationGenePairs = tableMutationGenePairs;
    dataSupTables.lstGenesStrongSupport = lstGenesStrongSupport;
    dataSupTables.lstColGenes = lstColGenes;
    dataSupTables.lstColGenesBlood = lstColGenesBlood;
    dataSupTables.lstColGenesDLBCL = lstColGenesDLBCL;
    dataSupTables.lstColTissues = lstColTissues;
    dataSupTables.lstColsMutationGenePairs = lstColsMutationGenePairs;
    %%
    toc
    %myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'dataSupTables', 'tableGencodeGenes', 'tableMutations_candidate');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'dataSupTables', 'tableGencodeGenes', 'tableMutations_candidate');
end
%%
% if (false)
%     tableMutationGenePairs.iSample = tableMutations_candidate.iSample(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.tissue = tableTissues_data1.tissuePrint(tableMutationGenePairs.iTissue);
%     tableMutationGenePairs.iPattern = tableMutations_candidate.iPattern(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.VAF = tableMutations_candidate.VAF(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.isPlusStrand = tableMutations_candidate.isPlusStrand(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.ref = tableMutations_candidate.ref(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.alt = tableMutations_candidate.alt(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.context = tableMutations_candidate.context(tableMutationGenePairs.iMutationCandidate);
%     tableMutationGenePairs.context3 = tableMutationGenePairs.context;
%     isOK = tableMutationGenePairs.iPattern>0 & tableMutationGenePairs.iPattern<97;
%     tableMutationGenePairs.context3(isOK) = cellfun(@(x) x(8:10), tableMutationGenePairs.context(isOK), 'UniformOutput', false);
% 
%     tmp1 = tableMutationGenePairs(tableMutationGenePairs.iTissue>1 & tableMutationGenePairs.isHighCADD & ~tableMutationGenePairs.isExcluded,:);
%     tmp2 = grpstats(tmp1(:,{'iSample', 'tissue', 'geneSymbol_thisGene'}), {'iSample', 'tissue', 'geneSymbol_thisGene'});
%     tmp2(tmp2.GroupCount>1,:)
%     writetable(tmp2(tmp2.GroupCount>1,:), 'multipleMutationsPerSamplePerGene.xlsx');
%     tmp1a = tableMutations_candidate(tableMutations_candidate.iSample == 38 & tableMutations_candidate.iTissue==4 & tableMutations_candidate.chrNumeric==20,:);
%     cellfun(@(x) x(8:10), tmp1a.context(1:end-1), 'UniformOutput', false)
%     tmp1a = tableMutations_candidate(tableMutations_candidate.iSample == 185 & tableMutations_candidate.iTissue==3 & tableMutations_candidate.chrNumeric>0,:);
%     cellfun(@(x) x(8:10), tmp1a.context, 'UniformOutput', false)
% 
%     %
%     tmp1.id = strrep(strcat(num2str(tmp1.iSample), '_', tmp1.tissue, '_', tmp1.geneSymbol_thisGene), ' ', '');
%     tmp3 = tmp1(ismember(tmp1.id, tmp2.Row(tmp2.GroupCount>1)),:);
%     fig = createMaximisedFigure(1);
%     histogram(tmp3.iPattern, 1:97);
%     %%
% end
