function [tableGencodeGenes, tableTissuesWithPancancer, sResults, tableMutations_candidate, matGeneGencodeIsCandidateMut, tableTissues] = loadData1()
%% Loads the main data (and runs the analysis if not precomputed).

saveFileData = 'save/data/data1.mat';
if (~exist(saveFileData, 'file'))
    tic
    %% 
    fprintf('Computing %s...\n', saveFileData);
    [tableTissues, sProperties] = loadParameters;
    runAgain = sProperties.runAgain; tailDirection = sProperties.tailDirection; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; mutTypeName = sProperties.mutTypeName; nGencodeGenes = sProperties.nGencodeGenes;
    %%
    nTissues = size(tableTissues, 1);
    tableTissues.nSamplesWGS = NaN*ones(nTissues, 1);
    tableTissues.nSamplesWGSandRNA = NaN*ones(nTissues, 1);
    tableTissues.pFisherCDG = NaN*ones(nTissues, 1);
    tableTissues.pFisherCDG_text = cell(nTissues, 1); tableTissues.pFisherCDG_text(:) = {''};
    tableTissues.enrichmentCDG = NaN*ones(nTissues, 1);
    tableTissues.PCAWG_projects = cell(nTissues, 1); tableTissues.PCAWG_projects(:) = {''};
    tableTissues.nUniqueEnhancers = NaN*ones(nTissues, 1);
    tableTissues.nEnhancerGenePairs = NaN*ones(nTissues, 1);
    tableTissues.nUsedGenes = NaN*ones(nTissues, 1);
    tableTissues.nDriverUpregulatedGenes = NaN*ones(nTissues, 1);
    tableTissues.nDriverDownregulatedGenes = NaN*ones(nTissues, 1);
    tableTissues.nCandidates = NaN*ones(nTissues, 1);
    tableTissues.nExpressedCDGs = NaN*ones(nTissues, 1);
    tableTissues.nExpressedONCOGENEs = NaN*ones(nTissues, 1);
    tableTissues.nExpressedTSGs = NaN*ones(nTissues, 1);
    %% These can be later removed:
    vecGencodeGeneTissue_isDriver = false(nGencodeGenes, 1);
    vecGencodeGeneTissue_isOncogene = false(nGencodeGenes, 1);
    vecGencodeGeneTissue_isTSG = false(nGencodeGenes, 1);
    %%
    matGencodeGeneTissue_isUsed = false(nGencodeGenes, nTissues);
    matGencodeGeneTissue_isCandidate = false(nGencodeGenes, nTissues);
    matGencodeGeneTissue_isCandidate_pM_only = false(nGencodeGenes, nTissues);
    matGencodeGeneTissue_isUpregulated = false(nGencodeGenes, nTissues);
    %%
    [~, tableGencodeGenes] = loadGenes(false);
    tableGencodeGenes.nMutations = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.nMutationsHighCADD = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.nMutSamples = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.nMutSamplesHighCADD = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.nMutSamplesHighCADD_hasRNA = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.iTissue = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.tissuePrint = cell(nGencodeGenes, 1); tableGencodeGenes.tissuePrint(:) = {''};
    tableGencodeGenes.pM = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.pE = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.pCombined = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.qCombined = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.sizeEffectM = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.sizeEffectE = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.isUP = false(nGencodeGenes, 1);
    %%
    tableTrinucleotides = readtable(sProperties.TABLE_TRINUCLEOTIDES); % 'data/tableTriNucl96.txt'
    %% Prognostic genes from The Human Protein Atlas
    tablePrognosticGenes = readtable(sProperties.GENES_SURVIVAL_PROTEIN_ATLAS);
    tablePrognosticGenes = grpstats(tablePrognosticGenes(:,{'GeneName', 'Cancer', 'prognostic_Favorable', 'unprognostic_Favorable', 'prognostic_Unfavorable', 'unprognostic_Unfavorable'}), {'GeneName', 'Cancer'}, 'min');
    lstTypesPrognostic = {'min_prognostic_Favorable', 'min_unprognostic_Favorable', 'min_prognostic_Unfavorable', 'min_unprognostic_Unfavorable'};
    [tablePrognosticGenes.minPValue, tablePrognosticGenes.iTypePrognostic] = min(table2array(tablePrognosticGenes(:,lstTypesPrognostic)), [], 2);
    tablePrognosticGenes.iTypePrognostic(isnan(tablePrognosticGenes.minPValue)) = 5;
    tablePrognosticGenes.isFavorable = tablePrognosticGenes.iTypePrognostic < 3;
    %%
    tableMutations_candidate = table();
    matGeneGencodeIsCandidateMut = [];
    sResults = cell(nTissues, 1);
    %%
    tic
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        tissueNameSV = tableTissues.tissueSV{iTissue};
        biosampleABC = tableTissues.biosampleABC{iTissue};
        levelOutputArguments = 3;
        %%
        [tableGenesNasserExpressed, tableGenes_pValues, ~, tableSamples, ~, ... % levelOutputArguments = 1 (Important: do not use tableGencodeGenes from here!)
            ~, ~, tableMutations, matMutationsEnhancers, matUniqueEnhancersGenes, matExpressionGenesSamples, ... % levelOutputArguments = 2
            matGenesSamplesNMut_SNVs_highCADD, ~, matCNV_genesSamples, ~, ~, tableUniqueEnhancers, ...
            ~, ~, ~, ~, ~, ~, tableUE_annotations_hyperUE, matGenesSamplesNMut_SNVs] = ... % levelOutputArguments = 3
        computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
        %% isCandidate vs isDriver
        isDriver = tableGenesNasserExpressed.isDriver;
        pM = tableGenes_pValues.(['p',xTestName,'_SNVs_highCADD']);        
        pE = tableGenes_pValues.(['p',yTestName,'_SNVs_highCADD']);
        sizeEffectE = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]);
        sizeEffectM = tableGenes_pValues.(['e',xTestName,'_',mutTypeName]);
        pCombined = combinePValues_EBM(pM,pE); %combinePValues(pM,pE,'Brown');  % Fisher
        qCombined = mafdr(pCombined, 'BHFDR', true); % ALTERNATIVELY: [~, ~, ~, qCombined] = fdr_bh(pCombined); % OLD:         
        
    
        tableGenesNasserExpressed.isUP = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]) > 0;
        tableGenesNasserExpressed.pM = pM;
        tableGenesNasserExpressed.pE = pE;
        tableGenesNasserExpressed.pCombined = pCombined;
        tableGenesNasserExpressed.qCombined = qCombined;

        P_cutoff = 0.05;
        Q_cutoff = 0.15;
        isP_M = pM < P_cutoff;
        isP_E = pE < P_cutoff;
        isCandidate = isP_M & isP_E & qCombined < Q_cutoff;
        isONCOGENE = contains(tableGenesNasserExpressed.role_CGC, 'oncogene');
        isTSG = contains(tableGenesNasserExpressed.role_CGC, 'TSG');
        isONCOGENE_notTSG = isONCOGENE & ~isTSG;
        isTSG_notONCOGENE = isTSG & ~isONCOGENE;

        sResults{iTissue}.tissuePrint = tableTissues.tissuePrint{iTissue};
        sResults{iTissue}.biosampleABC = biosampleABC;
        sResults{iTissue}.pM = pM;
        sResults{iTissue}.pE = pE;
        sResults{iTissue}.pCombined = pCombined;
        sResults{iTissue}.qCombined = qCombined;
        sResults{iTissue}.sizeEffectE = sizeEffectE;
        sResults{iTissue}.sizeEffectM = sizeEffectM;
        sResults{iTissue}.isUP = tableGenesNasserExpressed.isUP;
        sResults{iTissue}.isCandidate = isCandidate;
        sResults{iTissue}.isDriver = isDriver;
        sResults{iTissue}.isONCOGENE = isONCOGENE;
        sResults{iTissue}.isTSG = isTSG;
        sResults{iTissue}.geneName = tableGenesNasserExpressed.geneName;
        %
        [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isDriver, tailDirection, false);
        tableTissues.nSamplesWGS(iTissue) = sum(~tableSamples.isExcluded);
        tableTissues.nSamplesWGSandRNA(iTissue) = sum(~tableSamples.isExcluded & tableSamples.has_RNA);
        tableTissues.pFisherCDG(iTissue) = p;
        tableTissues.pFisherCDG_text{iTissue} = getPValueAsText(p);
        tableTissues.enrichmentCDG(iTissue) = enrichment;
        tableTissues.nObserved_candidateDrivers_CDGs(iTissue) = nObserved;
        tableTissues.nExpected_candidateDrivers_CDGs(iTissue) = nExpected;
        %% UP-DOWN vs ONCO-TSG
        isUP = isCandidate & tableGenesNasserExpressed.isUP;
        isDOWN = isCandidate & ~tableGenesNasserExpressed.isUP;
        %% onlyP_M (isP_M)
        [p, enrichment, nObserved, nExpected] = myFisherTest(isP_M, isDriver, tailDirection, false);
        tableTissues.onlyP_M_pFisherCDG(iTissue) = p;
        tableTissues.onlyP_M_pFisherCDG_text{iTissue} = getPValueAsText(p);
        tableTissues.onlyP_M_enrichmentCDG(iTissue) = enrichment;
        tableTissues.onlyP_M_nObserved_candidateDrivers_CDGs(iTissue) = nObserved;
        tableTissues.onlyP_M_nExpected_candidateDrivers_CDGs(iTissue) = nExpected;
        %%
        tableMutations.iSamplePCAWG = tableSamples.iSamplePCAWG(tableMutations.iSample);
        tableMutations.iTissue = iTissue*ones(size(tableMutations, 1), 1);
        %% tableMutations_candidate
        isUE_candidate = sum(matUniqueEnhancersGenes(:,isCandidate), 2)>0;
        isUE_candidateUP = sum(matUniqueEnhancersGenes(:,isCandidate & isUP), 2)>0;
        tableMutations.nUE = full(sum(matMutationsEnhancers, 2));
        tableMutations.iMut = (1:size(tableMutations, 1))';
        clear tableMutations_candidateOneTissue
        tableMutations.isCandidateDriver = false(size(tableMutations, 1), 1); tableMutations.isCandidateDriver(full(sum(matMutationsEnhancers(:,isUE_candidate), 2)>0)) = true;
        tableMutations.isCandidateDriverUP = false(size(tableMutations, 1), 1); tableMutations.isCandidateDriverUP(full(sum(matMutationsEnhancers(:,isUE_candidateUP), 2)>0)) = true;
        tableMutations.isCandidateDriverDOWN = tableMutations.isCandidateDriver & ~tableMutations.isCandidateDriverUP;
        tableMutations_candidateOneTissue = tableMutations(tableMutations.isCandidateDriver,:);
        tableMutations_candidateOneTissue.iTissue = iTissue*ones(size(tableMutations_candidateOneTissue, 1), 1);
        tableMutations_candidateOneTissue.iSamplePCAWG = tableSamples.iSamplePCAWG(tableMutations_candidateOneTissue.iSample);
        tableMutations_candidateOneTissue.candidateGenes = cell(size(tableMutations_candidateOneTissue, 1),1); tableMutations_candidateOneTissue.candidateGenes(:) = {''};
        tableMutations_candidateOneTissue.expressionBelowExpectation = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
        tableMutations_candidateOneTissue.expressionMedianWT = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
        tableMutations_candidateOneTissue.expressionThisMut = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
        tableMutations_candidateOneTissue.CNVThisMut = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
        matIsCandidateMutGene = false(size(tableMutations_candidateOneTissue, 1), size(tableGenesNasserExpressed, 1));
        for jMut = 1:size(tableMutations_candidateOneTissue, 1)
            iMut = tableMutations_candidateOneTissue.iMut(jMut);
            isUE = full(matMutationsEnhancers(iMut,:)==1);
            isGene = isCandidate & (sum(matUniqueEnhancersGenes(isUE,:), 1)>0)';
            matIsCandidateMutGene(jMut,isGene) = true;
            tableMutations_candidateOneTissue.candidateGenes{jMut} = strjoin(tableGenesNasserExpressed.geneName(isGene)');
        end

        for jMut = 1:size(tableMutations_candidateOneTissue, 1)
            for iGene = find(matIsCandidateMutGene(jMut,:))
                tableMutations_candidateOneTissue.expressionThisMut(jMut) = matExpressionGenesSamples(iGene, tableMutations_candidateOneTissue.iSample(jMut));
            end
        end

        for iGene = find(isCandidate)'
            isMut = matIsCandidateMutGene(:,iGene);
            isWT = ~tableSamples.isExcluded; isWT(tableMutations_candidateOneTissue.iSample(isMut)) = false;
            tableMutations_candidateOneTissue.expressionMedianWT(isMut) = median(matExpressionGenesSamples(iGene,isWT), 'omitnan');
        end
        [~,permMuts] = sortrows(tableMutations_candidateOneTissue,'candidateGenes','ascend');
        matIsCandidateMutGene = matIsCandidateMutGene(permMuts,:);
        tableMutations_candidateOneTissue = tableMutations_candidateOneTissue(permMuts,:);
        tableMutations_candidate = [tableMutations_candidate; tableMutations_candidateOneTissue];
        matGeneGencodeIsCandidateMut_oneTissue = false(nGencodeGenes, size(tableMutations_candidateOneTissue, 1));
        matGeneGencodeIsCandidateMut_oneTissue(tableGenesNasserExpressed.iGencode, :) = matIsCandidateMutGene';
        matGeneGencodeIsCandidateMut = [matGeneGencodeIsCandidateMut, matGeneGencodeIsCandidateMut_oneTissue];  
        %%
        tableTissues.pFisherCDG_up_oncogene(iTissue) = myFisherTest(isUP, isONCOGENE_notTSG, tailDirection, false);
        tableTissues.pFisherCDG_down_TSG(iTissue) = myFisherTest(isDOWN, isTSG_notONCOGENE, tailDirection, false);
        tableTissues.pFisherCDG_up_TSG(iTissue) = myFisherTest(isUP, isTSG_notONCOGENE, tailDirection, false);
        tableTissues.pFisherCDG_down_oncogene(iTissue) = myFisherTest(isDOWN, isONCOGENE_notTSG, tailDirection, false);
        %
        tableTissues.PCAWG_projects{iTissue} = strjoin(unique(tableSamples.projectName), '|');
        tableTissues.nUniqueEnhancers(iTissue) = size(matUniqueEnhancersGenes, 1);
        tableTissues.nEnhancerGenePairs(iTissue) = sum(matUniqueEnhancersGenes(:));
        tableTissues.nUsedGenes(iTissue) = size(matUniqueEnhancersGenes, 2);
        tableTissues.nDriverUpregulatedGenes(iTissue) = sum(isCandidate & isUP);
        tableTissues.nDriverDownregulatedGenes(iTissue) = sum(isCandidate & isDOWN);
        tableTissues.nCandidates(iTissue) = sum(isCandidate);
        tableTissues.nExpressedCDGs(iTissue) = sum(isDriver);
        tableTissues.nExpressedONCOGENEs(iTissue) = sum(isONCOGENE);
        tableTissues.nExpressedTSGs(iTissue) = sum(isTSG);
        %        
        lstGencodeCandidates = tableGenesNasserExpressed.iGencode(isCandidate);
        tableGencodeGenes.nMutations(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs(isCandidate, ~tableSamples.isExcluded), 2);
        tableGencodeGenes.nMutationsHighCADD(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs_highCADD(isCandidate, ~tableSamples.isExcluded), 2);
        tableGencodeGenes.nMutSamples(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs(isCandidate, ~tableSamples.isExcluded)>0, 2);
        tableGencodeGenes.nMutSamplesHighCADD(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs_highCADD(isCandidate, ~tableSamples.isExcluded)>0, 2);
        tableGencodeGenes.nMutSamplesHighCADD_hasRNA(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs_highCADD(isCandidate, tableSamples.has_RNA & ~tableSamples.isExcluded)>0, 2);
        %
        tableGencodeGenes.iTissue(lstGencodeCandidates) = iTissue;
        tableGencodeGenes.tissuePrint(lstGencodeCandidates) = tableTissues.tissuePrint(iTissue);
        tableGencodeGenes.pM(lstGencodeCandidates) = pM(isCandidate);
        tableGencodeGenes.pE(lstGencodeCandidates) = pE(isCandidate);
        tableGencodeGenes.pCombined(lstGencodeCandidates) = pCombined(isCandidate);
        tableGencodeGenes.qCombined(lstGencodeCandidates) = qCombined(isCandidate);
        tableGencodeGenes.sizeEffectM(lstGencodeCandidates) = sizeEffectM(isCandidate);
        tableGencodeGenes.sizeEffectE(lstGencodeCandidates) = sizeEffectE(isCandidate);
        tableGencodeGenes.isUP(lstGencodeCandidates) = isUP(isCandidate);
        %
        matGencodeGeneTissue_isUsed(tableGenesNasserExpressed.iGencode, iTissue) = true;
        matGencodeGeneTissue_isCandidate(tableGenesNasserExpressed.iGencode(isCandidate), iTissue) = true;
        matGencodeGeneTissue_isCandidate_pM_only(tableGenesNasserExpressed.iGencode(pM<P_cutoff), iTissue) = true;
        matGencodeGeneTissue_isUpregulated(tableGenesNasserExpressed.iGencode(tableGenesNasserExpressed.isUP), iTissue) = true;
        %% These can be later removed (now we have them for control purposes only)
        vecGencodeGeneTissue_isDriver(tableGenesNasserExpressed.iGencode(tableGenesNasserExpressed.isDriver)) = true;
        vecGencodeGeneTissue_isOncogene(tableGenesNasserExpressed.iGencode(isONCOGENE)) = true;
        vecGencodeGeneTissue_isTSG(tableGenesNasserExpressed.iGencode(isTSG)) = true;
        %% Save additional information for all driver-upregulated and driver-downregulated genes
        for iGene = find(isCandidate)'
            geneName = tableGenesNasserExpressed.geneName{iGene};
            saveForOneGeneVisualisation(tissueName, biosampleABC, geneName, pM(iGene), pE(iGene), qCombined(iGene), tableSamples, matCNV_genesSamples, matExpressionGenesSamples, matGenesSamplesNMut_SNVs_highCADD, ...
                tableMutations, matMutationsEnhancers, iGene, tableGencodeGenes, tableGenesNasserExpressed, matUniqueEnhancersGenes, tableUniqueEnhancers, tableUE_annotations_hyperUE, tableTrinucleotides);
        end
    end
    toc
    %save('workspace.mat', '-v7.3');
    %%
    isDriver = tableGencodeGenes.isDriver;
    isONCOGENE = tableGencodeGenes.isONCOGENE;
    isTSG = tableGencodeGenes.isTSG;

    tableTissuesWithPancancer = tableTissues([1:end,end,end],:);
    % We initialise the values of the newly added two rows as empty:
    colNames = tableTissuesWithPancancer.Properties.VariableNames;
    classPerColumn = cellfun(@(x) class(tableTissuesWithPancancer.(x)), colNames, 'UniformOutput', false);
    for iColumn = find(strcmp(classPerColumn, 'double'))
        tableTissuesWithPancancer.(colNames{iColumn})(end-1:end) = NaN;
    end
    for iColumn = find(strcmp(classPerColumn, 'cell'))
        tableTissuesWithPancancer.(colNames{iColumn})(end-1:end) = {''};
    end
    nRows = size(tableTissuesWithPancancer, 1);
    tableTissuesWithPancancer.observed_up_oncogene = NaN*ones(nRows, 1);
    tableTissuesWithPancancer.expected_up_oncogene = NaN*ones(nRows, 1);
    tableTissuesWithPancancer.observed_down_TSG = NaN*ones(nRows, 1);
    tableTissuesWithPancancer.expected_down_TSG = NaN*ones(nRows, 1);

    for iType = 1:2
        iRow = nTissues + iType;
        if (iType == 1)
            isRelTissue = ~contains(tableTissues.tissue, 'blood');
            tableTissuesWithPancancer.tissue{iRow} = 'pancancerWoBlood';
            tableTissuesWithPancancer.tissuePrint{iRow} = 'Pan-cancer wo Blood';
        else
            isRelTissue = true(nTissues, 1);
            tableTissuesWithPancancer.tissue{iRow} = 'pancancer';
            tableTissuesWithPancancer.tissuePrint{iRow} = 'Pan-cancer';
        end
        tableTissuesWithPancancer.tissueSV{iRow} = 'NA';
        tableTissuesWithPancancer.iTissue(iRow) = iRow;
        tableTissuesWithPancancer.nSamplesWGS(iRow) = sum(tableTissuesWithPancancer.nSamplesWGS(isRelTissue));
        tableTissuesWithPancancer.nSamplesWGSandRNA(iRow) = sum(tableTissuesWithPancancer.nSamplesWGSandRNA(isRelTissue));     
        %
        isUsedGene = sum(matGencodeGeneTissue_isUsed(:,isRelTissue), 2)>0;             % Used in at least one tissue
        isCandidate = sum(matGencodeGeneTissue_isCandidate(:,isRelTissue), 2)>0;       % Candidate in at least one tissue
        isCandidate_pM_only = sum(matGencodeGeneTissue_isCandidate_pM_only(:,isRelTissue), 2)>0;       % Candidate in at least one tissue - pM only

        isUP = sum(matGencodeGeneTissue_isCandidate(:,isRelTissue) & matGencodeGeneTissue_isUpregulated(:,isRelTissue), 2)>0;      % Upregulated-candidate in at least one tissue (can be both upregulated & downregulated in fact)
        isDOWN = sum(matGencodeGeneTissue_isCandidate(:,isRelTissue) & ~matGencodeGeneTissue_isUpregulated(:,isRelTissue), 2)>0;   % Downregulated-candidate in at least one tissue (can be both upregulated & downregulated in fact)

        if (~isequal(isDriver(isUsedGene), vecGencodeGeneTissue_isDriver(isUsedGene))), error('Wrong driver genes.'); end
        if (~isequal(isONCOGENE(isUsedGene), vecGencodeGeneTissue_isOncogene(isUsedGene))), error('Wrong oncogenes.'); end
        if (~isequal(isTSG(isUsedGene), vecGencodeGeneTissue_isTSG(isUsedGene))), error('Wrong TSG.'); end

        tableTissuesWithPancancer.nUsedGenes(iRow) = sum(isUsedGene);
        tableTissuesWithPancancer.nExpressedCDGs(iRow) = sum(isDriver);
        tableTissuesWithPancancer.nExpressedONCOGENEs(iRow) = sum(isONCOGENE);
        tableTissuesWithPancancer.nExpressedTSGs(iRow) = sum(isTSG);
        tableTissuesWithPancancer.nDriverUpregulatedGenes(iRow) = sum(isCandidate & isUP);
        tableTissuesWithPancancer.nDriverDownregulatedGenes(iRow) = sum(isCandidate & isDOWN);
        tableTissuesWithPancancer.nCandidates(iRow) = sum(isCandidate);
        %
        [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate(isUsedGene), isDriver(isUsedGene), tailDirection, false);
        tableTissuesWithPancancer.pFisherCDG(iRow) = p;
        tableTissuesWithPancancer.pFisherCDG_text{iRow} = getPValueAsText(p);
        tableTissuesWithPancancer.enrichmentCDG(iRow) = enrichment;
        tableTissuesWithPancancer.nObserved_candidateDrivers_CDGs(iRow) = nObserved;
        tableTissuesWithPancancer.nExpected_candidateDrivers_CDGs(iRow) = nExpected;
        %
        isOK = isUsedGene & ~(isONCOGENE&isTSG) & ~(isCandidate&~isUP);
        [p, enrichment, nObserved, nExpected] = myFisherTest(isUP(isOK), isONCOGENE(isOK) & ~isTSG(isOK), tailDirection, false)
        [p, enrichment, nObserved, nExpected] = myFisherTest(isUP(isUsedGene), isONCOGENE(isUsedGene) & ~isTSG(isUsedGene), tailDirection, false);
        tableTissuesWithPancancer.pFisherCDG_up_oncogene(iRow) = p;
        tableTissuesWithPancancer.observed_up_oncogene(iRow) = nObserved;
        tableTissuesWithPancancer.expected_up_oncogene(iRow) = nExpected;
        isOK = isUsedGene & ~(isONCOGENE&isTSG) & ~(isCandidate&~isDOWN);
        [p, enrichment, nObserved, nExpected] = myFisherTest(isDOWN(isOK), ~isONCOGENE(isOK) & isTSG(isOK), tailDirection, false)
        [p, enrichment, nObserved, nExpected] = myFisherTest(isDOWN(isUsedGene), isTSG(isUsedGene) & ~isONCOGENE(isUsedGene), tailDirection, false);
        tableTissuesWithPancancer.pFisherCDG_down_TSG(iRow) = p;
        tableTissuesWithPancancer.observed_down_TSG(iRow) = nObserved;
        tableTissuesWithPancancer.expected_down_TSG(iRow) = nExpected;
        %
        [p, enrichment, nObserved, nExpected] = myFisherTest(isUP(isUsedGene), isTSG(isUsedGene) & ~isONCOGENE(isUsedGene), tailDirection, false);
        tableTissuesWithPancancer.pFisherCDG_up_TSG(iRow) = p;
        [p, enrichment, nObserved, nExpected] = myFisherTest(isDOWN(isUsedGene), isONCOGENE(isUsedGene) & ~isTSG(isUsedGene), tailDirection, false);
        tableTissuesWithPancancer.pFisherCDG_down_oncogene(iRow) = p;
        %
        [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate_pM_only(isUsedGene), isDriver(isUsedGene), tailDirection, false);
        tableTissuesWithPancancer.onlyP_M_pFisherCDG(iRow) = p;
        tableTissuesWithPancancer.onlyP_M_pFisherCDG_text{iRow} = getPValueAsText(p);
        tableTissuesWithPancancer.onlyP_M_enrichmentCDG(iRow) = enrichment;
        tableTissuesWithPancancer.onlyP_M_nObserved_candidateDrivers_CDGs(iRow) = nObserved;
        tableTissuesWithPancancer.onlyP_M_nExpected_candidateDrivers_CDGs(iRow) = nExpected;
        %
        fprintf('%s: %d candidate genes upregulated, %d genes downregulated\n', tableTissuesWithPancancer.tissuePrint{iRow}, sum(isUP), sum(isDOWN));
        fprintf('%s: %d candidate driver genes upregulated, %d driver genes downregulated\n', tableTissuesWithPancancer.tissuePrint{iRow}, sum(isUP & isDriver), sum(isDOWN & isDriver));
        isOK = isUP & isONCOGENE;     fprintf('%s: %d genes upregulated oncogenes: %s\n', tableTissuesWithPancancer.tissuePrint{iRow}, sum(isOK), strjoin(tableGencodeGenes.geneSymbol(isOK)));
        isOK = isDOWN & isTSG;        fprintf('%s: %d genes downregulated TSGs: %s\n', tableTissuesWithPancancer.tissuePrint{iRow}, sum(isOK), strjoin(tableGencodeGenes.geneSymbol(isOK)));
    end
    tableTissuesWithPancancer.nCandidates_CDG_observed = tableTissuesWithPancancer.nObserved_candidateDrivers_CDGs;
    tableTissuesWithPancancer.nCandidates_CDG_expected = tableTissuesWithPancancer.nExpected_candidateDrivers_CDGs;
    %
    tableGencodeGenes.isUsedGene = isUsedGene;
    tableGencodeGenes.isCandidate = isCandidate;
    %%
    isNotBlood = ~contains(tableTissues.tissue, 'blood');
    tableGencodeGenes.isCandidateSolid = sum(matGencodeGeneTissue_isCandidate(:, isNotBlood), 2)>0;
    tableGencodeGenes.isUsedSolid = sum(matGencodeGeneTissue_isUsed(:, isNotBlood), 2)>0;
    %
    tableGencodeGenes.nTissues_candidate = sum(matGencodeGeneTissue_isCandidate, 2);
    if (max(tableGencodeGenes.nTissues_candidate)>1)
        tableGencodeGenes(tableGencodeGenes.nTissues_candidate>1,:)
        fprintf('%d gene (%s) is candidate in >1 tissue.\n', sum(tableGencodeGenes.nTissues_candidate>1), strjoin(tableGencodeGenes.geneSymbol(tableGencodeGenes.nTissues_candidate>1)));
        error('Some genes are candidates in more than one tissue! The code below needs to be adjusted.'); 
    end
    %%
    [tablePrognosticGenes.isInGencode, tablePrognosticGenes.indexGencode] = ismember(tablePrognosticGenes.GeneName, tableGencodeGenes.geneSymbol);
    tableGencodeGenes.pValuePrognostic_tissueMatched = NaN*ones(nGencodeGenes, 1);
    tableGencodeGenes.isFavorable_tissueMatched = NaN*ones(nGencodeGenes, 1);
    
    for iGG = find(tableGencodeGenes.isCandidate)'
        geneSymbol = tableGencodeGenes.geneSymbol{iGG};
        for iTissue = find(matGencodeGeneTissue_isCandidate(iGG,:))
            iRow = find(strcmp(tablePrognosticGenes.Cancer, tableTissues.tissuePrognostic{iTissue}) & strcmp(tablePrognosticGenes.GeneName, geneSymbol));
            if (~isempty(iRow) && ~isnan(tablePrognosticGenes.minPValue(iRow)))
                tableGencodeGenes.pValuePrognostic_tissueMatched(iGG) = tablePrognosticGenes.minPValue(iRow);
                tableGencodeGenes.isFavorable_tissueMatched(iGG) = tablePrognosticGenes.isFavorable(iRow);
            end
        end
    end
    %%
    tableGencodeGenes.typePrognostic_tissueMatched(tableGencodeGenes.isFavorable_tissueMatched==1 & tableGencodeGenes.pValuePrognostic_tissueMatched < 0.05) = {'F*'};
    tableGencodeGenes.typePrognostic_tissueMatched(tableGencodeGenes.isFavorable_tissueMatched==1 & tableGencodeGenes.pValuePrognostic_tissueMatched < 0.01) = {'F**'};
    tableGencodeGenes.typePrognostic_tissueMatched(tableGencodeGenes.isFavorable_tissueMatched==1 & tableGencodeGenes.pValuePrognostic_tissueMatched < 0.001) = {'F***'};
    tableGencodeGenes.typePrognostic_tissueMatched(tableGencodeGenes.isFavorable_tissueMatched==0 & tableGencodeGenes.pValuePrognostic_tissueMatched < 0.05) = {'U*'};
    tableGencodeGenes.typePrognostic_tissueMatched(tableGencodeGenes.isFavorable_tissueMatched==0 & tableGencodeGenes.pValuePrognostic_tissueMatched < 0.01) = {'U**'};
    tableGencodeGenes.typePrognostic_tissueMatched(tableGencodeGenes.isFavorable_tissueMatched==0 & tableGencodeGenes.pValuePrognostic_tissueMatched < 0.001) = {'U***'};
    %%
    toc
    %myPrintMemory
    %%
    isOK = sum(matGeneGencodeIsCandidateMut, 2)>0; % This matrix is too large. We save only the relevant genes (and indicate this in the tableGencodeGenes table)
    tableGencodeGenes.isIn_matGeneGencodeIsCandidateMut = isOK;
    matGeneGencodeIsCandidateMut = matGeneGencodeIsCandidateMut(isOK,:)==1; % To convert it to boolean
    %%
    if (ismember('QBiC_minPValue_increased', tableMutations_candidate.Properties.VariableNames))
        tableMutations_candidate.QBiC_minPValue_increased = [];
        tableMutations_candidate.QBiC_minPValue_decreased = [];
        tableMutations_candidate.QBiC_TF_increased = [];
        tableMutations_candidate.QBiC_TF_decreased = [];
        tableMutations_candidate.forQBiC_numeric = [];
        tableMutations_candidate.forQBiC = [];
        tableMutations_candidate.indexQBiC = [];
        tableMutations_candidate.driverGene = [];
        tableMutations_candidate.expressionBelowExpectation = [];
        tableMutations_candidate.CNVThisMut = [];
    end
    createDir(fileparts(saveFileData));
    save(saveFileData, 'tableGencodeGenes', 'tableTissuesWithPancancer', 'sResults', 'tableMutations_candidate', 'tableTissues', 'matGeneGencodeIsCandidateMut');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableGencodeGenes', 'tableTissuesWithPancancer', 'sResults', 'tableMutations_candidate', 'tableTissues', 'matGeneGencodeIsCandidateMut');
end
%