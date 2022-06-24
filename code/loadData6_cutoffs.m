function dataCutoffs = loadData6_cutoffs()

saveFileData = 'save/data/data6_cutoffs.mat'; 
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
    %%
    lstCutoffsPM = [0.02:0.005:0.1, 1.01];
    lstCutoffsPE = [0.02:0.005:0.1, 1.01];
    nCutoffPM = length(lstCutoffsPM);
    nCutoffPE = length(lstCutoffsPE);
    %%
    matGencodeGeneTissue_isCandidate_pM_pE = false(nGencodeGenes, nTissues, nCutoffPE, nCutoffPM);
    matGencodeGeneTissue_isCandidate_pM_pE_FDR = false(nGencodeGenes, nTissues, nCutoffPE, nCutoffPM);
    matGencodeGeneTissue_isUsed = false(nGencodeGenes, nTissues, 1);
    matGencodeGeneTissue_isDriver = false(nGencodeGenes, 1);
    %%
    tic
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        tissueNameSV = tableTissues.tissueSV{iTissue};
        biosampleABC = tableTissues.biosampleABC{iTissue};
        fprintf('\n=================== %s ===================\n', tissueName);
        levelOutputArguments = 1; % We require only the minimal number of output arguments, to speed everything up.
        [tableGenesNasserExpressed, tableGenes_pValues] = computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
        %
        pM = tableGenes_pValues.(['p',xTestName,'_',mutTypeName]); 
        pE = tableGenes_pValues.(['p',yTestName,'_',mutTypeName]); 
        pCombined = combinePValues_EBM(pM,pE);
        cutoffQ = 0.15;
        qCombined = mafdr(pCombined, 'BHFDR', true);
        pM(isnan(pM)) = 1;
        pE(isnan(pE)) = 1;
        for iCutoffPM = 1:nCutoffPM
            cutoffPM = lstCutoffsPM(iCutoffPM);
            fprintf('%d. %.2g...\n', iCutoffPM, cutoffPM);
            for iCutoffPE = 1:nCutoffPE
                cutoffPE = lstCutoffsPE(iCutoffPE);
                %
                isCandidate = pE < cutoffPE & pM < cutoffPM;
                matGencodeGeneTissue_isCandidate_pM_pE(tableGenesNasserExpressed.iGencode(isCandidate), iTissue, iCutoffPM, iCutoffPE) = true;
                isCandidate = pE < cutoffPE & pM < cutoffPM & qCombined < cutoffQ;
                matGencodeGeneTissue_isCandidate_pM_pE_FDR(tableGenesNasserExpressed.iGencode(isCandidate), iTissue, iCutoffPM, iCutoffPE) = true;
            end
        end
        matGencodeGeneTissue_isUsed(tableGenesNasserExpressed.iGencode, iTissue) = true;
        matGencodeGeneTissue_isDriver(tableGenesNasserExpressed.iGencode(tableGenesNasserExpressed.isDriver)) = true;
    end
    toc
    %%
    sResPanCancer = struct();
    sResPanCancer_FDR = struct();
    lstPrintNames = [tableTissues.tissuePrint', {'Pan-cancer Solid', 'Pan-cancer'}];
    nTypes = length(lstPrintNames);
    %%
    for iType = 1:nTypes
        if (iType <= nTissues)
            isRelTissue = iType;
        elseif (iType == nTissues + 1)
            isRelTissue = ~contains(tableTissues.tissue, 'blood');
        else
            isRelTissue = true(nTissues, 1);
        end
        isUsedGene = sum(matGencodeGeneTissue_isUsed(:,isRelTissue), 2)>0;                                                  % Used in at least one tissue
        isDriver = matGencodeGeneTissue_isDriver;
        for iCutoffPM = 1:nCutoffPM
            for iCutoffPE = 1:nCutoffPE
                isCandidate = sum(matGencodeGeneTissue_isCandidate_pM_pE(:,isRelTissue, iCutoffPM, iCutoffPE), 2)>0;        % Candidate in at least one tissue
                [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate(isUsedGene), isDriver(isUsedGene), tailDirection, false);
                sResPanCancer.pFisherCDG(iType, iCutoffPM, iCutoffPE) = p;
                sResPanCancer.pFisherCDG_text{iType, iCutoffPM, iCutoffPE} = getPValueAsTextShort(p);
                sResPanCancer.enrichmentCDG(iType, iCutoffPM, iCutoffPE) = enrichment;
                sResPanCancer.nObserved_candidateDrivers_CDGs(iType, iCutoffPM, iCutoffPE) = nObserved;
                sResPanCancer.nExpected_candidateDrivers_CDGs(iType, iCutoffPM, iCutoffPE) = nExpected;
                sResPanCancer.nCandidateDrivers(iType, iCutoffPM, iCutoffPE) = sum(isCandidate);
                %
                isCandidate = sum(matGencodeGeneTissue_isCandidate_pM_pE_FDR(:,isRelTissue, iCutoffPM, iCutoffPE), 2)>0;    % Candidate in at least one tissue
                [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate(isUsedGene), isDriver(isUsedGene), tailDirection, false);
                sResPanCancer_FDR.pFisherCDG(iType, iCutoffPM, iCutoffPE) = p;
                sResPanCancer_FDR.pFisherCDG_text{iType, iCutoffPM, iCutoffPE} = getPValueAsTextShort(p);
                sResPanCancer_FDR.enrichmentCDG(iType, iCutoffPM, iCutoffPE) = enrichment;
                sResPanCancer_FDR.nObserved_candidateDrivers_CDGs(iType, iCutoffPM, iCutoffPE) = nObserved;
                sResPanCancer_FDR.nExpected_candidateDrivers_CDGs(iType, iCutoffPM, iCutoffPE) = nExpected;
                sResPanCancer_FDR.nCandidateDrivers(iType, iCutoffPM, iCutoffPE) = sum(isCandidate);
            end
        end
    end
    %%
    toc
    %myPrintMemory
    dataCutoffs.sResPanCancer = sResPanCancer;
    dataCutoffs.sResPanCancer_FDR = sResPanCancer_FDR;
    dataCutoffs.tableTissues = tableTissues;
    dataCutoffs.lstCutoffsPM = lstCutoffsPM;
    dataCutoffs.lstCutoffsPE = lstCutoffsPE;
    createDir(fileparts(saveFileData));
    save(saveFileData, 'dataCutoffs');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'dataCutoffs');
end
