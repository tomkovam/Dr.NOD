function [sResCrossCADD, sResPanCancerCrossCADD, lstMinCADD_PHRED, tableTissues] = loadData2_crossCADD()
%% Loads the cross-CADD analysis data (and runs the analysis if not precomputed).

saveFileData = 'save/data/data2_crossCADD.mat';
if (~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
    [tableTissues, sProperties] = loadParameters;
    runAgain = sProperties.runAgain; tailDirection = sProperties.tailDirection; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; mutTypeName = sProperties.mutTypeName; nGencodeGenes = sProperties.nGencodeGenes;
    %%
    lstMinCADD_PHRED = 0:2:18; 
    nBinCADD = length(lstMinCADD_PHRED);
    nTissues = size(tableTissues, 1);
    sResCrossCADD.pFisherCDG = NaN*ones(nTissues, nBinCADD);
    sResCrossCADD.pFisherCDG_text = cell(nTissues, nBinCADD); sResCrossCADD.pFisherCDG_text(:) = {''};
    sResCrossCADD.enrichmentCDG = NaN*ones(nTissues, nBinCADD);
    sResCrossCADD.nObserved_candidateDrivers_CDGs = NaN*ones(nTissues, nBinCADD);
    sResCrossCADD.nExpected_candidateDrivers_CDGs = NaN*ones(nTissues, nBinCADD);
    sResCrossCADD.nCandidateDrivers = NaN*ones(nTissues, nBinCADD);
    %%
    matGencodeGeneTissue_isUsed = false(nGencodeGenes, nTissues, nBinCADD);
    matGencodeGeneTissue_isCandidate = false(nGencodeGenes, nTissues, nBinCADD);
    matGencodeGeneTissue_isDriver = false(nGencodeGenes, 1);
    %%
    tic
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        tissueNameSV = tableTissues.tissueSV{iTissue};
        biosampleABC = tableTissues.biosampleABC{iTissue};
        for iBinCADD = 1:nBinCADD
            minCADD_PHRED = lstMinCADD_PHRED(iBinCADD);
            fprintf('\n=================== %s %d ===================\n', tissueName, minCADD_PHRED);
            sProperties.minCADD_PHRED = minCADD_PHRED;
            levelOutputArguments = 1; % We require only the minimal number of output arguments, to speed everything up.
            [tableGenesNasserExpressed, tableGenes_pValues] = ...
                computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
            %
            pM = tableGenes_pValues.(['p',xTestName,'_',mutTypeName]);
            pE = tableGenes_pValues.(['p',yTestName,'_',mutTypeName]);
            pCombined = combinePValues_EBM(pM,pE);
            qCombined = mafdr(pCombined, 'BHFDR', true);

            P_cutoff = 0.05;
            Q_cutoff = 0.15;
            isCandidate = pE < P_cutoff & pM < P_cutoff & qCombined < Q_cutoff;
            isDriver = tableGenesNasserExpressed.isDriver;
            [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isDriver, tailDirection, false);
            sResCrossCADD.pFisherCDG(iTissue,iBinCADD) = p;
            sResCrossCADD.pFisherCDG_text{iTissue,iBinCADD} = getPValueAsText(p);
            sResCrossCADD.enrichmentCDG(iTissue,iBinCADD) = enrichment;
            sResCrossCADD.nObserved_candidateDrivers_CDGs(iTissue,iBinCADD) = nObserved;
            sResCrossCADD.nExpected_candidateDrivers_CDGs(iTissue,iBinCADD) = nExpected;
            sResCrossCADD.nCandidateDrivers(iTissue,iBinCADD) = sum(isCandidate);
            %
            isONCOGENE = contains(tableGenesNasserExpressed.role_CGC, 'oncogene') & ~contains(tableGenesNasserExpressed.role_CGC, 'TSG'); % tableGenesNasserExpressed.isDriver & tableGenesNasserExpressed.isOncogene;
            isTSG = contains(tableGenesNasserExpressed.role_CGC, 'TSG') & ~contains(tableGenesNasserExpressed.role_CGC, 'oncogene'); % tableGenesNasserExpressed.isDriver & tableGenesNasserExpressed.isTSG;
            tableGenesNasserExpressed.isOncogene = isONCOGENE;
            tableGenesNasserExpressed.isTSG = isTSG;
            tableGenesNasserExpressed.isUP = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]) > 0;
            matGencodeGeneTissue_isUsed(tableGenesNasserExpressed.iGencode, iTissue, iBinCADD) = true;
            matGencodeGeneTissue_isCandidate(tableGenesNasserExpressed.iGencode(isCandidate), iTissue, iBinCADD) = true;
            matGencodeGeneTissue_isDriver(tableGenesNasserExpressed.iGencode(tableGenesNasserExpressed.isDriver)) = true;
        end
    end
    toc
    %%
    sResPanCancerCrossCADD = struct();
    for iBinCADD = 1:nBinCADD
        for iType = 1:3
            if (iType == 1)
                isRelTissue = contains(tableTissues.tissue, 'blood');
                tissueName = 'blood';
            elseif (iType == 2)
                isRelTissue = ~contains(tableTissues.tissue, 'blood');
                tissueName = 'Pan-cancer Solid';
            else
                isRelTissue = true(nTissues, 1);
                tissueName = 'Pan-caner';
            end
            isUsedGene = sum(matGencodeGeneTissue_isUsed(:,isRelTissue, iBinCADD), 2)>0;             % Used in at least one tissue
            isCandidate = sum(matGencodeGeneTissue_isCandidate(:,isRelTissue, iBinCADD), 2)>0;       % Candidate in at least one tissue
            isDriver = matGencodeGeneTissue_isDriver;            
            %
            [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate(isUsedGene), isDriver(isUsedGene), tailDirection, false);
            sResPanCancerCrossCADD.pFisherCDG(iType,iBinCADD) = p;
            sResPanCancerCrossCADD.pFisherCDG_text{iType,iBinCADD} = getPValueAsTextShort(p);
            sResPanCancerCrossCADD.enrichmentCDG(iType,iBinCADD) = enrichment;
            sResPanCancerCrossCADD.nObserved_candidateDrivers_CDGs(iType,iBinCADD) = nObserved;
            sResPanCancerCrossCADD.nExpected_candidateDrivers_CDGs(iType,iBinCADD) = nExpected;
            sResPanCancerCrossCADD.nCandidateDrivers(iType,iBinCADD) = sum(isCandidate);
            sResPanCancerCrossCADD.tissueName{iType} = tissueName;
        end
    end
    %%
    toc
    %myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'sResCrossCADD', 'sResPanCancerCrossCADD', 'lstMinCADD_PHRED', 'tableTissues');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'sResCrossCADD', 'sResPanCancerCrossCADD', 'lstMinCADD_PHRED', 'tableTissues');
end