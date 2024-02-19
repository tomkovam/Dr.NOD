function [sResCrossCADD, sResPanCancerCrossCADD, lstMinCADD_PHRED, tableTissues] = loadData2_crossCADD(sProperties, tableTissues)
%% Loads the cross-CADD analysis data (and runs the analysis if not precomputed).

saveFileData = [sProperties.DIRECTORY_SAVE, '/main/data2_crossCADD.mat'];
if (~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
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
    sResCrossCADD.stats = cell(nTissues, nBinCADD);
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
            levelOutputArguments = 3; % We require only the minimal number of output arguments, to speed everything up.
            [tableGenesNasserExpressed, tableGenes_pValues, stats] = ...
                computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
            %
            %             pM = tableGenes_pValues.(['p',xTestName,'_',mutTypeName]);
            %             pE = tableGenes_pValues.(['p',yTestName,'_',mutTypeName]);
            %             pCombined = combinePValues_EBM(pM,pE);
            %             qCombined = mafdr(pCombined, 'BHFDR', true);
            %
            %             P_cutoff = 0.05;
            %             Q_cutoff = 0.15;
            %             isCandidate = pE < P_cutoff & pM < P_cutoff & qCombined < Q_cutoff;
            %             isDriver = tableGenesNasserExpressed.isDriver;
            %
            %             %
            %             tableGenesNasserExpressed.isUP = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]) > 0;
            %             isONCOGENE_notTSG = contains(tableGenesNasserExpressed.role_CGC, 'oncogene') & ~contains(tableGenesNasserExpressed.role_CGC, 'TSG'); % tableGenesNasserExpressed.isDriver & tableGenesNasserExpressed.isOncogene;
            %             isTSG_notONCOGENE = contains(tableGenesNasserExpressed.role_CGC, 'TSG') & ~contains(tableGenesNasserExpressed.role_CGC, 'oncogene'); % tableGenesNasserExpressed.isDriver & tableGenesNasserExpressed.isTSG;
            %
            [isCandidate, isDriver, pM, pE, tableGenesNasserExpressed, isONCOGENE_notTSG, isTSG_notONCOGENE] = computeCandidateDrivers(tableGenesNasserExpressed, tableGenes_pValues, sProperties, xTestName, yTestName, mutTypeName);
            %
            tableGenesNasserExpressed.isOncogene = isONCOGENE_notTSG;
            tableGenesNasserExpressed.isTSG = isTSG_notONCOGENE;
            matGencodeGeneTissue_isUsed(tableGenesNasserExpressed.iGencode, iTissue, iBinCADD) = true;
            matGencodeGeneTissue_isCandidate(tableGenesNasserExpressed.iGencode(isCandidate), iTissue, iBinCADD) = true;
            matGencodeGeneTissue_isDriver(tableGenesNasserExpressed.iGencode(tableGenesNasserExpressed.isDriver)) = true;

            %             stats.tableGenesNasserExpressed = tableGenesNasserExpressed;
            %             stats.tableGenes_pValues = tableGenes_pValues;
            %             stats.pM = pM;
            %             stats.pE = pE;
            %             stats.isCandidate = isCandidate;
            %             stats.isDriver = isDriver;

            [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isDriver, tailDirection, false);
            sResCrossCADD.pFisherCDG(iTissue,iBinCADD) = p;
            sResCrossCADD.pFisherCDG_text{iTissue,iBinCADD} = getPValueAsText(p);
            sResCrossCADD.enrichmentCDG(iTissue,iBinCADD) = enrichment;
            sResCrossCADD.nObserved_candidateDrivers_CDGs(iTissue,iBinCADD) = nObserved;
            sResCrossCADD.nExpected_candidateDrivers_CDGs(iTissue,iBinCADD) = nExpected;
            sResCrossCADD.nCandidateDrivers(iTissue,iBinCADD) = sum(isCandidate);
            %sResCrossCADD.stats{iTissue,iBinCADD} = stats;
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