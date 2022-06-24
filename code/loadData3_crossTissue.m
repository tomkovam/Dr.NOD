function [tableABC, sResTissues, tableTissues] = loadData3_crossTissue()

saveFileData = 'save/data/data3_crossTissue.mat';
if (~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
    [tableTissues, sProperties] = loadParameters;
    runAgain = sProperties.runAgain; tailDirection = sProperties.tailDirection; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; 
    %%
    nTissues = size(tableTissues, 1);
    tableABC = tableTissues(:,{'iTissue', 'biosampleABC'});
    nABC = size(tableABC, 1);
    %%
    sResTissues.pFisherCDG = NaN*ones(nTissues, nABC);
    sResTissues.pFisherCDG_text = cell(nTissues, nABC); sResTissues.pFisherCDG_text(:) = {''};
    sResTissues.enrichmentCDG = NaN*ones(nTissues, nABC);
    sResTissues.nObserved_candidateDrivers_CDGs = NaN*ones(nTissues, nABC);
    sResTissues.nExpected_candidateDrivers_CDGs = NaN*ones(nTissues, nABC);
    sResTissues.nCandidateDrivers = NaN*ones(nTissues, nABC);
    sResTissues.tissueName = cell(nTissues, 1);
    %%
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        tissueNameSV = tableTissues.tissueSV{iTissue};
        for iABC =  1:nABC
            biosampleABC = tableABC.biosampleABC{iABC};
            levelOutputArguments = 1; % We require only the minimal number of output arguments, to speed everything up.
            [tableGenesNasserExpressed, tableGenes_pValues] = ...
                computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
            pM = tableGenes_pValues.(['p',xTestName,'_SNVs_highCADD']);
            pE = tableGenes_pValues.(['p',yTestName,'_SNVs_highCADD']);
            pCombined = combinePValues_EBM(pM,pE);
            qCombined = mafdr(pCombined, 'BHFDR', true);

            P_cutoff = 0.05;
            Q_cutoff = 0.15;
            isCandidate = pE < P_cutoff & pM < P_cutoff & qCombined < Q_cutoff;

            isDriver = tableGenesNasserExpressed.isDriver;
            [p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isDriver, tailDirection, false);
            sResTissues.pFisherCDG(iTissue, iABC) = p;
            sResTissues.pFisherCDG_text{iTissue, iABC} = getPValueAsText(p);
            sResTissues.enrichmentCDG(iTissue, iABC) = enrichment;
            sResTissues.nObserved_candidateDrivers_CDGs(iTissue, iABC) = nObserved;
            sResTissues.nExpected_candidateDrivers_CDGs(iTissue, iABC) = nExpected;
            sResTissues.nCandidateDrivers(iTissue, iABC) = sum(isCandidate);
        end
        sResTissues.tissueName{iTissue} = tissueName;
    end
    %%
    toc
    %myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'tableABC', 'sResTissues', 'tableTissues');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableABC', 'sResTissues', 'tableTissues');
end