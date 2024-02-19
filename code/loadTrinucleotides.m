function cTrinucleotides = loadTrinucleotides(sProperties, tableTissues)

saveFileData = [sProperties.DIRECTORY_SAVE, '/main/data8_trinucleotides.mat'];
if (~exist(saveFileData, 'file'))
    tic
    %%
    fprintf('Computing %s...\n', saveFileData);
    runAgain = sProperties.runAgain; tailDirection = sProperties.tailDirection; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; mutTypeName = sProperties.mutTypeName; nGencodeGenes = sProperties.nGencodeGenes;
    
    nTissues = size(tableTissues, 1);
    cTrinucleotides = cell(nTissues, 1);
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        tissueNameSV = tableTissues.tissueSV{iTissue};
        biosampleABC = tableTissues.biosampleABC{iTissue};
        levelOutputArguments = 3;
        %%
        [tableGenesNasserExpressed, tableGenes_pValues, stats, tableSamples, tableGencodeGenes, ... % levelOutputArguments = 1
            tableGenes_pValues_hyperUE, stats_hyperUE, tableMutations, matMutationsEnhancers, matUniqueEnhancersGenes, matExpressionGenesSamples, ... % levelOutputArguments = 2
            matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_INDEL, matCNV_genesSamples, tableGenes_annotations, tableGenes_mean_trinucleotides, tableUniqueEnhancers, ...
            tableDriverMutations, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL, tableUE_annotations, tableUE_mean_trinucleotdies, matUESamplesIsMut_SNVs_highCADD_hyperUE, tableUE_annotations_hyperUE, matGenesSamplesNMut_SNVs] = ... % levelOutputArguments = 3
            computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
        [~, ~, ~, ~, tableGenesNasserExpressed] = computeCandidateDrivers(tableGenesNasserExpressed, tableGenes_pValues, sProperties, xTestName, yTestName, mutTypeName);
        cTrinucleotides{iTissue}.tableGenesNasserExpressed = tableGenesNasserExpressed;
        cTrinucleotides{iTissue}.tableGenes_pValues = tableGenes_pValues;
        cTrinucleotides{iTissue}.tableSamples = tableSamples;
        cTrinucleotides{iTissue}.tableGencodeGenes = tableGencodeGenes;
        cTrinucleotides{iTissue}.tableGenes_mean_trinucleotides = tableGenes_mean_trinucleotides; % tableGenes_mean_trinucleotdies = array2table(matGenes_sum_trinucleotides./tableGenes_annotations.nPositionsInEnhancers);
        cTrinucleotides{iTissue}.tableGenes_annotations = tableGenes_annotations;
    end
    %%
    toc
    %myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'cTrinucleotides');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'cTrinucleotides');
end
