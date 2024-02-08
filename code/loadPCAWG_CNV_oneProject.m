function [matCNV_genesSamples, sample_has_CNV] = loadPCAWG_CNV_oneProject(runAgain, suffix, tissueName, biosampleABC, tableGenesNasserExpressed, tableGencodeGenes, tableSamples, sProperties, doSave)
%% Loads copy number variation (CNV) data

saveFileData = ['save/PCAWG/PCAWG_CNV_',suffix,'.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    %% 
    fprintf('Computing %s...\n', saveFileData);    
    %% First, we load the input file
    inputFileName = [sProperties.PCAWG_CNV_DIR, 'gene_CNV_calls_perTissue/CNV_',tissueName,'.txt'];   % PCAWG/CNV/
    if (~exist(inputFileName, 'file') && exist([inputFileName, '.gz'], 'file'))
        gunzip([inputFileName, '.gz']);
    end
    try
        tableCNV_oneTissue = readtable(inputFileName, 'ReadVariableNames', false, 'delimiter', '\t');
        tableCNV_oneTissue.Properties.VariableNames = {'geneIndex', 'geneName', 'aliquot_id', 'CNV_total', 'CNV_major', 'CNV_minor'};
        tableSamplesWithCNV = readtable([sProperties.PCAWG_CNV_DIR, 'CNV_samples_perTissue/samplesWithCNV_',tissueName,'.txt'], 'ReadVariableNames', false);  % PCAWG/CNV/
        tableSamplesWithCNV = tableSamplesWithCNV.Var1;
        fprintf('%s %s: %s CNV rows\n', tissueName, biosampleABC, num2sepNumStr(size(tableCNV_oneTissue, 1)));
    catch
        fprintf('File %s not found.\n', inputFileName);
        return
    end
    %% Next we map the genes to tableGenesNasser genes.
    tableGencodeGenes.iGeneNasserExpressed = NaN*ones(size(tableGencodeGenes, 1), 1);
    tableGencodeGenes.iGeneNasserExpressed(tableGenesNasserExpressed.iGencode) = tableGenesNasserExpressed.iGene;
    if (~isequal(tableGenesNasserExpressed.iGene, (1:size(tableGenesNasserExpressed, 1))')), error('ERROR: sample indices do not match.'); end
    if (~isequal(tableGencodeGenes.geneSymbol(tableCNV_oneTissue.geneIndex), tableCNV_oneTissue.geneName)), error('Gene names do not match.\n'); end
    tableCNV_oneTissue.indexNasser = tableGencodeGenes.iGeneNasserExpressed(tableCNV_oneTissue.geneIndex); 
    %% Next we map samples using the aliquot_id (WGS).
    sample_has_CNV = ismember(tableSamples.aliquotID_WGS, tableSamplesWithCNV);
    [~, tableCNV_oneTissue.iSample] = ismember(tableCNV_oneTissue.aliquot_id, tableSamples.aliquotID_WGS);
    fprintf('%.1f%% rows have gene and %.1f%% rows have sample, %d (%.1f%%) samples have CNV.\n', 100*mean(tableCNV_oneTissue.indexNasser>0), 100*mean(tableCNV_oneTissue.iSample>0), sum(sample_has_CNV), 100*mean(sample_has_CNV));
    %% We keep only rows with known genes and samples and check that there are no duplicate rows
    tableCNV_oneTissue = tableCNV_oneTissue(tableCNV_oneTissue.indexNasser>0 & tableCNV_oneTissue.iSample>0,:);
    if (size(tableCNV_oneTissue, 1) ~= size(unique(tableCNV_oneTissue(:,{'indexNasser', 'iSample'})), 1)), error('Some duplicate rows. This needs to be explored.\n'); end
    %% Finally, we reshape the CNV calls into matCNV_genesSamples matrix, which has CNV_total calls per gene x sample (always the call which covers at least half of the gene).
    nGenes = size(tableGenesNasserExpressed, 1);
    nSamples = size(tableSamples, 1);
    ind = sub2ind([nGenes, nSamples],tableCNV_oneTissue.indexNasser,tableCNV_oneTissue.iSample);
    tmp = NaN*ones(nGenes, nSamples); tmp(ind) = tableCNV_oneTissue.iSample;
    tmp2 = nanmax(tmp);
    tmp3 = NaN*(1:nSamples); tmp3(unique(tableCNV_oneTissue.iSample)) = unique(tableCNV_oneTissue.iSample);
    if (~isequal(tmp2(~isnan(tmp2)), tmp3(~isnan(tmp3))) || sum(tmp2<tmp3 | tmp2>tmp3)>0), error('Reshaping is suspicious/did not work.\n'); end
    matCNV_genesSamples = 2*ones(nGenes, nSamples); 
    matCNV_genesSamples(tableGenesNasserExpressed.iGencode==0,:) = NaN;  % For these genes, the CNV calls are not available
    matCNV_genesSamples(:,~sample_has_CNV) = NaN;             % For these samples, the CNV calls are not available
    matCNV_genesSamples(ind) = tableCNV_oneTissue.CNV_total;
    %% We save the results into saveFileData.
    %myPrintMemory
    if (doSave)
        createDir(fileparts(saveFileData));
        save(saveFileData, 'matCNV_genesSamples', 'sample_has_CNV');
    end
else
    load(saveFileData, 'matCNV_genesSamples', 'sample_has_CNV');
end
