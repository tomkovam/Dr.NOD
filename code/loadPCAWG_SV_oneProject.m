function [matSV_genesSamples_minDistance, matSV_genesSamples_nameSVs, sample_has_SV] = loadPCAWG_SV_oneProject(runAgain, suffix, tissueName, biosampleABC, tableGenesNasserExpressed, tableGencodeGenes, tableSamples, sProperties, doSave)
%% Loads structural (SV) data

saveFileData = [sProperties.DIRECTORY_SAVE, '/PCAWG/PCAWG_SV_',suffix,'.mat'];
if (runAgain || ~exist(saveFileData, 'file'))
    %% 
    fprintf('Computing %s...\n', saveFileData);    
    %% First, we load the input file
    inputFileName = [sProperties.PCAWG_SV_DIR, 'gene_SV_calls_perTissue/SV_',tissueName,'.txt'];   % PCAWG/SV/
    try
        tableSV_oneTissue = readtable(inputFileName, 'ReadVariableNames', false, 'delimiter', '\t');
        tableSV_oneTissue.Properties.VariableNames = {'geneIndex', 'geneName', 'SV_name', 'aliquot_id', 'distance'};
        tableSamplesWithSV = readtable([sProperties.PCAWG_SV_DIR, 'SV_samples_perTissue/samplesWithSV_',tissueName,'.txt'], 'ReadVariableNames', false); % PCAWG/SV/
        tableSamplesWithSV = tableSamplesWithSV.Var1;
        fprintf('%s %s: %s SV rows, %d samples\n', tissueName, biosampleABC, num2sepNumStr(size(tableSV_oneTissue, 1)));
    catch
        fprintf('File %s not found.\n', inputFileName);
        return
    end
    %% Next we map the genes to tableGenesNasser genes.
    tableGencodeGenes.iGeneNasserExpressed = NaN*ones(size(tableGencodeGenes, 1), 1);
    tableGencodeGenes.iGeneNasserExpressed(tableGenesNasserExpressed.iGencode) = tableGenesNasserExpressed.iGene;
    if (~isequal(tableGenesNasserExpressed.iGene, (1:size(tableGenesNasserExpressed, 1))')), error('ERROR: sample indices do not match.'); end
    if (~isequal(tableGencodeGenes.geneSymbol(tableSV_oneTissue.geneIndex), tableSV_oneTissue.geneName)), error('Gene names do not match.\n'); end
    tableSV_oneTissue.indexNasser = tableGencodeGenes.iGeneNasserExpressed(tableSV_oneTissue.geneIndex); 
    %% Next we map samples using the aliquot_id (WGS).
    tableSamples.has_SV = ismember(tableSamples.aliquotID_WGS, tableSamplesWithSV);
    [~, tableSV_oneTissue.iSample] = ismember(tableSV_oneTissue.aliquot_id, tableSamples.aliquotID_WGS);
    fprintf('%.1f%% rows have gene and %.1f%% rows have sample, %d (%.1f%%) samples have SV.\n', 100*mean(tableSV_oneTissue.indexNasser>0), 100*mean(tableSV_oneTissue.iSample>0), sum(tableSamples.has_SV), 100*mean(tableSamples.has_SV));
    %% We keep only rows with known genes and samples and check that there are no duplicate rows
    tableSV_oneTissue = tableSV_oneTissue(tableSV_oneTissue.indexNasser>0 & tableSV_oneTissue.iSample>0,:);
    nGenes = size(tableGenesNasserExpressed, 1);
    nSamples = size(tableSamples, 1);
    %% Compute minimal distance per gene and sample
    grp_tableSV_oneTissue = grpstats(tableSV_oneTissue(:,{'indexNasser', 'iSample', 'distance'}), {'indexNasser', 'iSample'}, 'min');
    grp_tableSV_oneTissue.SV_name = cell(size(grp_tableSV_oneTissue, 1), 1); grp_tableSV_oneTissue.SV_name(:) = {''};
    for iRow = find(~isnan(grp_tableSV_oneTissue.min_distance))'
        indexNasser = grp_tableSV_oneTissue.indexNasser(iRow);
        iSample = grp_tableSV_oneTissue.iSample(iRow);
        isOK = tableSV_oneTissue.indexNasser == indexNasser & tableSV_oneTissue.iSample == iSample;
        grp_tableSV_oneTissue.SV_name{iRow} = strjoin(tableSV_oneTissue.SV_name(isOK), '|');        
    end
    %% Finally, we reshape the SV calls into matSV_genesSamples matrix, which has SV_total calls per gene x sample (always the call which covers at least half of the gene).    
    ind = sub2ind([nGenes, nSamples],grp_tableSV_oneTissue.indexNasser,grp_tableSV_oneTissue.iSample);
    tmp = NaN*ones(nGenes, nSamples); tmp(ind) = grp_tableSV_oneTissue.iSample;
    tmp2 = max(tmp, [], 'omitnan');
    tmp3 = NaN*(1:nSamples); tmp3(unique(grp_tableSV_oneTissue.iSample)) = unique(grp_tableSV_oneTissue.iSample);
    if (~isequal(tmp2(~isnan(tmp2)), tmp3(~isnan(tmp3))) || sum(tmp2<tmp3 | tmp2>tmp3)>0), error('Reshaping is suspicious/did not work.\n'); end
    matSV_genesSamples_minDistance = NaN*ones(nGenes, nSamples); 
    matSV_genesSamples_minDistance(tableGenesNasserExpressed.iGencode==0,:) = NaN;  % For these genes, the SV calls are not available
    matSV_genesSamples_minDistance(:,~tableSamples.has_SV) = NaN;             % For these samples, the SV calls are not available
    matSV_genesSamples_minDistance(ind) = grp_tableSV_oneTissue.min_distance;
    
    matSV_genesSamples_nameSVs = cell(nGenes, nSamples); matSV_genesSamples_nameSVs(:) = {''};
    matSV_genesSamples_nameSVs(ind) = grp_tableSV_oneTissue.SV_name;
    sample_has_SV = tableSamples.has_SV;
    %% We save the results into saveFileData.
    %myPrintMemory
    if (doSave)
        createDir(fileparts(saveFileData));
        save(saveFileData, 'matSV_genesSamples_minDistance', 'matSV_genesSamples_nameSVs', 'sample_has_SV');
    end
else
    load(saveFileData, 'matSV_genesSamples_minDistance', 'matSV_genesSamples_nameSVs', 'sample_has_SV');
end
