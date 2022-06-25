function dataDepMap = loadData5_DepMap(tableTissues_data1, sResults, tableGencodeGenes, sProperties)
%% Loads the DepMap analysis data (and runs the analysis if not precomputed).

rng(1);
nPermutations = 10000;

if (~exist('datasetDepMap', 'var'))
    datasetDepMap = 'Achilles_gene_dependency';
end
%%
saveFileData = ['save/data/data5_DepMap_',datasetDepMap,'.mat'];
if (~exist(saveFileData, 'file')) 
    tic
    fprintf('Computing %s...\n', saveFileData);
    %%
    tableDepMap = readtable([sProperties.DEPMAP_DIR, datasetDepMap,'.csv']);
    lstCellLines = tableDepMap.DepMap_ID;
    lstGenesDepMap = cellfun(@(x) x(1:strfind(x, '_')-1), tableDepMap.Properties.VariableNames(2:end)', 'UniformOutput', false); % SYMBOL (INDEX) --> we take just the SYMBOL
    matDepMat = table2array(tableDepMap(:,2:end))';
    clear tableDepMap
    %%
    tableDepMapCellLines = readtable(sProperties.DEPMAP_SAMPLES); % 'data/DepMap/sample_info.csv'
    [isOK, index] = ismember(lstCellLines, tableDepMapCellLines.DepMap_ID);
    if (sum(~isOK)>0), error('%d cell lines not found'); end
    tableCellLines = tableDepMapCellLines(index,:); % {'primary_disease', 'lineage'}
    [~, tableCellLines.iTissue] = ismember(tableCellLines.primary_disease, tableTissues_data1.tissueDepMap);
    clear tableDepMapCellLines lstCellLines
    %%
    matCCLE_Expression = readmatrix(sProperties.DEPMAP_CCLE_EXPRESSION); % 'data/DepMap/CCLE_expression.csv'
    matCCLE_Expression = (matCCLE_Expression(:,2:end)');
    tableCellLinesExpression = readtable(sProperties.DEPMAP_CCLE_EXPRESSION); % 'data/DepMap/CCLE_expression.csv'
    lstCCLE_CellLines = tableCellLinesExpression.Var1;
    lstCCLE_GeneSymbols = cellfun(@(x) x(1:strfind(x, '_')-1), tableCellLinesExpression.Properties.VariableNames(2:end)', 'UniformOutput', false); % SYMBOL (INDEX) --> we take just the SYMBOL
    clear tableCellLinesExpression
    %%
    isOK = ismember(lstGenesDepMap, lstCCLE_GeneSymbols);
    lstGenesDepMap = lstGenesDepMap(isOK);
    matDepMat = matDepMat(isOK,:);
    isNormalCellLine = strcmp(tableCellLines.primary_disease, 'Non-Cancerous'); 
    isOK = ismember(tableCellLines.DepMap_ID, lstCCLE_CellLines) & ~isNormalCellLine; % We include only cancer cell-lines
    tableCellLines = tableCellLines(isOK,:);
    matDepMat = matDepMat(:,isOK);
    %%
    [~, index] = ismember(lstGenesDepMap, lstCCLE_GeneSymbols);
    lstCCLE_GeneSymbols = lstCCLE_GeneSymbols(index);
    matCCLE_Expression = matCCLE_Expression(index,:);
    [~, index] = ismember(tableCellLines.DepMap_ID, lstCCLE_CellLines); 
    lstCCLE_CellLines = lstCCLE_CellLines(index);
    matCCLE_Expression = matCCLE_Expression(:,index);    
    if (~isequal(lstCCLE_CellLines, tableCellLines.DepMap_ID)), error('Cell lines do not match.'); end
    if (~isequal(lstCCLE_GeneSymbols, lstGenesDepMap)), error('Genes do not match.'); end
    %%
    matIsExpressed = matCCLE_Expression>1; % At least TPM = 1, i.e., log2(1+1)=1
    tmp = matDepMat(~matIsExpressed);
    mean(~isnan(tmp(:)))
    matDepMat(~matIsExpressed) = NaN; % Values in non-expressed cell-lines are removed.
    %%
    nTissues = size(tableTissues_data1, 1);
    nGenes = length(lstGenesDepMap);
    vDepMapGenes_average = mean(matDepMat, 2, 'omitnan');                   % mean dependency across all cell-lines per each gene
    matDepMatGenesTissues_average = NaN*ones(nGenes, nTissues);             % mean dependency score per this tissue
    vDepMapGenes_pAbove50PerGene = 100*sum(matDepMat>0.5, 2)./sum(~isnan(matDepMat), 2);   % percentage of cell-lines with dependency score above 0.5 (50%)
    matDepMatGenesTissues_pAbove50 = NaN*ones(nGenes, nTissues);            % percentage of cell-lines from this tissue with dependency score above 0.5 (50%)
    vDepMapGenes_pAbove50PerGene_tissueMatched = NaN*vDepMapGenes_pAbove50PerGene;
    %
    nTissues = size(tableTissues_data1, 1);
    lstUsedGenes = tableGencodeGenes.geneSymbol(tableGencodeGenes.isUsedGene);
    isUsedGene = ismember(lstGenesDepMap, lstUsedGenes);
    nGenesDepMap = length(lstGenesDepMap);

    nRows = nTissues;
    tableTissuesWithPancancer_DepMap = table();
    tableTissuesWithPancancer_DepMap.tissuePrint = cell(nRows, 1); tableTissuesWithPancancer_DepMap.tissuePrint = tableTissues_data1.tissuePrint;
    tableTissuesWithPancancer_DepMap.n_control = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.n_driverUP = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.nGenesEssential_control = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.nGenesEssential_driverUP = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.pGenesEssential_driverUP = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.pGenesEssential_nonDriver = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.pFisher = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.enrichment = NaN*ones(nRows, 1);
    tableTissuesWithPancancer_DepMap.crossTissue_pPermutation = NaN*ones(nTissues, 1);

    matEnrichmentIterationTissue_crossTissue = NaN*ones(nPermutations, nTissues);
    matEssentialIterationTissue_crossTissue = NaN*ones(nPermutations, nTissues);
    for iTissue = 1:nTissues
        isCellLine = tableCellLines.iTissue == iTissue; %strcmp(tableCellLines.primary_disease, tableTissues.tissueDepMap{iTissue});
        fprintf('%s: %d cell lines\n', tableTissues_data1.tissueDepMap{iTissue}, sum(isCellLine));
        matDepMatGenesTissues_average(:,iTissue) = mean(matDepMat(:,isCellLine), 2, 'omitnan');
        matDepMatGenesTissues_pAbove50(:,iTissue) = 100*sum(matDepMat(:,isCellLine)>0.5, 2)./sum(~isnan(matDepMat(:,isCellLine)), 2);
        nMatchedCellLines = sum(isCellLine);
        %
        isGeneGroup1 = isUsedGene & ~ismember(lstGenesDepMap, sResults{iTissue}.geneName(sResults{iTissue}.isCandidate));
        isGeneGroup2 = isUsedGene & ismember(lstGenesDepMap, sResults{iTissue}.geneName(sResults{iTissue}.isCandidate & sResults{iTissue}.isUP));

        groups = NaN*ones(nGenesDepMap, 1);
        groups(isGeneGroup1) = 1;
        groups(isGeneGroup2) = 2;
        sResults{iTissue}.groups = groups;
        
        vDepMapGenes_pAbove50PerGene_tissueMatched(groups==2) = matDepMatGenesTissues_pAbove50(groups==2,iTissue);
        yValues = matDepMatGenesTissues_pAbove50(:,iTissue);
        isOK = ~isnan(groups) & ~isnan(yValues);
        tableTissuesWithPancancer_DepMap.n_control(iTissue) = sum(groups==1 & isOK);
        tableTissuesWithPancancer_DepMap.n_driverUP(iTissue) = sum(groups==2 & isOK);
        tableTissuesWithPancancer_DepMap.nGenesEssential_control(iTissue) = sum(yValues(groups==1 & isOK)>0);
        tableTissuesWithPancancer_DepMap.nGenesEssential_driverUP(iTissue) = sum(yValues(groups==2 & isOK)>0);
        tableTissuesWithPancancer_DepMap.pGenesEssential_driverUP(iTissue) = 100*mean(yValues(groups==2 & isOK)>0);
        tableTissuesWithPancancer_DepMap.pGenesEssential_nonDriver(iTissue) = 100*mean(yValues(groups==1 & isOK)>0);
        [p, enrichment] = myFisherTest(yValues(isOK)>0, groups(isOK)==2, 'both', false);
        tableTissuesWithPancancer_DepMap.pFisher(iTissue) = p;
        tableTissuesWithPancancer_DepMap.enrichment(iTissue) = enrichment;

        lstCellLines = find(~isCellLine);
        for iPermutation = 1:nPermutations
            perm = randperm(length(lstCellLines)); tmp_lstCellLines = lstCellLines(perm(1:nMatchedCellLines));

            yValues = 100*sum(matDepMat(:,tmp_lstCellLines)>0.5, 2);%./sum(~isnan(matDepMat(:,tmp_lstCellLines)), 2);   % Number of cell lines where the gene is essential
            isOK = ~isnan(groups) & ~isnan(yValues);
            tmp_groups = groups;

            isGroup1 = yValues(isOK)>0;
            isGroup2 = tmp_groups(isOK)==2;
            nExpected = sum(isGroup1) * mean(isGroup2); 
            nObserved = sum(isGroup1 & isGroup2);
            enrichment = nObserved/nExpected;
            matEssentialIterationTissue_crossTissue(iPermutation, iTissue) = nObserved;
            matEnrichmentIterationTissue_crossTissue(iPermutation, iTissue) = enrichment;
        end
        p1 = mean(matEnrichmentIterationTissue_crossTissue(:,iTissue)>=tableTissuesWithPancancer_DepMap.enrichment(iTissue));
        p2 = mean(matEnrichmentIterationTissue_crossTissue(:,iTissue)<=tableTissuesWithPancancer_DepMap.enrichment(iTissue));
        p = 2*min(p1, p2);
        tableTissuesWithPancancer_DepMap.crossTissue_pPermutation(iTissue) = p;
        fprintf('%s: %.1g larger\n', tableTissuesWithPancancer_DepMap.tissuePrint{iTissue}, p1);
        fprintf('%s: %.1g smaller\n', tableTissuesWithPancancer_DepMap.tissuePrint{iTissue}, p2);
    end
    %%
    dataDepMap.lstGenesDepMap = lstGenesDepMap;
    dataDepMap.vDepMapGenes_average = vDepMapGenes_average;
    dataDepMap.vDepMapGenes_pAbove50PerGene = vDepMapGenes_pAbove50PerGene;
    dataDepMap.vDepMapGenes_pAbove50PerGene_tissueMatched = vDepMapGenes_pAbove50PerGene_tissueMatched;
    dataDepMap.tableCellLines = tableCellLines;
    dataDepMap.tableTissuesWithPancancer_DepMap = tableTissuesWithPancancer_DepMap;
    dataDepMap.sResults = sResults;
    dataDepMap.matDepMatGenesTissues_average = matDepMatGenesTissues_average;
    dataDepMap.matDepMatGenesTissues_pAbove50 = matDepMatGenesTissues_pAbove50;
    dataDepMap.matEnrichmentIterationTissue_crossTissue = matEnrichmentIterationTissue_crossTissue;
    %%
    tableGencodeGenes.isCandidateUP_solid = (tableGencodeGenes.isCandidateSolid & tableGencodeGenes.isUP);

    isGeneGroup1 = isUsedGene & ~ismember(dataDepMap.lstGenesDepMap, tableGencodeGenes.geneSymbol(tableGencodeGenes.isCandidateSolid));
    isGeneGroup2 = isUsedGene & ismember(dataDepMap.lstGenesDepMap, tableGencodeGenes.geneSymbol(tableGencodeGenes.isCandidateUP_solid));

    groups = NaN*ones(nGenesDepMap, 1);
    groups(isGeneGroup1) = 1;
    groups(isGeneGroup2) = 2;
    dataDepMap.groups = groups;
    %%
    dataDepMap.isCandidateGene = ismember(dataDepMap.lstGenesDepMap, tableGencodeGenes.geneSymbol(tableGencodeGenes.isCandidate));
    dataDepMap.matDepMatCandidatesCellLines = matDepMat(dataDepMap.isCandidateGene,:);
    %%
    toc
    %myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'dataDepMap');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'dataDepMap');
end