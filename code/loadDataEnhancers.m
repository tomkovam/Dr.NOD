function [tableEnhancers, tableGenesNasserExpressed, tableUniqueEnhancers, tableUniqueEnhancers_regions, matUniqueEnhancersGenes, nUniqueEnhancers] = loadDataEnhancers(biosampleABC, enhancerAnalysis, tableGenesNasser, verbose, sProperties)
% Here we load the ABC enhancer data of the given type (enhancerAnalysis, such as non-coding regions only), create a list of unique enhancers (just
% the regions, without the TargetGene information) and their discontinuous regions (when coding parts are subracted), create enhancer-gene matrix etc.

% enhancerAnalysis = 'Slop250bpNoncoding';
% enhancerBiosample = 'bipolar_neuron_from_iPSC-ENCODE';

fprintf('Computing loadDataEnhancers %s...\n', biosampleABC);

if (~ismember(enhancerAnalysis, {'All', 'Noncoding', 'Slop250bpAll', 'Slop250bpNoncoding'}))
    error('Allowed values for enhancerAnalysis are: All|Noncoding|Slop250bpAll|Slop250bpNoncoding, while %s was used.\n', enhancerAnalysis);
end
%% Load tableEnhancers
tableEnhancers = readtable(sProperties.ABC_ENHANCER_GENE_MAPS, 'delimiter', '\t'); % ['data/Nasser2021/',enhancerAnalysis,'Predictions.together.selCols.txt']
tableEnhancers.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'name', 'class', 'activity_base', 'TargetGene', 'TargetGeneTSS', 'TargetGeneIsExpressed', 'distance', 'isSelfPromoter', 'hic_contact_pl_scaled_adj', 'ABC.Score.Numerator', 'ABC.Score', 'biosample'};
lstBiosamplesABC = unique(tableEnhancers.biosample);

if (sum(strcmp(lstBiosamplesABC, biosampleABC)==1))
    biosampleABC_full = biosampleABC;
else
    isOK = contains(lstBiosamplesABC, biosampleABC);
    if (sum(isOK)==1)
        biosampleABC_full = lstBiosamplesABC{isOK};
    else
        error('biosample %d not found', biosampleABC);
    end
end
isOK = strcmp(tableEnhancers.biosample, biosampleABC_full);
if (sum(isOK)==0)
    error('biosample %d not found', biosampleABC_full);
end
tableEnhancers = tableEnhancers(isOK,:); % We keep only our biosample
nEnhancers = size(tableEnhancers, 1);
%% Adjust some columns
tableEnhancers.TargetGeneIsExpressed = []; % This is always true, so we don't need it;
if (~isscalar(tableEnhancers.isSelfPromoter(1)))
    tableEnhancers.isSelfPromoter = strcmp(tableEnhancers.isSelfPromoter, 'True');
end
tableEnhancers.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableEnhancers.chr);
tableEnhancers.chrNumeric(strcmp(tableEnhancers.chr, 'chrX')) = 23;
tableEnhancers.chrNumeric(strcmp(tableEnhancers.chr, 'chrY')) = 24;
tableEnhancers.length = tableEnhancers.pos1-tableEnhancers.pos0;
%% Some of the tableEnhancers.TargetGene may not be in tableGenesNasser as we have removed genes that are not in GENCODE (as we would not have RNA-seq data for them anyway).
isOK = ismember(unique(tableEnhancers.TargetGene), tableGenesNasser.geneName);
fprintf('%d TargetGene (%.1f%%) were not in GENCODE and we will ignore them therefore.\n', sum(~isOK), mean(~isOK));
isOK = ismember(tableEnhancers.TargetGene, tableGenesNasser.geneName);
fprintf('%d enhancers (%.1f%%) have a gene that was not in GENCODE and we will ignore them therefore.\n', sum(~isOK), mean(~isOK));
tableEnhancers = tableEnhancers(isOK,:);
%% Filter for expressed genes with enhancers
lstGenesExpressedWithEnhancers = unique(tableEnhancers.TargetGene);
tableGenesNasser.isExpressed = ismember(tableGenesNasser.geneName, lstGenesExpressedWithEnhancers);
if (verbose)
    fprintf('%d (%.1f%%) genes are expressed in %s.\n', sum(tableGenesNasser.isExpressed), 100*mean(tableGenesNasser.isExpressed), biosampleABC);
end
tableGenesNasserExpressed = tableGenesNasser(tableGenesNasser.isExpressed, :); 
nGenes = size(tableGenesNasserExpressed, 1);
tableGenesNasserExpressed.iGene = (1:nGenes)';
[~, tableEnhancers.iGene] = ismember(tableEnhancers.TargetGene, tableGenesNasserExpressed.geneName);
tableGenesNasserExpressed.nEnhancers = histcounts(tableEnhancers.iGene, 1:nGenes+1)';
%% Annotate with geneNameGencode
tableEnhancers.geneNameGencode = cell(size(tableEnhancers, 1), 1); tableEnhancers.geneNameGencode(:) = {''};
tableEnhancers.geneNameGencode(tableEnhancers.iGene>0) = tableGenesNasserExpressed.geneNameGencode(tableEnhancers.iGene(tableEnhancers.iGene>0));
tableEnhancers.hasGencode = ~cellfun(@isempty, tableEnhancers.geneNameGencode);
lstClasses = unique(tableEnhancers.class); %nClasses = length(lstClasses);
[~, tableEnhancers.iClass] = ismember(tableEnhancers.class, lstClasses);
%% Due to removal of non-coding regions, some enhancers have been split into multiple regions (rows). We will keep these in tableUniqueEnhancers_regions.
tableUniqueEnhancers_regions = unique(tableEnhancers(:,{'chr', 'pos0', 'pos1', 'name', 'class', 'activity_base', 'chrNumeric', 'length', 'iClass'}));
nUniqueEnhancers_regions = size(tableUniqueEnhancers_regions, 1);
%% The list of truly unique enhancers (with unique ID name) will be stored in tableUniqueEnhancers. We will link all together.
tableUniqueEnhancers = unique(tableEnhancers(:,{'chr', 'name', 'class', 'activity_base', 'chrNumeric', 'iClass'}));
nUniqueEnhancers = size(tableUniqueEnhancers, 1);
if (verbose)
    fprintf('%d enhancer-gene pairs, %d unique enhancer regions, %d unique enhancers (with unique ID name).\n', nEnhancers, nUniqueEnhancers_regions, nUniqueEnhancers);
end
tableUniqueEnhancers.iUniqueEnhancers = (1:nUniqueEnhancers)';
[~, tableEnhancers.iUniqueEnhancer] = ismember(tableEnhancers.name, tableUniqueEnhancers.name);
[~, tableUniqueEnhancers_regions.iUniqueEnhancer] = ismember(tableUniqueEnhancers_regions.name, tableUniqueEnhancers.name);
%% Number of genes regulated by this iUniqueEnhancer
tmp = unique(tableEnhancers(:,{'iUniqueEnhancer', 'TargetGene'}));
tableUniqueEnhancers.nGenes = histcounts(tmp.iUniqueEnhancer, 1:nUniqueEnhancers+1)'; clear tmp
%% Sum of lengths of tableUniqueEnhancers_regions
tmp = grpstats(tableUniqueEnhancers_regions(:, {'iUniqueEnhancer', 'length'}), 'iUniqueEnhancer', 'sum');
tableUniqueEnhancers.nPositions = NaN*ones(nUniqueEnhancers, 1);
tableUniqueEnhancers.nRegions = NaN*ones(nUniqueEnhancers, 1);
tableUniqueEnhancers.nPositions(tmp.iUniqueEnhancer) = tmp.sum_length; 
tableUniqueEnhancers.nRegions(tmp.iUniqueEnhancer) = tmp.GroupCount;clear tmp
%%
matUniqueEnhancersGenes = false(nUniqueEnhancers, nGenes);
for iRow = 1:size(tableEnhancers, 1)
    matUniqueEnhancersGenes(tableEnhancers.iUniqueEnhancer(iRow), tableEnhancers.iGene(iRow)) = true;
end


