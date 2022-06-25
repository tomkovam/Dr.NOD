function [tableGenesNasser, tableGencodeGenes, tableChrSizes, tableCGC, tableDriverGenesPCAWG] = loadGenes(runAgain, sProperties, tableDriverMutationsPCAWG)
%% Here we load all the relevant gene lists, including various cancer driver gene (CDG) definitions, and map them together.
% Later on, we will work with tableGenesNasser, which is the list of all target genes of ABC enhancers in the Nasser2021 data. 
% We will further restrict this to only genes expressed in the given tissue and with at least one non-coding enhancer.
% The tableGencodeGenes is the list of all genes for which there is PCAWG expression calls.

if (~exist('tableDriverMutationsPCAWG', 'var'))
    [~, ~, ~, ~, ~, ~, tableDriverMutationsPCAWG] = loadDataPCAWG();
end

saveFileData = 'save/genes.mat';
if (runAgain || ~exist(saveFileData, 'file'))
    tic
    fprintf('Computing %s...\n', saveFileData);
        
    tableGenesNasser = readtable(sProperties.GENES_ABC, 'delimiter', '\t'); % [DIR_DATA_ORIG, 'genes/Nasser2021_genesTSS.formatted.txt']
    tableGenesNasser.Properties.VariableNames = {'nEnhancersAllBiosamples', 'chr', 'geneName', 'TSS'};
    tableGenesNasser = tableGenesNasser(:,{'chr', 'geneName', 'TSS', 'nEnhancersAllBiosamples'});
    tableGenesNasser.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableGenesNasser.chr);
    tableGenesNasser.chrNumeric(strcmp(tableGenesNasser.chr, 'chrX')) = 23;
    tableGenesNasser.chrNumeric(strcmp(tableGenesNasser.chr, 'chrY')) = 24;
    nGenesNasser = size(tableGenesNasser, 1);
    %%
    tableGencodeGenes = readtable(sProperties.GENES_GENCODE); % [DIR_DATA_ORIG, 'genes/sortedGenes.allGenes.hg19v19.txt']
    tableGencodeGenes.Properties.VariableNames = {'chromosome', 'pos0', 'pos1', 'geneSymbol', 'geneType1', 'strand', 'geneNameGencode', 'geneType2'};
    tableGencodeGenes.geneNameGencodeFirstPart = cellfun(@(x) x(1:min(strfind(x, '.'))-1), tableGencodeGenes.geneNameGencode, 'UniformOutput', false);
    %%
    [tableGenesNasser.isInGencode, tableGenesNasser.iGencode] = ismember(tableGenesNasser.geneName, tableGencodeGenes.geneSymbol);
    fprintf('%d (%.1f%%) Nasser genes are not in tableGencodeGenes.\n', sum(~tableGenesNasser.isInGencode), 100*mean(~tableGenesNasser.isInGencode));
    %% 1927 NASSER genes not found in GENCODE --> we could try matching them by distance to TSS, but not doing that for now.
    tableGenesNasser.geneType = cell(nGenesNasser, 1); tableGenesNasser.geneType(:) = {''};
    tableGenesNasser.geneNameGencode = cell(nGenesNasser, 1); tableGenesNasser.geneNameGencode(:) = {''};
    tableGenesNasser.strand_GENCODE = cell(nGenesNasser, 1); tableGenesNasser.strand_GENCODE(:) = {''};
    tableGenesNasser.chr_GENCODE = cell(nGenesNasser, 1); tableGenesNasser.chr_GENCODE(:) = {''};
    tableGenesNasser.pos0_GENCODE = NaN*ones(nGenesNasser, 1);
    tableGenesNasser.pos1_GENCODE = NaN*ones(nGenesNasser, 1);
    isOK = tableGenesNasser.iGencode>0;
    tableGenesNasser.geneType(isOK) = tableGencodeGenes.geneType2(tableGenesNasser.iGencode(isOK));
    tableGenesNasser.geneNameGencode(isOK) = tableGencodeGenes.geneNameGencode(tableGenesNasser.iGencode(isOK));
    tableGenesNasser.strand_GENCODE(isOK) = tableGencodeGenes.strand(tableGenesNasser.iGencode(isOK));
    tableGenesNasser.chr_GENCODE(isOK) = tableGencodeGenes.chromosome(tableGenesNasser.iGencode(isOK));
    tableGenesNasser.pos0_GENCODE(isOK) = tableGencodeGenes.pos0(tableGenesNasser.iGencode(isOK));
    tableGenesNasser.pos1_GENCODE(isOK) = tableGencodeGenes.pos1(tableGenesNasser.iGencode(isOK));
    %% Some have matching name but different chromosome --> 28 RNA and pseudogene genes --> these are wrong
    isWRONG = ~strcmp(tableGenesNasser.chr, tableGenesNasser.chr_GENCODE) & tableGenesNasser.isInGencode;
    fprintf('%d genes are on different chromosomes in NASSER vs GENCODE --> we will ignore these\n', sum(isWRONG)); % 28
    tableGenesNasser(isWRONG, :)
    tableGenesNasser.iGencode(isWRONG) = 0;
    tableGenesNasser.geneNameGencode(isWRONG) = {''};
    %% Some values of distanceTTS_gencode are very high - these are suspicious and may be wrong and we will ignore them
    isOK = strcmp(tableGenesNasser.strand_GENCODE, '+');
    tableGenesNasser.distanceTTS_gencode(isOK) = tableGenesNasser.TSS(isOK)-tableGenesNasser.pos0_GENCODE(isOK);
    tableGenesNasser.distanceTTS_gencode(~isOK) = tableGenesNasser.TSS(~isOK)-tableGenesNasser.pos1_GENCODE(~isOK);
    tableGenesNasser.geneNameGencodeIncDistantTSS = tableGenesNasser.geneNameGencode;
    isSUSPICIOUS = tableGenesNasser.distanceTTS_gencode>1e4;
    fprintf('%d genes have distance larger than 1e4in NASSER vs GENCODE --> we will ignore these\n', sum(isSUSPICIOUS)); % 59CGC
    tmp = tableGenesNasser(isSUSPICIOUS, :);
    tableGenesNasser.iGencode(isSUSPICIOUS) = 0;
    tableGenesNasser.geneNameGencode(isSUSPICIOUS) = {''};
    fprintf('TOTAL %d genes found in GENCODE. We remove %d genes that were NOT found in GENCODE (as we do not have expression data for them).\n', sum(tableGenesNasser.iGencode>0), sum(tableGenesNasser.iGencode==0));
    tableGenesNasser = tableGenesNasser(tableGenesNasser.iGencode>0,:);
    nGenesNasser = size(tableGenesNasser, 1);
    %% CGC DRIVER GENES
    tableCGC = readtable(sProperties.GENES_CGC); % [DIR_DATA_ORIG, 'genes/CGC_Census_allTue Aug 10 20_28_04 2021.txt']
    [tableCGC.isInNasser, tableCGC.iNasser] = ismember(tableCGC.GeneSymbol, tableGenesNasser.geneName);
    isOK = tableCGC.isInNasser;
    fprintf('%d (%.1f%%) CGC genes found in NASSER by symbol (%d NOT).\n', sum(isOK), 100*mean(isOK), sum(~isOK));
    % tableCGC(~isOK,:)
    tableGenesNasser.isDriver_CGC = false(nGenesNasser, 1); tableGenesNasser.isDriver_CGC(tableCGC.iNasser(tableCGC.iNasser>0)) = true;
    tableGenesNasser.tier_CGC = NaN*ones(nGenesNasser, 1); tableGenesNasser.isDriver_CGC(tableCGC.iNasser(tableCGC.iNasser>0)) = true;
    tableGenesNasser.role_CGC = cell(nGenesNasser, 1); tableGenesNasser.role_CGC(:) = {''}; tableGenesNasser.role_CGC(tableCGC.iNasser(tableCGC.iNasser>0)) = tableCGC.RoleInCancer(tableCGC.iNasser>0);
    tableGenesNasser.cancerTypes_CGC = cell(nGenesNasser, 1); tableGenesNasser.cancerTypes_CGC(:) = {''}; tableGenesNasser.cancerTypes_CGC(tableCGC.iNasser(tableCGC.iNasser>0)) = tableCGC.TumourTypes_Somatic_(tableCGC.iNasser>0);
    %% PCAWG DRIVER-MUTATIONS
    %tableDriverMutationsPCAWG = readtable([DIR_DATA_ORIG, 'PCAWG/TableS3_panorama_driver_mutations_ICGC_samples.controlled.txt']);
    % size(unique(tableDriverMutations(:, {'gene', 'driver', 'driver_statement'}))) % 958
    % size(unique(tableDriverMutations(:, {'gene', 'driver'}))) % 878
    % size(unique(tableDriverMutations(:, {'gene'}))) % 573

    tableDriverMutations_genes = unique(tableDriverMutationsPCAWG(:, {'gene'}));
    tableDriverMutations_genes = sortrows(tableDriverMutations_genes,'gene','ascend');
    [tableDriverMutations_genes.isInNasser, tableDriverMutations_genes.iNasser] = ismember(tableDriverMutations_genes.gene, tableGenesNasser.geneName);
    isOK = tableDriverMutations_genes.isInNasser;
    fprintf('%d (%.1f%%) PCAWG DRIVER-MUTATIONS genes found in NASSER by symbol (%d NOT).\n', sum(isOK), 100*mean(isOK), sum(~isOK));
    % tableDriverMutations_genes.gene(~isOK,:) % Of the protein coding, these were not mapped:
    %     {'APITD1'                                                    }
    %     {'ARID1B'                                                    }
    %     {'C17orf70'                                                  }
    %     {'ERBB2IP'                                                   }
    %     {'MRE11A'                                                    }
    %     {'MUC16'                                                     }
    %     {'PARK2'                                                     }
    %     {'PVRL4'                                                     }
    %     {'SDHA'                                                      }
    %     {'STRA13'                                                    }
    %     {'WHSC1'                                                     }
    tableGenesNasser.isDriver_PCAWG_MUTATIONS = false(nGenesNasser, 1); tableGenesNasser.isDriver_PCAWG_MUTATIONS(tableDriverMutations_genes.iNasser(tableDriverMutations_genes.iNasser>0)) = true;
    tableDriverMutations_coding = tableDriverMutationsPCAWG(contains(tableDriverMutationsPCAWG.category, 'coding'),:);
    tableGenesNasser.isDriver_PCAWG_MUTATIONS_CODING = ismember(tableGenesNasser.geneName, tableDriverMutations_coding.gene);
    %% PCAWG DRIVER GENES
    tableDriverGenesPCAWG = readtable(sProperties.GENES_DRIVERS_PCAWG_S1, 'Sheet', 'TableS1_compendium_mutational_d'); % [DIR_DATA_ORIG, 'genes/PCAWG_TableS1_compendium_mutational_drivers_formatted.xlsx']
    
    [tableDriverGenesPCAWG.isInNasser, tableDriverGenesPCAWG.iNasser] = ismember(tableDriverGenesPCAWG.Gene, tableGenesNasser.geneName);
    isOK = tableDriverGenesPCAWG.isInNasser;
    fprintf('%d (%.1f%%) PCAWG DRIVER-GENES genes found in NASSER by symbol (%d NOT).\n', sum(isOK), 100*mean(isOK), sum(~isOK));
    % tableDriverGenesPCAWG.Gene(~isOK,:) % Of the protein coding, these were not mapped:
    %     {'ARID1B'       }
    %     {'ERBB2IP'      }
    %     {'HIC2'         }
    %     {'MLLT4'        }
    %     {'MRE11A'       }
    %     {'POM121'       }
    %     {'RQCD1'        }
    %     {'SDHA'         }
    %     {'WHSC1'        }
    %     {'WHSC1L1'      }
    %     {'RMRP'         }
    %     {'RMRP'         }
    isOK = tableDriverGenesPCAWG.iNasser>0;
    lstOK = tableDriverGenesPCAWG.iNasser(isOK);
    tableGenesNasser.isDriver_PCAWG_GENES = false(nGenesNasser, 1); tableGenesNasser.isDriver_PCAWG_GENES(lstOK) = true;
    
    tableGenesNasser.driverGenePCAWG_Ensembl = cell(nGenesNasser, 1); tableGenesNasser.driverGenePCAWG_Ensembl(:) = {''};
    tableGenesNasser.driverGenePCAWG_Element_type = cell(nGenesNasser, 1); tableGenesNasser.driverGenePCAWG_Element_type(:) = {''};
    tableGenesNasser.driverGenePCAWG_Category = cell(nGenesNasser, 1); tableGenesNasser.driverGenePCAWG_Category(:) = {''};
    tableGenesNasser.driverGenePCAWG_MoF = cell(nGenesNasser, 1); tableGenesNasser.driverGenePCAWG_MoF(:) = {''};
    tableGenesNasser.driverGenePCAWG_Tissue = cell(nGenesNasser, 1); tableGenesNasser.driverGenePCAWG_Tissue(:) = {''};
    tableGenesNasser.driverGenePCAWG_Element = cell(nGenesNasser, 1); tableGenesNasser.driverGenePCAWG_Element(:) = {''};
    
    tableGenesNasser.driverGenePCAWG_Ensembl(lstOK) = tableDriverGenesPCAWG.Ensembl(isOK);
    tableGenesNasser.driverGenePCAWG_Element_type(lstOK) = tableDriverGenesPCAWG.Element_type(isOK);
    tableGenesNasser.driverGenePCAWG_Category(lstOK) = tableDriverGenesPCAWG.Category(isOK);
    tableGenesNasser.driverGenePCAWG_MoF(lstOK) = tableDriverGenesPCAWG.MoF(isOK);
    tableGenesNasser.driverGenePCAWG_Tissue(lstOK) = tableDriverGenesPCAWG.Tissue(isOK);
    tableGenesNasser.driverGenePCAWG_Element(lstOK) = tableDriverGenesPCAWG.Tissue(isOK);

    tableGenesNasser.role_PCAWG = tableGenesNasser.driverGenePCAWG_MoF;

    tableDriverGenesPCAWG_coding = tableDriverGenesPCAWG(strcmp(tableDriverGenesPCAWG.Element_type, 'cds'),:);
    tableGenesNasser.isDriver_PCAWG_GENES_CODING = ismember(tableGenesNasser.geneName, tableDriverGenesPCAWG_coding.Gene);
    %%
    tableGenesNasser.geneNameGencodeFirstPart = cellfun(@(x) x(1:min(strfind(x, '.'))-1), tableGenesNasser.geneNameGencode, 'UniformOutput', false);
    tableGenesNasser.iGene = (1:nGenesNasser)';
    %% PROGNOSTIC GENES
    tablePrognosticGenes = readtable(sProperties.GENES_SURVIVAL_PROTEIN_ATLAS); % [DIR_DATA_ORIG, 'genes/proteinAtlas_pathology.txt']
    [tablePrognosticGenes.isInTableNasser, tablePrognosticGenes.indexTableNasser] = ismember(tablePrognosticGenes.GeneName, tableGenesNasser.geneName);
    groupedPrognosticAnyTissue = grpstats(tablePrognosticGenes(:, {'indexTableNasser', 'prognostic_Favorable', 'unprognostic_Favorable', 'prognostic_Unfavorable', 'unprognostic_Unfavorable'}), 'indexTableNasser', 'min');
    groupedPrognosticAnyTissue.isPrognostic = groupedPrognosticAnyTissue.min_prognostic_Favorable<1 | groupedPrognosticAnyTissue.min_prognostic_Unfavorable<1;
    groupedPrognosticAnyTissue = groupedPrognosticAnyTissue(groupedPrognosticAnyTissue.indexTableNasser>0,:);
    
    tableGenesNasser.hasPrognosticData = false(nGenesNasser, 1);
    tableGenesNasser.hasPrognosticData(groupedPrognosticAnyTissue.indexTableNasser) = true;
    tableGenesNasser.isPrognosticAnyTissue = false(nGenesNasser, 1);
    tableGenesNasser.isPrognosticAnyTissue(groupedPrognosticAnyTissue.indexTableNasser) = groupedPrognosticAnyTissue.isPrognostic;
    fprintf('%d genes have prognostic data, %d genes are prognostic in any tissue.\n', sum(tableGenesNasser.hasPrognosticData), sum(tableGenesNasser.isPrognosticAnyTissue));
    %%
    tableGenesNasser.isDriver = tableGenesNasser.isDriver_CGC | tableGenesNasser.isDriver_PCAWG_GENES_CODING | tableGenesNasser.isDriver_PCAWG_MUTATIONS_CODING;
    tableGenesNasser.isOncogene = contains(tableGenesNasser.role_CGC, 'oncogene');
    tableGenesNasser.isTSG = contains(tableGenesNasser.role_CGC, 'TSG');
    %%
    tableChrSizes = readtable(sProperties.GENOME_SIZE, 'delimiter', '\t'); % [DIR_DATA_ORIG, 'genomeSize_hg19.main.txt']
    tableChrSizes.Properties.VariableNames = {'chr', 'nPositions'};
    tableChrSizes.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableChrSizes.chr);
    tableChrSizes.chrNumeric(strcmp(tableChrSizes.chr, 'chrX')) = 23;
    tableChrSizes.chrNumeric(strcmp(tableChrSizes.chr, 'chrY')) = 24;
    tableChrSizes = sortrows(tableChrSizes,'chrNumeric','ascend');
    %% Annotate GENCODE genes
    nGencodeGenes = size(tableGencodeGenes, 1);
    tableGencodeGenes.isDriver = false(nGencodeGenes, 1);
    tableGencodeGenes.isInCGC = false(nGencodeGenes, 1); 
    tableGencodeGenes.role_CGC = cell(nGencodeGenes, 1); tableGencodeGenes.role_CGC(:) = {''}; 
    [isOK, index] = ismember(tableCGC.GeneSymbol, tableGencodeGenes.geneSymbol); % 707 (97.8%) CGC genes found in tableGencodeGenes by symbol (16 NOT).
    tableCGC = tableCGC(isOK,:); index = index(isOK);
    tableGencodeGenes.isInCGC(index) = true;
    tableGencodeGenes.role_CGC(index) = tableCGC.RoleInCancer;
    tableGencodeGenes.isONCOGENE = contains(tableGencodeGenes.role_CGC, 'oncogene');
    tableGencodeGenes.isTSG = contains(tableGencodeGenes.role_CGC, 'TSG');
    tableGencodeGenes.isDriver = tableGencodeGenes.isInCGC;
    tableGencodeGenes.isDriver(tableGenesNasser.iGencode) = tableGenesNasser.isDriver;
    %%
    toc
    myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'tableGenesNasser', 'tableGencodeGenes', 'tableChrSizes', 'tableCGC', 'tableDriverGenesPCAWG'); %, 'tableIntOGen');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableGenesNasser', 'tableGencodeGenes', 'tableChrSizes', 'tableCGC', 'tableDriverGenesPCAWG'); %, 'tableIntOGen');
end