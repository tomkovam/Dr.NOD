function [tableTissues, sProperties] = loadParameters()

sProperties = readPropertiesFile('inputParameters.properties');
sProperties.datasetDepMap = 'Achilles_gene_dependency';
sProperties.minCADD_PHRED = 10;
sProperties.mutTypeName = 'SNVs_highCADD';
sProperties.enhancerAnalysis = 'Slop250bpNoncoding';
sProperties.exclusionType = 'excludePOLE_MSI';
sProperties.name_scoreM='M_fullModel';  % Alternative (unused) simpler options: {'M_basic', 'M_CADDnormalized', 'M_fullModel'};
sProperties.name_scoreE='E_poisson';    % Alternative (unused) simpler options: {'E_normalLog2', 'E_poisson'};
sProperties.nGencodeGenes = 57820;
sProperties.runAgain = false; 
sProperties.doSave = true; 
sProperties.tailDirection = 'both'; 

tableTissues = readtable(sProperties.TABLE_TISSUES); % 'data/tableTissues.xlsx'
