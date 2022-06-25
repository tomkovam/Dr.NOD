% scriptRunAnalysis
clear; clc; close all; addpath(genpath('code/')); rng(1);
imagesPath = 'results/mainFigures/'; createDir(imagesPath);
%% Loading data
[tableTissues, sProperties] = loadParameters;
[tableGencodeGenes, tableTissuesWithPancancer, sResults, tableMutations_candidate, matGeneGencodeIsCandidateMut, tableTissues_data1] = loadData1();
[sResCrossCADD, sResPanCancerCrossCADD, lstMinCADD_PHRED, tableTissues_data2] = loadData2_crossCADD();
[tableABC, sResCrossTissues, tableTissues_data3] = loadData3_crossTissue();
[tableMutations_candidate, tableTissues_data4, tableTissuesWithPancancer_data4] = loadData4_TFBS(tableMutations_candidate);
[dataTFBS, tableMutations_candidate] = loadMotifs(tableMutations_candidate, tableTissues_data4, sProperties);
dataDepMap = loadData5_DepMap(tableTissues_data1, sResults, tableGencodeGenes, sProperties);
dataCutoffs = loadData6_cutoffs();
[dataSupTables, tableGencodeGenes, tableMutations_candidate] = loadData7_annotatedMutationGenes(tableGencodeGenes, tableMutations_candidate, matGeneGencodeIsCandidateMut, dataDepMap, sProperties);
sColours = getColours();
%% Checking that the tissue tables match
if (~isequal(tableTissues.tissue, tableTissues_data1.tissue)), error('Tissues data1 do not match.'); end
if (~isequal(tableTissues.tissue, tableTissues_data2.tissue)), error('Tissues data2 do not match.'); end
if (~isequal(tableTissues.tissue, tableTissues_data3.tissue)), error('Tissues data3 do not match.'); end
if (~isequal(tableTissues.tissue, tableTissues_data4.tissue)), error('Tissues data4 do not match.'); end
if (~isequal(tableTissues.tissuePrint(2:end), dataDepMap.tableTissuesWithPancancer_DepMap.tissuePrint(2:end))), error('Tissues data5 do not match.'); end
if (~isequal(tableTissues.tissue, dataCutoffs.tableTissues.tissue)), error('Tissues data6 do not match.'); end
%% Supplementary tables
writetable(tableTissuesWithPancancer(:,dataSupTables.lstColTissues), [imagesPath, 'SupplementaryTable1.xlsx']);
writetable(dataSupTables.tableGencodeGenesCandidates(dataSupTables.tableGencodeGenesCandidates.isCandidateSolid,dataSupTables.lstColGenes), [imagesPath, 'SupplementaryTable2.xlsx'], 'sheet', 'solid'); 
writetable(dataSupTables.tableGencodeGenesCandidates(~dataSupTables.tableGencodeGenesCandidates.isCandidateSolid,dataSupTables.lstColGenesBlood), [imagesPath, 'SupplementaryTable2.xlsx'], 'sheet', 'blood');
writetable(dataSupTables.tableMutationGenePairs(dataSupTables.tableMutationGenePairs.iTissue > 1 & dataSupTables.tableMutationGenePairs.isHighCADD & ~dataSupTables.tableMutationGenePairs.isExcluded,dataSupTables.lstColsMutationGenePairs), [imagesPath, 'SupplementaryTable4.xlsx']);
%% Main and supplementary figures
plotFigure1(imagesPath, sColours, sResults, dataSupTables.lstGenesStrongSupport);
plotFigure2(imagesPath, sColours, tableTissuesWithPancancer, tableTissues_data3, tableABC, sResCrossTissues, sResCrossCADD, sResPanCancerCrossCADD, tableMutations_candidate, lstMinCADD_PHRED, sProperties);
plotFigure3(imagesPath, sColours, tableMutations_candidate, tableTissues_data1, sResults, sProperties);
plotFigure4(imagesPath, sColours, tableGencodeGenes, tableTissues_data1, dataDepMap, sProperties);
plotFigure5(imagesPath, sColours, tableTissuesWithPancancer_data4, tableTissues_data4, dataTFBS, tableMutations_candidate, tableTissues_data1);
plotFigure6(dataSupTables.tableMutationGenePairs, imagesPath, sColours, tableTissues_data1);
plotSupFigure1(imagesPath, dataCutoffs);
plotSupFigure2(imagesPath, sColours, dataDepMap);
%% List of used matlab code/function files
% [fList,pList] = matlab.codetools.requiredFilesAndProducts('scriptPlotFigures.m');
% writetable(cell2table(fList'), 'usedFiles2.txt');
