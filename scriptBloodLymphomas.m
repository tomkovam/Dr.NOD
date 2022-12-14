% scriptBloodLymphomas.m
clear; clc; close all; addpath(genpath('code/')); rng(1);
imagesPath = 'results/mainFiguresBlood/'; createDir(imagesPath);
%% A more detailed analysis of the DLBCL (lymphoma) subground of blood cancers
%%
[tableTissues, sProperties] = loadParameters;
runAgain = sProperties.runAgain; tailDirection = sProperties.tailDirection; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; mutTypeName = sProperties.mutTypeName; nGencodeGenes = sProperties.nGencodeGenes;
%%
nTissues = size(tableTissues, 1);
iTissue = 1;
tissueName = tableTissues.tissue{iTissue};
tissueNameSV = tableTissues.tissueSV{iTissue};
biosampleABC = tableTissues.biosampleABC{iTissue};
exclusionType = 'exclude_cll';
sProperties.exclusionType = exclusionType;
%%
 levelOutputArguments = 3;
 [tableGenesNasserExpressed, tableGenes_pValues, stats, tableSamples, tableGencodeGenes, ... % levelOutputArguments = 1
    tableGenes_pValues_hyperUE, stats_hyperUE, tableMutations, matMutationsEnhancers, matUniqueEnhancersGenes, matExpressionGenesSamples, ... % levelOutputArguments = 2
    matGenesSamplesNMut_SNVs_highCADD, matGenesSamplesNMut_INDEL, matCNV_genesSamples, tableGenes_annotations, tableGenes_mean_trinucleotides, tableUniqueEnhancers, ...
    tableDriverMutations, matUESamplesIsMut_SNVs_highCADD, matUESamplesIsMut_INDEL, tableUE_annotations, tableUE_mean_trinucleotdies, matUESamplesIsMut_SNVs_highCADD_hyperUE, tableUE_annotations_hyperUE, matGenesSamplesNMut_SNVs] = ... % levelOutputArguments = 3
    computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
 %%
 tableGencodeGenes.isGeneIG = startsWith(tableGencodeGenes.geneType2, 'IG_');
 tableGencodeGenes.chrNumeric = cellfun(@(x) str2double(x(4:end)), tableGencodeGenes.chromosome);
 tableGencodeGenes.chrNumeric(strcmp(tableGencodeGenes.chromosome, 'chrX')) = 23;
 tableGencodeGenes.chrNumeric(strcmp(tableGencodeGenes.chromosome, 'chrY')) = 24;
 maxDistance = 1e3; % 10 kb
 tableMutations.isCloseToGeneIG = false(size(tableMutations, 1), 1);
 for iGeneIG = find(tableGencodeGenes.isGeneIG)'
     isCloseMutation = (tableMutations.chrNumeric == tableGencodeGenes.chrNumeric(iGeneIG)) & ...
         ((tableMutations.pos1 - tableGencodeGenes.pos0(iGeneIG)) <= maxDistance) & ...
         ((tableGencodeGenes.pos1(iGeneIG) - tableMutations.pos0) <= maxDistance);
     tableMutations.isCloseToGeneIG(isCloseMutation) = true;
 end
 sum(tableMutations.isCloseToGeneIG)
 mean(tableMutations.isCloseToGeneIG) % 35% (3467) candidate driver mutations happen in up to 10 kb from an IG gene | 6% (9063) mutations happen in up to 10 kbp from an IG gene
 %%
 tableGenesNasserExpressed.isMutationCloseToGeneIG = (sum(matUniqueEnhancersGenes(tableMutations.iUniqueEnhancer(tableMutations.isCloseToGeneIG & tableMutations.iUniqueEnhancer>0), :), 1)>0)';
 sum(tableGenesNasserExpressed.isMutationCloseToGeneIG)
 mean(tableGenesNasserExpressed.isMutationCloseToGeneIG) % 0.22% (26) genes have mutated regulatory regions close to an IG gene - of them,  were annotated as candidates (CRIP1, MTA1, PPM1F, PRAMENP, TOP3B)
%  tableGenesNasserExpressed.geneName(isCandidate & tableGenesNasserExpressed.isMutationCloseToGeneIG)
 tableGenes_pValues.(['p',xTestName,'_',mutTypeName])(tableGenesNasserExpressed.isMutationCloseToGeneIG) = NaN; % Genes with mutated regulatory regions close to IG genes will be excluded from this analysis
%% isCandidate vs isDriver
[isCandidate, isDriver, pM, pE, tableGenesNasserExpressed, isONCOGENE_notTSG, isTSG_notONCOGENE, isONCOGENE, isTSG, sizeEffectE, sizeEffectM, pCombined, qCombined, isP_M, P_cutoff] = computeCandidateDrivers(tableGenesNasserExpressed, tableGenes_pValues, sProperties, xTestName, yTestName, mutTypeName);
sResults{iTissue}.tissuePrint = tableTissues.tissuePrint{iTissue};
sResults{iTissue}.biosampleABC = biosampleABC;
sResults{iTissue}.pM = pM;
sResults{iTissue}.pE = pE;
sResults{iTissue}.pCombined = pCombined;
sResults{iTissue}.qCombined = qCombined;
sResults{iTissue}.sizeEffectE = sizeEffectE;
sResults{iTissue}.sizeEffectM = sizeEffectM;
sResults{iTissue}.isUP = tableGenesNasserExpressed.isUP;
sResults{iTissue}.isCandidate = isCandidate;
sResults{iTissue}.isDriver = isDriver;
sResults{iTissue}.isONCOGENE = isONCOGENE;
sResults{iTissue}.isTSG = isTSG;
sResults{iTissue}.geneName = tableGenesNasserExpressed.geneName;
[p, enrichment, nObserved, nExpected] = myFisherTest(isCandidate, isDriver, tailDirection, true);
tableTissues.nSamplesWGS(iTissue) = sum(~tableSamples.isExcluded);
tableTissues.nSamplesWGSandRNA(iTissue) = sum(~tableSamples.isExcluded & tableSamples.has_RNA);
tableTissues.pFisherCDG(iTissue) = p;
tableTissues.pFisherCDG_text{iTissue} = getPValueAsText(p);
tableTissues.enrichmentCDG(iTissue) = enrichment;
tableTissues.nObserved_candidateDrivers_CDGs(iTissue) = nObserved;
tableTissues.nExpected_candidateDrivers_CDGs(iTissue) = nExpected;
%%

isUP = isCandidate & tableGenesNasserExpressed.isUP;
isDOWN = isCandidate & ~tableGenesNasserExpressed.isUP;
lstGencodeCandidates = tableGenesNasserExpressed.iGencode(isCandidate);
%%
tableGencodeGenes.nMutations = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.nMutationsHighCADD = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.nMutSamples = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.nMutSamplesHighCADD = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.nMutSamplesHighCADD_hasRNA = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.iTissue = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.tissuePrint = cell(nGencodeGenes, 1); tableGencodeGenes.tissuePrint(:) = {''};
tableGencodeGenes.pM = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.pE = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.pCombined = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.qCombined = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.sizeEffectM = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.sizeEffectE = NaN*ones(nGencodeGenes, 1);
tableGencodeGenes.isUP = false(nGencodeGenes, 1);
tableGencodeGenes.isCandidate = false(nGencodeGenes, 1);
tableGencodeGenes.isDriver = false(nGencodeGenes, 1);
tableGencodeGenes.isONCOGENE = false(nGencodeGenes, 1);
tableGencodeGenes.isTSG = false(nGencodeGenes, 1);
%%
tableGencodeGenes.nMutations(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs(isCandidate, ~tableSamples.isExcluded), 2);
tableGencodeGenes.nMutationsHighCADD(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs_highCADD(isCandidate, ~tableSamples.isExcluded), 2);
tableGencodeGenes.nMutSamples(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs(isCandidate, ~tableSamples.isExcluded)>0, 2);
tableGencodeGenes.nMutSamplesHighCADD(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs_highCADD(isCandidate, ~tableSamples.isExcluded)>0, 2);
tableGencodeGenes.nMutSamplesHighCADD_hasRNA(lstGencodeCandidates) = sum(matGenesSamplesNMut_SNVs_highCADD(isCandidate, tableSamples.has_RNA & ~tableSamples.isExcluded)>0, 2);
tableGencodeGenes.iTissue(lstGencodeCandidates) = iTissue;
tableGencodeGenes.tissuePrint(lstGencodeCandidates) = tableTissues.tissuePrint(iTissue);
tableGencodeGenes.pM(lstGencodeCandidates) = pM(isCandidate);
tableGencodeGenes.pE(lstGencodeCandidates) = pE(isCandidate);
tableGencodeGenes.pCombined(lstGencodeCandidates) = pCombined(isCandidate);
tableGencodeGenes.qCombined(lstGencodeCandidates) = qCombined(isCandidate);
tableGencodeGenes.sizeEffectM(lstGencodeCandidates) = sizeEffectM(isCandidate);
tableGencodeGenes.sizeEffectE(lstGencodeCandidates) = sizeEffectE(isCandidate);
tableGencodeGenes.isUP(lstGencodeCandidates) = isUP(isCandidate);
tableGencodeGenes.isCandidate(lstGencodeCandidates) = true;
tableGencodeGenes.isDriver(lstGencodeCandidates) = isDriver(isCandidate);
tableGencodeGenes.isONCOGENE(lstGencodeCandidates) = isONCOGENE(isCandidate);
tableGencodeGenes.isTSG(lstGencodeCandidates) = isTSG(isCandidate);
%%
[tableMutations_FunSeq2] = loadTableMutations_FunSeq2(runAgain, tissueName, biosampleABC, tableMutations, tableSamples, sProperties);
if (~isequal(tableMutations.pos1, tableMutations_FunSeq2.pos1)), error('Pos1 do not match.'); end
tableMutations = [tableMutations, tableMutations_FunSeq2(:,2:end)];
%% tableMutations_candidate
isUE_candidate = sum(matUniqueEnhancersGenes(:,isCandidate), 2)>0;
isUE_candidateUP = sum(matUniqueEnhancersGenes(:,isCandidate & isUP), 2)>0;
tableMutations.nUE = full(sum(matMutationsEnhancers, 2));
tableMutations.iMut = (1:size(tableMutations, 1))';
clear tableMutations_candidateOneTissue
tableMutations.isCandidateDriver = false(size(tableMutations, 1), 1); tableMutations.isCandidateDriver(full(sum(matMutationsEnhancers(:,isUE_candidate), 2)>0)) = true;
tableMutations.isCandidateDriverUP = false(size(tableMutations, 1), 1); tableMutations.isCandidateDriverUP(full(sum(matMutationsEnhancers(:,isUE_candidateUP), 2)>0)) = true;
tableMutations.isCandidateDriverDOWN = tableMutations.isCandidateDriver & ~tableMutations.isCandidateDriverUP;
tableMutations_candidateOneTissue = tableMutations(tableMutations.isCandidateDriver,:);
tableMutations_candidateOneTissue.iTissue = iTissue*ones(size(tableMutations_candidateOneTissue, 1), 1);
tableMutations_candidateOneTissue.iSamplePCAWG = tableSamples.iSamplePCAWG(tableMutations_candidateOneTissue.iSample);
tableMutations_candidateOneTissue.candidateGenes = cell(size(tableMutations_candidateOneTissue, 1),1); tableMutations_candidateOneTissue.candidateGenes(:) = {''};
% tableMutations_candidateOneTissue.expressionBelowExpectation = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
tableMutations_candidateOneTissue.expressionMedianWT = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
tableMutations_candidateOneTissue.expressionMedianWT_orLowCADD = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
tableMutations_candidateOneTissue.expressionThisMut = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
% tableMutations_candidateOneTissue.CNVThisMut = NaN*ones(size(tableMutations_candidateOneTissue, 1), 1);
matIsCandidateMutGene = false(size(tableMutations_candidateOneTissue, 1), size(tableGenesNasserExpressed, 1));
for jMut = 1:size(tableMutations_candidateOneTissue, 1)
    iMut = tableMutations_candidateOneTissue.iMut(jMut);
    isUE = full(matMutationsEnhancers(iMut,:)==1);
    isGene = isCandidate & (sum(matUniqueEnhancersGenes(isUE,:), 1)>0)';
    matIsCandidateMutGene(jMut,isGene) = true;
    tableMutations_candidateOneTissue.candidateGenes{jMut} = strjoin(tableGenesNasserExpressed.geneName(isGene)');
end

for jMut = 1:size(tableMutations_candidateOneTissue, 1)
    for iGene = find(matIsCandidateMutGene(jMut,:))
        tableMutations_candidateOneTissue.expressionThisMut(jMut) = matExpressionGenesSamples(iGene, tableMutations_candidateOneTissue.iSample(jMut));
    end
end

for iGene = find(isCandidate)'
    % WT
    isMut = matIsCandidateMutGene(:,iGene);
    isWT = ~tableSamples.isExcluded; isWT(tableMutations_candidateOneTissue.iSample(isMut)) = false;
    tableMutations_candidateOneTissue.expressionMedianWT(isMut) = median(matExpressionGenesSamples(iGene,isWT), 'omitnan');
    % WT or low-CADD
    isMut = matIsCandidateMutGene(:,iGene) & tableMutations_candidateOneTissue.isHighCADD;
    isWT = ~tableSamples.isExcluded; isWT(tableMutations_candidateOneTissue.iSample(isMut)) = false;
    tableMutations_candidateOneTissue.expressionMedianWT_orLowCADD(isMut) = median(matExpressionGenesSamples(iGene,isWT), 'omitnan');
end
[~,permMuts] = sortrows(tableMutations_candidateOneTissue,'candidateGenes','ascend');
matIsCandidateMutGene = matIsCandidateMutGene(permMuts,:);
tableMutations_candidateOneTissue = tableMutations_candidateOneTissue(permMuts,:);
%%
matGeneGencodeIsCandidateMut_oneTissue = false(nGencodeGenes, size(tableMutations_candidateOneTissue, 1));
matGeneGencodeIsCandidateMut_oneTissue(tableGenesNasserExpressed.iGencode, :) = matIsCandidateMutGene';
matGeneGencodeIsCandidateMut = matGeneGencodeIsCandidateMut_oneTissue;
%%
tableMutations_candidate = tableMutations_candidateOneTissue;
tableMutations_candidate_FunSeq2 = loadTableMutations_FunSeq2_context50(runAgain, tissueName, biosampleABC, tableMutations_candidate, tableSamples, sProperties);
% tableMutations_candidate = [tableMutations_candidate, tableMutations_candidate_FunSeq2(:,{'FunSeq2_motifAnalysis', 'context50bp'})];
tableMutations_candidate.FunSeq2_motif_analysis = tableMutations_candidate_FunSeq2.FunSeq2_motifAnalysis;
tableMutations_candidate.context50bp = tableMutations_candidate_FunSeq2.context50bp;
% tableMutations_candidate.chr = tableMutations_candidate_FunSeq2.chr;
tableMutations_candidate.isMOTIFG = tableMutations_candidate.FunSeq2_isMOTIFG;
tableMutations_candidate.isMOTIFBR = tableMutations_candidate.FunSeq2_isMOTIFBR;
tableMutations_candidateFullTable = tableMutations_candidate;
isOK = ~tableMutations_candidate.isIndel & ~tableMutations_candidate.isExcluded & tableMutations_candidate.iTissue==1;
tableMutations_candidate = tableMutations_candidate(isOK,:);
[dataTFBS, tableMutations_candidate] = loadMotifs(tableMutations_candidate, nTissues, sProperties);
tableMotifs = dataTFBS.tableMotifs;
tableMutations_candidate.tissuePrint = tableTissues.tissuePrint(tableMutations_candidate.iTissue);
tableMutations_candidate.chr = strrep(cellstr(num2str(tableMutations_candidate.chrNumeric,'chr%d')), ' ', '');
tableMutations_candidate.pos1_text = num2sepNumStr(tableMutations_candidate.pos1);
tableMutations_candidate.iRow = (1:size(tableMutations_candidate, 1))';
tableMutations_candidate.expressionFC = tableMutations_candidate.expressionThisMut./tableMutations_candidate.expressionMedianWT;
%% Manual edit - in one of the top hits, the first TF (motif) is not a negative regulator, and the break is not as strong. Therefore, we use the second TF (motif) in the list.
% It would be more elegant to do this automatically. However, we don't need it for the other examples. Therefore, it is done manually.
iRow = 193;
tableMutations_candidate.FunSeq2_motif_analysis{iRow} = strrep(tableMutations_candidate.FunSeq2_motif_analysis{iRow}, 'CTCF,DHS,EBF1,EGR1,ELF1,FOXA1,FOXA2,TAF1,ZEB1#FOXJ2_1#60988218#60988236#-#12#0.341463#0.658537,', '');
tableMutations_candidate.MOTIFBR{iRow} = regexp(tableMutations_candidate.FunSeq2_motif_analysis{iRow}, 'MOTIFBR=.+', 'match', 'once');
tmp2 = strsplit(tableMutations_candidate.MOTIFBR{iRow}, {'=', '#'}); % , ','
tableMutations_candidate.MOTIFBR_TFs{iRow} = tmp2{2};
tableMutations_candidate.MOTIFBR_motifName{iRow} = tmp2{3};
tmp3 = regexp(tmp2{3}, '.*_', 'match', 'once');
tableMutations_candidate.MOTIFBR_motifNamePrefix{iRow} = tmp3(1:end-1);
tableMutations_candidate.MOTIFBR_motifStart(iRow) = str2double(tmp2{4});
tableMutations_candidate.MOTIFBR_motifEnd(iRow) = str2double(tmp2{5});
tableMutations_candidate.MOTIFBR_isMinusStrand(iRow) = strcmp(tmp2{6}, '-');
tableMutations_candidate.MOTIFBR_positionMutation(iRow) = str2double(tmp2{7});
tableMutations_candidate.MOTIFBR_scoreAlt(iRow) = str2double(tmp2{8});               % alternative allele frequency in PWM
tmp3 = strsplit(tmp2{9}, ',');
tableMutations_candidate.MOTIFBR_scoreRef(iRow) = str2double(tmp3{1});               % reference allele frequency in PWM
%% Save additional information for all driver-upregulated and driver-downregulated genes
tableTrinucleotides = readtable(sProperties.TABLE_TRINUCLEOTIDES); % 'data/tableTriNucl96.txt'
for iGene = find(isCandidate | ismember(tableGenesNasserExpressed.geneName, {'BCL7A', 'BCL6', 'CXCR4', 'BCL2', 'HIST1H2BG', 'SGK1', 'EBF1', 'IRF1'}))'
    geneName = tableGenesNasserExpressed.geneName{iGene};
    saveForOneGeneVisualisation(tissueName, biosampleABC, geneName, pM(iGene), pE(iGene), qCombined(iGene), tableSamples, matCNV_genesSamples, matExpressionGenesSamples, matGenesSamplesNMut_SNVs_highCADD, ...
        tableMutations, matMutationsEnhancers, iGene, tableGencodeGenes, tableGenesNasserExpressed, matUniqueEnhancersGenes, tableUniqueEnhancers, tableUE_annotations_hyperUE, tableTrinucleotides, exclusionType);
end
%%
sColours = getColours();
%%
lstColsBR = {'iRow', 'candidateGenes', 'expressionFC', 'CADD_PHRED', 'VAF', 'qtlVAF', 'iSample', 'pos1_text', 'MOTIFBR_motifName', 'MOTIFBR_scoreAlt', 'MOTIFBR_scoreRef', 'nMotifPrefix_MOTIFBR', 'isMOTIFBR_positive', 'isMOTIFBR_negative', 'isMOTIFBR_activator', 'isMOTIFBR_repressor', 'MOTIFBR_TFs', 'FunSeq2_motif_analysis', 'isMOTIFBR_negative'};
lstColsG = {'iRow', 'candidateGenes', 'expressionFC', 'CADD_PHRED', 'VAF', 'qtlVAF', 'iSample', 'pos1_text', 'MOTIFG_motifName', 'MOTIFG_scoreAlt', 'MOTIFG_scoreRef', 'nMotifPrefix_MOTIFG', 'isMOTIFG_positive', 'isMOTIFG_negative', 'isMOTIFG_activator', 'isMOTIFG_repressor', 'FunSeq2_motif_analysis', 'isMOTIFG_positive'};
%%
tMC = tableMutations_candidate;
tMC = sortrows(tMC,'pos1','ascend');
isStrongMOTIFG = tMC.isMOTIFG_positive & tMC.isHighCADD & tMC.MOTIFG_scoreDiff>1; %  & tMC.MOTIFG_scoreAlt<.1 & tMC.MOTIFG_scoreRef>.1
isStrongMOTIFBR = tMC.isMOTIFBR_negative & tMC.MOTIFBR_scoreAlt<.1 & tMC.MOTIFBR_scoreRef>.1 & tMC.isHighCADD;
tMC_selectionG = tMC(isStrongMOTIFG & tMC.isMOTIFG_activator & ~tMC.isMOTIFG_repressor,lstColsG)
tMC_selectionG_2 = tMC(isStrongMOTIFG & tMC.isMOTIFG_activator,lstColsG)
tMC_selectionBR = tMC(isStrongMOTIFBR & tMC.isMOTIFBR_repressor & ~tMC.isMOTIFBR_activator,lstColsBR)
tMC_selectionBR_2 = tMC(isStrongMOTIFBR & tMC.isMOTIFBR_repressor,lstColsBR)
% Nothing interesting driver-DOWN-regulated
% tMC(tMC.isCandidateDriverDOWN & tMC.expressionFC<1 & (tMC.isMOTIFG | tMC.isMOTIFBR),:)
% tMC(tMC.isCandidateDriverDOWN & tMC.isMOTIFBR, lstColsBR)
% tMC(tMC.isCandidateDriverDOWN & tMC.isMOTIFG, lstColsG)
%%
tMC(tMC.pos1 == 60988031,lstColsG) % BCL2 driver-upregulating mutation/position in 2 samples... but only 2 lead to motif gain - SOX10, which is activator and not repressor

tMC(tMC.pos1 == 128748995,lstColsBR) % MYC driver-upregulating mutation in 6 samples --> breaks binding of RFX2, RFX3, RFX5 (RFX2 is an activator, others seem to be both positive and negative regulators, especially RFX3)
iGene = find(strcmp(tableGenesNasserExpressed.geneName, 'MYC'));
matExpressionGenesSamples(iGene, tMC.iSample(tMC.pos1 == 128748995))
iGene = find(strcmp(tableGenesNasserExpressed.geneName, 'CASC11'));
matExpressionGenesSamples(iGene, tMC.iSample(tMC.pos1 == 128748995))
%%
fig = createMaximisedFigure(4, [0 0 15 15]); hold on; fontSize = 12;
showMoreLabels = true; doPlotXLabel = true; doPlotYLabel = true; doPlotLegend = true;
plotTissueScatter(sColours, sResults, iTissue, fontSize, showMoreLabels, doPlotXLabel, doPlotYLabel, doPlotLegend);
set(gca, 'YTick', 0:4:16);
title('DLBCL');
mySaveAs(fig, imagesPath, 'SupFig_DLBCL_scatter', false, true);
savefig([imagesPath, 'SupFig_DLBCL_genomicView.fig']);
% %%
dataDepMap = loadData5_DepMap; % The general data, just as dummy (we will not save those results into the supplementary table)
%%
isOK = sum(matGeneGencodeIsCandidateMut, 2)>0; % This matrix is too large. We save only the relevant genes (and indicate this in the tableGencodeGenes table)
tableGencodeGenes.isIn_matGeneGencodeIsCandidateMut = isOK;
matGeneGencodeIsCandidateMut = matGeneGencodeIsCandidateMut(isOK,:)==1; % To convert it to boolean
%%
listChromosomes = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'}';
tableMutations_candidateFullTable.chr = listChromosomes(tableMutations_candidateFullTable.chrNumeric);
tableMutations_candidateFullTable.isInFunSeq2 = tableMutations_candidateFullTable.FunSeq2_isAnnotated;
[dataSupTables, tableGencodeGenes, tableMutations_candidateFullTable] = loadData7_annotatedMutationGenes(tableGencodeGenes, tableMutations_candidateFullTable, matGeneGencodeIsCandidateMut, dataDepMap, sProperties); 
writetable(dataSupTables.tableGencodeGenesCandidates(~dataSupTables.tableGencodeGenesCandidates.isCandidateSolid,dataSupTables.lstColGenesDLBCL), [imagesPath, 'SupplementaryTable2.xlsx'], 'sheet', 'DLBCL');
%%
tableMutations_candidate.pos1_text_simple = num2sepNumStr(tableMutations_candidate.pos1);
lstCols = {'chr', 'chrNumeric', 'pos1_text_simple', 'ref', 'alt', 'candidateGenes', 'isHighCADD', 'isMOTIFBR_negative', 'isMOTIFG_positive', 'expressionFC', 'expressionThisMut'};
%%
isSelection = tableMutations_candidate.isHighCADD==1 & ... %     (tableMutations_candidate.FunSeq2_isMOTIFBR+tableMutations_candidate.FunSeq2_isMOTIFG)>0 & ...
    (tableMutations_candidate.isMOTIFBR_negative+tableMutations_candidate.isMOTIFG_positive)>0;
tablePositionsAll = grpstats(tableMutations_candidate(isSelection,{'chr', 'chrNumeric', 'pos1', 'ref', 'candidateGenes', 'CADD_PHRED', 'expressionFC', 'isMOTIFBR_negative', 'isMOTIFG_positive'}), ... % , 'isMOTIFBR_repressor', 'isMOTIFG_activator', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'isNearTSS_250bp', 'expressionMedianWT', 'expressionMedianWT_orLowCADD', 'expressionThisMut', 'isHighCADD', 'VAF', 'qtlVAF', 'FunSeq2_isAnnotated', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'isMOTIFBR_positive', 'isMOTIFBR_negative', 'isMOTIFG_positive', 'isMOTIFG_negative', 'isMOTIFG_repressor'
    {'chr', 'chrNumeric', 'pos1', 'ref', 'candidateGenes'}, 'mean'); % , 'alt'
tablePositionsAll.posText = strcat(tablePositionsAll.chr, ':', num2sepNumStr(tablePositionsAll.pos1));
tablePositionsAll = sortrows(tablePositionsAll,{'GroupCount', 'mean_expressionFC'},'descend');
%%
cTopPositions = cell(2, 1);
for iDirection = 1:2
    if (iDirection == 1) % top TFBS gain of positive regulator
        tablePositions = tablePositionsAll(tablePositionsAll.mean_expressionFC > 3 & tablePositionsAll.GroupCount > 1 & tablePositionsAll.mean_isMOTIFG_positive==1,:);
    else                 % top TFBS break of negative regulator
        tablePositions = tablePositionsAll(tablePositionsAll.mean_expressionFC > 3 & tablePositionsAll.GroupCount > 2 & tablePositionsAll.mean_isMOTIFBR_negative==1,:);
    end

    nSelPositions = size(tablePositions, 1);
    %tablePositions.ref = cell(nSelPositions, 1); tablePositions.ref(:) = {''};
    tablePositions.alt = cell(nSelPositions, 1); tablePositions.alt(:) = {''};
    tablePositions.breakNegativeRegulatorTFBS = cell(nSelPositions, 1); tablePositions.breakNegativeRegulatorTFBS(:) = {''};
    tablePositions.gainPositiveRegulatorTFBS = cell(nSelPositions, 1); tablePositions.gainPositiveRegulatorTFBS(:) = {''};
    tablePositions.expressionFC_text = cell(nSelPositions, 1); tablePositions.expressionFC_text(:) = {''};
    tablePositions.iRow = NaN*ones(nSelPositions, 1);
    tableMutations_candidate.isSelPosition = false(size(tableMutations_candidate, 1), 1);
    for iSelPosition = 1:nSelPositions
        isOK = isSelection & tableMutations_candidate.chrNumeric == tablePositions.chrNumeric(iSelPosition) & tableMutations_candidate.pos1 == tablePositions.pos1(iSelPosition);
        %tablePositions.ref{iSelPosition} = strjoin(unique(tableMutations_candidate.ref(isOK)), ',');
        tablePositions.alt{iSelPosition} = strjoin(unique(tableMutations_candidate.alt(isOK)), ',');
        if (tablePositions.mean_isMOTIFBR_negative(iSelPosition)==1)
            tablePositions.breakNegativeRegulatorTFBS{iSelPosition} = strjoin(dataTFBS.tableMotifPrefix.motifPrefix(dataTFBS.tableMotifPrefix.isNegativeRegulator' & sum(dataTFBS.isMutMotifPrefix_MOTIFBR(isOK,:))>0), ',');
        end
        if (tablePositions.mean_isMOTIFG_positive(iSelPosition)==1)
            tablePositions.gainPositiveRegulatorTFBS{iSelPosition} = strjoin(dataTFBS.tableMotifPrefix.motifPrefix(dataTFBS.tableMotifPrefix.isPositiveRegulator' & sum(dataTFBS.isMutMotifPrefix_MOTIFG(isOK,:))>0)', ',');
        end
        tablePositions.context{iSelPosition} = strjoin(unique(tableTrinucleotides.patternName(tableMutations_candidate.iPattern(isOK))), ',');
        tmp = sortrows(tableMutations_candidate(isOK,:), 'expressionThisMut', 'descend');
        tablePositions.iRow(iSelPosition) = tmp.iRow(1);
        tablePositions.expressionFC_text{iSelPosition} = sprintf('%.1f \\pm %.1f', mean(tableMutations_candidate.expressionFC(isOK)), std(tableMutations_candidate.expressionFC(isOK)));
        tableMutations_candidate.isSelPosition = false(size(tableMutations_candidate, 1), 1);
    end
    tablePositions.geneName = strrep(strrep(tablePositions.candidateGenes, 'CASC11 ', ''), 'HIST1H2AM', ''); % In the double-hits, we take the one with stronger expression effect and literature evidence
    cTopPositions{iDirection} = tablePositions;
end
tableMutations_candidate.candidateGenes_toPlot = strrep(strrep(tableMutations_candidate.candidateGenes, 'CASC11 ', ''), 'HIST1H2AM', ''); % In the double-hits, we take the one with stronger expression effect and literature evidence
%%
tmp = grpstats(tableMutations_candidate(isSelection & tableMutations_candidate.isMOTIFG_positive,{'chr', 'chrNumeric', 'pos1', 'ref', 'candidateGenes', 'CADD_PHRED', 'expressionFC', 'isMOTIFBR_negative', 'isMOTIFG_positive', 'iRow'}), ... % , 'isMOTIFBR_repressor', 'isMOTIFG_activator', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'isNearTSS_250bp', 'expressionMedianWT', 'expressionMedianWT_orLowCADD', 'expressionThisMut', 'isHighCADD', 'VAF', 'qtlVAF', 'FunSeq2_isAnnotated', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'isMOTIFBR_positive', 'isMOTIFBR_negative', 'isMOTIFG_positive', 'isMOTIFG_negative', 'isMOTIFG_repressor'
    {'chr', 'chrNumeric', 'pos1', 'ref', 'candidateGenes'}, 'mean'); % , 'alt'
tmp = sortrows(tmp,{'GroupCount', 'candidateGenes', 'mean_expressionFC'},'descend');
tmp = tmp(tmp.mean_expressionFC > 3,:);
%%
tmp = tableMutations_candidate(isSelection & tableMutations_candidate.isMOTIFBR_negative,:);
tablePositions = cTopPositions{2};
tableMutations_candidate.FunSeq2_motif_analysis(tablePositions.iRow)
% NR3C1_known13#60986415#60986437#+#5#0.000000#0.916667
% FOXO3_1#60988219#60988233#+#11#0.090909#0.863636
% RFX5_known2#128748989#128749007#+#6#0.031250#0.531250
% REST_disc9#60986978#60986992#+#2#0.000000#0.594326
% BHLHE40_disc2#60986945#60986962#+#15#0.000000#0.605485
% HDAC2_disc6#60986656#60986671#+#13#0.133333
% NR3C1_known13#60986415#60986437#+#8#0.000000#1.000000
% FOXO3_1#60988219#60988233#+#6#0.000000#1.000000
% HIC1_4#131825801#131825819#+#14#0.002039
%%
lstBases = {'A', 'C', 'G', 'T'};
sColours.A = [18,151,72]/256; %[0,1,0];
sColours.C = [32,98,156]/256; %[0,0,1];
sColours.G = [250,177,48]/256; %[16,151,72]/256; %[1,1,0];
sColours.T = [216,39,59]/256; %[1,0,0];

tableMutations_candidate.isTopHitRegion = false(size(tableMutations_candidate, 1), 1);
fig = createMaximisedFigure(6, [0 0 30 25]);
nR = 6; nC = 3; iS = 1; xS = 0.8; yS = 0.5; xB = 0.05; yB = 0.05; xM = -0.02; yM = -0.02;


lst_iS = [1,2,3, 7,8,9, 13,14,15];
for iSelPos = 1:size(tablePositions, 1)
    iRow = tablePositions.iRow(iSelPos);
    isOK = tableMutations_candidate.isMOTIFBR_negative & ...
        strcmp(tableMutations_candidate.candidateGenes_toPlot, tableMutations_candidate.candidateGenes_toPlot{iRow}) & ...
        strcmp(tableMutations_candidate.MOTIFBR_motifName, tableMutations_candidate.MOTIFBR_motifName{iRow}) & ...
        tableMutations_candidate.MOTIFBR_motifStart == tableMutations_candidate.MOTIFBR_motifStart(iRow) & ...
        tableMutations_candidate.MOTIFBR_motifEnd == tableMutations_candidate.MOTIFBR_motifEnd(iRow);
    %tableMutations_candidate(isOK,lstCols)
    tableMutations_candidate.isTopHitRegion(isOK) = true;
    
    % Bottom part: motif logo
    myGeneralSubplot(nR,nC,lst_iS(iSelPos)+3,xS,yS+.45,xB,yB,xM,yM); hold on;
    geneName = tableMutations_candidate.candidateGenes_toPlot{iRow};
    plotTFBS_logos_region(tableMutations_candidate, tableMotifs, iRow, false);
    xLimVal = get(gca, 'XLim');
    xlabel(sprintf('%s:%s-%s', tableMutations_candidate.chr{iRow}, num2sepNumStr(tableMutations_candidate.MOTIFBR_motifStart(iRow)), num2sepNumStr(tableMutations_candidate.MOTIFBR_motifEnd(iRow))),...
        'Color', .5*[1,1,1], 'VerticalAlignment', 'bottom');
    drawnow;

    % Top part: mutations by position (x-axis) and expression (y-axis) and alt base (colour)
    myGeneralSubplot(nR,nC,lst_iS(iSelPos),xS,yS,xB,yB,xM,yM); hold on;
    xValues = tableMutations_candidate.pos1 - tableMutations_candidate.MOTIFBR_motifStart(iRow);
    yValues = tableMutations_candidate.expressionFC;
    for iAlt = 4:-1:1
        isOK2 = isOK & strcmp(tableMutations_candidate.alt, lstBases{iAlt});
        colour = sColours.(lstBases{iAlt});
        plot(xValues(isOK2), yValues(isOK2), 'h', 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
    end
    xlim(xLimVal);
    ylim([0, 1.05*max(yValues(isOK))]);
    set(gca, 'XTick', [], 'XColor', 'none', 'FontSize', 10); % , 'TickLength', [0 0]
    ylabel('Expr. FC');
    title(sprintf('{\\it%s}\n%s break', geneName, strrep(tableMutations_candidate.MOTIFBR_motifName{iRow}, '_', '\_')));
end
mySaveAs(fig, imagesPath, 'SupFig_DLBCL_TFBS_breakRegions', true, true);
savefig([imagesPath, 'SupFig_DLBCL_TFBS_breakRegions.fig']);
%%
nSamples = size(tableSamples, 1);
lstGenes =  {'BCL2', 'MYC', 'IRF1', 'SGK1', 'PIM1', 'EBF1'}; % tableGenesNasserExpressed.geneName(tableGenesNasserExpressed.isCandidate); %
for jGene = 1:length(lstGenes)
    geneName = lstGenes{jGene};
    iGene = find(strcmp(tableGenesNasserExpressed.geneName, geneName));
    expressionPerSample = matExpressionGenesSamples(iGene, :)';
    isMutatedSample = matGenesSamplesNMut_SNVs_highCADD(iGene, :)'>0;
    lstSamplesInTopHitRegion = unique(tableMutations_candidate.iSample(tableMutations_candidate.isTopHitRegion & strcmp(tableMutations_candidate.candidateGenes_toPlot, geneName)));
    lstSamplesInSelection = unique(tableMutations_candidate.iSample(isSelection & strcmp(tableMutations_candidate.candidateGenes_toPlot, geneName)));
    lstSamplesHighCADD = unique(tableMutations_candidate.iSample(tableMutations_candidate.isHighCADD & strcmp(tableMutations_candidate.candidateGenes_toPlot, geneName)));
    isUpregulated = expressionPerSample >= 3*median(expressionPerSample(~isMutatedSample & ~tableSamples.isExcluded), 'omitnan');
    isInTopHitRegion = false(nSamples, 1);      isInTopHitRegion(lstSamplesInTopHitRegion) = true;
    isInSelection = false(nSamples, 1);         isInSelection(lstSamplesInSelection) = true;
    isHighCADD = false(nSamples, 1);         isHighCADD(lstSamplesHighCADD) = true;
    fprintf('%s upregulated samples\n\tall: %d, high-CADD-mutated: %d, high-CADD-mutated: %d, mutated-and-TFBS: %d, topHitRegion-mutated: %d\n', geneName, ...
        sum(isUpregulated & ~tableSamples.isExcluded), ...
        sum(isMutatedSample & isUpregulated & ~tableSamples.isExcluded), ...
        sum(isHighCADD & isUpregulated & ~tableSamples.isExcluded), ...
        sum(isInSelection & isUpregulated & ~tableSamples.isExcluded), ...
        sum(isInTopHitRegion & isUpregulated & ~tableSamples.isExcluded));
end
%%
fontSize = 12; fontSizeAnnotation = 14;
cListRows = cell(2,1);

cListRows{1} = [1179, 1073, 1120, ... % GAIN-UP MYC (not that convincing gain: 1235, 1150, 1031 | not that much overexpressed: 1222, 1251)
    1146, 1055, 1054, 1018, ... MYC
    1841, 1858, ... % PIM (not that much overexpressed: ?)
    1405, 1412, 1399, ... % EBF1
    2263, ... % SGK1 (not that much overexpressed: 8582, 8794) 2224, 
    1614, ... % IRF1 1577 (not that much overexpressed: 7898)
    937, ... % BCL2
    ];

% Not used in the end:
cListRows{2} = [...
    95, 675, 870, ... % BREAK-UP BCL2-NR3C1 (also 356, 614, 938)
    937, 303, 159, ... % BREAK-UP BCL2-'EGR1_disc4' (E2F1) (also 539, 201, 595)
    936, 863, 651, ... % BREAK-UP BCL2 'PAX4_4' (also , 418)
    294, 378, 247, ... % BREAK-UP BCL2 'REST_disc9' (also , 707)
    1581, 1606, 1592, ... % BREAK-UP IRF1 'HIC1_4'
    ]; 
cListRows{3} = cTopPositions{1}.iRow'; % top TFBS gain of positive regulator
cListRows{4} = cTopPositions{2}.iRow'; % top TFBS break of negative regulator

suffixes = {'gain', 'break', 'gainTopPositions', 'breakTopPositions'};
%%
for iDirection = 4%1:4
    nC = 9; xS = 0.85; yS = 0.7; xB = 0.1; yB = 0.05; xM = -0.03; yM = 0.012; iS = 1;
    nR = ceil(length(cListRows{iDirection})/(nC/3));
    if (iDirection > 2)
        nR = nR + 2;
        iS = 2*nC + 1;
    end
    fig = createMaximisedFigure(iDirection, [0 0 35 nR*5]); jRow = 1;
    for iRow = cListRows{iDirection}
        iTissue = tableMutations_candidate.iTissue(iRow);
        geneName = tableMutations_candidate.candidateGenes{iRow};
        if (contains(geneName, ' '))
            if (contains(geneName, 'MYC'))
                geneName = 'MYC';
            elseif (contains(geneName, 'HIST1H3J'))
                geneName = 'HIST1H3J';
            else
                tmp2 = strsplit(geneName, ' ');
                geneName = tmp2{end};
            end
        end
        myGeneralSubplot(nR,nC,iS,.6+xS,yS,xB,yB,xM,yM); hold on; iS = iS + 2;
        yAltRelative1 = plotTFBS_logos(tableMutations_candidate, tableMotifs, iRow, mod(iDirection, 2)==1, geneName);
        axPos1 = get(gca, 'Position');

        myGeneralSubplot(nR,nC,iS,.5,yS,xB,yB,xM,yM); hold on; iS = iS + 1;
        isOK = tableMutations_candidate.chrNumeric == tableMutations_candidate.chrNumeric(iRow) & tableMutations_candidate.pos1 == tableMutations_candidate.pos1(iRow) & strcmp(tableMutations_candidate.alt, tableMutations_candidate.alt{iRow});
        [xAltRelative2, yAltRelative2] = plotGene_boxplot_forLogos(tableTissues.tissue{iTissue}, biosampleABC, geneName, sColours, tableMutations_candidate.iSample(isOK), exclusionType);
        [yAltRelative2, iMax] = max(yAltRelative2);
        xAltRelative2 = xAltRelative2(iMax);
        if (iDirection < 3)
            title(sprintf('s%d', tableMutations_candidate.iSample(iRow)), 'Color', .5*[1,1,1]);
        end
        axPos2 = get(gca, 'Position');

        xa = [axPos1(1) + axPos1(3), axPos2(1) + axPos2(3)*xAltRelative2 - 0.005];
        ya = [axPos1(2) + axPos1(4)*yAltRelative1, axPos1(2) + axPos1(4)*yAltRelative2];
        annotation('arrow',xa,ya, 'Color', sColours.mutated, 'LineWidth', 1.5);

        if (mod(jRow, 3) == 1)
            if (iDirection == 1)
                str = 'TFBS gain';
            elseif (iDirection == 2)
                str = [tableMutations_candidate.MOTIFBR_motifNamePrefix{iRow}, ' TFBS break'];
            elseif (iDirection == 3)
                str = 'TFBS gain';
            else
                str = 'TFBS break';
            end
            try
                dim = [.1, axPos2(2), .05, axPos2(4)];
                annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
            catch
                % In Matlab<2022a, class TextBox does not have property Rotation.
            end
        end
        jRow = jRow + 1;
    end

    if (iDirection > 2)
        tablePositions = cTopPositions{iDirection - 2};
        axes('Position', [.02, .65, .92, .3]); hold on;
        lstXVal = 1:10; lstXVal(1) = 1.5;
        xlim(lstXVal([1,end])); ylim([-.5, nSelPositions]); axis off
        lstColNames = {'', 'Position', 'Base change', 'Samples', 'Target gene', 'Direction', {'Expression', 'fold-change'}, {'CADD', 'PHRED'}, {'Negative regulator', 'TFBS break'}, {'Positive regulator', 'TFBS gain'}};
        for iSelPosition = 1:nSelPositions
            iCol = 1;
            text(lstXVal(iCol), iSelPosition, sprintf('%d', iSelPosition), 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, tablePositions.posText{iSelPosition}, 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, strcat(tablePositions.ref{iSelPosition}, '>', tablePositions.alt{iSelPosition}), 'HorizontalAlignment', 'center');    iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, sprintf('%d', tablePositions.GroupCount(iSelPosition)), 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, tablePositions.geneName{iSelPosition}, 'FontAngle', 'italic', 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, 'upregulation', 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, tablePositions.expressionFC_text{iSelPosition}, 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, sprintf('%.1f', tablePositions.mean_CADD_PHRED(iSelPosition)), 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, tablePositions.breakNegativeRegulatorTFBS{iSelPosition}, 'HorizontalAlignment', 'center'); iCol = iCol + 1;
            text(lstXVal(iCol), iSelPosition, tablePositions.gainPositiveRegulatorTFBS{iSelPosition}, 'HorizontalAlignment', 'center'); iCol = iCol + 1;
        end
        for iCol = 1:length(lstColNames)
            text(lstXVal(iCol), -.5, lstColNames{iCol}, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        set(gca, 'YDir', 'reverse');

        fontSizeLetters = 26;
        dim = [.005 .975 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
        dim = [.005 .61 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
    end

    mySaveAs(fig, imagesPath, ['SupFig_DLBCL_TFBS_', suffixes{iDirection}], true, true);
    savefig([imagesPath, ['SupFig_DLBCL_TFBS_', suffixes{iDirection}, '.fig']]);
end
%%
tMC_toPlot = tableMutations_candidate(unique([find(tableMutations_candidate.isTopHitRegion)', cListRows{1}, cListRows{2}, cListRows{3}, cListRows{4}]),:);

fig = createMaximisedFigure(5, [0 0 30 30]); 
nR = 6; nC = 4; iS = 1; xS = 0.7; yS = 0.5; xB = 0.05; yB = 0.05; xM = -0.05; yM = -0.01;

lstGenes = {'BCL2', 'MYC', 'PIM1', 'SGK1', 'EBF1', 'IRF1'}; % HIST1H2BG
for iExample = 1:length(lstGenes)
    geneName = lstGenes{iExample};

    myGeneralSubplot(nR,nC,iS,2+xS,yS,xB,yB,xM,yM); iS = iS + 3;
    plotGene_genomicView(tissueName, biosampleABC, geneName, sColours, true, false, sProperties, exclusionType);
    
    isOK = contains(tMC_toPlot.candidateGenes, geneName);
    plot(tMC_toPlot.pos1(isOK), zeros(sum(isOK), 1), 'or', 'MarkerFaceColor', (1+[1,0,0])/2);

    myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1;
    plotGene_boxplot(tissueName, biosampleABC, geneName, sColours, false, exclusionType); drawnow;
end
mySaveAs(fig, imagesPath, 'SupFig_DLBCL_TFBS_genomicView', true, true);
savefig([imagesPath, 'SupFig_DLBCL_TFBS_genomicView.fig']);
