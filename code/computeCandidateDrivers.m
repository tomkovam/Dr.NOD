function [isCandidate, isDriver, pM, pE, tableGenesNasserExpressed, isONCOGENE_notTSG, isTSG_notONCOGENE, isONCOGENE, isTSG, sizeEffectE, sizeEffectM, pCombined, qCombined, isP_M, P_cutoff, isP_E, Q_cutoff] = ...
    computeCandidateDrivers(tableGenesNasserExpressed, tableGenes_pValues, sProperties, xTestName, yTestName, mutTypeName)

isDriver = tableGenesNasserExpressed.isDriver;
pM = tableGenes_pValues.(['p',xTestName,'_',mutTypeName]); 
pE = tableGenes_pValues.(['p',yTestName,'_',mutTypeName]);
sizeEffectM = tableGenes_pValues.(['e',xTestName,'_',mutTypeName]); % fold-change observed/expected (stats.observed_mf./stats.expected_mf) | values > 1 mean enrichment of mutations, values < 1 mean depletion of mutations
sizeEffectE = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]); % beta value of MUT predictor in the Poisson GLM of expression | values > 0 mean upregulation in mutated samples, values < 0 mean downregulation in mutated samples
pCombined = combinePValues_EBM(pM,pE); %combinePValues(pM,pE,'Brown');  % Fisher
qCombined = mafdr(pCombined, 'BHFDR', true); % ALTERNATIVELY: [~, ~, ~, qCombined] = fdr_bh(pCombined); % OLD:

pM(sizeEffectM<1) = 1; % Here, we are only interested in positive selection (enrichment of mutations), not negative selection (depletion of mutations)


P_cutoff = sProperties.P_cutoff; %0.05;
Q_cutoff = sProperties.Q_cutoff; % 0.15;

isP_M = pM < P_cutoff;
isP_E = pE < P_cutoff;
isCandidate = isP_M & isP_E & qCombined < Q_cutoff;


isONCOGENE = contains(tableGenesNasserExpressed.role_CGC, 'oncogene');
isTSG = contains(tableGenesNasserExpressed.role_CGC, 'TSG');
isONCOGENE_notTSG = isONCOGENE & ~isTSG;
isTSG_notONCOGENE = isTSG & ~isONCOGENE;


tableGenesNasserExpressed.isUP = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]) > 0;
tableGenesNasserExpressed.pM = pM;
tableGenesNasserExpressed.pE = pE;
tableGenesNasserExpressed.pCombined = pCombined;
tableGenesNasserExpressed.qCombined = qCombined;
tableGenesNasserExpressed.isCandidate = isCandidate;
