function [tableMutations_candidate, tableTissues, tableTissuesWithPancancer] = loadData4_TFBS_2(tableMutations_candidate)
%% Loads the TFBS analysis data (and runs the analysis if not precomputed).

saveFileData = 'save/main/data4_TFBS_2.mat';
if (~exist(saveFileData, 'file'))
    tic
    %% 
    fprintf('Computing %s...\n', saveFileData);
    [tableTissues, sProperties] = loadParameters;
    runAgain = sProperties.runAgain; xTestName = sProperties.name_scoreM; yTestName = sProperties.name_scoreE; mutTypeName = sProperties.mutTypeName; nGencodeGenes = sProperties.nGencodeGenes;
    %%
    nTissues = size(tableTissues, 1);
    tableAllMutations_TFBS = table();
    %%
    tic
    for iTissue = 1:nTissues
        tissueName = tableTissues.tissue{iTissue};
        tissueNameSV = tableTissues.tissueSV{iTissue};
        %iABC =  tableTissues.iABC(iTissue);
        %biosampleABC = tableABC.biosampleABC{iABC};
        biosampleABC = tableTissues.biosampleABC{iTissue};
        levelOutputArguments = 2; % We require only the minimal number of output arguments, to speed everything up.
        %%
        [~, tableGenes_pValues, ~, tableSamples, ~, ~, ~, tableMutations, matMutationsEnhancers, matUniqueEnhancersGenes] = ...
            computeMainAnalysis(runAgain, levelOutputArguments, tissueName, biosampleABC, sProperties, tissueNameSV);
        tableMutations.iSamplePCAWG = tableSamples.iSamplePCAWG(tableMutations.iSample);
        tableMutations.iTissue = iTissue*ones(size(tableMutations, 1), 1);
        %% isCandidate vs isDriver        
        pM = tableGenes_pValues.(['p',xTestName,'_SNVs_highCADD']);
        pE = tableGenes_pValues.(['p',yTestName,'_SNVs_highCADD']);      

        P_cutoff = 0.05;
        Q_cutoff = 0.15;
        isP_M = pM < P_cutoff;
        isP_E = pE < P_cutoff;
        pCombined = combinePValues_EBM(pM,pE); 
        qCombined = mafdr(pCombined, 'BHFDR', true);        
        isCandidate = isP_M & isP_E & qCombined < Q_cutoff;
        isUP = tableGenes_pValues.(['e',yTestName,'_',mutTypeName]) > 0;
        %% Annotate with FunSeq2 predictions of TFBS change
        tableMutations_FunSeq2 = loadTableMutations_FunSeq2(false, tissueName, biosampleABC, tableMutations, tableSamples, sProperties); 
        %% Annotate with candidate driver mutations
        isUE_candidate = sum(matUniqueEnhancersGenes(:,isCandidate), 2)>0;
        isUE_candidateUP = sum(matUniqueEnhancersGenes(:,isCandidate & isUP), 2)>0;
        tableMutations.isCandidateDriver = false(size(tableMutations, 1), 1); tableMutations.isCandidateDriver(full(sum(matMutationsEnhancers(:,isUE_candidate), 2)>0)) = true;
        tableMutations.isCandidateDriverUP = false(size(tableMutations, 1), 1); tableMutations.isCandidateDriverUP(full(sum(matMutationsEnhancers(:,isUE_candidateUP), 2)>0)) = true;
        %% Pool mutations across tissues
        tableAllMutations_TFBS = [tableAllMutations_TFBS; ...
            [tableMutations(:,{'iTissue', 'isIndel', 'isExcluded', 'iSamplePCAWG', 'isCandidateDriver', 'isCandidateDriverUP', 'CADD_PHRED'}), ...
            tableMutations_FunSeq2(:,{'FunSeq2_isAnnotated', 'FunSeq2_isMOTIFBR', 'FunSeq2_isMOTIFG', 'FunSeq2_gerp', 'FunSeq2_noncodingScore', 'isInABCEnhancer', 'isNearTSS_250bp'})]];        
    end
    %% The table below was obtained through the offline FunSeq2 tool
    tableFunSeq2_driver = readtable(sProperties.FUNSEQ2_ANNOTATED_CANDIDATE_MUTATIONS, 'delimiter', {'\t', ';'}); % 'data/FunSeq2/annotatedPCAWG/driverVariants_allNonBloodIncluded.intersectedFunSeq2.context50bp.txt'
    tableFunSeq2_driver.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'ref', 'alt', 'FunSeq2_chr', 'FunSeq2_pos0', 'FunSeq2_pos1', 'FunSeq2_ref', 'FunSeq2_alt', 'FunSeq2_sample', ...
        'FunSeq2_gerp', 'FunSeq2_cds', 'FunSeq2_variant.annotation_cds', 'FunSeq2_network_hub', 'FunSeq2_gene_under_negative_selection', 'FunSeq2_ENCODE_annotated', 'FunSeq2_hot_region', ...
        'FunSeq2_motif_analysis', 'FunSeq2_sensitive', 'FunSeq2_ultra_sensitive', 'FunSeq2_ultra_conserved', 'FunSeq2_target_gene', 'FunSeq2_repeat', 'FunSeq2_coding_score', ...
        'FunSeq2_noncoding_score', 'FunSeq2_recurrence_within_samples', 'FunSeq2_recurrence_database', 'context50bp'};
    tableFunSeq2_driver.isMOTIFBR = contains(tableFunSeq2_driver.FunSeq2_motif_analysis, 'MOTIFBR');
    tableFunSeq2_driver.isMOTIFG = contains(tableFunSeq2_driver.FunSeq2_motif_analysis, 'MOTIFG');
    % Next, we keep only rows where the alt column matches the FunSeq2_alt (such as G and A|C|G)
    isOK = arrayfun(@(iRow) strcmp(tableFunSeq2_driver.ref{iRow}, tableFunSeq2_driver.FunSeq2_ref{iRow}) & contains(tableFunSeq2_driver.FunSeq2_alt{iRow}, tableFunSeq2_driver.alt{iRow}), (1:size(tableFunSeq2_driver, 1))', 'UniformOutput', true);
    tableFunSeq2_driver = tableFunSeq2_driver(isOK,:);
    %
    tableFunSeq2_driver.mutID = strcat(tableFunSeq2_driver.chr, '_', strrep(cellstr(num2str(tableFunSeq2_driver.pos1)), ' ', ''), '_', tableFunSeq2_driver.ref, '_', tableFunSeq2_driver.alt);
    
    tableMutations_candidate.chr = strrep(cellstr(num2str(tableMutations_candidate.chrNumeric,'chr%d')), ' ', '');
    tableMutations_candidate.chr(tableMutations_candidate.chrNumeric == 23) = {'chrX'};
    tableMutations_candidate.chr(tableMutations_candidate.chrNumeric == 24) = {'chrY'};
    tableMutations_candidate.mutID = strcat(tableMutations_candidate.chr, '_', strrep(cellstr(num2str(tableMutations_candidate.pos1)), ' ', ''), '_', tableMutations_candidate.ref, '_', tableMutations_candidate.alt);
    [tableMutations_candidate.isInFunSeq2, tableMutations_candidate.indexFunSeq2] = ismember(tableMutations_candidate.mutID, tableFunSeq2_driver.mutID);
    fprintf('%.1f%% driver variants have FunSeq2 annotation.\n', 100*mean(tableMutations_candidate.isInFunSeq2(tableMutations_candidate.iTissue>1 & ~tableMutations_candidate.isExcluded & ~tableMutations_candidate.isIndel)));
    tableMutations_candidate.isMOTIFBR = NaN*ones(size(tableMutations_candidate, 1), 1);
    tableMutations_candidate.isMOTIFG = NaN*ones(size(tableMutations_candidate, 1), 1);
    tableMutations_candidate.FunSeq2_motif_analysis = cell(size(tableMutations_candidate, 1), 1); tableMutations_candidate.FunSeq2_motif_analysis(:) = {''};
    tableMutations_candidate.FunSeq2_ENCODE_annotated = cell(size(tableMutations_candidate, 1), 1); tableMutations_candidate.FunSeq2_ENCODE_annotated(:) = {''};
    tableMutations_candidate.FunSeq2_target_gene = cell(size(tableMutations_candidate, 1), 1); tableMutations_candidate.FunSeq2_target_gene(:) = {''};
    tableMutations_candidate.context50bp = cell(size(tableMutations_candidate, 1), 1); tableMutations_candidate.context50bp(:) = {''};
    isOK = tableMutations_candidate.isInFunSeq2;
    tableMutations_candidate.isMOTIFBR(isOK) = 0+tableFunSeq2_driver.isMOTIFBR(tableMutations_candidate.indexFunSeq2(isOK));
    tableMutations_candidate.isMOTIFG(isOK) = 0+tableFunSeq2_driver.isMOTIFG(tableMutations_candidate.indexFunSeq2(isOK));
    tableMutations_candidate.isFunSeq2_motifChanged = tableMutations_candidate.isMOTIFBR + tableMutations_candidate.isMOTIFG;
    tableMutations_candidate.FunSeq2_motif_analysis(isOK) = tableFunSeq2_driver.FunSeq2_motif_analysis(tableMutations_candidate.indexFunSeq2(isOK));
    tableMutations_candidate.FunSeq2_ENCODE_annotated(isOK) = tableFunSeq2_driver.FunSeq2_ENCODE_annotated(tableMutations_candidate.indexFunSeq2(isOK));
    tableMutations_candidate.FunSeq2_target_gene(isOK) = tableFunSeq2_driver.FunSeq2_target_gene(tableMutations_candidate.indexFunSeq2(isOK));
    tableMutations_candidate.context50bp(isOK) = tableFunSeq2_driver.context50bp(tableMutations_candidate.indexFunSeq2(isOK));
    %% Checking that the middle 17bp of the context50bp do match the context field:
    tmp1 = tableMutations_candidate(tableMutations_candidate.isInFunSeq2,:);
    tmp1.context8bp = cellfun(@(x) upper(x(43:59)), tmp1.context50bp, 'UniformOutput', false);
    if (~isequal(tmp1.context, tmp1.context8bp)), error('Contexts do not match'); end
    %%
    save('workspace_TFBS.mat', '-v7.3');
    %%
    tableTissuesWithPancancer = table();
    nTissues = size(tableTissues, 1);
    tableTissuesWithPancancer.tissuePrint = cell(nTissues, 1);
    tableTissuesWithPancancer.nControlMutations = NaN*ones(nTissues, 1);
    tableTissuesWithPancancer.nCandidateDriverMutations = NaN*ones(nTissues, 1);
    lstTypes = {'MOTIF_any', 'MOTIFBR', 'MOTIFG'};
    for iType = 1:3
        motifType = lstTypes{iType};
        tableTissuesWithPancancer.(['nControlMutations_motifChange_',motifType]) = NaN*ones(nTissues, 1);
        tableTissuesWithPancancer.(['nCandidateDriverMutations_motifChange_',motifType]) = NaN*ones(nTissues, 1);
        tableTissuesWithPancancer.(['pControlMutations_motifChange_',motifType]) = NaN*ones(nTissues, 1);
        tableTissuesWithPancancer.(['pCandidateDriverMutations_motifChange_',motifType]) = NaN*ones(nTissues, 1);
        tableTissuesWithPancancer.(['enrichment_motifChange_',motifType]) = NaN*ones(nTissues, 1);
        tableTissuesWithPancancer.(['pValue_motifChange_',motifType]) = NaN*ones(nTissues, 1);
    end
    
    isGeneralOK = ~tableAllMutations_TFBS.isExcluded & tableAllMutations_TFBS.FunSeq2_isAnnotated & tableAllMutations_TFBS.isInABCEnhancer; % only SNVs are FunSeq2_isAnnotated
    
    tableAllMutations_TFBS.FunSeq2_isMOTIF_any = tableAllMutations_TFBS.FunSeq2_isMOTIFBR | tableAllMutations_TFBS.FunSeq2_isMOTIFG;
    sMatByCADD = struct();
    lstMinCADD_PHRED = [0:2:18,Inf]; 
    %thresholdCADD = 10;
    nBinsCADD = length(lstMinCADD_PHRED)-1;
    for iType = 1:3
        motifType = lstTypes{iType};
        sMatByCADD.(['mat_',motifType]) = NaN*ones(nTissues, nBinsCADD); % pCandidateDriverMutations_motifChange_motifType
        sMatByCADD.(['binary_',motifType]) = NaN*ones(nTissues, 2); % pCandidateDriverMutations_motifChange_motifType
    end

    for jTissue = 1:nTissues
        if (jTissue == 1)
            isOK = tableAllMutations_TFBS.iTissue > 1 & isGeneralOK;
            tissuePrint = 'Pan-cancer Solid';
        else
            iTissue = jTissue;
            isOK = tableAllMutations_TFBS.iTissue==iTissue & isGeneralOK;
            tissuePrint = tableTissues.tissuePrint{iTissue};
        end
        isCandidateDriver = tableAllMutations_TFBS.isCandidateDriver(isOK);
        tableTissuesWithPancancer.tissuePrint{jTissue} = tissuePrint;
        tableTissuesWithPancancer.nControlMutations(jTissue) = sum(~isCandidateDriver);
        tableTissuesWithPancancer.nCandidateDriverMutations(jTissue) = sum(isCandidateDriver);
        for iType = 1:3
            motifType = lstTypes{iType};
            isMotifChange = tableAllMutations_TFBS.(['FunSeq2_is', motifType])(isOK);
            tableTissuesWithPancancer.(['nControlMutations_motifChange_',motifType])(jTissue) = sum(isMotifChange(~isCandidateDriver));
            tableTissuesWithPancancer.(['nCandidateDriverMutations_motifChange_',motifType])(jTissue) = sum(isMotifChange(isCandidateDriver));
            tableTissuesWithPancancer.(['pControlMutations_motifChange_',motifType])(jTissue) = 100*mean(isMotifChange(~isCandidateDriver));
            tableTissuesWithPancancer.(['pCandidateDriverMutations_motifChange_',motifType])(jTissue) = 100*mean(isMotifChange(isCandidateDriver));
            tableTissuesWithPancancer.(['enrichment_motifChange_',motifType])(jTissue) = mean(isMotifChange(isCandidateDriver))/mean(isMotifChange(~isCandidateDriver));
            tableTissuesWithPancancer.(['pValue_motifChange_',motifType])(jTissue) = myFisherTest(isMotifChange, isCandidateDriver, 'both', false);
            for iBin = 1:nBinsCADD
                sMatByCADD.(['mat_',motifType])(jTissue,iBin) = 100*mean(isMotifChange(isCandidateDriver & tableAllMutations_TFBS.CADD_PHRED(isOK)>=lstMinCADD_PHRED(iBin) & tableAllMutations_TFBS.CADD_PHRED(isOK)<lstMinCADD_PHRED(iBin+1)));
            end
            sMatByCADD.(['binary_',motifType])(jTissue,1) = 100*mean(isMotifChange(isCandidateDriver & tableAllMutations_TFBS.CADD_PHRED(isOK)<=4));
            sMatByCADD.(['binary_',motifType])(jTissue,2) = 100*mean(isMotifChange(isCandidateDriver & tableAllMutations_TFBS.CADD_PHRED(isOK)>=10));
        end
    end
    %%
    %     fig = createMaximisedFigure(1); hold on; xValues = lstMinCADD_PHRED(1:end-1);
    %     subplot(1,3,1); plot(xValues, sMatByCADD.mat_MOTIFG(1,:)', 'o', 'MarkerFaceColor','c'); lsline; title('Break or gain'); set(gca, 'FontSize', 12); ylabel('TFBS event (%)'); xlabel('CADD PHRED');
    %     subplot(1,3,2); plot(xValues, sMatByCADD.mat_MOTIFBR(1,:)', 'o', 'MarkerFaceColor','c'); lsline; title('Break'); set(gca, 'FontSize', 12); ylabel('TFBS event (%)'); xlabel('CADD PHRED');
    %     subplot(1,3,3); plot(xValues, sMatByCADD.mat_MOTIF_any(1,:)', 'o', 'MarkerFaceColor','c'); lsline; title('Gain'); set(gca, 'FontSize', 12); ylabel('TFBS event (%)'); xlabel('CADD PHRED');
    %     %lsline; grid on; grid minor;
    %     mySaveAs(fig, imagesPath, 'Fig_TFBS_CADD_bins');
    %     %%
    %     fig = createMaximisedFigure(2);
    %     hB = bar([sMatByCADD.binary_MOTIF_any(1,:); sMatByCADD.binary_MOTIFBR(1,:); sMatByCADD.binary_MOTIFG(1,:)], 'EdgeColor', 'none');
    %     hB(1).FaceColor = [0    0.4470    0.7410];
    %     hB(2).FaceColor = [0.4940    0.1840    0.5560];
    %     legend({'CADD <= 4', 'CADD >= 10'});
    %     ylabel('TFBS event (%)'); set(gca, 'XTickLabel', {'break or gain', 'break', 'gain'}, 'TickLength', [0, 0], 'FontSize', 16);
    %     mySaveAs(fig, imagesPath, 'Fig_TFBS_CADD_binary');
    %% TODELETE
    %%% computed with: isGeneralOK = ~tableAllMutations_TFBS.isExcluded & tableAllMutations_TFBS.FunSeq2_isAnnotated & tableAllMutations_TFBS.isInABCEnhancer & tableAllMutations_TFBS.CADD_PHRED>=10; %  & ~tableAllMutations_TFBS.isIndel
    %     fig = createMaximisedFigure(2);
    %     subplot(1,3,1); bar([tableTissuesWithPancancer_data4.pControlMutations_motifChange_MOTIF_any, tableTissuesWithPancancer_data4.pCandidateDriverMutations_motifChange_MOTIF_any,...
    %         tableTissuesWithPancancer_highCADD.pControlMutations_motifChange_MOTIF_any, tableTissuesWithPancancer_highCADD.pCandidateDriverMutations_motifChange_MOTIF_any]); ylabel('ANY');
    %     subplot(1,3,2); bar([tableTissuesWithPancancer_data4.pControlMutations_motifChange_MOTIFBR, tableTissuesWithPancancer_data4.pCandidateDriverMutations_motifChange_MOTIFBR,...
    %         tableTissuesWithPancancer_highCADD.pControlMutations_motifChange_MOTIFBR, tableTissuesWithPancancer_highCADD.pCandidateDriverMutations_motifChange_MOTIFBR]); ylabel('BREAK');
    %     subplot(1,3,3); bar([tableTissuesWithPancancer_data4.pControlMutations_motifChange_MOTIFG, tableTissuesWithPancancer_data4.pCandidateDriverMutations_motifChange_MOTIFG,...
    %         tableTissuesWithPancancer_highCADD.pControlMutations_motifChange_MOTIFG, tableTissuesWithPancancer_highCADD.pCandidateDriverMutations_motifChange_MOTIFG]); ylabel('GAIN');
    %     legend({'control all SNVs', 'candidates all SNVs', 'control CADD>=10', 'candidates CADD>=10'}, 'Location', 'NorthWest');
    %     mySaveAs(fig, imagesPath, 'Fig_TFBS_CADD_10_ANY_BR_G');
    %%
    toc
    %myPrintMemory
    createDir(fileparts(saveFileData));
    save(saveFileData, 'tableMutations_candidate', 'tableTissues', 'tableTissuesWithPancancer', 'sMatByCADD');
else
    fprintf('Loading data from %s...\n', saveFileData);
    load(saveFileData, 'tableMutations_candidate', 'tableTissues', 'tableTissuesWithPancancer', 'sMatByCADD');
end
