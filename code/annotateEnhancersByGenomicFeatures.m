function [tableGenes_annotations, tableGenes_mean_trinucleotdies, nSamplesInFlanks, tableUE_annotations, tableUE_mean_trinucleotdies] = ...
    annotateEnhancersByGenomicFeatures(runAgain, suffix, minCADD_PHRED, biosampleABC, enhancerAnalysis, tableSamples, tableUniqueEnhancers, tableGenesNasserExpressed, matUniqueEnhancersGenes, sProperties)

fileNameBMM = ['save/annotatedEnhancers_', suffix, '_', num2str(minCADD_PHRED),'.mat'];
if (~runAgain && exist(fileNameBMM, 'file'))
    fprintf('Loading %s...\n', fileNameBMM);
    load(fileNameBMM, 'tableGenes_annotations', 'tableGenes_mean_trinucleotdies', 'nSamplesInFlanks', 'tableUE_annotations', 'tableUE_mean_trinucleotdies'); 
else
    t1 = tic;
    fprintf('annotateEnhancersByGenomicFeatures %s...\n', fileNameBMM);
    tableGenesNasserExpressed.nPositionsInEnhancers = matUniqueEnhancersGenes'*tableUniqueEnhancers.nPositions; % Number of positions in all unique enhancers regulating the given gene (I have checked it gives the same results as the previous slower implementation)
    tableGenes_annotations = tableGenesNasserExpressed(:,{'geneName', 'nPositionsInEnhancers'});
    tableUE_annotations = tableUniqueEnhancers(:,{'name', 'nPositions'});
    %nUE = size(tableUE_annotations, 1);
    nGenes = size(tableGenes_annotations, 1);
    %% We read the genomic features for each non-coding enhancer region and the trinucleotide counts in them
    % This file was created in nasserPreparations_step2.sh
    inFile = sProperties.ABC_ANNOTATED_PREDICTORS; % [DIR_DATA_ORIG, 'AnnotatedEnhancers\',enhancerAnalysis,'Predictions.together.6cols.replicationTiming.gc.blacklisted.txt'];
    tableAnnotatedEnhancers = readtable(inFile);
    tableAnnotatedEnhancers.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'name', 'baseActivity', 'biosample', 'replicationTiming', 'GC', 'blacklisted'};
    %% We load annotations of the theoretical number of CADD mutations with PHRED >= minCADD_PHRED
    % This file was created in mixMT\singer2\ mapCADDtoNasser2021.sh and mapCADD_bin_toAllABC.sh
    % OLD: inFile = [DIR_DATA_ORIG, 'AnnotatedEnhancers\AnnotatedByCADD\',enhancerAnalysis,'Predictions.together.6cols.CADD.geq',num2str(minCADD_PHRED),'.txt'];
    inFile = [sProperties.CADD_DIRECTORY,enhancerAnalysis,'Predictions.together.6cols.CADD.geq',num2str(minCADD_PHRED),'.txt'];
    tableUE_annotatedCADD = readtable(inFile);
    tableUE_annotatedCADD.Properties.VariableNames = {'chr', 'pos0', 'pos1', 'name', 'baseActivity', 'biosample', 'nMutations_PHRED_geqCUTOFF'};
    if (~isequal(tableUE_annotatedCADD(:,1:6), tableAnnotatedEnhancers(:,1:6))), error('ERROR: Columns 1-6 in tableUE_annotatedCADD and tableAnnotatedEnhancers do not match.'); end
    tableAnnotatedEnhancers.nMutations_PHRED_geqCUTOFF = tableUE_annotatedCADD.nMutations_PHRED_geqCUTOFF;
    if (~isequal(size(unique(tableAnnotatedEnhancers(:,1:6))), size(tableAnnotatedEnhancers(:,1:6)))), error('ERROR: duplicate rows in tableAnnotatedEnhancers(:,1:6).'); end
    %% We load trinucleotide counts
    inFile = sProperties.ABC_ANNOTATED_PREDICTORS_TRINUCLEOTIDES; %[DIR_DATA_ORIG, 'AnnotatedEnhancers\',enhancerAnalysis,'Predictions.together.trinucleotides.txt'];
    tableAnnotatedEnhancers_trinucleotides = readtable(inFile);
    tableAnnotatedEnhancers_trinucleotides = tableAnnotatedEnhancers_trinucleotides(:,2:end);        % 2:end because this file starts with an empty column atm
    iColStartTrinucleotides = size(tableAnnotatedEnhancers, 2) + 1;
    if (~isequal(size(tableAnnotatedEnhancers, 1), size(tableAnnotatedEnhancers_trinucleotides, 1))), error('ERROR: Number of rows in tableAnnotatedEnhancers_trinucleotides and tableAnnotatedEnhancers do not match.'); end
    tableAnnotatedEnhancers = [tableAnnotatedEnhancers, tableAnnotatedEnhancers_trinucleotides];
    lstTrinucleotides = tableAnnotatedEnhancers_trinucleotides.Properties.VariableNames;     
    %% We find the correct enhancer name (of which enhancersCT should be a substring)
    tmp = unique(tableAnnotatedEnhancers.biosample);
    if (sum(ismember(biosampleABC, tmp))==1)
        biosampleABC_full = biosampleABC;
    elseif (sum(ismember(biosampleABC, tmp))==0 && sum(contains(tmp, biosampleABC))==1)
        biosampleABC_full = tmp{contains(tmp, biosampleABC)};
        fprintf('INFO: biosample %s not found in AnnotatedEnhancers, however %s was found.\n', biosampleABC, biosampleABC_full);
    else
        error('INFO: biosample %s NOT found in AnnotatedEnhancers.\n', biosampleABC);
    end
    %% We select only rows belonging to the enhancersCT_full biosample
    isOK = strcmp(tableAnnotatedEnhancers.biosample, biosampleABC_full);
    tableAnnotatedEnhancers_1 = unique(tableAnnotatedEnhancers(isOK,:));
    tableAnnotatedEnhancers_2 = tableAnnotatedEnhancers(isOK,:);
    if (~isequal(tableAnnotatedEnhancers_1, tableAnnotatedEnhancers_2)), error('ERROR: duplicate rows in tableAnnotatedEnhancers'); end
    tableAnnotatedEnhancers = tableAnnotatedEnhancers(isOK,:);
    %% Check that trinucleotide counts are as expected
    tmp1 = sum(table2array(tableAnnotatedEnhancers(:,iColStartTrinucleotides:end)), 2);
    tmp2 = tableAnnotatedEnhancers.pos1 - tableAnnotatedEnhancers.pos0;
    %     sum(tmp1 ~= tmp2) % The sum of trinucleotides does not match length of the region in only two regions |
    if (100*mean(tmp1 ~= tmp2)>1), error('ERROR: the sum of trinucleotides does not match the number of positions in %d (%.1f%%, which is >1%%) unique enhancers.', sum(tmp1 ~= tmp2), 100*mean(tmp1 ~= tmp2)); end
    %%  We link these regions to the unique enhancers and compute average/total value per enhancer
    [tableAnnotatedEnhancers.isInCurrentTissue, tableAnnotatedEnhancers.indexUniqueEnhancer] = ismember(tableAnnotatedEnhancers.name, tableUE_annotations.name);
    % Some enhancers may not be found - for genes that do not have expression data (we removed them from tableUniqueEnhancers and tableGenesNasserExpressed)
    % if (mean(tableAnnotatedEnhancers.isInCurrentTissue)<1), warning('WARNING: %d unique enhancers not found.\n', sum(~tableAnnotatedEnhancers.isInCurrentTissue)); end 

    tmp = grpstats(tableAnnotatedEnhancers(tableAnnotatedEnhancers.isInCurrentTissue,{'indexUniqueEnhancer', 'nMutations_PHRED_geqCUTOFF'}), 'indexUniqueEnhancer', 'sum');
    tableUE_annotations.nTheoreticalMutations_PHRED_geqCUTOFF = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.nTheoreticalMutations_PHRED_geqCUTOFF(tmp.indexUniqueEnhancer) = tmp.sum_nMutations_PHRED_geqCUTOFF;

    tmp = grpstats(tableAnnotatedEnhancers(tableAnnotatedEnhancers.isInCurrentTissue,[{'indexUniqueEnhancer', 'replicationTiming', 'GC', 'blacklisted', 'baseActivity'}, lstTrinucleotides]), 'indexUniqueEnhancer', 'mean');
    % We don't need the mean_baseActivity in the end, as we already have tableUniqueEnhancers.activity_base (and I checked that it matches)
    tableUE_annotations.mean_baseActivity = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.mean_baseActivity(tmp.indexUniqueEnhancer) = tmp.mean_baseActivity;
    tableUE_annotations.mean_replicationTiming = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.mean_replicationTiming(tmp.indexUniqueEnhancer) = tmp.mean_replicationTiming;
    tableUE_annotations.mean_GC = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.mean_GC(tmp.indexUniqueEnhancer) = tmp.mean_GC;
    tableUE_annotations.mean_blacklisted = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.mean_blacklisted(tmp.indexUniqueEnhancer) = tmp.mean_blacklisted;
    tmp = grpstats(tableAnnotatedEnhancers(tableAnnotatedEnhancers.isInCurrentTissue,[{'indexUniqueEnhancer'}, lstTrinucleotides]), 'indexUniqueEnhancer', 'sum');
    matUniqueEnhancers_sum_trinucleotides = NaN*ones(size(tableUE_annotations, 1), length(lstTrinucleotides));
    matUniqueEnhancers_sum_trinucleotides(tmp.indexUniqueEnhancer,:) = table2array(tmp(:,strcat('sum_', lstTrinucleotides)));
    tableUE_mean_trinucleotdies = array2table(matUniqueEnhancers_sum_trinucleotides./tableUE_annotations.nPositions);
    tableUE_mean_trinucleotdies.Properties.VariableNames = lstTrinucleotides;
    %% We read the flanking mutations around the enhancers - we can read the Slop250bpNoncoding version, as we then take +- 50kb and remove non-coding and enhancers anyway, so it is fine.
    inFile = sProperties.PCAWG_PROJECT_CODES; % [DIR_DATA_NEW, 'PROJECT_CODES.txt']
    tableProjectCodes = readtable(inFile, 'ReadVariableNames', false, 'Delimiter', '\t');
    inFile = [sProperties.ABC_ANNOTATED_PREDICTORS_FLANKING_MF_PREFIX,biosampleABC_full,'.slop50kbp.SNVs.uniq.txt']; % [DIR_DATA_ORIG, 'AnnotatedEnhancers/AnnotatedBySNVs/Slop250bpNoncodingPredictions.',biosampleABC_full,'.slop50kbp.SNVs.uniq.txt']
    tableEnhancersFlankingSNVs = readtable(inFile, 'ReadVariableNames', false, 'Delimiter', '\t');
    tableEnhancersFlankingSNVs.Properties.VariableNames = [{'chr', 'pos0', 'pos1', 'name', 'baseActivity', 'biosample'}, tableProjectCodes.Var1'];
    if (~isequal(size(unique(tableEnhancersFlankingSNVs)), size(tableEnhancersFlankingSNVs))), error('tableEnhancersFlankingSNVs has duplicate rows'); end
    %% We compute a sum over all cancer projects belonging to this tissue
    tableSamples.projectNameSave = tableSamples.project;
    isToBeShortened = count(tableSamples.project, '_')>1;
    tableSamples.projectNameSave(isToBeShortened) = cellfun(@(x) x(1:max(strfind(x, '_'))-1), tableSamples.project(isToBeShortened), 'UniformOutput', false);
    lstProjects = unique(tableSamples.projectNameSave)';
    %
    nSamplesInFlanks = 0; % This is the number of samples used for computation of the flanking-regions mutation frequency (this may be higher than nSamples = rows in tableSamples, as those filter out multiple samples from one donor).
    tableEnhancersFlankingSNVs.nMutatedSamplesInFlanks = 0*tableEnhancersFlankingSNVs.pos1;
    for iProject = 1:length(lstProjects)
        projectName = ['PCAWG_', lstProjects{iProject}];
        if (ismember(projectName, tableEnhancersFlankingSNVs.Properties.VariableNames))
            tableEnhancersFlankingSNVs.nMutatedSamplesInFlanks = tableEnhancersFlankingSNVs.nMutatedSamplesInFlanks + tableEnhancersFlankingSNVs.(projectName);
            inFile = [sProperties.PCAWG_SAMPLES_DIR,'samples.',projectName,'.WGS.txt']; % [DIR_DATA_NEW, 'out04_samples/samples.',projectName,'.WGS.txt']
            tableSamplesOneProjectOne = readtable(inFile, 'ReadVariableNames', false, 'Delimiter', '\t');
            nSamplesInFlanks = nSamplesInFlanks + size(tableSamplesOneProjectOne, 1);
        else
            warning('WARNING: %s not found as a column of tableEnhancersFlankingSNVs', projectName);
        end
    end
    %% Then we link these regions to the unique enhancers and compute the sum of mutations in any non-coding regions flanking the enhancers
    [tableEnhancersFlankingSNVs.isInCurrentTissue, tableEnhancersFlankingSNVs.indexUniqueEnhancer] = ismember(tableEnhancersFlankingSNVs.name, tableUE_annotations.name); % NOTE: some tableEnhancersFlankingSNVs.isInCurrentTissue=0, as those are enhancers fully in coding regions
    tableEnhancersFlankingSNVs.nFlankingPositions = tableEnhancersFlankingSNVs.pos1 - tableEnhancersFlankingSNVs.pos0;
    tmp = grpstats(tableEnhancersFlankingSNVs(tableEnhancersFlankingSNVs.isInCurrentTissue,{'indexUniqueEnhancer', 'nMutatedSamplesInFlanks', 'nFlankingPositions'}), 'indexUniqueEnhancer', 'sum');
    tableUE_annotations.nFlankingRegions = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.nFlankingRegions(tmp.indexUniqueEnhancer) = tmp.GroupCount; 
    tableUE_annotations.sum_nFlankingPositions = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.sum_nFlankingPositions(tmp.indexUniqueEnhancer) = tmp.sum_nFlankingPositions; 
    tableUE_annotations.sum_nMutatedSamplesInFlanks = NaN*tableUE_annotations.nPositions;
    tableUE_annotations.sum_nMutatedSamplesInFlanks(tmp.indexUniqueEnhancer) = tmp.sum_nMutatedSamplesInFlanks; % If the same sample has mutations in multiple non-coding regions flanking one enhancer, they will count as >1 (however, within each region, we will count always just one)
    tableUE_annotations.mfInFlanks = tableUE_annotations.sum_nMutatedSamplesInFlanks./(tableUE_annotations.sum_nFlankingPositions * nSamplesInFlanks);
    %% Average value per gene
    matGenes_sum_trinucleotides = NaN*ones(nGenes, length(lstTrinucleotides));
    tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF = 0*tableGenes_annotations.nPositionsInEnhancers;
    tableGenes_annotations.mean_baseActivity = NaN*tableGenes_annotations.nPositionsInEnhancers;
    tableGenes_annotations.mean_replicationTiming = NaN*tableGenes_annotations.nPositionsInEnhancers;
    tableGenes_annotations.mean_GC = NaN*tableGenes_annotations.nPositionsInEnhancers;
    tableGenes_annotations.mean_mfInFlanks = NaN*tableGenes_annotations.nPositionsInEnhancers;
    for iGene = 1:nGenes
        matGenes_sum_trinucleotides(iGene,:) = sum(matUniqueEnhancers_sum_trinucleotides(matUniqueEnhancersGenes(:,iGene), :), 1);
        tableGenes_annotations.nTheoreticalMutations_PHRED_geqCUTOFF(iGene) = sum(tableUE_annotations.nTheoreticalMutations_PHRED_geqCUTOFF(matUniqueEnhancersGenes(:,iGene)));
        tableGenes_annotations.mean_baseActivity(iGene) = mean(tableUE_annotations.mean_baseActivity(matUniqueEnhancersGenes(:,iGene)));
        tableGenes_annotations.mean_replicationTiming(iGene) = mean(tableUE_annotations.mean_replicationTiming(matUniqueEnhancersGenes(:,iGene)));
        tableGenes_annotations.mean_GC(iGene) = mean(tableUE_annotations.mean_GC(matUniqueEnhancersGenes(:,iGene)));
        tableGenes_annotations.mean_mfInFlanks(iGene) = mean(tableUE_annotations.mfInFlanks(matUniqueEnhancersGenes(:,iGene)));
    end
    tableGenes_mean_trinucleotdies = array2table(matGenes_sum_trinucleotides./tableGenes_annotations.nPositionsInEnhancers);
    tableGenes_mean_trinucleotdies.Properties.VariableNames = lstTrinucleotides;
    %% Finally, we save the resulting file
    save(fileNameBMM, 'tableGenes_annotations', 'tableGenes_mean_trinucleotdies', 'nSamplesInFlanks', 'tableUE_annotations', 'tableUE_mean_trinucleotdies');
end
%%
if (exist('tableGenesNasserExpressed', 'var'))
    if (~isequal(tableGenesNasserExpressed.geneName, tableGenes_annotations.geneName)), error('ERROR gene names do not match.'); end
end
if (exist('tableUniqueEnhancers', 'var'))
    if (~isequal(tableUniqueEnhancers.name, tableUE_annotations.name)), error('ERROR UE names do not match.'); end
end

