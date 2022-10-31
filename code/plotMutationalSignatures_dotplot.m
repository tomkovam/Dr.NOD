function plotMutationalSignatures_dotplot(tableMutations_candidate, tableTissues, sProperties, doNormaliseTrinucleotides, cTrinucleotides)

nTissues = size(tableTissues, 1);

tableTrinucleotides = readtable(sProperties.TABLE_TRINUCLEOTIDES); % 'data/tableTriNucl96.txt'
nPatterns = size(tableTrinucleotides, 1);
tableSignaturesCOSMIC = readtable(sProperties.COSMIC_SIGNATURES_SBS); % 'data/COSMIC_v3.2_SBS_GRCh37.txt'
tableSignaturesCOSMIC.refBase = cellfun(@(x) x(3), tableSignaturesCOSMIC.Type, 'UniformOutput', false);             % C from 'G[C>G]G'
tableSignaturesCOSMIC.altBase = cellfun(@(x) x(5), tableSignaturesCOSMIC.Type, 'UniformOutput', false);             % G from 'G[C>G]G'
tableSignaturesCOSMIC.trinucleotide = cellfun(@(x) x([1,3,7]), tableSignaturesCOSMIC.Type, 'UniformOutput', false); % GCG from 'G[C>G]G'
tableSignaturesCOSMIC = sortrows(tableSignaturesCOSMIC,{'refBase', 'altBase', 'trinucleotide'},'ascend');
tableSignaturesCOSMIC.patternName = strcat(tableSignaturesCOSMIC.trinucleotide, '>', tableSignaturesCOSMIC.altBase);
if (~isequal(tableSignaturesCOSMIC.patternName, tableTrinucleotides.patternName)), error('WRONG ORDER'); end
tableSignaturesCOSMIC = tableSignaturesCOSMIC(:,2:end-4);
lstSignatures = tableSignaturesCOSMIC.Properties.VariableNames;
nSignatures = length(lstSignatures);
%% The COSMIC signatures are scaled to the whole-genome trinucleotide counts. 
if (doNormaliseTrinucleotides)
    vecTrinucleotidesWGS = tableTrinucleotides.genome';
    for iSignature = 1:nSignatures
        tableSignaturesCOSMIC.(lstSignatures{iSignature}) = tableSignaturesCOSMIC.(lstSignatures{iSignature})./vecTrinucleotidesWGS';
        tableSignaturesCOSMIC.(lstSignatures{iSignature}) = tableSignaturesCOSMIC.(lstSignatures{iSignature})/sum(tableSignaturesCOSMIC.(lstSignatures{iSignature}));
    end
    [~, indexTrinucleotide] = ismember(tableTrinucleotides.type, cTrinucleotides{1}.tableGenes_mean_trinucleotides.Properties.VariableNames);
end
%%
matSimilarityTissuesSignatures = zeros(nTissues, nSignatures);

for iTissue = 1:nTissues
    isOK = tableMutations_candidate.iTissue == iTissue & ~tableMutations_candidate.isExcluded & tableMutations_candidate.isHighCADD;
    catalogue = histcounts(tableMutations_candidate.iPattern(isOK), 1:nPatterns+1);
    if (doNormaliseTrinucleotides)
        isGene = cTrinucleotides{iTissue}.tableGenesNasserExpressed.isCandidate;
        matGenesTrinucleotides = table2array(cTrinucleotides{iTissue}.tableGenes_mean_trinucleotides(isGene,:)).*cTrinucleotides{iTissue}.tableGenes_annotations.nPositionsInEnhancers(isGene); % tableGenes_mean_trinucleotdies = array2table(matGenes_sum_trinucleotides./tableGenes_annotations.nPositionsInEnhancers);
        matGenesTrinucleotides = matGenesTrinucleotides(:, indexTrinucleotide); % Rescale from 32 to 96 channels
        vecTrinucleotidesCandidateDrivers = sum(matGenesTrinucleotides, 1);
        %vecTrinucleotideFrequencyCandidateDrivers = vecTrinucleotidesCandidateDrivers/sum(vecTrinucleotidesCandidateDrivers); % Frequency of each trinucleotide in all regulatory regions of all regulatory driver candidates in this tissue
        catalogueRescaled = (catalogue./vecTrinucleotidesCandidateDrivers); % Mutation frequency per trinucleotide (normalised by trinucelotide counts)
        catalogueRescaled = catalogueRescaled/sum(catalogueRescaled);
        matSimilarityTissuesSignatures(iTissue,:) = computeSimilarityWithSignatures(tableSignaturesCOSMIC, catalogueRescaled');
    else
        matSimilarityTissuesSignatures(iTissue,:) = computeSimilarityWithSignatures(tableSignaturesCOSMIC, catalogue');
    end
end

isSignAbove50 = max(matSimilarityTissuesSignatures, [], 1)>0.5;
matSimilarityTissuesSignatures_above50 = matSimilarityTissuesSignatures(:,isSignAbove50);
nSignAbove50 = sum(isSignAbove50);

minValue = min(matSimilarityTissuesSignatures_above50(:));
maxValue = max(matSimilarityTissuesSignatures_above50(:));


hold on; cmap = flipud(lbmap(101, 'RedBlue')); fontSizeLarge = 12; fontSizeSmall = 8;

iLeg = 1; hLeg = NaN*ones(10,1);  legValues = cell(10,1); 
for value = linspace(maxValue, minValue, 10)
    valueColour = 1+round(100*(value-minValue)/(maxValue-minValue));
    colour = cmap(valueColour,:);
    hLeg(iLeg) = plot(-1, -1, 'o', 'MarkerSize', 1+floor(valueColour/5), 'Color', colour, 'MarkerFaceColor', (1+colour)/2);
    legValues{iLeg} = sprintf('%.1f', value);
    iLeg = iLeg + 1;
end

for iTissue = 1:nTissues
    for jSignature = 1:nSignAbove50
        value = matSimilarityTissuesSignatures_above50(iTissue, jSignature);
        valueColour = 1+round(100*(value-minValue)/(maxValue-minValue));
        colour = cmap(valueColour,:);
        plot(jSignature, iTissue, 'o', 'MarkerSize', 1+floor(valueColour/5), 'Color', colour, 'MarkerFaceColor', (1+colour)/2); 
        if (value > 0.75)
            text(jSignature, iTissue, sprintf('%.1f', value), 'Rotation', 0, 'HorizontalAlignment','center', 'FontSize', fontSizeSmall);
        end
    end
    drawnow;
end

legend(hLeg, legValues, 'Location', 'EastOutside'); ylim(0.5+[0,nTissues]);

xlabel('SBS mutational signatures'); xlim(0.5+[0,nSignAbove50]);
set(gca, 'YDir', 'reverse', 'YTick', 1:nTissues, 'YTickLabel', tableTissues.tissuePrint, 'XTick', 1:nSignAbove50, 'XTickLabel', strrep(lstSignatures(isSignAbove50), 'SBS', ''), 'XTickLabelRotation', 0, 'TickLength', [0 0], 'FontSize', fontSizeLarge);
