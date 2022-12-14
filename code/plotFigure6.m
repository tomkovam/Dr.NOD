function plotFigure6(imagesPath, sColours, tableTissues_data1, dataSupTables, tableMutations_candidate)
%% Here, we take the stringent QC together, and annotate the target genes accordingly, and then plot a list of target genes in each tissue and our confidence in them
%%
tmp = dataSupTables.tableGencodeGenesCandidates(:,{'geneSymbol', 'sizeEffectE', 'pE', 'sizeEffectM', 'pM', 'qCombined', 'tissuePrint', 'isDriver', 'isONCOGENE', 'isTSG', 'isUP', 'iTissue', 'literatureEvidenceOncogene', 'literatureEvidenceTSG'});
tmp.abs_sizeEffectE = abs(tmp.sizeEffectE);
tmp = sortrows(tmp,'abs_sizeEffectE','descend');
tmp.literatureEvidence = max([tmp.literatureEvidenceOncogene, tmp.literatureEvidenceTSG], [], 2);
%%
fprintf('In fact, the absolute size effects are larger in blood (median %.1f) than in solid cancers (median %.1f), and the top %d targets with the largest size effects are also in blood (including MYC).\n', ...
    median(tmp.abs_sizeEffectE(tmp.iTissue==1)), median(tmp.abs_sizeEffectE(tmp.iTissue>1)), find(tmp.iTissue>1, 1, 'first')-1);
%%
tmp2 = unique(tableMutations_candidate.candidateGenes(tableMutations_candidate.iTissue==1));
tmp3 = tmp2(contains(tmp2, ' '));
tmp = sortrows(tmp,'pE','ascend'); % qCombined
lstGenesPotentialFalsePositives_group2_blood = [];
for iRow = 1:length(tmp3)
    lstGenes = strsplit(tmp3{iRow}, ' ');
    tmp4 = tmp(ismember(tmp.geneSymbol, lstGenes),:)
    tmp4.geneSymbol(tmp4.pE>tmp4.pE(1))
    lstGenesPotentialFalsePositives_group2_blood = [lstGenesPotentialFalsePositives_group2_blood; tmp4.geneSymbol(tmp4.pE>tmp4.pE(1))];
end
lstGenesPotentialFalsePositives_group2_blood = unique(lstGenesPotentialFalsePositives_group2_blood);
lstGenesPotentialFalsePositives_group1_blood = {'ZNF876P', 'WEE1', 'AICDA', 'BCAT1', 'C12orf77', 'PPM1F', 'TOP3B', 'PRAMENP'}; % computed in scriptBloodLymphomas.m
%%
fig = createMaximisedFigure(2); hold on;
isOK = tmp.iTissue>1;
yValues = tmp.abs_sizeEffectE(isOK); 
xValues = tmp.literatureEvidence(isOK); xValues(xValues>2) = 2;
labels = tmp.geneSymbol(isOK);
boxchart(xValues, yValues);
p = ranksum(yValues(xValues==0), yValues(xValues==2));
yVal = 1.8;
plot([0,2], yVal*[1,1], '.-k', 'LineWidth', 2);
plot(0*[1,1], yVal + [-.02,0], '-k', 'LineWidth', 2);
plot(2*[1,1], yVal + [-.02,0], '-k', 'LineWidth', 2);
text(1, yVal, sprintf('{\\itp = %s}', getPValueAsText(p)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16);
set(gca, 'XTick', 0:2, 'XTickLabel', {'no evidence', 'evidence level 1', 'evidence level 2-4'}, 'FontSize', 16);
yVal = 0.1;
for iValue = [0,1,2]
    text(iValue, yVal, sprintf('(%d genes)', sum(xValues==iValue)), 'HorizontalAlignment', 'center', 'FontSize', 14);
end
ylabel('Expression size effect');
mySaveAs(fig, imagesPath, 'SupFig_expressionSizeEffect');

cutoff_value = quantile(yValues(xValues==0), 1/3);
%%
tmp.potentialFalsePositive = 0*tmp.sizeEffectE;
tmp.potentialFalsePositive(ismember(tmp.geneSymbol, [{'HCG15', 'CPOX', 'CLTC'}, lstGenesPotentialFalsePositives_group1_blood])) = 1; % Genes in flanking regions have also regulatory mutation count higher than expected by the background mutagenesis model.
tmp.potentialFalsePositive(ismember(tmp.geneSymbol, [{'ZFP62', 'CCNB1IP1', 'ALOXE3'}, lstGenesPotentialFalsePositives_group2_blood'])) = 2; % Genes that share regulatory driver mutations (always the gene with the highest scoreE is kept).
tmp.potentialFalsePositive(tmp.abs_sizeEffectE <= cutoff_value) = 3;
tmp = sortrows(tmp,'geneSymbol','ascend'); % qCombined

tmp(tmp.potentialFalsePositive==3 & tmp.iTissue==1, {'geneSymbol', 'tissuePrint'})
%%
tmp(tmp.potentialFalsePositive==3 & tmp.iTissue>1, {'geneSymbol', 'tissuePrint'})
%%
isOK = tmp.iTissue>1 & tmp.isUP & tmp.potentialFalsePositive >= 0;
fprintf('Solid cancer driver-upregulated targets with oncogenic evidence: %.0f %% (%d/%d)\n', 100*mean(tmp.literatureEvidenceOncogene(isOK)>0), sum(tmp.literatureEvidenceOncogene(isOK)>0), sum(isOK));
isOK = tmp.iTissue>1 & tmp.isUP & tmp.potentialFalsePositive == 0;
fprintf('Solid cancer driver-upregulated targets with oncogenic evidence: %.0f %% (%d/%d)\n', 100*mean(tmp.literatureEvidenceOncogene(isOK)>0), sum(tmp.literatureEvidenceOncogene(isOK)>0), sum(isOK));
% %% Solid - in the terms of the percentage of CDGs, the third type of filtering helps in driver-downregulated, but doesn't make a difference in driver-upregulated.
% isOK = tmp.iTissue>1 & tmp.potentialFalsePositive >= 0;% & ~tmp.isUP;
% fprintf('Solid cancer CDGs: %.0f %% (%d/%d)\n', 100*mean(tmp.isDriver(isOK)), sum(tmp.isDriver(isOK)>0), sum(isOK));
% isOK = tmp.iTissue>1 & tmp.potentialFalsePositive == 0;% & ~tmp.isUP;
% fprintf('Solid cancer CDGs: %.0f %% (%d/%d)\n', 100*mean(tmp.isDriver(isOK)), sum(tmp.isDriver(isOK)>0), sum(isOK));
%% Blood - in the terms of the percentage of CDGs, the third type of filtering helps in driver-downregulated, but doesn't make a difference in driver-upregulated.
isOK = tmp.iTissue==1 & tmp.potentialFalsePositive >= 0;% & ~tmp.isUP;
fprintf('Blood cancer CDGs: %.0f %% (%d/%d)\n', 100*mean(tmp.isDriver(isOK)), sum(tmp.isDriver(isOK)>0), sum(isOK));
isOK = tmp.iTissue==1 & tmp.potentialFalsePositive == 0;% & ~tmp.isUP;
fprintf('Blood cancer CDGs: %.0f %% (%d/%d)\n', 100*mean(tmp.isDriver(isOK)), sum(tmp.isDriver(isOK)>0), sum(isOK));
%%
tmp(tmp.iTissue==1 & tmp.potentialFalsePositive>0 & tmp.isDriver,:)
tmp(tmp.potentialFalsePositive>0 & tmp.iTissue==1, :)
%% Fig 5
nTissues = size(tableTissues_data1, 1);
fontSize = 16;
fig = createMaximisedFigure(3); axes('Position', [.02, .05, .95, .88]); hold on;
for iDirection = 1:2
    for iTissue = 2:nTissues
        if (iDirection == 1)
            lstGenes = find(tmp.iTissue == iTissue & tmp.isUP);
            yVal = 0;
        else
            lstGenes = find(tmp.iTissue == iTissue & ~tmp.isUP);
            yVal = 15.5;
        end
        text(iTissue, yVal, tableTissues_data1.tissuePrint{iTissue}, 'FontSize', fontSize);
        for jGene = 1:length(lstGenes)
            iGene = lstGenes(jGene);
            if (tmp.literatureEvidenceOncogene(iGene)>2)
                colour = sColours.ONCOGENE;
            elseif (tmp.literatureEvidenceOncogene(iGene)>0)
                colour = (1+sColours.ONCOGENE)/2;
            elseif (tmp.literatureEvidenceTSG(iGene)>0)
                colour = (1+sColours.TSG)/2;
            else
                colour = .5*[1,1,1];
            end
            text(iTissue, yVal + jGene, tmp.geneSymbol{iGene}, 'FontAngle', 'italic', 'Color', colour, 'FontSize', fontSize);
            if (tmp.potentialFalsePositive(iGene)>0)
                text(iTissue-.08, yVal + jGene, '?');
            end
        end
    end
end
xlim([1.5, nTissues+.7]); ylim([0, 17]);
text(1.5, 0, 'a', 'FontSize', 26);
text(1.5, 15.5, 'b', 'FontSize', 26);
% LEGEND
xVal = 7; yVal = 11; yShift = .7;
text(xVal, yVal, 'Oncogene (strong evidence)', 'Color', sColours.ONCOGENE, 'FontSize', fontSize); yVal= yVal + yShift;
text(xVal, yVal, 'Oncogene (weak evidence)', 'Color', (1+sColours.ONCOGENE)/2, 'FontSize', fontSize); yVal= yVal + yShift;
text(xVal, yVal, 'Tumour-suppressor gene', 'Color', (1+sColours.TSG)/2, 'FontSize', fontSize); yVal= yVal + yShift;
text(xVal, yVal, '? = low confidence', 'FontSize', fontSize); yVal= yVal + .5; % possible false positive
annotation('rectangle',[.73, .22, .26, .17]);
%
set(gca, 'YDir', 'reverse');
axis off
mySaveAs(fig, imagesPath, 'Fig5', true, true);