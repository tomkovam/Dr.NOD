function plotFigure5(imagesPath, sColours, tableTissuesWithPancancer_data4, tableTissues_data4, dataTFBS, tableMutations_candidate, tableTissues_data1)


tableMotifs = dataTFBS.tableMotifs;
tableMutations_candidate.tissuePrint = tableTissues_data4.tissuePrint(tableMutations_candidate.iTissue);

%%
isOK = ~tableMutations_candidate.isIndel & ~tableMutations_candidate.isExcluded & tableMutations_candidate.iTissue>1;
tableMutations_candidate = tableMutations_candidate(isOK,:);
% tableMutations_candidate = tableMutations_candidate(tableMutations_candidate.isOK,:);
lstTypes = {'MOTIF_any', 'MOTIFBR', 'MOTIFG'};
lstTypesPrint = {'break or gain', 'break', 'gain'};
%%
fig = createMaximisedFigure(5, [0 0 30 30]);
fontSize = 12;
nR = 4; nC = 3; xS = 0.85; yS = 0.8; xB = 0.1; yB = 0.05; xM = -0.03; yM = 0.012;
%
for iType = 1:3
    myGeneralSubplot(nR,nC,iType,xS,yS,xB,yB,xM,yM); hold on;
    plotTFBS_barPlot(tableTissuesWithPancancer_data4, lstTypes{iType}, lstTypesPrint{iType});
end
% axPos1 = get(gca, 'Position');
%
%         nR = 4; nC = 5; iS = 1; xS = 0.7; yS = 0.8; xB = 0.1; yB = 0.05; xM = -0.03; yM = 0.05;

nR = 6; nC = 9; iS = 2*nC + 1; yS = 0.7;  xM = -0.03;

cListRows = cell(2,1);
cListRows{1} = [143, 89, 336, 218, 119, 283]; % GAIN-UP | 119 TRIM41 breast (nice and ok upregulation) | ID3: 331 (nice but no expression) | 62 (nice but not as impressive upregulation) | 42 MRRF breast (nice but not well knonw) | 68 PRKACA (alt in motif) | 314 BCAR1 ovary (alt in motif) | 90 SLC20A1 breast (alt in motif) | 144 IER3 CRC (alt in motif)
cListRows{2} = [103, 147, 283, 301, 68, 267]; % BREAK-UP | 314, 147, 263, 103, 145 (nice but outside context) | 263 CCND1 not great match, not great upregulation | 145 IKBKB good match, not huge upregulation
% cListRows{2} = [6, 164, 185, 345, 3, 226]; 
% cListRows{2} = [68, 257, 229, 288, 99]; 
% Good: 68 PRKACA, 229 ACD, 6 CDON, 345 ZNF37BP ovary, 
% BREAK-DOWN: 267, 269 (both CLTC lung)
% GAINs in regulatory regions of upregulated genes:
%     motifPrefix    nMOTIFBR_UP    nMOTIFG_UP    nMOTIFBR_DOWN    nMOTIFG_DOWN    isPositiveRegulator    isNegativeRegulator    isPositiveRegulator_prefix    isNegativeRegulator_prefix
%     ___________    ___________    __________    _____________    ____________    ___________________    ___________________    __________________________    __________________________
%     {'ARID3A'}          0             1               0               1                 true                   false                     true                          false           
%     {'ETS'   }          4             3               0               0                 false                  false                     true                          true            
%     {'GATA'  }          2             5               0               0                 false                  false                     true                          true            
%     {'HNF4'  }          3             2               1               0                 false                  false                     true                          true            
%     {'LHX3'  }          0             1               0               0                 true                   false                     true                          false           
%     {'NFATC2'}          0             1               0               0                 true                   true                      true                          true            
%     {'RAD21' }          3             1               0               0                 false                  false                     false                         false           
% BREAKs in regulatory regions of upregulated genes:
%     motifPrefix    nMOTIFBR_UP    nMOTIFG_UP    nMOTIFBR_DOWN    nMOTIFG_DOWN    isPositiveRegulator    isNegativeRegulator    isPositiveRegulator_prefix    isNegativeRegulator_prefix
%     ___________    ___________    __________    _____________    ____________    ___________________    ___________________    __________________________    __________________________
%     {'E2F'   }          4             0               0               0                 false                  false                     true                          true            
%     {'HDAC2' }          5             0               0               0                 true                   true                      true                          true            
%     {'IRF'   }          4             1               1               0                 false                  false                     true                          true            
%     {'MYC'   }          2             0               0               0                 true                   true                      true                          true            
% BREAKs in regulatory regions of downregulated genes:
%     motifPrefix    nMOTIFBR_UP    nMOTIFG_UP    nMOTIFBR_DOWN    nMOTIFG_DOWN    isPositiveRegulator    isNegativeRegulator    isPositiveRegulator_prefix    isNegativeRegulator_prefix
%     ___________    ___________    __________    _____________    ____________    ___________________    ___________________    __________________________    __________________________
%     {'ZNF143'}          1             0               1               0                 true                   false                     true                          false           


for iDirection = 1:2
    for iRow = cListRows{iDirection}
        iTissue = tableMutations_candidate.iTissue(iRow);
        geneName = tableMutations_candidate.candidateGenes{iRow};
        if (contains(geneName, ' '))
            if (contains(geneName, 'PARP2'))
                geneName = 'PARP2';
            elseif (contains(geneName, 'TRIM41'))
                geneName = 'TRIM41';
            else
                tmp2 = strsplit(geneName, ' ');
                geneName = tmp2{end};
            end
        end
        myGeneralSubplot(nR,nC,iS,.6+xS,yS,xB,yB,xM,yM); hold on; iS = iS + 2;
        yAltRelative1 = plotTFBS_logos(tableMutations_candidate, tableMotifs, iRow, iDirection==1, geneName);
        axPos1 = get(gca, 'Position');

        myGeneralSubplot(nR,nC,iS,.5,yS,xB,yB,xM,yM); hold on; iS = iS + 1;
        [xAltRelative2, yAltRelative2] = plotGene_boxplot_forLogos(tableTissues_data4.tissue{iTissue}, tableTissues_data1.biosampleABC{iTissue}, geneName, sColours, tableMutations_candidate.iSample(iRow));
        axPos2 = get(gca, 'Position');

        xa = [axPos1(1) + axPos1(3), axPos2(1) + axPos2(3)*xAltRelative2 - 0.005];
        ya = [axPos1(2) + axPos1(4)*yAltRelative1, axPos1(2) + axPos1(4)*yAltRelative2];
        annotation('arrow',xa,ya, 'Color', sColours.mutated, 'LineWidth', 1.5)
    end
end
% 68, 314, 52, 144

fontSizeAnnotation = 14;
dim = [.07 .55 .05 .05]; str = 'TFBS gain'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
dim = [.07 .40 .05 .05]; str = 'TFBS gain'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
dim = [.07 .23 .05 .05]; str = 'TFBS break'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
dim = [.07 .08 .05 .05]; str = 'TFBS break'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');

fontSizeLetters = 26;
dim = [.005 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.36 .99 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.67 .99 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.005 .66 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.005 .35 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');

mySaveAs(fig, imagesPath, 'Fig5', false, true);
%%
    function yAltRelative = plotTFBS_logos(tableMutations_candidate, tableMotifs, iRow, isMOTIFG, geneName)
        %lstCols = {'candidateGenes', 'MOTIFG_scoreDiff', 'MOTIFG_motifNamePrefix', 'tissuePrint', 'expressionMedianWT', 'expressionThisMut', 'gene', 'MOTIFG_motifName', 'MOTIFG_scoreAlt', 'MOTIFG_scoreRef', 'isHighCADD', 'VAF', 'qtlVAF'};
        %tableMutations_candidate(iRow,lstCols)
        %% 19
        if (isMOTIFG)
            motifName = tableMutations_candidate.MOTIFG_motifName{iRow};
            motifStart = tableMutations_candidate.MOTIFG_motifStart(iRow);
            motifEnd = tableMutations_candidate.MOTIFG_motifEnd(iRow);
            positionMutation = tableMutations_candidate.MOTIFG_positionMutation(iRow);
            isMinusStrand = tableMutations_candidate.MOTIFG_isMinusStrand(iRow);
        else
            % 49
            motifName = tableMutations_candidate.MOTIFBR_motifName{iRow};
            motifStart = tableMutations_candidate.MOTIFBR_motifStart(iRow);
            motifEnd = tableMutations_candidate.MOTIFBR_motifEnd(iRow);
            positionMutation = tableMutations_candidate.MOTIFBR_positionMutation(iRow);
            isMinusStrand = tableMutations_candidate.MOTIFBR_isMinusStrand(iRow);
        end
        %
        context_ref = upper(tableMutations_candidate.context50bp{iRow});
        lengthContext = length(context_ref);
        iPosMutation = ceil(lengthContext/2);
        if (~strcmp(context_ref(iPosMutation), tableMutations_candidate.ref{iRow})), error('Ref bases do not match.'); end
        context_alt = context_ref; context_alt(iPosMutation) = tableMutations_candidate.alt{iRow};
        x1 = motifStart-tableMutations_candidate.pos1(iRow)+1;    % iPosMutation + x1 >= 1 --> x1 >= 1-iPosMutation
        x2 = motifEnd-tableMutations_candidate.pos1(iRow);        % iPosMutation + x2 <= lengthContext --> x2 <= lengthContext-iPosMutation = iPosMutation-1
        if ((x1 < 1-iPosMutation) || (x2 > iPosMutation-1))
            error('Motif outside our context');
        end
        % context_ref(iPosMutation + (x1:x2))
        sequence_ref = context_ref(iPosMutation + (x1:x2));
        sequence_alt = context_alt(iPosMutation + (x1:x2));
        % positionMutation = x2 + 1;
        if (isMinusStrand)
            sequence_ref = seqrcomplement(sequence_ref);
            sequence_alt = seqrcomplement(sequence_alt);
        end
        %

        % strsplit(tableMotifs.motifName{1877}, '→')

        iMotif = find(contains(tableMotifs.motifName, ['>', motifName]));
        if (isempty(iMotif))
            error('Motif not found.')
        elseif (length(iMotif)>1)
            warning('Multiple motif matches')
            iMotif = iMotif(1);
        end

        % iMotif = find(strcmp(cellfun(@(x) x(2:length(motifName)), tableMotifs.motifName), ['>', motifName]));
        % tmp2 = cellfun(@(x) x(2:min([end, length(motifName)+1])), tableMotifs.motifName, 'UniformOutput', false);
        % iMotif = find(strcmp(tmp2, motifName));
        % if (length(iMotif)~=1)
        %     error('Motif not found.')
        % end
        
        yAltRelative = plotMotif(tableMotifs.motifMatrix{iMotif}', 10, 2, sequence_ref, sequence_alt, motifName, positionMutation);
        title(sprintf('%s in %s\n%s:%s', geneName, tableMutations_candidate.tissuePrint{iRow}, ...
            tableMutations_candidate.chr{iRow}, num2sepNumStr(tableMutations_candidate.pos1(iRow))));
        %         title(sprintf('%s %s\n%s:%s %s>%s', tableMutations_candidate.candidateGenes{iRow}, tableMutations_candidate.tissuePrint{iRow}, ...
        %             tableMutations_candidate.chr{iRow}, num2sepNumStr(tableMutations_candidate.pos1(iRow)), tableMutations_candidate.ref{iRow}, tableMutations_candidate.alt{iRow}));
    end
%%
    function plotTFBS_barPlot(tableTissuesWithPancancer, motifType, motifTypePrint)
        matValues = [tableTissuesWithPancancer.(['pControlMutations_motifChange_',motifType]), tableTissuesWithPancancer.(['pCandidateDriverMutations_motifChange_',motifType])];
        xValues = (1:size(tableTissuesWithPancancer, 1))';
        yValues = max(matValues, [], 2);
        hB = bar(matValues, 'EdgeColor', 'flat', 'FaceColor', 'flat');
        hB(1).CData = sColours.gray;
        hB(2).CData = sColours.darkRed;
        %text(xValues, 3 + yValues, tableTissuesWithPancancer.pFisherCDG_text, 'HorizontalAlignment', 'center', 'FontSize', 14);
        text(xValues, 3 + yValues, strcat(num2str(tableTissuesWithPancancer.(['enrichment_motifChange_',motifType]), '%.1fx')), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues, 5 + yValues, arrayfun(@getPValueStarsAsText, tableTissuesWithPancancer.(['pValue_motifChange_',motifType]), 'UniformOutput', false), 'HorizontalAlignment', 'center', 'FontSize', fontSize-2, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        % text(xValues, 2 + yValues, strcat({'FC: '}, num2str(tableTissuesWithPancancer.enrichmentCDG, '%-.1f')), 'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!
        %text(xValues, 1 + yValues, strcat({'n = '}, num2str(tableTissuesWithPancancer.nSamplesWGSandRNA, '%-d')), 'HorizontalAlignment', 'center', 'FontSize', fontSize-8, 'Color', .5*[1,1,1]); % {'n = '} insetad of 'n = ' will keep the space in there!

        maxVal = max(5 + yValues); yGap = maxVal/20;
        yVal = 0.2*yGap + maxVal;        ylim([0, yVal]);
        yVal1 = 1.5*yGap + maxVal;
        yVal2 = 3*yGap + maxVal;
        yVal3 = 4.5*yGap + maxVal;
        
        text(xValues+.3, yVal1+0*xValues, arrayfun(@num2sepNumStr, round(tableTissuesWithPancancer.nControlMutations/1e3), 'UniformOutput', false), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.gray); % {'n = '} insetad of 'n = ' will keep the space in there!
        text(xValues+.3, yVal2+0*xValues, num2str(tableTissuesWithPancancer.nCandidateDriverMutations, '%d'), 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color',  sColours.darkRed); % {'n = '} insetad of 'n = ' will keep the space in there!
            
        if (strcmp(motifType, 'MOTIF_any'))
            text(0, yVal1, 'Control \times 10^3', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', sColours.gray); 
            text(0, yVal2, 'Cand. driver', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color',  sColours.darkRed);
            text(0, yVal3, 'Mutations', 'HorizontalAlignment', 'right', 'FontSize', fontSize-2, 'Color', 'k'); 
        end

        set(gca, 'XTick', xValues, 'XTickLabel', strrep(tableTissuesWithPancancer.tissuePrint, 'wo Blood', 'Solid'), 'XTickLabelRotation', 45, 'FontSize', fontSize, 'TickLength', [0 0]);
        ylabel(['TFBS ', motifTypePrint, ' (%)']);
        legend({'Control', 'Driver'}, 'Location', 'NorthWest', 'FontSize', fontSize-2); legend boxoff 
        box off; xlim([0, xValues(end)+1]);
    end
%%

end