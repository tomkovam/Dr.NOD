function plotMutationalSignature(vector, textToShow, tableTrinucleotides) % imagesPath, saveNameFile, 

% vector = tableSignaturesCOSMIC.SBS84;
% textToShow = 'SBS84';
% saveNameFile = 'Fig3_signatureSBS84';

vector = 100*vector/sum(vector);

colours.CtoA = [82,195,242]/255; colours.lightCtoA = (colours.CtoA+1)/2;
colours.CtoG = [35,31,32]/255; colours.lightCtoG = (colours.CtoG+1)/2;
colours.CtoT = [230,34,36]/255; colours.lightCtoT = (colours.CtoT+1)/2;
colours.TtoA = 0.8*[203,202,200]/255; colours.lightTtoA = (colours.TtoA+1)/2;
colours.TtoC = [151,213,76]/255; colours.lightTtoC = (colours.TtoC+1)/2;
colours.TtoG = [237,192,195]/255; colours.lightTtoG = (colours.TtoG+1)/2;
colours.cmap = [colours.CtoA; colours.CtoG; colours.CtoT; colours.TtoA; colours.TtoC; colours.TtoG];
colours.gray = .5*[1,1,1];


% fig = createMaximisedFigure(8, [0 0 35 10]);  axes('Position', [0.05, 0.1, 0.95, 0.75]); 

hold on; 
hB = bar(vector, .5,'EdgeColor','flat', 'FaceColor','flat');
yLimVal = [0, ceil(max(vector))]; % [0,20]; 
yVal = -yLimVal(2)/70;
for iBar = 1:96
    colour = colours.cmap(ceil(iBar/16),:);
    hB.CData(iBar,:) = colour;
    textValue = sprintf('\\color[rgb]{0.5,0.5,0.5}%s\\color[rgb]{%f,%f,%f}%s\\color[rgb]{0.5,0.5,0.5}%s', tableTrinucleotides.type{iBar}(1), colour(1), colour(2), colour(3), tableTrinucleotides.type{iBar}(2), tableTrinucleotides.type{iBar}(3)); % colours.gray(1), colours.gray(2), 
    text(iBar, yVal*1, textValue, 'Rotation', 90, 'HorizontalAlignment','right', 'FontName', 'FixedWidth', 'FontWeight', 'bold', 'FontSize', 4); % Consolas
%     text(iBar, yVal*3, tableTrinucleotides.type{iBar}(1), 'Color', colours.gray, 'Rotation', 90, 'HorizontalAlignment','right', 'FontName', 'LucidaSansRegular');
%     text(iBar, yVal*2, tableTrinucleotides.type{iBar}(2), 'Color', colour, 'Rotation', 90, 'HorizontalAlignment','right', 'FontName', 'LucidaSansRegular');
%     text(iBar, yVal*1, tableTrinucleotides.type{iBar}(3), 'Color', colours.gray, 'Rotation', 90, 'HorizontalAlignment','right', 'FontName', 'LucidaSansRegular');
end

maxValY = yLimVal(2); 
epsilon = 0.5; x = epsilon; w = 16-2*epsilon; %h = maxValPos/7; 
y0 = maxValY-yVal; h = -5*yVal;
listTypes = {'CtoA','CtoG','CtoT','TtoA','TtoC','TtoG'}; listTypesPrint = {'C>A','C>G','C>T','T>A','T>C','T>G'}; listTypesPrintRC = {'G>T','G>C','G>A','A>T','A>G','A>C'};
for iType = 1:length(listTypes)
    x2=(iType-1)*16+x; rectangle('Position',[x2,y0,w,h],'FaceColor',colours.([listTypes{iType}]),'EdgeColor','none');
    text(x2+w/2,y0+h, listTypesPrint{iType},'Color',0*[1,1,1],'HorizontalAlignment', 'center','FontSize',16,'FontWeight','bold', 'VerticalAlignment','bottom');
end

text(1, maxValY+yVal, textToShow,'Color',0*[1,1,1],'HorizontalAlignment', 'left','FontSize',18,'FontWeight','bold', 'VerticalAlignment','top');

ylim(yLimVal); set(gca,'Clipping', 'off');
set(gca, 'XColor', 'none', 'YColor', colours.gray, 'XTick', [],'FontWeight','bold'); grid on;
ax = gca;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 12;
ytickformat(ax, 'percentage'); % 'YTickLabel', strcat(get(gca, 'YTickLabel'), '%'), 

% mySaveAs(fig, imagesPath, saveNameFile, false, false);
%% TODELETE from MDS project
% signatureMutationFrequency = [nanmean(mutationalCatalogsMutFreq, 2), zeros(96,1)];
% fig = createMaximisedFigure(2); axes('Position', [0.07,0.03,0.91,0.9]); hold on; n96 = n182/2; zeroVector = zeros(n96,1);
% maxValPos = max(signatureMutationFrequency(1:n96)); maxValNeg = max(signatureMutationFrequency(n96+1:end));
% maxValBoth = max([maxValPos, maxValNeg]);
% for i182 = 1:n182
%     if (i182>n96)
%         xVal = i182-n96; mult = -1; yText = -maxValNeg; alignment = 'right'; alignmentOther = 'left'; otherIndex = (i182-n96);
%     else
%         xVal = i182; mult = 1; yText = maxValPos; alignment = 'left'; alignmentOther = 'right'; otherIndex = (i182+n96);
%     end
%     zeroVectorTmp = zeroVector; zeroVectorTmp(xVal) = mult*signatureMutationFrequency(i182); bar(zeroVectorTmp, 0.5, 'FaceColor', tableSignCatalog.colour(i182,:), 'EdgeColor', 'none');
%     if (signatureMutationFrequency(i182)>maxValBoth/7)
%         text(xVal, mult*signatureMutationFrequency(i182), sprintf(' {\\bf%s} ', tableSignCatalog.contexts{i182}), 'Rotation', 90, 'HorizontalAlignment', alignment, 'FontSize', 12, 'Color', tableSignCatalog.colour(i182,:)); %[' ', num2str(tableSignCatalog.positions(i182)), tableSignCatalog.contexts{i182}, ' ']
%     elseif (signatureMutationFrequency(i182)>maxValBoth/5)
%         text(xVal, mult*signatureMutationFrequency(i182), sprintf(' {%s} ', tableSignCatalog.contexts{i182}), 'Rotation', 90, 'HorizontalAlignment', alignment, 'FontSize', 8, 'Color', tableSignCatalog.colour(i182,:)); %[' ', num2str(tableSignCatalog.positions(i182)), tableSignCatalog.contexts{i182}, ' ']
%         %             if (signatureMutationFrequency(otherIndex)==0)
%         %                 text(xVal, 0, sprintf(' {\\bf%s}', tableSignCatalog.contexts{otherIndex}), 'Rotation', 90, 'HorizontalAlignment', alignmentOther, 'FontSize', 8, 'Color', tableSignCatalog.colour(i182,:)); %[' ', tableSignCatalog.contexts{otherIndex}, ' ']
%         %             end
%     end
% end
% xlim([0.5,n96+.5]);
% set(gca, 'XTick', []); grid on; drawnow;
% %%% WORKS ONLY IN NEWER MATLAB
% ax = gca;
% ax.YAxis.MinorTick = 'on';
% ax.YMinorGrid = 'on';
% ax.MinorGridLineStyle = ':';
% ax.MinorGridColor = 0.4*[1,1,1];
% %%% WORKS ONLY IN NEWER MATLAB-END
%     set(gca, 'FontSize', fontSize);
% ylabel('Mutation frequency', 'FontSize', 20); 
% epsilon = 0.5; x = 0.5+epsilon; y1 = maxValPos*1.25; w = 16-2*epsilon; h = maxValPos/7; y2 = y1+2*h; fontSizeType = 18;
% y0 = maxValPos+h; y1 = -maxValNeg-h*2;
% listTypes = {'CtoA','CtoG','CtoT','TtoA','TtoC','TtoG'}; listTypesPrint = {'C>A','C>G','C>T','T>A','T>C','T>G'}; listTypesPrintRC = {'G>T','G>C','G>A','A>T','A>G','A>C'};
% for iType = 1:length(listTypes)
%     x2=(iType-1)*16+x; rectangle('Position',[x2,y0,w,h],'FaceColor',colours.([listTypes{iType}]),'EdgeColor','none');
%     text(x2+w/2,y0+h/2, listTypesPrint{iType},'Color',1*[1,1,1],'HorizontalAlignment', 'center','FontSize',fontSizeType,'FontWeight','bold');
%     x2=(iType-1)*16+x; rectangle('Position',[x2,y1,w,h],'FaceColor',colours.([listTypes{iType}]),'EdgeColor','none');
%     text(x2+w/2,y1+h/2, listTypesPrintRC{iType},'Color',1*[1,1,1],'HorizontalAlignment', 'center','FontSize',fontSizeType,'FontWeight','bold');
% end
% ylim([y1, y0+h]);
% set(gca,'Clipping', 'off');
% title({'{\itPOLE}-mutated and MMR-deficient cancer patients'}, 'FontSize', 24);
% mySaveAs(fig, imagesPathBasic, ['cancerPatientSignature_', expName, '_', sampleName, '.png']);