function plotSupFigure1(imagesPath, dataCutoffs)

sResPanCancer = dataCutoffs.sResPanCancer;
% sResPanCancer_FDR = dataCutoffs.sResPanCancer_FDR;
tableTissues = dataCutoffs.tableTissues;
lstCutoffsPM = dataCutoffs.lstCutoffsPM;
lstCutoffsPE = dataCutoffs.lstCutoffsPE;
% lstCutoffsFDR = dataCutoffs.lstCutoffsFDR;
%%
nCutoffPM = length(lstCutoffsPM);
nCutoffPE = length(lstCutoffsPE);
%%
lstCutoffsPMText = cellstr(num2str(lstCutoffsPM', '%.3f')); lstCutoffsPMText{end} = 'all';
lstCutoffsPEText = cellstr(num2str(lstCutoffsPE', '%.3f')); lstCutoffsPEText{end} = 'all';
% lstCutoffsFDRText = cellstr(num2str(lstCutoffsFDR', '%.2f')); lstCutoffsFDRText{end} = 'all';
lstPrintNames = [tableTissues.tissuePrint', {'Pan-cancer Solid', 'Pan-cancer'}];
nTypes = length(lstPrintNames);
%%
fig = createMaximisedFigure(3, [0 0 35 25]); 
cmap = flipud(lbmap(600,'RedBlue')); cmap(290:310,:) = 1; colormap([.5*[1,1,1]; cmap]);
cmap = flipud(lbmap(600,'RedBlue')); cmap(290:310,:) = 1; colormap([.5*[1,1,1]; cmap]); %colormap([.5*[1,1,1]; cmap(1:end/2-1,:); ones(100,3); cmap(end/2+1:end,:)]);
colormap([.5*[1,1,1]; flipud(lbmap(600,'RedBlue'))]);
nR = round(sqrt(nTypes)); nC = ceil(nTypes/nR); xS = 0.8; yS = 0.75; xB = 0.05; yB = 0.08; xM = -0.02; yM = -0.02;
for iType = 1:nTypes
    myGeneralSubplot(nR,nC,iType,xS,yS,xB,yB,xM,yM);
    matToPlot = log2(squeeze(sResPanCancer.enrichmentCDG(iType,:,:)))';
    maxVal = max(abs(matToPlot(~isinf(matToPlot(:)))));
    imagesc(matToPlot); hBar = colorbar; title(hBar, 'log_2 enr.'); caxis([-maxVal, maxVal]); % caxis(3*[-1,1]);
    set(gca, 'YTick', 1:nCutoffPE, 'YTickLabel', lstCutoffsPEText, 'XTick', 1:nCutoffPM, 'XTickLabel', lstCutoffsPMText, 'TickLength', [0 0], 'XTickLabelRotation', 45, 'FontSize', 7);
    xlabel('{\itp_M} cut-off') % if (iType > nTypes-nC), xlabel('{\itp_M} cut-off'); end
    ylabel('{\itp_E} cut-off');
    title(lstPrintNames{iType}); drawnow;
end
mySaveAs(fig, imagesPath, 'ExtDataFig1.png', false, true);
savefig([imagesPath, 'ExtDataFig1.fig']);
