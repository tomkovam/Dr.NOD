function fig = createMaximisedFigure(figNumber, PaperPositionValue)
% createMaximisedFigure(figNumber)

if (~exist('PaperPositionValue', 'var'))
    PaperPositionValue = [0 0 30 20];
end

fig = figure(figNumber); clf(fig);  set(gcf,'PaperUnits','centimeters','PaperPosition',PaperPositionValue,'units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'Position', get(0,'Screensize'));

% fig = figure(1); clf(fig);  set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20],'units','normalized','outerposition',[0 0 1 1]);
% nR = 5; nC = 5; iS = 1; xS = 0.8; yS = 0.8; xB = 0.05; yB = 0.05; xM = 0.01; yM = 0.01;
% myGeneralSubplot(nR,nC,iS,xS,yS,xB,yB,xM,yM); iS = iS + 1; hold on;