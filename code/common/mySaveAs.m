function mySaveAs(fig, myPath, pictureName, saveToEps, saveDPNG300)
% mySaveAs(fig, myPath, pictureName) saves the fig with pictureName into path as png

createDir(myPath);                                              % Checkes, that the directory exists.

if (nargin > 4 && saveDPNG300)
    print(fig, [myPath, pictureName], '-dpng', '-r400');     % Saves the image.
else
    saveas(fig, [myPath, pictureName], 'png');     % Saves the image.
end


% export_fig([myPath, pictureName], '-m2', '-nocrop')

if (nargin > 3 && saveToEps)
    print(fig, '-depsc2', '-painters', [myPath, pictureName,'.eps']) %THIS IS THE BEST!
end
% saveas(fig, [myPath, pictureName], 'pdf');     % Saves the image.

% plot2svg([myPath, pictureName,'.svg'], fig);   % svg

% fprintf('Now trying myaa...\n');
% myaa(8, 'standard', fig);
% print(fig, [myPath, pictureName], '-dpng', '-r300');     % Saves the image.

% fprintf('Now trying export_fig...\n');
% export_fig([myPath, pictureName], '-m2', '-nocrop')
% % export_fig [myPath, pictureName] -m3