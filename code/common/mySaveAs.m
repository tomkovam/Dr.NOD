function mySaveAs(fig, myPath, pictureName, saveToEps, saveDPNG300)
% mySaveAs(fig, myPath, pictureName) saves the fig with pictureName into path as png

createDir(myPath);                                              % Checkes, that the directory exists.

if (nargin > 4 && saveDPNG300)
    try
        print(fig, [myPath, pictureName], '-dpng', '-r400');     % Saves the image.
    catch % If the above does not work, we use standard plotting
        saveas(fig, [myPath, pictureName], 'png');
    end
else
    saveas(fig, [myPath, pictureName], 'png');     % Saves the image.
end


if (nargin > 3 && saveToEps)
    print(fig, '-depsc2', '-painters', [myPath, pictureName,'.eps']) %THIS IS THE BEST!
end
