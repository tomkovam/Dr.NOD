function outputText = insertColourIntoText(inputText, colour)

outputText = sprintf('\\color[rgb]{%f,%f,%f}{%s}', colour(1), colour(2), colour(3), inputText);
