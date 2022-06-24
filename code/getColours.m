function sColours = getColours()

sColours.candidate = .7*[1,1,1];
sColours.candidateEdge = .5*[1,1,1];
sColours.cadidateCDG = ([1,0,0]+1)/2;
sColours.cadidateCDGEdge = ([1,0,0]+1)/3;
sColours.other = .9*[1,1,1];
sColours.otherEdge = .8*[1,1,1];
sColours.otherCDG = ([1,0,0]+2)/3;
sColours.otherCDGEdge = ([1,0,0]+2)/4;
sColours.gray = .5*[1,1,1];
sColours.darkRed = [.65,0,0];
sColours.darkCyan = [0,.5,.5];
sColours.lightRed = (1+sColours.darkRed)/2;
sColours.TSG = [0 .5 .5];
sColours.ONCOGENE = [.65,0,0];
sColours.fusion = [.5 0 .5];
sColours.nonCDG = .7*[1,1,1];
sColours.mutated = sColours.darkRed;
sColours.WT = sColours.gray;
sColours.highExpression = [227,170,224]/256; %.2*[1,1,1];
sColours.lowExpression = [64,119,190]/256; %.7*[1,1,1];
sColours.crossTissueColour = [0    0.4470    0.7410];

% https://www.google.com/imgres?imgurl=https%3A%2F%2Fwww.itakeyou.co.uk%2Fwp-content%2Fuploads%2F2020%2F08%2Fcolor-hex-3.jpg&imgrefurl=https%3A%2F%2Fwww.itakeyou.co.uk%2Femerald-green-and-salmon-colour-scheme%2F&tbnid=U_PpsYTopoCZKM&vet=12ahUKEwi27M7lvZT4AhU0JH0KHdnaBwIQMygTegUIARC1Ag..i&docid=6X6SDAQmsAHOVM&w=757&h=1050&q=colour%20pairs%20teal&ved=2ahUKEwi27M7lvZT4AhU0JH0KHdnaBwIQMygTegUIARC1Ag
sColours.closeMutation = [0,112,116]/256; %[3,76,83]/256; %'#034C52'; 
sColours.distantMutation = [243,140,121]/256; %'#F48D79';
sColours.distanceBackground = (1 + .5*[1,1,1] + sColours.closeMutation)/3;
sColours.distanceBackground = (3 + sColours.distantMutation + sColours.closeMutation)/5;