function plotSupFigureSignaturesNormalised(imagesPath, tableMutations_candidate, tableTissues_data3, sProperties, cTrinucleotides)

fig = createMaximisedFigure(6, [0 0 30 10]); 

axes('Position', [.1, .2, .85, .7]);
plotMutationalSignatures_dotplot(tableMutations_candidate, tableTissues_data3, sProperties, true, cTrinucleotides);

mySaveAs(fig, imagesPath, 'SupFig_signaturesNormalised.png', false, true);
savefig([imagesPath, 'SupFig_signaturesNormalised.fig']);
