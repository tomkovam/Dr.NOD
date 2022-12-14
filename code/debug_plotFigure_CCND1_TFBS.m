function plotFigure_CCND1_TFBS(imagesPath, sColours, tableTissuesWithPancancer_data4, tableTissues_data4, dataTFBS, tableMutations_candidate, tableTissues_data1)


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
% for iType = 1:3
%     myGeneralSubplot(nR,nC,iType,xS,yS,xB,yB,xM,yM); hold on;
%     plotTFBS_barPlot(tableTissuesWithPancancer_data4, lstTypes{iType}, lstTypesPrint{iType});
% end
% axPos1 = get(gca, 'Position');
%
%         nR = 4; nC = 5; iS = 1; xS = 0.7; yS = 0.8; xB = 0.1; yB = 0.05; xM = -0.03; yM = 0.05;

nR = 6; nC = 9; iS = 2*nC + 1; yS = 0.7;  xM = -0.03;

cListRows = cell(2,1);
% cListRows{1} = [143, 89, 336, 218, 119, 283]; % GAIN-UP | 119 TRIM41 breast (nice and ok upregulation) | ID3: 331 (nice but no expression) | 62 (nice but not as impressive upregulation) | 42 MRRF breast (nice but not well knonw) | 68 PRKACA (alt in motif) | 314 BCAR1 ovary (alt in motif) | 90 SLC20A1 breast (alt in motif) | 144 IER3 CRC (alt in motif)
% cListRows{2} = [103, 147, 283, 301, 68, 267]; % BREAK-UP | 314, 147, 263, 103, 145 (nice but outside context) | 263 CCND1 not great match, not great upregulation | 145 IKBKB good match, not huge upregulation

cListRows{1} = [143, 336, 301, 89, 218, 119]; % GAIN-UP | 62-PPM1D | 119 TRIM41 breast (nice and ok upregulation) | ID3: 331 (nice but no expression) | 62 (nice but not as impressive upregulation) | 42 MRRF breast (nice but not well knonw) | 68 PRKACA (alt in motif) | 314 BCAR1 ovary (alt in motif) | 90 SLC20A1 breast (alt in motif) | 144 IER3 CRC (alt in motif)
cListRows{2} = [256, 262, 263]; % BREAK-UP 


find(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))
tableMutations_candidate(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'),:)
tableMutations_candidate(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'),{'mutID', 'expressionThisMut', 'FunSeq2_motif_analysis'})
 tableMutations_candidate.mutID(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))
%     {'chr11_69064853_G_A'}
%     {'chr11_69453653_C_T'}
%     {'chr11_69451432_G_T'}
tableMutations_candidate.expressionThisMut(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))
%   159.4191
%    41.5170
%    90.1008

tableMutations_candidate.distance_thisGene(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))
%       391019
%         2219
%         4440

tableMutations_candidate.FunSeq2_motif_analysis(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))

%     {['MOTIFBR=' ...
%         'BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A#HIC1_4#69064838#69064856#-#4#0.000343#0.989722,' ...
%         'BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A#TFAP2_known13#69064842#69064855#-#3#0.000000#1.000000,' ...
%         'BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A#RAD21_disc8#69064843#69064857#-#5#0.000000#0.519947,' ...
%         'BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A#AP1_disc10#69064847#69064857#+#6#0.000000#1.000000']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 }
%     {['MOTIFBR=' ...
%         'DHS,EGR1,STAT1#ZIC2_2#69453645#69453660#-#8#0.071535#0.631760,' ...
%         'DHS,EGR1,STAT1#ZIC3_2#69453645#69453660#-#8#0.059681#0.665195,' ...
%         'DHS,EGR1,STAT1#TFAP2_known13#69453646#69453659#-#7#0.125887#0.292553,' ...
%         'DHS,EGR1,STAT1#ZIC1_2#69453646#69453660#-#8#0.045273#0.607615']                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      }
%     {['MOTIFBR=' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#ATF3_disc2#69451415#69451436#+#17#0.005405#0.978379,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known1#69451424#69451438#+#8#0.000000#1.000000,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known1#69451424#69451438#-#7#0.000000#1.000000,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MXI1_disc2#69451425#69451435#+#7#0.000000#1.000000,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known15#69451425#69451435#+#7#0.000000#0.941176,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_disc2#69451425#69451436#-#5#0.000000#0.995726,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known16#69451425#69451436#+#7#0.000000#0.952381,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#CREB3L1_2#69451425#69451437#-#6#0.000297#0.997329,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known14#69451425#69451437#+#7#0.000000#1.000000,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MNT_1#69451426#69451436#+#6#0.000000#0.978448,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYCN_2#69451426#69451436#+#6#0.002283#0.924658,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known18#69451426#69451436#+#6#0.008811#0.898678,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known21#69451426#69451436#+#6#0.002616#0.946716,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_known6#69451426#69451436#-#5#0.000000#1.000000,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MYC_disc2#69451426#69451437#+#6#0.000000#0.816239,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MXI1_disc2#69451427#69451437#-#6#0.000000#0.967442,' ...
%         'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#NFE2_disc3#69451428#69451438#+#4#0.000000#0.689013']}



tableMutations_candidate.FunSeq2_ENCODE_annotated(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))
%     {'DHS(MCV-87|chr11:69064800-69064950),Enhancer(chmm/segway|chr11:69064400-69065600),TFM(BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A|AP1_disc10|chr11:69064847-69064857),TFM(BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A|HIC1_4|chr11:69064838-69064856),TFM(BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A|RAD21_disc8|chr11:69064843-69064857),TFM(BCL11A,BCL3,CEBPB,DHS,E2F1,EBF1,EGR1,ELF1,ELK4,EP300,ESR1,GATA2,HDAC2,MAX,MEF2A,MEF2C,MYC,NFKB1,NR3C1,PAX5,POU2F2,RXRA,SMC3,STAT3,TBP,TCF12,TCF7L2,TRIM28,USF1,YY1,ZBTB7A|TFAP2_known13|chr11:69064842-69064855),TFP(BATF|chr11:69064314-69066631),TFP(BCL11A|chr11:69064458-69066514),TFP(BCL3|chr11:69064219-69065696),TFP(CCNT2|chr11:69064456-69066408),TFP(CEBPB|chr11:69063894-69065534),TFP(CEBPB|chr11:69063942-69065118),TFP(CHD2|chr11:69064563-69066401),TFP(E2F1|chr11:69064720-69065276),TFP(E2F6|chr11:69062514-69065530),TFP(EBF1|chr11:69064485-69066959),TFP(EBF1|chr11:69064567-69065738),TFP(EGR1|chr11:69064290-69066760),TFP(ELF1|chr11:69064165-69067098),TFP(ELF1|chr11:69064516-69065547),TFP(ELK4|chr11:69064261-69065087),TFP(EP300|chr11:69063969-69065374),TFP(EP300|chr11:69064497-69065271),TFP(EP300|chr11:69064518-69066303),TFP(EP300|chr11:69064519-69066332),TFP(EP300|chr11:69064681-69065287),TFP(ESR1|chr11:69064519-69065286),TFP(ESR1|chr11:69064543-69065239),TFP(FOSL1|chr11:69064752-69066374),TFP(FOSL2|chr11:69064819-69065304),TFP(FOS|chr11:69064276-69066476),TFP(FOS|chr11:69064465-69066394),TFP(FOS|chr11:69064474-69065240),TFP(GATA2|chr11:69064241-69065317),TFP(HDAC2|chr11:69063892-69065531),TFP(HDAC2|chr11:69064117-69065413),TFP(HMGN3|chr11:69064698-69065224),TFP(IRF4|chr11:69064390-69066745),TFP(IRF4|chr11:69064395-69066342),TFP(JUNB|chr11:69064197-69066326),TFP(JUND|chr11:69061203-69066745),TFP(JUND|chr11:69064017-69065362),TFP(JUND|chr11:69064739-69065256),TFP(JUN|chr11:69064092-69065449),TFP(JUN|chr11:69064479-69065267),TFP(JUN|chr11:69064539-69065284),TFP(JUN|chr11:69064637-69065280),TFP(JUN|chr11:69064664-69065392),TFP(MAFK|chr11:69064568-69065381),TFP(MAX|chr11:69060764-69066412),TFP(MAX|chr11:69064074-69066371),TFP(MAX|chr11:69064109-69065261),TFP(MEF2A|chr11:69064227-69066471),TFP(MEF2C|chr11:69064171-69065584),TFP(MXI1|chr11:69064100-69066334),TFP(MYC|chr11:69064021-69066353),TFP(MYC|chr11:69064079-69065543),TFP(MYC|chr11:69064092-69065208),TFP(MYC|chr11:69064097-69065223),TFP(MYC|chr11:69064252-69065251),TFP(MYC|chr11:69064548-69066315),TFP(MYC|chr11:69064745-69065632),TFP(NFKB1|chr11:69064113-69066538),TFP(NFKB1|chr11:69064209-69066423),TFP(NFKB1|chr11:69064230-69066443),TFP(NFKB1|chr11:69064267-69066367),TFP(NFKB1|chr11:69064279-69066382),TFP(NFKB1|chr11:69064282-69066381),TFP(NFKB1|chr11:69064392-69066562),TFP(NFKB1|chr11:69064420-69066352),TFP(NFKB1|chr11:69064440-69066378),TFP(NFKB1|chr11:69064457-69066405),TFP(NFKB1|chr11:69064475-69065626),TFP(NR3C1|chr11:69064482-69065262),TFP(NR3C1|chr11:69064664-69066349),TFP(PAX5|chr11:69064116-69066977),TFP(PAX5|chr11:69064265-69066510),TFP(PAX5|chr11:69064439-69066404),TFP(PAX5|chr11:69064496-69065524),TFP(PBX3|chr11:69064552-69066144),TFP(POU2F2|chr11:69064366-69066362),TFP(POU2F2|chr11:69064448-69066362),TFP(RAD21|chr11:69064009-69065387),TFP(RAD21|chr11:69064344-69066428),TFP(RAD21|chr11:69064417-69065277),TFP(RAD21|chr11:69064628-69065173),TFP(RFX5|chr11:69064011-69065255),TFP(RXRA|chr11:69064291-69065209),TFP(SETDB1|chr11:69064360-69065029),TFP(SIN3A|chr11:69064066-69066570),TFP(SMARCA4|chr11:69063920-69066402),TFP(SMARCB1|chr11:69064129-69066475),TFP(SMARCC1|chr11:69064093-69068174),TFP(SMARCC2|chr11:69064108-69066100),TFP(SMC3|chr11:69063978-69066516),TFP(SMC3|chr11:69064222-69066414),TFP(SMC3|chr11:69064586-69065191),TFP(SP1|chr11:69064242-69066456),TFP(SPI1|chr11:69064563-69066404),TFP(SPI1|chr11:69064579-69065607),TFP(SPI1|chr11:69064615-69065503),TFP(SRF|chr11:69064440-69066708),TFP(STAT3|chr11:69064520-69065448),TFP(STAT3|chr11:69064671-69065472),TFP(STAT3|chr11:69064692-69065476),TFP(STAT3|chr11:69064693-69065284),TFP(STAT3|chr11:69064725-69065473),TFP(TAF1|chr11:69063981-69065490),TFP(TAF1|chr11:69064108-69066412),TFP(TAF1|chr11:69064512-69065262),TFP(TBP|chr11:69064263-69065341),TFP(TCF12|chr11:69064471-69066384),TFP(TCF4|chr11:69064192-69065317),TFP(TCF4|chr11:69064692-69065161),TFP(TFAP2A|chr11:69064055-69065229),TFP(TFAP2C|chr11:69064047-69066565),TFP(USF1|chr11:69064687-69065082),TFP(YY1|chr11:69063532-69066783),TFP(YY1|chr11:69064071-69066327),TFP(ZBTB33|chr11:69064278-69065529),TFP(ZBTB7A|chr11:69064089-69065457),TFP(ZEB1|chr11:69064261-69066694)'}
%     {'DHS(MCV-13|chr11:69453585-69453735),Enhancer(chmm/segway|chr11:69453200-69454000),TFM(DHS,EGR1,STAT1|TFAP2_known13|chr11:69453646-69453659),TFM(DHS,EGR1,STAT1|ZIC1_2|chr11:69453646-69453660),TFM(DHS,EGR1,STAT1|ZIC2_2|chr11:69453645-69453660),TFM(DHS,EGR1,STAT1|ZIC3_2|chr11:69453645-69453660),TFP(CEBPB|chr11:69453585-69454047),TFP(E2F1|chr11:69452798-69454307),TFP(ELF1|chr11:69453651-69454078),TFP(EP300|chr11:69452257-69454126),TFP(EP300|chr11:69453040-69454513),TFP(EP300|chr11:69453244-69454365),TFP(EP300|chr11:69453601-69454211),TFP(FAM48A|chr11:69453470-69453958),TFP(HDAC2|chr11:69452633-69454468),TFP(JUND|chr11:69452804-69454403),TFP(RFX5|chr11:69452751-69454465),TFP(RFX5|chr11:69453143-69454303),TFP(RFX5|chr11:69453217-69456773),TFP(SMARCC1|chr11:69453473-69454118),TFP(SMC3|chr11:69453610-69454187),TFP(STAT1|chr11:69453270-69454190),TFP(TAF1|chr11:69452305-69454210),TFP(TAF7|chr11:69452547-69454479),TFP(TBP|chr11:69452289-69454199),TFP(TCF4|chr11:69451874-69454359),TFP(YY1|chr11:69452776-69454396)'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      }
%     {'DHS(MCV-96|chr11:69451380-69451530),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|ATF3_disc2|chr11:69451415-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|CREB3L1_2|chr11:69451425-69451437),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MNT_1|chr11:69451426-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MXI1_disc2|chr11:69451425-69451435),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MXI1_disc2|chr11:69451427-69451437),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYCN_2|chr11:69451426-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_disc2|chr11:69451425-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_disc2|chr11:69451426-69451437),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known14|chr11:69451425-69451437),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known15|chr11:69451425-69451435),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known16|chr11:69451425-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known18|chr11:69451426-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known1|chr11:69451424-69451438),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known21|chr11:69451426-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|MYC_known6|chr11:69451426-69451436),TFM(DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2|NFE2_disc3|chr11:69451428-69451438),TFP(E2F1|chr11:69451092-69451819),TFP(E2F1|chr11:69451200-69451666),TFP(E2F1|chr11:69451248-69451703),TFP(EP300|chr11:69451011-69451757),TFP(MAX|chr11:69451008-69451618),TFP(MAX|chr11:69451126-69451679),TFP(MAX|chr11:69451196-69451650),TFP(MYC|chr11:69451060-69451615),TFP(MYC|chr11:69451115-69451654),TFP(MYC|chr11:69451164-69451600),TFP(RFX5|chr11:69450973-69451601),TFP(SMC3|chr11:69450987-69451644),TFP(TCF4|chr11:69451147-69451639),TFP(USF1|chr11:69451180-69451638),TFP(USF1|chr11:69451226-69451611),TFP(USF1|chr11:69451237-69451628),TFP(USF2|chr11:69451160-69451596)'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     }


tableMutations_candidate.FunSeq2_target_gene(tableMutations_candidate.isMOTIFBR_negative & contains(tableMutations_candidate.candidateGenes, 'CCND1'))
%     {'MYEOV(Intron&Promoter)'                         }
%     {'CCND1(Promoter)[DNA_repair][actionable][cancer]'}
%     {'CCND1(Promoter)[DNA_repair][actionable][cancer]'}


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


%%
iRow = 263;
% tableMutations_candidate.FunSeq2_motif_analysis{iRow} = strrep(tableMutations_candidate.FunSeq2_motif_analysis{iRow}, 'CTCF,DHS,EBF1,EGR1,ELF1,FOXA1,FOXA2,TAF1,ZEB1#FOXJ2_1#60988218#60988236#-#12#0.341463#0.658537,', '');
tableMutations_candidate.MOTIFBR{iRow} = 'DHS,E2F1,MAX,MXI1,MYC,NR3C1,RFX5,TCF12,TCF7L2,USF1,USF2#MXI1_disc2#69451425#69451435#+#7#0.000000#1.000000'; %regexp(tableMutations_candidate.FunSeq2_motif_analysis{iRow}, 'MOTIFBR=.+', 'match', 'once');
tmp2 = strsplit(tableMutations_candidate.MOTIFBR{iRow}, {'=', '#'}); % , ','
tableMutations_candidate.MOTIFBR_TFs{iRow} = tmp2{2};
tableMutations_candidate.MOTIFBR_motifName{iRow} = tmp2{3};
tmp3 = regexp(tmp2{3}, '.*_', 'match', 'once');
tableMutations_candidate.MOTIFBR_motifNamePrefix{iRow} = tmp3(1:end-1);
tableMutations_candidate.MOTIFBR_motifStart(iRow) = str2double(tmp2{4});
tableMutations_candidate.MOTIFBR_motifEnd(iRow) = str2double(tmp2{5});
tableMutations_candidate.MOTIFBR_isMinusStrand(iRow) = strcmp(tmp2{6}, '-');
tableMutations_candidate.MOTIFBR_positionMutation(iRow) = str2double(tmp2{7});
tableMutations_candidate.MOTIFBR_scoreAlt(iRow) = str2double(tmp2{8});               % alternative allele frequency in PWM
tmp3 = strsplit(tmp2{9}, ',');
tableMutations_candidate.MOTIFBR_scoreRef(iRow) = str2double(tmp3{1});               % reference allele frequency in PWM
%%
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

try
    fontSizeAnnotation = 14;
    dim = [.07 .55 .05 .05]; str = 'TFBS gain'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
    dim = [.07 .40 .05 .05]; str = 'TFBS gain'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
    dim = [.07 .23 .05 .05]; str = 'TFBS break'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
    dim = [.07 .08 .05 .05]; str = 'TFBS break'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeAnnotation, 'EdgeColor','none', 'Rotation', 90, 'HorizontalAlignment', 'center');
catch
    % In Matlab<2022a, class TextBox does not have property Rotation.
end


fontSizeLetters = 26;
dim = [.005 .99 .01 .01]; str = 'a'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.36 .99 .01 .01]; str = 'b'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.67 .99 .01 .01]; str = 'c'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.005 .66 .01 .01]; str = 'd'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');
dim = [.005 .35 .01 .01]; str = 'e'; annotation('textbox',dim,'String',str, 'FontSize', fontSizeLetters, 'EdgeColor','none', 'FontWeight','bold');

mySaveAs(fig, imagesPath, 'FigCCND1_TFBS', true, true);
savefig([imagesPath, 'FigCCND1_TFBS.fig']);
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