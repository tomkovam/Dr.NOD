#!/usr/bin/env bash
#SBATCH -J test
#SBATCH -o test."%j".out
#SBATCH -e test."%j".err
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH --time=7-23:00:00
#SBATCH -p production

# prepareAnnotatedEnhancers_step1

DIR_NASSER=data/Nasser2021/
INFILE=$DIR_NASSER/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz

SORTED_PROTEIN_CODING_EXOME=data/genes/GENCODE/mergedSorted.protein_coding.exons.hg19v19.bed
SORTED_PROTEIN_CODING_CDS=data/genes/GENCODE/mergedSorted.protein_coding.CDS.hg19v19.bed
INFILE_REPLICATION_TIMING=data/AnnotatedEnhancers/per_base_territories_20kb_line_numbers.bed.txt
INFILE_REFERENCE_FASTA=/share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref.fa
INFILE_CHR_SIZES=data/hg19.chrom.sizes
INFILE_BLACKLISTED_REGIONS=data/blacklisted/blacklisted.bed.gz

FILE_ALL_PREDICTIONS=$DIR_NASSER/AllPredictions.together.selCols.sorted.txt
OUTDIR_ANNOTATIONS=data/AnnotatedEnhancers2; [ -d $OUTDIR_ANNOTATIONS ] || mkdir -p $OUTDIR_ANNOTATIONS

if [ ! -f $FILE_ALL_PREDICTIONS ]; then
	echo -n "#" > $FILE_ALL_PREDICTIONS # This is to start the header line with a dash
	gunzip -c $INFILE | cut -f 1-8,11-13,19-21,24 | sort -k1,1 -k2,2n >> $FILE_ALL_PREDICTIONS
	sh printFile.sh $FILE_ALL_PREDICTIONS
fi

SLOP_TEXT="Slop250bpNoncodingPredictions"

FILE_SLOPPED=$OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.selCols.txt

if [ $SLOP_TEXT = "Slop250bpNoncodingPredictions" ]; then
	bedtools slop -i $FILE_ALL_PREDICTIONS -g $INFILE_CHR_SIZES -b 250 | bedtools subtract -a stdin -b $SORTED_PROTEIN_CODING_CDS > $FILE_SLOPPED
elif [ $SLOP_TEXT = "NoncodingPredictions" ]; then
	bedtools subtract -a $FILE_ALL_PREDICTIONS -b $SORTED_PROTEIN_CODING_CDS > $FILE_SLOPPED
elif [ $SLOP_TEXT = "Slop250bpAllPredictions" ]; then
	bedtools slop -i $FILE_ALL_PREDICTIONS -g $INFILE_CHR_SIZES -b 250 > $FILE_SLOPPED
fi
sh printFile.sh $FILE_SLOPPED

cut -f 1-4,6,15 $FILE_SLOPPED | sort -k 1,1 -k2,2n | uniq > $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.txt	# {'chr', 'pos0', 'pos1', 'name', 'baseActivity', 'biosample'} 
sh printFile.sh $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.txt

bedtools map -a $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.txt -b $INFILE_REPLICATION_TIMING -c 9 -o mean -null NaN > $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.txt	
sh printFile.sh $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.txt

NCOLS=`head -n 1 $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.txt | awk '{print NF}'`
COL_GC=$(($NCOLS + 2))
bedtools nuc -fi $INFILE_REFERENCE_FASTA -bed $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.txt | cut -f 1-$NCOLS,$COL_GC | tail -n +2 > $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.txt
sh printFile.sh $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.txt 					# {'chr', 'pos0', 'pos1', 'name', 'baseActivity', 'biosample', 'replicationTiming', 'GC', 'blacklisted'}


############ BLACKLISTED-START ############ preparation of $INFILE_BLACKLISTED_REGIONS
# 551  mkdir blacklisted
# 552  cd blacklisted/
# 553  wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/Anshul_Hg19UltraHighSignalArtifactRegions.bed.gz
# 554  wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/Duke_Hg19SignalRepeatArtifactRegions.bed.gz
# 555  wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz
# gunzip -c data/blacklisted/Duke_Hg19SignalRepeatArtifactRegions.bed.gz data/blacklisted/Duke_Hg19SignalRepeatArtifactRegions.bed.gz data/blacklisted/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz | cut -f 1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin | gzip -c > $INFILE_BLACKLISTED_REGIONS
# The following file is almost identical to our wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz (two regions are defined slightly differently)
# wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz # from https://www.encodeproject.org/annotations/ENCSR636HFF/  used in https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008749
# cut -f 1-3 $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.txt | head > tmp.bed
# awk '{printf "%s\t%d\t%d\n", $1, $2-1, $3+1}' tmp.bed | bedtools getfasta -fi $INFILE_REFERENCE_FASTA -bed stdin -fo stdout -tab 2>&1 | cut -f 2 | paste tmp.bed - | grep -v "WARNING" > tmp2.txt


bedtools map -a $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.txt -b $INFILE_BLACKLISTED_REGIONS -c 3 -o count -null 0 > $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.blacklisted.txt	
sh printFile.sh $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.blacklisted.txt


############ TRINUCLEOTIDES ############
sh computeTrinucleotides.sh $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.replicationTiming.gc.blacklisted.txt $INFILE_REFERENCE_FASTA > $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.trinucleotides.txt
sh printFile.sh $OUTDIR_ANNOTATIONS/$SLOP_TEXT.together.trinucleotides.txt

################################################################################
# 1	chr
# 2	start
# 3	end
# 4	name
# 5	class
# 6	activity_base
# 7	TargetGene
# 8	TargetGeneTSS
# 9	TargetGeneExpression
# 10	TargetGenePromoterActivityQuantile
# 11	TargetGeneIsExpressed
# 12	distance
# 13	isSelfPromoter
# 14	hic_contact
# 15	powerlaw_contact
# 16	powerlaw_contact_reference
# 17	hic_contact_pl_scaled
# 18	hic_pseudocount
# 19	hic_contact_pl_scaled_adj
# 20	ABC.Score.Numerator
# 21	ABC.Score
# 22	powerlaw.Score.Numerator
# 23	powerlaw.Score
# 24	CellType

# cut -f 1-3 data/AnnotatedEnhancers2/Slop250bpAllPredictions.together.selCols.txt | bedtools merge -i stdin > data/AnnotatedEnhancers2/Slop250bpAllPredictions.together.merged.bed
