#!/usr/bin/env bash
#SBATCH -J test
#SBATCH -o test."%j".out
#SBATCH -e test."%j".err
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH --time=7-23:00:00
#SBATCH -p production

# prepareAnnotatedEnhancers_step2 after step1

DIR_NASSER=data/Nasser2021/
INFILE=$DIR_NASSER/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz

SORTED_PROTEIN_CODING_EXOME=data/genes/GENCODE/mergedSorted.protein_coding.exons.hg19v19.bed
SORTED_PROTEIN_CODING_CDS=data/genes/GENCODE/mergedSorted.protein_coding.CDS.hg19v19.bed
INFILE_CHR_SIZES=data/hg19.chrom.sizes
INFILE_PCAWG_PROJECTS=data/PCAWG2/PROJECT_CODES.txt
INDIR_PCAWG_SNVs=data/PCAWG2/out02_annotating

SLOP_TEXT="Slop250bpNoncodingPredictions"

INDIR_ANNOTATIONS=data/AnnotatedEnhancers
OUTDIR_ANNOTATIONS2=data/AnnotatedEnhancers/AnnotatedBySNVs; [ -d $OUTDIR_ANNOTATIONS2 ] || mkdir -p $OUTDIR_ANNOTATIONS2

INFILE=$INDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.txt # from prepareAnnotatedEnhancers_step1
SLOPPED_FILE=$INDIR_ANNOTATIONS/$SLOP_TEXT.together.6cols.slop50kbp.bed

# GOAL: normalize for the mutation frequency in the surrounding region. How this normalization could work... 
# We extent each unique enhancer by +- 50 kbp (100 kbp in total) and then subtract all coding space and all enhancer space. 
# Then we use bedtools map to count number of different samples with a mutation in the given enhancer. 
# Then we use this as one of the features in our model. 

bedtools slop -i $INFILE -g $INFILE_CHR_SIZES -b 50000 | bedtools subtract -a stdin -b $SORTED_PROTEIN_CODING_CDS | bedtools subtract -a stdin -b $INFILE | sort -k1,1 -k2,2n -k4,4 -k6,6 | uniq > $SLOPPED_FILE
sh printFile.sh $SLOPPED_FILE 

OUTFILE=$OUTDIR_ANNOTATIONS2/$SLOP_TEXT.together.6cols.slop50kbp.SNVs.txt
printf "chr\tpos0\tpos1\tname\tactivity\tbiosample\n" > $OUTFILE
cat $SLOPPED_FILE >> $OUTFILE
sh printFile.sh $OUTFILE
while read PROJECT_CODE OTHER; do
	INFILE_SNVs=$INDIR_PCAWG_SNVs/substitutions.$PROJECT_CODE.WGS.annotated.bedlike.gz.withContext.txt
	TMP_OUTFILE_COUNTED_SNVs=$OUTDIR_ANNOTATIONS2/$SLOP_TEXT.together.6cols.slop50kbp.SNVs.$PROJECT_CODE.txt
	printf "$PROJECT_CODE\n" > $TMP_OUTFILE_COUNTED_SNVs
	echo "bedtools map -a $SLOPPED_FILE -b $INFILE_SNVs -c 4 -o count_distinct -null 0 >> $TMP_OUTFILE_COUNTED_SNVs ..."
	bedtools map -a $SLOPPED_FILE -b $INFILE_SNVs -c 4 -o count_distinct -null 0 | cut -f 7 >> $TMP_OUTFILE_COUNTED_SNVs
	paste $OUTFILE $TMP_OUTFILE_COUNTED_SNVs > $OUTFILE$$
	mv $OUTFILE$$ $OUTFILE
	sh printFile.sh $OUTFILE
	rm $TMP_OUTFILE_COUNTED_SNVs
done < $INFILE_PCAWG_PROJECTS

uniq $OUTDIR_ANNOTATIONS2/$SLOP_TEXT.together.6cols.slop50kbp.SNVs.txt > $OUTDIR_ANNOTATIONS2/$SLOP_TEXT.together.6cols.slop50kbp.SNVs.uniq.txt

awk -v OUT_PREFIX="$OUTDIR_ANNOTATIONS2/$SLOP_TEXT" '{
	chromosome=$1; pos0=$2; pos1=$3; name=$4; baseActivity=$5; biosample=$6;
	print $0 > OUT_PREFIX"."biosample".slop50kbp.SNVs.uniq.txt"
}' $OUTDIR_ANNOTATIONS2/$SLOP_TEXT.together.6cols.slop50kbp.SNVs.uniq.txt
