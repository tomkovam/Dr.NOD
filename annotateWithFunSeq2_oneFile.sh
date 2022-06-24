#!/bin/bash

#######################################################################
# SGE RUN ALL SUBMISSION SCRIPTs annotateWithFunSeq2_oneFile.sh       #
#######################################################################
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -M marketa.tomkova@ndm.ox.ac.uk
#$ -m as
#$ -l h_rt=680:00:00
#$ -l h_vmem=20G
#$ -pe smp 1

# Example usage:
# qsub annotateWithFunSeq2_oneFile.sh /mnt/lustre/users/mtomkova/singer2/data/out02_annotating/substitutions.ICGC_WXS_PBCA_DE.WXS.annotated.bedlike.gz.withContext.txt ICGC_WXS_PBCA_DE data/ICGC_WXS_PBCA_DE.annotated.FunSeq2.bed
# qsub annotateWithFunSeq2_oneFile.sh /mnt/lustre/users/mtomkova/singer2/data/out02_annotating/substitutions.PCAWG_GBM_US.WGS.annotated.bedlike.gz.withContext.txt PCAWG_GBM_US data/PCAWG_GBM_US.annotated.FunSeq2.bed

module load apps/bedtools/2.27.0/gcc-4.8.5

NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script annotateWithFunSeq2_oneFile.sh ..."
set -o xtrace


INBEDFILE=$1		# /mnt/lustre/users/mtomkova/singer2/data/out02_annotating/substitutions.ICGC_WXS_PBCA_DE.WXS.annotated.bedlike.gz.withContext.txt
PROJECT_NAME=$2 	# ICGC_WXS_PBCA_DE

# which perl

DIR_DATA=data2
INTERMEDIATE_DIR=$DIR_DATA/$PROJECT_NAME

[ -d $INTERMEDIATE_DIR ] || mkdir -p $INTERMEDIATE_DIR

INTERMEDIATE_BED=$DIR_DATA/$PROJECT_NAME.bed
OUTFILE1=$DIR_DATA/$PROJECT_NAME.annotated.FunSeq2.bed
OUTFILE2=$DIR_DATA/$PROJECT_NAME.motif.FunSeq2.bed
OUTFILE3=$DIR_DATA/intersected.$PROJECT_NAME.motif.FunSeq2.bed.txt
OUTFILE4=$DIR_DATA/intersectedB.$PROJECT_NAME.motif.FunSeq2.bed.txt

if [ ! -f $INTERMEDIATE_BED ]; then
	cut -f 1-3,6,7,14 $INBEDFILE | awk '{ m["A"]="T";m["C"]="G";m["G"]="C";m["T"]="A";
		if ($6 == "-") {
			$4 = m[$4]
			$5 = m[$5]
		}
		printf "%s\t%d\t%d\t%s\t%s\n", $1, $2, $3, $4, $5
	}' > $INTERMEDIATE_BED
	printFile $INTERMEDIATE_BED
fi

./funseq2.sh -f $INTERMEDIATE_BED -inf BED -nc -outf BED -o $INTERMEDIATE_DIR

cp $INTERMEDIATE_DIR/Output.BED $OUTFILE1
printFile $OUTFILE1

cat $OUTFILE1 | tr ';' '\t' | cut -f 1-5,14,20 > $OUTFILE2
printFile $OUTFILE2

# bedtools intersect -a data/hg19_score.funseq216.sorted.selectedColumns.bed.gz -b $INTERMEDIATE_BED -sorted > $OUTFILE3

# bedtools intersect -sorted -a $INTERMEDIATE_BED -b data/hg19_score.funseq216.sorted.selectedColumns.bed.gz -wa -wb > $OUTFILE4

set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script annotateWithFunSeq2_oneFile.sh FINISHED."
echo ""