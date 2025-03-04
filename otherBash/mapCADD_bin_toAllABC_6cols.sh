#!/bin/bash

#########################################################
# SGE RUN ALL SUBMISSION SCRIPTs mapCADD_bin_toAllABC.sh    #
#########################################################
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -M marketa.tomkova@ndm.ox.ac.uk
#$ -m as
#$ -l h_rt=680:00:00
#$ -l h_vmem=20G
#$ -pe smp 1

module load apps/bedtools/2.27.0/gcc-4.8.5

NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script mapCADD_bin_toAllABC.sh ..."
set -o xtrace

MIN_PHRED=$1
SLOP_TEXT=$2 # "Slop250bpNoncodingPredictions" "NoncodingPredictions"

INFILE_CADD=/mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz
INFILE_ABC=/mnt/lustre/users/mtomkova/singer2/data/Nasser2021/$SLOP_TEXT.together.6cols.txt
OUTFILE_CADD_THRESHOLDED=/mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.geq${MIN_PHRED}.bed.gz
OUTFILE_ABC_MAPPED=/mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions_CADD/$SLOP_TEXT.together.6cols.CADD.geq${MIN_PHRED}.txt

if [ ! -f $OUTFILE_CADD_THRESHOLDED ]; then
	gunzip -c $INFILE_CADD | awk -v MIN_PHRED=$MIN_PHRED '{if ($5>=MIN_PHRED) {$5=1} else {$5=0} printf "%s\t%d\t%d\t%s\t%s\n", $1, $2, $3, $4, $5 }' | gzip -c > $OUTFILE_CADD_THRESHOLDED
	printFile $OUTFILE_CADD_THRESHOLDED
fi

# if [ ! -f $OUTFILE_ABC_MAPPED ]; then
	bedtools map -a $INFILE_ABC -b $OUTFILE_CADD_THRESHOLDED -c 5 -o sum -null NaN > $OUTFILE_ABC_MAPPED
	printFile $OUTFILE_ABC_MAPPED
# fi

# /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/lstABC.txt
# B_cell-ENCODE
# bipolar_neuron_from_iPSC-ENCODE
# body_of_pancreas-ENCODE
# breast_epithelium-ENCODE
# breast-CONSENSUS
# CD19-positive_B_cell-Roadmap
# epithelial_cell_of_prostate-ENCODE
# fibroblast_of_dermis-Roadmap
# hepatocyte-ENCODE
# liver-ENCODE
# mammary_epithelial_cell-Roadmap
# osteoblast-ENCODE
# ovary-Roadmap
# pancreas-Roadmap
# PC-9-ENCODE
# sigmoid_colon-ENCODE
# stomach-Roadmap
# transverse_colon-ENCODE
# uterus-ENCODE


set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script mapCADD_bin_toAllABC.sh FINISHED."
echo ""