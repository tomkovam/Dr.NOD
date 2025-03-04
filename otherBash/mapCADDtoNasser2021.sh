#!/bin/bash

#########################################################
# SGE RUN ALL SUBMISSION SCRIPTs mapCADDtoNasser2021.sh #
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
echo "${NOW}| Script mapCADDtoNasser2021.sh ..."
# set -o xtrace

# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz | cut -f 1-4 | sort -k1,1 -k2,2n | uniq > /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.bed

# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.tsv.gz | awk '{printf "chr%s\t%d\t%d\t%s\t%s\n", $1, $2-1, $2, $4, $6}' | tail -n +2 | gzip -c > /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz 

# 1       13841   T       A       -0.243645       0.536
# 1       13841   T       C       -0.232541       0.578
# 1       13843   T       C       0.065279        3.680
# 1       13843   T       G       0.033086        3.145
# 1       13844   G       A       0.124706        4.669
# 1       13844   G       C       0.069969        3.758
# 1       13844   G       T       0.071092        3.777
# 1       13845   C       A       0.018412        2.908
# 1       13845   C       G       0.013048        2.824
# 1       13845   C       T       0.087330        4.050

# chr1    10000   10001   A       4.066
# chr1    10000   10001   C       4.353

# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz | head
# chr##   -1      0       (c)     of
# chr#Chrom       -1      0       Alt     PHRED
# chr1    10000   10001   A       4.066
# chr1    10000   10001   C       4.353
# chr1    10000   10001   G       3.939
# chr1    10001   10002   C       3.934

# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.withHeader.bed.gz | tail -n +3 | gzip -c > /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz 

# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz | head

# bedtools map -a tmp2.bed -b /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz -c 5 -o mean

# bedtools map -a /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.bed -b /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz -c 5 -o mean > /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.CADD.mean.bed

# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.bed.gz | awk '{if ($5>=10) {$5=1} else {$5=0} printf "%s\t%d\t%d\t%s\t%s\n", $1, $2, $3, $4, $5 }' | gzip -c > /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.geq10.bed.gz

# bedtools map -a /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.bed -b /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.geq10.bed.gz -c 5 -o count > /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.CADD.geq10count.bed

# bedtools map -a /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.bed -b /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.geq10.bed.gz -c 5 -o sum > /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/AllPredictions.4cols.CADD.geq10sum.bed


# gunzip -c /mnt/lustre/users/mtomkova/singer2/data/CADD/CADD_whole_genome_SNVs.GRCh37.geq10.bed.gz | head -n 10000 > tmp.txt
# sort -k1,1 -k2,2n /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions/Slop250bpNoncodingPredictions.B_cell-ENCODE.selCols.txt | bedtools map -a stdin -b tmp.txt -c 5 -o sum -null NaN | head


#####################################################################
# while read ABC_NAME; do
	# qsub -N mapOne mapOneABC_toCADD.sh "$ABC_NAME"
# done < /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/lstABC_v2.txt

############### Create a file with all ABC biosamples ###############
# OUTFILE=/mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions.together.5cols.txt
# echo -n "" > $OUTFILE$$

# while read ABC_NAME; do
	# cut -f 1-4 "/mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions/Slop250bpNoncodingPredictions.${ABC_NAME}.selCols.txt" | sed -e "s/$/\t${ABC_NAME}/" >> $OUTFILE$$
# done < /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/lstABC.txt

# sort -k 1,1 -k2,2n $OUTFILE$$ | uniq > $OUTFILE
# rm $OUTFILE$$

# printFile $OUTFILE
#####################################################################


for MIN_PHRED in 1 3 5 7 9 11 13 15 17 19 21; do #0 2 6 10 14 18 22 26 30 34 38; do # 4 8 12 16 20
	#qsub -N mapOne_$MIN_PHRED mapCADD_bin_toAllABC.sh $MIN_PHRED "Slop250bpNoncodingPredictions"
	qsub -N mapOne_$MIN_PHRED mapCADD_bin_toAllABC_6cols.sh $MIN_PHRED "Slop250bpNoncodingPredictions" # "Slop250bpNoncodingPredictions" "NoncodingPredictions" Slop250bpAllPredictions
done 




# "/mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions/Slop250bpNoncodingPredictions.${ABC_NAME}.selCols.txt"
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


# OCI-LY7
# GM12878
# Karpas-422
# H1_Derived_Neuronal_Progenitor_Cultured_Cells-Roadmap
# MCF-7
# MDA-MB-231
# HCT116
# HT29
# LoVo
# fibroblast_of_lung-Roadmap
# A549
# Panc1
# keratinocyte-Roadmap
# HepG2
# stomach_fetal-Roadmap


# cat /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions/Slop250bpNoncodingPredictions.* | cut -f 1-4 | wc -l
# cat /mnt/lustre/users/mtomkova/singer2/data/Nasser2021/Slop250bpNoncodingPredictions/Slop250bpNoncodingPredictions.* | cut -f 1-4 | sort | uniq | wc -l

# set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script mapCADDtoNasser2021.sh FINISHED."
echo ""