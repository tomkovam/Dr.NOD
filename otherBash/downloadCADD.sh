#!/bin/bash

#######################################################
# SGE RUN ALL SUBMISSION SCRIPTs downloadCADD.sh      #
#######################################################
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -M marketa.tomkova@ndm.ox.ac.uk
#$ -m as
#$ -l h_rt=680:00:00
#$ -l h_vmem=30G
#$ -pe smp 1

NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script downloadCADD.sh ..."
set -o xtrace

######## EXAMPLES ########
# CADD_whole_genome_SNVs.GRCh37.tsv.gz
## CADD GRCh37-v1.6 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2019. All rights reserved.
#Chrom  Pos     Ref     Alt     RawScore        PHRED
# 1       10001   T       A       0.088260        4.066
# 1       10001   T       C       0.105475        4.353


# CADD_whole_genome_indels.GRCh37.tsv.gz
## CADD GRCh37-v1.6 (c) University of Washington, Hudson-Alpha Institute for Biotechnology and Berlin Institute of Health 2013-2020. All rights reserved.
#Chrom  Pos     Ref     Alt     RawScore        PHRED
# 1       10001   T       TC      -0.115832       1.266
# 1       10009   A       AC      -0.115586       1.268


######## DATA ########
DIR_CADD=data/CADD; [ -d $DIR_CADD ] || mkdir -p $DIR_CADD
DIR_CADD_SPLIT=data/CADD/splitUnsorted; [ -d $DIR_CADD_SPLIT ] || mkdir -p $DIR_CADD_SPLIT
DIR_CADD_SPLIT_SORTED=data/CADD/splitSorted; [ -d $DIR_CADD_SPLIT_SORTED ] || mkdir -p $DIR_CADD_SPLIT_SORTED
# cd $DIR_CADD

######## DOWNLOADING ########
# wget -O "CADD_whole_genome_SNVs.GRCh37.tsv.gz" -q  https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz

# wget -O "CADD_whole_genome_SNVs.GRCh37.tsv.gz.tbi" https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi

# wget -O "CADD_whole_genome_indels.GRCh37.tsv.gz" https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz

# wget -O "CADD_whole_genome_indels.GRCh37.tsv.gz.tbi" https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi

# gunzip -c CADD_whole_genome_SNVs.GRCh37.tsv.gz | awk '{printf "chr%s_%d_%s_%s\t%f\t%f\n", $1, $2, $3, $4, $5, $6}' > CADD_whole_genome_SNVs.GRCh37.reformatted.txt

# cd ../../


######## SPLITTING ########
# split --verbose -d -a2 -l 290000000 data/CADD/CADD_whole_genome_SNVs.GRCh37.reformatted.txt data/CADD/splitUnsorted/CADD_whole_genome_SNVs.GRCh37.reformatted.


######## SORTING ########
# for indexName in {00..29}; do 
	# qsub -N s${indexName} sortCADD.sh $indexName
# done


######## MERGING ########
sort -m -k1,1 -T /mnt/lustre/users/mtomkova/singer2/tmp/ data/CADD/splitSorted/CADD_whole_genome_SNVs.GRCh37.reformatted.{00..29} > data/CADD/CADD_whole_genome_SNVs.GRCh37.reformatted2.sorted.txt
head data/CADD/CADD_whole_genome_SNVs.GRCh37.reformatted2.sorted.txt

######## ANNOTATING ########
# awk '{printf "%s_%d_%s_%s\t%s\n", $1, $3, $6, $7, $0}' /mnt/lustre/users/mtomkova/singer2/data/out01_parsing/substitutions.ICGC_WXS_BLCA_CN.WXS.bedlike | sort -k1,1 > tmp.txt
# join tmp.txt data/CADD/splitSorted/CADD_whole_genome_SNVs.GRCh37.reformatted.00 -a1 | tr " " "\t" | cut -f 2- | sort -k1,1 -k2,2n > tmp3.txt
# join --nocheck-order tmp.txt data/CADD/CADD_whole_genome_SNVs.GRCh37.reformatted.sorted.txt -a1 | tr " " "\t" | cut -f 2- | sort -k1,1 -k2,2n > tmp3full.txt

# join --nocheck-order tmp.txt data/CADD/CADD_whole_genome_SNVs.GRCh37.reformatted.sorted.txt -a1 | tr " " "\t" > tmp4full.unsorted.txt

set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script downloadCADD.sh FINISHED."
echo ""


######## SORTING OLD ########

# sort -T /mnt/lustre/users/mtomkova/singer2/tmp/ -S 50% --parallel=4 CADD_whole_genome_SNVs.GRCh37.reformatted.txt > CADD_whole_genome_SNVs.GRCh37.reformatted.sorted.txt
# sort -T /mnt/lustre/users/mtomkova/singer2/tmp/ -S 50% CADD_whole_genome_SNVs.GRCh37.reformatted.txt > CADD_whole_genome_SNVs.GRCh37.reformatted.sortedWithoutP.txt

# gunzip -c CADD_whole_genome_indels.GRCh37.tsv.gz | awk '{printf "chr%s_%d_%s_%s\t%f\t%f\n", $1, $2, $3, $4, $5, $6}' | sort | gzip > CADD_whole_genome_indels.GRCh37.sorted.bed.gz

# gunzip -c data/CADD/CADD_whole_genome_SNVs.GRCh37.tsv.gz | tail -n +3 | awk '{printf "chr%s_%d_%s_%s\t%f\t%f\n", $1, $2, $3, $4, $5, $6}' | head -n 10000000 | sort > tmp2.txt
# awk '{printf "%s_%d_%s_%s\n", $1, $3, $6, $7}' /mnt/lustre/users/mtomkova/singer2/data/out01_parsing/substitutions.ICGC_WXS_BLCA_CN.WXS.bedlike | sort > tmp.txt
# join tmp.txt tmp2.txt -a1 > tmp3.txt

# split --verbose -d -a2 -n30 CADD_whole_genome_SNVs.GRCh37.reformatted.txt splitUnsorted/CADD_whole_genome_SNVs.GRCh37.reformatted.