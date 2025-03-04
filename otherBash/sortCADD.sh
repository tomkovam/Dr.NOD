#!/bin/bash

#######################################################
# SGE RUN ALL SUBMISSION SCRIPTs sortCADD.sh      #
#######################################################
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -M marketa.tomkova@ndm.ox.ac.uk
#$ -m as
#$ -l h_rt=680:00:00
#$ -l h_vmem=20G
#$ -pe smp 1

NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script sortCADD.sh ..."
set -o xtrace

indexName=$1 # such as 00 (up to 29)

sort -k1,1 -T /mnt/lustre/users/mtomkova/singer2/tmp/ -S 50% data/CADD/splitUnsorted/CADD_whole_genome_SNVs.GRCh37.reformatted.${indexName} > data/CADD/splitSorted/CADD_whole_genome_SNVs.GRCh37.reformatted.${indexName}

set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script sortCADD.sh FINISHED."
echo ""