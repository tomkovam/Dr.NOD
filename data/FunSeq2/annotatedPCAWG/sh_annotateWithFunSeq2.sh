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

module load apps/bedtools/2.27.0/gcc-4.8.5

NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script annotateWithFunSeq2_oneFile.sh ..."
set -o xtrace

while read PCAWG_PROJECT; do
	sh annotateWithFunSeq2_oneFile_CRC.sh /mnt/lustre/users/mtomkova/singer2/data/out02_annotating/substitutions.$PCAWG_PROJECT.WGS.annotated.bedlike.gz.withContext.txt $PCAWG_PROJECT data/$PCAWG_PROJECT.annotated.FunSeq2.bed
	#qsub -N annotate_$PCAWG_PROJECT annotateWithFunSeq2_oneFile.sh /mnt/lustre/users/mtomkova/singer2/data/out02_annotating/substitutions.$PCAWG_PROJECT.WGS.annotated.bedlike.gz.withContext.txt $PCAWG_PROJECT data/$PCAWG_PROJECT.annotated.FunSeq2.bed
done < data/projectsPCAWG_enhCRC.txt


set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script annotateWithFunSeq2_oneFile.sh FINISHED."
echo ""