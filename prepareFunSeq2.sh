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
#$ -l h_vmem=40G
#$ -pe smp 1

module load apps/bedtools/2.27.0/gcc-4.8.5

NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script annotateWithFunSeq2_oneFile.sh ..."
set -o xtrace

#################################################################
### Prepare FunSeq2 annotation of all potential hg19 variants ### From: http://org.gersteinlab.funseq.s3-website-us-east-1.amazonaws.com/funseq2.1.2/hg19_NCscore_funseq216.tsv.bgz (website http://funseq2.gersteinlab.org/downloads)
#################################################################
[ -d tmp ] || mkdir tmp
gunzip -c data/hg19_score.funseq216.bed.gz | sort -k1,1 -k2,2n -k5,5 -T tmp -S 30G | gzip -c > data/hg19_score.funseq216.sorted.bed.gz
gunzip -c data/hg19_score.funseq216.sorted.bed.gz | tr ';' '\t' | cut -f 1-5,7,8,14,21 | gzip -c > data/hg19_score.funseq216.sorted.selectedColumns.bed.gz # chr, pos0, pos1, ref, alt, gerp, cds (Yes/No), motif, noncodingScore

####################################################################################################################
### Annotate the solid cancer candidate regulatory driver variants with FunSeq2 predictions of motif alterations ### 
####################################################################################################################

GENOME_REFERENCE_FILE=../../../hg19/new/hg19.fa
INFILE_CANDIDATES=data/driverVariants_allNonBloodIncluded.txt
OUTFILE_CANDIDATES_0=data/driverVariants_allNonBloodIncluded.sorted.bed
OUTFILE_CANDIDATES_1=data/driverVariants_allNonBloodIncluded.intersectedFunSeq2.txt
OUTFILE_CANDIDATES_2=data/driverVariants_allNonBloodIncluded.intersectedFunSeq2.context50bp.txt

sort -k1,1 -k2,2n $INFILE_CANDIDATES > $OUTFILE_CANDIDATES_0
printFile $OUTFILE_CANDIDATES_0
bedtools intersect -a $OUTFILE_CANDIDATES_0 -b data/hg19_score.funseq216.sorted.bed.gz -sorted -wa -wb > $OUTFILE_CANDIDATES_1
printFile $OUTFILE_CANDIDATES_1


awk '{printf "%s\t%d\t%d\n", $1, $2-50, $3+50}' $OUTFILE_CANDIDATES_1 | \
	bedtools getfasta -fi $GENOME_REFERENCE_FILE -bed stdin -fo stdout -tab 2>&1 | \
	cut -f 2 | paste $OUTFILE_CANDIDATES_1 - | grep -v "WARNING" > $OUTFILE_CANDIDATES_2 # Ignore rows outside the genome (lines with WARNING. chromosome (chrXYZ) was not found in the FASTA file. Skipping.)
printFile.sh $OUTFILE_CANDIDATES_2
nRows1=`cat $OUTFILE_CANDIDATES_1 | wc -l`
nRows2=`cat $OUTFILE_CANDIDATES_2 | wc -l`
if [ "$nRows1" != "$nRows2" ]; then 
	echo "WARNING: nRows1=${nRows1}, nRows2=${nRows2} ===============" 1>&2
	WARNING_DIFFERENCES=$OUTFILE_CANDIDATES_2.warning.txt
	cut -f 1-10 $OUTFILE_CANDIDATES_2 | diff $OUTFILE_CANDIDATES_1 - >> $WARNING_DIFFERENCES 
	echo "CHROMOSOMES ONLY IN INPUT FILE:" 1>&2
	grep "<" $WARNING_DIFFERENCES | cut -f 1 | sort -u
	echo "CHROMOSOMES ONLY IN OUTPUT FILE:" 1>&2
	grep ">" $WARNING_DIFFERENCES | cut -f 1 | sort -u
	echo "END OF WARNING ==========================================" 1>&2
	rm $WARNING_DIFFERENCES
fi


set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script annotateWithFunSeq2_oneFile.sh FINISHED."
echo ""
