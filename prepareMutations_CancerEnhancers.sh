#!/usr/bin/env bash
#SBATCH -J test
#SBATCH -o test."%j".out
#SBATCH -e test."%j".err
#SBATCH -n 1
#SBATCH --mem=8000
#SBATCH --time=7-23:00:00
#SBATCH -p production

# prepareMutations_CancerEnhancers.sh

DIR_NASSER=data/Nasser2021/

DIR_PCAWG=data/PCAWG2/
OUTDIR=$DIR_PCAWG/out03_inEnhancersPromoters_250bp/; [ -d $OUTDIR ] || mkdir $OUTDIR # withContext SNVs and basic INDELs that lie within 250bp of any enhancer or 200bp of any TSS

INFILE_PROJECTS=$DIR_PCAWG/PROJECT_CODES.txt

INFILE_ENHANCERS_GZ=data/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz
INFILE_REFERENCE=data/hg19.chrom.sizes

FILE_ENHANCERS_PROMOTERS=data/Nasser2021/MergedAllPredictions.slop250bpAndPromoters.bed 


runPrepareABC=true
runPCAWGMutations=true

if [ $runPrepareABC = true ]; then
	gunzip -c $INFILE_ENHANCERS_GZ | tail -n +2 | cut -f 1-3 | sort -k1,1 -k2,2n | uniq | bedtools merge -i stdin >  data/Nasser2021/MergedAllPredictions.bed
	bedtools slop -i data/Nasser2021/MergedAllPredictions.bed -g $INFILE_REFERENCE -b 250 | bedtools merge -i stdin > data/Nasser2021/MergedAllPredictions.slop250bp.bed	
	gunzip -c $DIR_NASSER/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz | cut -f 1,7,8 | tail -n +2 | sort -u > $DIR_NASSER/AllPredictions.TSS.txt
	awk '{chr=$1; TargetGene=$2; TargetGeneTSS=$3; printf "%s\t%d\t%d\t%s\n", chr, TargetGeneTSS-201, TargetGeneTSS+200, TargetGene}' $DIR_NASSER/AllPredictions.TSS.txt > $DIR_NASSER/AllPredictions.promoters.bed
	cut -f 1-3 $DIR_NASSER/AllPredictions.promoters.bed | cat $DIR_NASSER/MergedAllPredictions.slop250bp.bed - | sort -k1,1 -k2,2n | bedtools merge -i stdin > $FILE_ENHANCERS_PROMOTERS
	sh printBedfile.sh $FILE_ENHANCERS_PROMOTERS # 221616 lines, 355,017,579 bp (355 Mbp)
fi


if [ $runPCAWGMutations = true ]; then
	while read PROJECT OTHER; do
		SEQUENCING_TYPE=WGS
		gunzip -c $DIR_PCAWG/out01_parsing/indels.$PROJECT.$SEQUENCING_TYPE.bedlike.gz | bedtools intersect -a stdin -b $FILE_ENHANCERS_PROMOTERS > $OUTDIR/indels.$PROJECT.bedlike.txt
		bedtools intersect -a $DIR_PCAWG/out02_annotating/substitutions.$PROJECT.$SEQUENCING_TYPE.annotated.bedlike.gz.withContext.txt -b $FILE_ENHANCERS_PROMOTERS > $OUTDIR/substitutions.$PROJECT.annotated.txt
		sh printBedfile.sh $OUTDIR/substitutions.$PROJECT.annotated.txt
	done < $INFILE_PROJECTS 
fi
