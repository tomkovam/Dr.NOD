#!/bin/bash -l

##################################################
# SGE SUBMISSION SCRIPT run_mutAnnotate.sh       #
##################################################
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -j y
#$ -M marketa.tomkova@ndm.ox.ac.uk
#$ -m as
#$ -l h_rt=480:00:00
#$ -l h_vmem=10G
#$ -pe smp 1

###############################
# DESCRIPTION                 #
###############################
# Takes a bedlike file as input, sorts it and annotates with sequence context and gene names if in exome.
# Rewrites all mutations according to the mutated pyrimidine (+ if the pyrimidine is the reference base; - otherwise).

###############################
# CODE                        #
###############################
# set -o xtrace

INPUT_FILE_BEDLIKE=$1			# TAB-delimited: chr=$1; pos0=$2; pos1=$3; iSample=$4; sample=$5; fromBase=$6; toBase=$7; supportingReadsTumour=$8; coverageTumour=$9;
OUTPUT_FILE_ANNOTATED=$2		# TAB-delimited: chr=$1; pos0=$2; pos1=$3; iSample=$4; sample=$5; fromBase=$6; toBase=$7; supportingReadsTumour=$8; coverageTumour=$9; rawCADD=$10; phredCADD=$11; context=$12; iPattern=$13; strand=$14; geneSymbol=$15;
OUTPUT_FILE_SAMPLES=$3			# TAB-delimited: iSample=$1; sample=$2;
GENOME_REFERENCE_FILE=$4		# Reference genome (fasta file)
INPUT_FILE_EXONS=$5				# TAB-delimited: chr=$1; pos0=$2; pos1=$3; geneSymbol=$4 (coordinates of exones, sorted)
INPUT_FILE_CADD=data/CADD/CADD_whole_genome_SNVs.GRCh37.reformatted2.sorted.txt

# INPUT_FILE_BEDLIKE=data/out01/substitutions.TCGA_CRC_MSI.WGS.bedlike.gz
# OUTPUT_FILE_ANNOTATED=data/out02/substitutions.TCGA_CRC_MSI.WGS.annotated.bedlike.gz
# GENOME_REFERENCE_FILE=/mnt/lustre/users/mtomkova/hg19/hg19_full.fa
# INPUT_FILE_EXONS=/mnt/lustre/users/mtomkova/singer/genes/exons.hg19.sorted.bed

if [ -z "$INPUT_FILE_BEDLIKE" ] || [ -z "$OUTPUT_FILE_ANNOTATED" ] || [ -z "$OUTPUT_FILE_SAMPLES" ] || [ -z "$GENOME_REFERENCE_FILE" ] || [ -z "$INPUT_FILE_EXONS" ]; then
	echo "Not enough parameters." 1>&2
	exit 1
fi
NOW="$(date +'%d/%m/%Y %H:%M:%S')"; SECONDS=0
echo "${NOW}| Script run_mutAnnotate.sh $INPUT_FILE_BEDLIKE $OUTPUT_FILE_ANNOTATED $OUTPUT_FILE_SAMPLES $GENOME_REFERENCE_FILE $INPUT_FILE_EXONS ..."

if [ ! -f "$INPUT_FILE_BEDLIKE" ]; then
	echo "Input file $INPUT_FILE_BEDLIKE is missing." 1>&2; exit 1
fi
sh printFile.sh $INPUT_FILE_BEDLIKE

FILE_INTERMEDIATE0="${OUTPUT_FILE_ANNOTATED}.intermediate.0.bedlike"
FILE_INTERMEDIATE1="${OUTPUT_FILE_ANNOTATED}.intermediate.1.bedlike"
FILE_INTERMEDIATE2="${OUTPUT_FILE_ANNOTATED}.intermediate.2.bedlike"
FILE_INTERMEDIATE3="${OUTPUT_FILE_ANNOTATED}.intermediate.3.bedlike"
FILE_INTERMEDIATE4="${OUTPUT_FILE_ANNOTATED}.withContext.txt"
WARNING_DIFFERENCES="${OUTPUT_FILE_ANNOTATED}.warning.differences.txt"
WARNING_DIFFERENT_REFERENCE="${OUTPUT_FILE_ANNOTATED}.warning.differentReference.txt"
MUTATIONS_UNIQUE="${OUTPUT_FILE_ANNOTATED}.info.uniqueMutationTypes.txt"
[ -f $MUTATIONS_UNIQUE ] && rm $MUTATIONS_UNIQUE
[ -f $FILE_INTERMEDIATE1 ] && rm $FILE_INTERMEDIATE1
[ -f $FILE_INTERMEDIATE2 ] && rm $FILE_INTERMEDIATE2
[ -f $FILE_INTERMEDIATE3 ] && rm $FILE_INTERMEDIATE3
	
runStep0=true
runStep1=true
runStep2=true
runStep3=true
runStep4=true
runStep5=true
runStep6=true
################################################################################
######## STEP0: annotate with CADD, and reformat. ########
################################################################################
if [ "$runStep0" == true ]; then
	gunzip -c $INPUT_FILE_BEDLIKE | awk '{printf "%s_%d_%s_%s\t%s\n", $1, $3, $6, $7, $0}' | sort -k1,1 > tmpA$$.txt 	# create SNV id and sort according to this column (it looks like chr10_100015345_G_A)
	printFile tmpA$$.txt
	join --nocheck-order tmpA$$.txt $INPUT_FILE_CADD -a1 | tr " " "\t" > tmpB$$.txt										# annotate with CADD using this SNV id
	printFile tmpB$$.txt
	cut -f 2- tmpB$$.txt | sort -k1,1 -k2,2n > $FILE_INTERMEDIATE0														# remove the SNV id column and resort as bedtools requires
	printFile $FILE_INTERMEDIATE0
	rm tmpA$$.txt
	rm tmpB$$.txt
fi


################################################################################
######## STEP1: sort input file and annotate with +-1 sequence context. ########
################################################################################
if [ "$runStep1" == true ]; then
	awk '$1 == "chrX" || $1 == "chrY" || $1 ~ /^chr[0-9]+$/' $FILE_INTERMEDIATE0 > $FILE_INTERMEDIATE1
	sh printFile.sh $FILE_INTERMEDIATE1
	awk '{printf "%s\t%d\t%d\n", $1, $2-1, $3+1}' $FILE_INTERMEDIATE1 | \
		bedtools getfasta -fi $GENOME_REFERENCE_FILE -bed stdin -fo stdout -tab 2>&1 | \
		cut -f 2 | paste $FILE_INTERMEDIATE1 - | grep -v "WARNING" > $FILE_INTERMEDIATE2 # Ignore rows outside the genome (lines with WARNING. chromosome (chrXYZ) was not found in the FASTA file. Skipping.)
	sh printFile.sh $FILE_INTERMEDIATE2
	nRows1=`cat $FILE_INTERMEDIATE1 | wc -l`
	nRows2=`cat $FILE_INTERMEDIATE2 | wc -l`
	if [ "$nRows1" != "$nRows2" ]; then 
		echo "WARNING: nRows1=${nRows1}, nRows2=${nRows2} ===============" 1>&2
		cut -f 1-10 $FILE_INTERMEDIATE2 | diff $FILE_INTERMEDIATE1 - > $WARNING_DIFFERENCES #sh printFile.sh $WARNING_DIFFERENCES 1>&2
		echo "CHROMOSOMES ONLY IN INPUT FILE:" 1>&2
		grep "<" $WARNING_DIFFERENCES | cut -f 1 | sort -u
		echo "CHROMOSOMES ONLY IN OUTPUT FILE:" 1>&2
		grep ">" $WARNING_DIFFERENCES | cut -f 1 | sort -u
		echo "END OF WARNING ==========================================" 1>&2
		rm $WARNING_DIFFERENCES
	fi
fi
####################################################################################################################################
######## STEP2: sort input file ####################################################################################################
### First, check that fromBase corresponds to the reference. If not, complement both the fromBase and the toBase. 			########
### If it again does not match the reference, exit with error.																########
### Second, if the fromBase is not pyrimidine (C or T), reverse complement the fromBase, toBase, and the sequence context.	########
####################################################################################################################################
if [ "$runStep2" == true ]; then
	gawk -v WARNING_DIFFERENT_REFERENCE=$WARNING_DIFFERENT_REFERENCE ' # if baseMiddle !~ fromBase --> complement; if fromBase !~ [CT] --> reverseComplement
		function normalComplement(sequence) {
			p = ""
			for(i=1; i <= length(sequence); i++) { p = p m[substr(sequence, i, 1)] }
			return p
		}
		function reverseComplement(sequence) {
			p = ""
			for(i=length(sequence); i > 0; i--) { p = p m[substr(sequence, i, 1)] }
			return p
		}
		BEGIN{
			IGNORECASE = 1; m["A"]="T";m["C"]="G";m["G"]="C";m["T"]="A"; bases[1]="A"; bases[2]="C"; bases[3]="G"; bases[4]="T"; iPattern=1;
			for (fromBase=2; fromBase <= 4; fromBase += 2) { # C, T
				for (toBase=1; toBase <= 4; toBase++) {
				if (fromBase != toBase) {
						for (leftBase=1; leftBase <= 4; leftBase++) {
							for (rightBase=1; rightBase <= 4; rightBase++) {
								arrayPatterns["|"bases[leftBase]""bases[fromBase]""bases[rightBase]"to"bases[toBase]"|"] = iPattern
								iPattern++
							}
						}
					}
				}
			}
		}{ 
		chrNumber=$1; pos0=$2; pos1=$3; iSample=$4; sample=$5; fromBase=$6; toBase=$7; supportingReadsTumour=$8; coverageTumour=$9; rawCADD=$10; phredCADD=$11; context=toupper($12); fromBaseReference=substr(context, 2, 1);
		if (length(fromBase) == 1 && length(toBase) == 1 && fromBaseReference !~ "N" && fromBase !~ "N") {				# only substitutions
			if (fromBaseReference !~ fromBase) {fromBase = normalComplement(fromBase); toBase = normalComplement(toBase)}
			if (fromBaseReference !~ fromBase) {
				printf "ERROR: fromBaseReference=%s, fromBase=%s, NR=%d, line=%s, FILENAME=%s\n", fromBaseReference, fromBase, NR, $0, FILENAME > WARNING_DIFFERENT_REFERENCE
			} else {
				strand = "+";
				if (fromBase !~ /[CT]/) {strand = "-"; context = reverseComplement(context); fromBase = reverseComplement(fromBase); toBase = reverseComplement(toBase); }
				pattern = "|"context"to"toBase"|"
				if (pattern in arrayPatterns) {
					iPattern = arrayPatterns[pattern]
				} else {
					iPattern = 0
				}
				printf "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%f\t%f\t%s\t%d\t%s\n", chrNumber, pos0, pos1, iSample, sample, fromBase, toBase, supportingReadsTumour, coverageTumour, rawCADD, phredCADD, context, iPattern, strand
			}
		}
	}' $FILE_INTERMEDIATE2 > $FILE_INTERMEDIATE3
	sh printFile.sh $FILE_INTERMEDIATE3
	nRows3=`cat $FILE_INTERMEDIATE3 | wc -l`
	if [ "$nRows3" != "$nRows2" ]; then 
		echo "WARNING: nRows3=${nRows3}, nRows2=${nRows2}. ============" 1>&2; 
		sh printFile.sh $WARNING_DIFFERENT_REFERENCE 1>&2; 
		echo "END OF WARNING ==========================================" 1>&2
	else
		cut -f 6,7 $FILE_INTERMEDIATE2 | paste - $FILE_INTERMEDIATE3 | cut -f 1,2,8,9 | tr '[a-z]' '[A-Z]' | sort -u > $MUTATIONS_UNIQUE # with seq context: cut -f 6,7,10 tmp.bedlike | paste - tmp2.bedlike | cut -f 1,2,3,9,10,13 | tr '[a-z]' '[A-Z]' | sort -u
		sh printFullFile.sh $MUTATIONS_UNIQUE
		nRows4=`cat $MUTATIONS_UNIQUE | wc -l`
		if [ "$nRows4" -gt 12 ]; then echo "WARNING: nRows4=${nRows4}."; 1>&2; fi
	fi
fi
########################################################################
######## STEP3: annotate with gene names (only if inside exome) ########
########################################################################
if [ "$runStep3" == true ]; then
	bedtools map -a $FILE_INTERMEDIATE3 -b $INPUT_FILE_EXONS -c 4 -o distinct | gzip > $OUTPUT_FILE_ANNOTATED
	sh printFile.sh $OUTPUT_FILE_ANNOTATED
fi


################################################################################
######## STEP4: annotate with +-8 sequence context. ########
################################################################################
if [ "$runStep4" == true ]; then
	gunzip -c $OUTPUT_FILE_ANNOTATED > $FILE_INTERMEDIATE3
	awk '{printf "%s\t%d\t%d\n", $1, $2-8, $3+8}' $FILE_INTERMEDIATE3 | \
		bedtools getfasta -fi $GENOME_REFERENCE_FILE -bed stdin -fo stdout -tab 2>&1 | \
		cut -f 2 | paste $FILE_INTERMEDIATE3 - | grep -v "WARNING" > $FILE_INTERMEDIATE4 # Ignore rows outside the genome (lines with WARNING. chromosome (chrXYZ) was not found in the FASTA file. Skipping.)
	sh printFile.sh $FILE_INTERMEDIATE4
	nRows1=`cat $FILE_INTERMEDIATE3 | wc -l`
	nRows2=`cat $FILE_INTERMEDIATE4 | wc -l`
	if [ "$nRows1" != "$nRows2" ]; then 
		echo "WARNING: nRows1=${nRows1}, nRows2=${nRows2} ===============" 1>&2
		#cut -f 1-10 $FILE_INTERMEDIATE4 | diff $FILE_INTERMEDIATE3 - >> $WARNING_DIFFERENCES #sh printFile.sh $WARNING_DIFFERENCES 1>&2
		echo "CHROMOSOMES ONLY IN INPUT FILE:" 1>&2
		grep "<" $WARNING_DIFFERENCES | cut -f 1 | sort -u
		echo "CHROMOSOMES ONLY IN OUTPUT FILE:" 1>&2
		grep ">" $WARNING_DIFFERENCES | cut -f 1 | sort -u
		echo "END OF WARNING ==========================================" 1>&2
		rm $WARNING_DIFFERENCES
	fi
fi


##########################################
######## STEP4: print sample list ######## NEWLY, we do this already in the previous (parsing) step
##########################################
if [ "$runStep5" == true ]; then
	gunzip -c $OUTPUT_FILE_ANNOTATED | cut -f 4,5 | sort -u -k1,1n -k2,2 > $OUTPUT_FILE_SAMPLES
	sh printFile.sh $OUTPUT_FILE_SAMPLES
fi

#################################
######## STEP5: clean up ########
#################################
if [ "$runStep6" == true ]; then
	rm $MUTATIONS_UNIQUE
	rm $FILE_INTERMEDIATE0
	rm $FILE_INTERMEDIATE1
	rm $FILE_INTERMEDIATE2
	rm $FILE_INTERMEDIATE3
fi

# set +o xtrace
NOW="$(date +'%d/%m/%Y %H:%M:%S')"
DURATION=$SECONDS
echo "$(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds elapsed."
echo "${NOW}| Script run_mutAnnotate.sh FINISHED."
echo ""

