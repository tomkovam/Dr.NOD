#!/usr/bin/env bash
# computeTrinucleotides.sh $INPUT_BEDFILE $INPUT_REFERENCE_FASTA

INPUT_BEDFILE=$1
INPUT_REFERENCE_FASTA=$2 # /share/hormozdiarilab/Data/ReferenceGenomes/Hg19/hg19.ref.fa

if [ -z "$INPUT_BEDFILE" ] || [ -z "$INPUT_REFERENCE_FASTA" ]; then
	echo "Not enough parameters $0." 1>&2
	exit 1
fi

awk '{printf "%s\t%d\t%d\n", $1, $2-1, $3+1}' $INPUT_BEDFILE | bedtools getfasta -fi $INPUT_REFERENCE_FASTA -bed stdin -fo stdout -tab 2>&1 | cut -f 2 | gawk '
	function reverseComplement(sequence) {
		p = ""
		for(i=length(sequence); i > 0; i--) { p = p m[substr(sequence, i, 1)] }
		return p
	}
	BEGIN{
		m["A"]="T";m["C"]="G";m["G"]="C";m["T"]="A"; 
		bases[1]="A"; bases[2]="C"; bases[3]="G"; bases[4]="T"; 
		for (fromBase=1; fromBase <= 4; fromBase++) {						# First, we initialise arrayTrinucleotides with zeros
			for (leftBase=1; leftBase <= 4; leftBase++) {
				for (rightBase=1; rightBase <= 4; rightBase++) {
					trinucleotide = bases[leftBase]""bases[fromBase]""bases[rightBase]
					arrayTrinucleotides[trinucleotide] = 0
				}
			}
		}
		for (fromBase=2; fromBase <= 4; fromBase += 2) { # 2=C, 4=T			# Next we print the header (the names of the trinucleotides)
			for (leftBase=1; leftBase <= 4; leftBase++) {
				for (rightBase=1; rightBase <= 4; rightBase++) {
					trinucleotide = bases[leftBase]""bases[fromBase]""bases[rightBase]
					printf "\t%s", trinucleotide
				}
			}
		}
		printf "\n"
	}{ 
		sequence = toupper($1); nPositions = length(sequence);				# Next we read the input file - it has one column only, which represents the sequence of the given region - and we compute trinucleotides in this sequence
		for (iPosition = 2; iPosition < nPositions; iPosition++) { 			# We ignore the first and last positions (they are here just to provide context from the second and penultimate positions)
			trinucleotide = substr(sequence, iPosition-1, 3)
			fromBase = substr(trinucleotide, 2, 1)
			if (fromBase !~ /[CT]/) { 
				trinucleotide = reverseComplement(trinucleotide) 
			}
			arrayTrinucleotides[trinucleotide]++
			#printf "%s %d: %s %d\n", sequence, iPosition, trinucleotide, arrayTrinucleotides[trinucleotide]
		}
		for (fromBase=2; fromBase <= 4; fromBase += 2) { # 2=C, 4=T			# We print the final values of all 2*4*4=32 trinucleotides
			for (leftBase=1; leftBase <= 4; leftBase++) {
				for (rightBase=1; rightBase <= 4; rightBase++) {
					trinucleotide = bases[leftBase]""bases[fromBase]""bases[rightBase]
					printf "\t%d", arrayTrinucleotides[trinucleotide]
					arrayTrinucleotides[trinucleotide] = 0
				}
			}
		}
		printf "\n"
	}' #| paste $INPUT_BEDFILE - -d ''
