#!/usr/bin/bash

# Load modules and activate environment

source ~/env/RSeQC/bin/activate

# MAIN

REFBED=$1
INFILE=$2
OUTFILE=$3
OUTSTAT=$4
BAMID=$5

# Determine library type using RSeQC

echo ":: Inferring library type"

infer_experiment.py -r ${REFBED} -s 1000000 -i ${INFILE} > ${OUTFILE}

libtype=$(cat $OUTFILE | grep End | awk '{print $3}')

if [ -z "$libtype" ]
then
    echo "# Error: Lib type not found for $INFILE"
    exit 1
fi

# Determine strandedness (FROM RACHED)

echo ":: Inferring library strandedness"

stranded=$(
    cat $OUTFILE | awk -v FS=': '\
    'BEGIN { isPaired; i = 1; forward_s; reverse_s; } /This is .+End Data|Fraction of reads explained by/ {
    	
    	# Get the fraction of reads supporting each hypothesis
	    if (isPaired == "yes") {

		    # cDNA library is paired end
		    if ($1 == "Fraction of reads explained by \"1++,1--,2+-,2-+\"") {
			    forward_s = $2;
		    } else if ($1 == "Fraction of reads explained by \"1+-,1-+,2++,2--\"") {
			    reverse_s = $2;
	    	}

	    } else if (isPaired == "no") {

	    	# cDNA library is single end
		    if ($1 == "Fraction of reads explained by \"++,--\"") {
		    forward_s = $2;
		    } else if ($1 == "Fraction of reads explained by \"+-,-+\"") {
		    	reverse_s = $2;
			}
	    }                                                             
	    if (i == 1) {
	    	if ($0 == "This is PairEnd Data") { 
	    		isPaired = "yes" } 
	    	else if ($0 == "This is SingleEnd Data") { 
	    		isPaired = "no" 
	    	};
	    	i = 0;
	    }
    } END {
	    if (forward_s > 0.7) {
	    	printf("yes");
	    } else if (reverse_s > 0.7) {
	    	printf("reverse");
	    } else {
	    	printf("no");
	    }
    }'
)

statsline="$BAMID\t$libtype\t$stranded"
echo -e $statsline >> "$OUTSTAT"

deactivate