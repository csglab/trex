#!/usr/bin/bash

# Estimate read length

INFILE=$1
OUTLEN=$2
BAMID=$3

echo ":: Estimating read length"

rlen=$(samtools view ${INFILE} | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq)

if [ -z "$rlen" ]
then
	echo "# Error: could not estimate read length for $INFILE"
	exit 1
else
	echo -e "$BAMID\t$rlen" >> ${OUTLEN}
fi

