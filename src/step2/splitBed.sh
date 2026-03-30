#!/bin/sh

## Usage ./splitBed.sh input.bed

input=$1
outdir=$2

for chr in `cut -f 1 $input | sort | uniq`;
do
	echo $chr
	awk -v chr=$chr '{if ($1 == chr) print $_}' $input > $outdir/$chr.bed
	#awk -v chr='$chr' 'BEGIN {if ($1 == chr) print $_}' $input > $outdir/$chr.bed
done
