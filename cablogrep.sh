#!/bin/bash
help()
{
	echo "Cablogrep version 1 by Daniel Gerber"
	echo ""
	echo "Tool for inferring horse mitochondrial haplogroups (A-R + S + X) on reference NC_001640"
	echo ""
	echo "Available at Github under https://github.com/ArchGenIn/cablogrep"
	echo ""
	echo "Citation:"
	echo "Dicso et al. (2023) A genetic study on horsekeeping during the Bronze Age Carpathian Basin"
	echo ""
	echo "Options:"
	echo "	-p	path/to/fasta_files (must have .fasta extension)"
	echo "	-o	output name"
	echo "	-h	prints this message"
}
while getopts p:o:h opts
do
	case "${opts}" in
		p) path=$OPTARG;;
		o) oname=$OPTARG;;
		h) help; exit 2;;
		\?) echo "Unknown option: -$OPTARG" >&2; exit 1;;
		:) echo "Missing option for: -$OPTARG" >&2; exit 1;;
	esac
done
if (( $OPTIND == 1 )); then help; fi

cat "$path"/*fasta > "$path"/"$oname".cgall.fasta
nucmer Reference.fasta "$path"/"$oname".cgall.fasta && show-snps "$path"/"$oname".delta > "$path"/"$oname".snps
Rscript --vanilla cab_id.r "$path"/"$oname".snps ref_nodup_rename_SX.csv hg_map_v2.csv
