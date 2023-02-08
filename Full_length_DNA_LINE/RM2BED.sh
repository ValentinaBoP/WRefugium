#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ./RM2BED.sh -i <input.out> -o <output.bed>
# =========================================================
# Convert RepeatMasker output file into a BED file
# =========================================================
# Valentina Peona                               18 Jan 2019
# =========================================================
# Example:
# ./RM2BED.sh -i lycPyr.out -o lycPyr.bed
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
        *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

awk '{print $5, $6-1, $7, $10"__"$11, $2, $9}' OFS='\t' $INPUT | awk 'NR > 3' | awk '{ sub(/C$/, "-", $6) }1' OFS='\t' > $OUTPUT
