#!/bin/bash

# Check if an input file name is provided
if [ $# -ne 3 ]; then
    echo "Usage: ./plink_edward.sh <input_file> <num_pcs> <out_prefix>"
    exit 1
fi

# Assign the input file name to a variable
input_file="$1"
num_pcs="$2"
prefix="$3"

# Call plink command with the input file
echo Plink....
echo `time ./plink --vcf "$input_file" --pca "$num_pcs" --out "$prefix"`

# Call edward command with the input file
echo EDWARD....
echo `time edward -t v -i "$input_file" --pca -p "$prefix" -n "$num_pcs"`
