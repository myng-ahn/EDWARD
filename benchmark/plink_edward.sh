#!/bin/bash

# Check if an input file name is provided
if [ $# -eq 0 ]; then
    echo "Usage: ./plink_edward.sh <input_file>"
    exit 1
fi

# Assign the input file name to a variable
input_file="$1"

# Call plink command with the input file
echo Plink....
echo start: `date +%m-%d-%Y_%H-%M-%S` 
plink --vcf "$input_file" --pca --out "$input_file"
echo end: `date +%m-%d-%Y_%H-%M-%S`

# Call edward command with the input file
echo EDWARD....
echo start: `date +%m-%d-%Y_%H-%M-%S` 
edward -t v -i "$input_file" --pca -p "$input_file"
echo end: `date +%m-%d-%Y_%H-%M-%S`