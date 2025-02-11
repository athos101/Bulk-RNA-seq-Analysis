#!/bin/bash

if [ -z "$1" ]; then
  echo "Error: input file needed!"
  echo "Use: $0 file_with_sra_code.txt"
  exit 1
fi

input_file="$1"
output_file="ncbi_fastq_$(basename "$input_file" .txt).fastq"
> "$output_file"
while read -r srr_code; do
  echo "Downloading $srr_code..."
  prefetch "$srr_code"
  echo "Converting $srr_code to FASTQ..."
  fastq-dump --split-files --gzip "$srr_code" | cat >> "$output_file"
  rm -rf "$srr_code"

done < "$input_file"

echo "Fastq files concatenated in $output_file."