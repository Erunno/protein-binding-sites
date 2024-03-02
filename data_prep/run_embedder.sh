#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <esm|t5|bert] <input_directory> <output_directory>"
    exit 1
fi

EMBEDDER="$1" # T5, ESM and BERT possible
input_dir="$2"
output_dir="$3"

BASE_URL="http://proteinspector.projekty.ms.mff.cuni.cz:42013"


if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

for file in "$input_dir"/*.fasta; do
    filename=$(basename -- "$file")
    filename_no_ext="${filename%.*}"
    out_file="$output_dir/sequences_${EMBEDDER}_${filename_no_ext}.zip"

    echo processing file $filename
    curl \
     -X GET \
     -F "fasta=@$file" \
     "$BASE_URL?embedder=$EMBEDDER" \
     --output $out_file

    echo
done
