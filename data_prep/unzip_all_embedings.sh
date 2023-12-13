#!/bin/bash

input_folder="$1"
target_folder="$2"

# Check if input arguments are provided
if [ -z "$input_folder" ] || [ -z "$target_folder" ]; then
    echo "Usage: $0 <input_folder> <target_folder>"
    exit 1
fi

# Create the target folder if it doesn't exist
mkdir -p "$target_folder"

# Unzip all .zip files in the input folder
for file in "$input_folder"/*.zip; do
    echo unzip file $file
    unzip -o -d "$target_folder" "$file" | awk -v target="$target_folder" '
        BEGIN {
            FS=OFS=":"
        }
        {
            if ($1 == "  inflating") {
                filename = $2
                gsub(/^[ \t]+|[ \t]+$/, "", filename)
                if (filename in files) {
                    new_filename = "copy_of_" filename
                    print "Renaming duplicate:", filename, "->", new_filename > "/dev/stderr"
                    system("mv \"" target "/" filename "\" \"" target "/" new_filename "\"")
                } else {
                    files[filename] = 1
                }
            }
        }
    '
done
