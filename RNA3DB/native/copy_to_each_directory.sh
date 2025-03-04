#!/bin/bash

# Set the source directory where the .pdb files are located
source_directory="/Users/arthurdayne/rna/rna3db_pdb/RNA3DB/NATIVE"

# Change to the current directory
cd "$(pwd)" || exit

# Read the target file names from targetlist.txt
while IFS= read -r file_name; do
    # Construct the full source file path with the .pdb extension
    source_file="$source_directory/$file_name.pdb"
    
    # Construct the new file name with the _native suffix
    new_file_name="/Users/arthurdayne/rna/rna3db_analysis/rna3db_results/${file_name}/${file_name}_native.pdb"

    # Check if the source file exists
    if [ -f "$source_file" ]; then
        # Copy the file to the current directory with the new name
        cp "$source_file" "$new_file_name"
    else
        echo "Warning: Source file '$source_file' does not exist."
    fi
done < "targetlist.txt"
