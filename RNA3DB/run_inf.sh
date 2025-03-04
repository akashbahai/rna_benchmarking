#!/bin/bash

# Set the input directory containing the extracted PDB files
INPUT_DIR="extracted_pdb"
# Set the output directory for results
RESULTS_DIR="results"

# Create the results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

# Read the list of folder names from targetlist.txt
while IFS= read -r foldername; do
    # Construct paths for the native file and output CSV
    native_file="${INPUT_DIR}/${foldername}/${foldername}_native_cut.pdb"
    output_file="${RESULTS_DIR}/inf_${foldername}.csv"

    # Change to the directory containing the PDB files
    pdb_directory="${INPUT_DIR}/${foldername}"

    # Check if the native file exists
    if [ -f "$native_file" ]; then
        # Run the command for all PDB files in the directory
        rna_calc_inf.py -t "$native_file" "$pdb_directory"/*.pdb -o "$output_file"

        # Check for errors in command execution
        if [ $? -ne 0 ]; then
            echo "Error processing folder: $foldername"
        else
            echo "Successfully processed folder: $foldername"
        fi
    else
        echo "Warning: Native file not found for folder: $foldername"
    fi
done < "targetlist.txt"
