#!/bin/bash

# Loop through all directories
for dirname in */; do
    # Remove trailing slash to get the directory name
    dir_name="${dirname%/}"
    
    # Check if dirname_native_rpr.pdb exists
    native_pdb="${dir_name}_native_rpr.pdb"
    if [ -f ${dir_name}/"$native_pdb" ]; then
        # Loop through other *_rpr.pdb files
        for model_pdb in ${dir_name}/"${dir_name}"*_rpr.pdb; do
            if [ "$model_pdb" != "$native_pdb" ]; then
                fasta_file="${dir_name}.fasta"
                output_file="${dir_name}_native.txt"
                mycommand_output=$(python rna_comparison/native_contact_frac.py "${dir_name}/${native_pdb}" "${model_pdb}" "${dir_name}/${fasta_file}")
                echo "${native_pdb} ${model_pdb} ${mycommand_output}" >> "${dir_name}/${output_file}"
            fi
        done
    else
        echo "Warning: ${native_pdb} does not exist in ${dir_name}"
    fi
done
