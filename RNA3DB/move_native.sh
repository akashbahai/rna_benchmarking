#!/bin/bash

# Loop through all .pdb files in the current directory
for target in $(cat targetlist.txt); do
    
    
    full_loc="/Users/akash/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Macbook Pro backup/rna/rna/rna_targets/latest_analysis/new_targets"

    initial_loc="/Users/akash/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Macbook Pro backup/rna/rna/rna_targets/alltargets"
    
    
    # Move the file to the new directory with the desired name
    cp "${initial_loc}/${target}/${target}_native.pdb" "${full_loc}/${target}/"
    cp "${initial_loc}/${target}/${target}.fasta" "${full_loc}/${target}/"
done
