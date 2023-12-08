#!/bin/bash

# Loop through all .pdb files in the current directory
for target in $(cat targetlist.txt); do
    
    drna=${target}_3drna.pdb
    rnac=${target}_rnac.pdb
    trr=${target}_trr.pdb
    
    full_loc="/Users/akash/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Macbook Pro backup/rna/rna/rna_targets/latest_analysis/casp15"

    initial_loc="/Users/akash/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Macbook Pro backup/rna/rna/rna_targets/alltargets"
    
    
    # Move the file to the new directory with the desired name
    cp "${initial_loc}/${target}/${drna}" "${full_loc}/${target}/${drna}"
    cp "${initial_loc}/${target}/${rnac}" "${full_loc}/${target}/${rnac}"
    cp "${initial_loc}/${target}/${trr}" "${full_loc}/${target}/${trr}"
done
