#!/bin/bash

# Function to calculate RMSD using rnaalign tool
calculate_rmsd() {
    native=$1
    model=$2
    fasta=$3
    nativefrac_output=$(python /Users/akash/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/Macbook Pro backup/rna/rna/rna_targets/latest_analysis/rna_comparison/native_contact_frac.py $native $model $fasta)
    local nat_frac=$(echo "$nativefrac_output" | grep -oP 'Fnat: \K.*')
    echo "$nativefrac_output"
}

# CSV header
csv_header="Target,df,drfold,rf,rhofold, rnac, trr, 3drna"

# Output CSV file name
output_csv="nativefrac_comparison.csv"

work_dir=$(readlink -f "$(dirname "$1")")

# Create or clear the output CSV file
echo "$csv_header" > ${work_dir}/"$output_csv"

# Read the list of target directories from targetlist.txt
while IFS= read -r target_dir
do
    cd ${work_dir}
    # Check if the target directory exists
    if [ -d "$target_dir" ]; then
        # Prepare the line for the CSV file with the target name
        csv_line="$target_dir"

        # Native PDB file name pattern
        native_pdb_pattern="${target_dir}_native.pdb"

        # Array of model tool directories
        model_tools=("df" "drfold" "rf" "rhofold" "rnac" "trr" "3drna")

        # Iterate through each model tool directory and calculate RMSD
        for model_dir in "${model_tools[@]}"
        do
            # Find the PDB files in the target directory and model tool directory
            native_pdb=$(find "$target_dir" -type f -name "$native_pdb_pattern")
            model_pdb=$(find "$target_dir" -type f -name "${target_dir}_${model_dir}.pdb" | head -n 1)

            # Check if both native and model PDB files were found
            if [ -f "$native_pdb" ] && [ -f "$model_pdb" ]; then
                # Calculate RMSD using rnaalign and store the value
                #echo $native_pdb $model_pdb
                rmsd_value=$(calculate_rmsd ${work_dir}/"$native_pdb" ${work_dir}/"$model_pdb" ${work_dir}/${target_dir}/"${target_dir}.fasta" )
                csv_line="$csv_line,$rmsd_value"
            else
                # If either native or model PDB file is missing, set RMSD value as "NA"
                csv_line="$csv_line,NA"
            fi
        done

        # Append the line to the CSV file
        echo "$csv_line" >> "$output_csv"
    else
        echo "Target directory '$target_dir' does not exist."
    fi
done < $1
