#!/bin/bash

# Check if the input file exists
input_file=$1
if [ ! -e "$input_file" ]; then
    echo "Targetlist file '$input_file' not found."
    exit 1
fi

calculate_rmsd() {
    native=$1
    model=$2
    target=$3
    #echo $target
    rmsd_output=$(rna_calc_rmsd.py --target-selection=$target --model-selection=${target} -t $native $model)
    #echo $rmsd_output
    local rmsd=$(echo "$rmsd_output"| awk -F ',' 'END{print $2}')
    echo "$rmsd"
}

calculate_nativefrac() {
    native=$1
    model=$2
    fasta=$3
    nativefrac_output=$(python /Users/arthurdayne/rna/rna_targets/latest_analysis/rna_comparison/native_contact_frac.py $native $model $fasta)
    local nat_frac=$(echo "${nativefrac_output}" | grep -o 'Fnat: .*' | cut -d ' ' -f 2-)
    echo "${nat_frac}"
}

# CSV header
csv_header="Target,df,drfold,rf,rhofold, rnac, trr, 3drna"

# Output CSV file name
output_csv="rmsd_comparison.csv"
nativefrac_output_csv="nativefrac_comparison.csv"
nativefrac_cut="nativefrac_comparison_cut.csv"

work_dir=$(readlink -f "$(dirname "$1")")

# Create or clear the output CSV csvfile
echo "$csv_header" > "$output_csv"

#
echo "$csv_header" > "${nativefrac_output_csv}"

#
echo "$csv_header" > "${nativefrac_cut}"

mainextractdir="extracted_pdb"
mkdir $mainextractdir

# Loop through each directory listed in the input file
while IFS= read -r target || [[ -n "$target" ]]; do

    # Check if the directory exists
    if [ ! -d "$target" ]; then
        echo "Target '$target' not found. Skipping."
        continue
    fi

    csv_line="$target"
    nativefrac_csvline="$target"
    nativefrac_cutline="$target"

    # Check if the target pdb file exists in the directory
    target_pdb="${target}/${target}_3drna.pdb"
    if [ ! -e "$target_pdb" ]; then
        echo "Target pdb file '$target_pdb' not found in '$target'. Skipping."
        continue
    fi

    # Run your command (replace 'mycommand' with your actual command)
    echo "Correcting 3DRNA files!"

    python correct_residue.py ${target_pdb}

    echo "Correct DeepFoldRNA files!"

    python correct_residue_df.py $"${target}/${target}_df.pdb"

    # Check if the target_native.pdb file exists in the directory
    target_native_pdb="${target}/${target}_native.pdb"
    if [ ! -e "$target_native_pdb" ]; then
        echo "Target native pdb file '$target_native_pdb' not found in '$target'. Skipping."
        continue
    fi

    fastafile="${target}/${target}.fasta"
    #preprocess native pdb files
    python preprocessing_latest.py "${target_native_pdb}" "${fastafile}"

    processed_native_pdb="${target}/${target}_native_remapped.pdb"

    #make new directory to store processed files
    mkdir -p "${target}/processing"

    #take care of altlocations in the remapped pdb file
    pdb_selaltloc ${processed_native_pdb} > ${target}/native_altloc.pdb

    #fix missing atoms in the native pdb file
    pdbfixer ${target}/native_altloc.pdb --output ${processed_native_pdb} --add-atoms heavy
    rm ${target}/native_altloc.pdb

    final_native_file="${target}/processing/${target}_native_rpr.pdb"

    #make native pdb puzzle ready
    rna_pdb_tools.py --rpr --fix "${processed_native_pdb}" > "${final_native_file}"

    #reading the fragment file
    fragment_file="${target}/${target}_fragments.txt"
    fragment=$(head -n 1 ${fragment_file} | sed -e "s/'//g")
    toextract=$(head -n 1 ${fragment_file})

    #extract native file
    extractdir="${mainextractdir}/${target}"
    extracted_native="${extractdir}/${target}_native_cut.pdb"
    echo "$extractdir"
    mkdir -p "${extractdir}"
    rna_pdb_tools.py --extract $fragment ${final_native_file} > ${extractdir}/${target}_native_cut.pdb

    # Array of model tool directories
    model_tools=("df" "drfold" "rf" "rhofold" "rnac" "trr" "3drna")

    for model in "${model_tools[@]}"
        do
            model_pdb=$(find "$target" -type f -name "${target}_${model}.pdb" | head -n 1)
            echo $model_pdb
            # Check if model PDB files were found
            if [ -f "$model_pdb" ]; then
                #Set correct chain name
                chainoutput="${target}/processing/${target}_${model}_chainA.pdb"
                #echo "rna_pdb_tools --set-chain A ${model_pdb} > ${chainoutput}"
                #rna_pdb_tools --set-chain "A" "${model_pdb}" > "${chainoutput}"
                pdb_chain -A "${model_pdb}" > "${chainoutput}"

                final_model_file="${target}/processing/${target}_${model}_rpr.pdb"
                fastafile="${target}/${target}.fasta"
                # make model file puzzle-ready

                rna_pdb_tools.py --rpr "${chainoutput}" > "${final_model_file}"

                #Calculate RMSD
                echo "Calculating RMSD for ${fragment} for ${final_model_file} against ${final_native_file} "
                rmsd_value=$(calculate_rmsd ${final_native_file} ${final_model_file} ${fragment})
                echo ${rmsd_value}
                csv_line="$csv_line,$rmsd_value"

                #Calculate NAtivefraction
                echo "Calculating NATIVE_CONTACT_FRACTON for ${fragment} for ${final_model_file} against ${final_native_file} "
                ncf_value=$(calculate_nativefrac ${final_native_file} ${final_model_file} ${fastafile})
                echo ${ncf_value}
                nativefrac_csvline="${nativefrac_csvline},${ncf_value}"

                #remove chain_A files
                rm $chainoutput
                
                #extract files for each model based on fragments
                extracted_model="${extractdir}/${target}_${model}_cut.pdb"
                rna_pdb_tools.py --extract $fragment ${final_model_file} > ${extracted_model}


                #Calculate nativefraction based on cut file
                echo "Calculating NATIVE_CUT_CONTACT_FRACTON for ${fragment} for ${extracted_model} against ${extracted_native} "
                cutncf_value=$(calculate_nativefrac ${extracted_native} ${extracted_model} ${fastafile})
                echo ${cutncf_value}
                nativefrac_cutline="${nativefrac_cutline},${cutncf_value}"
            
            else
                # If either native or model PDB file is missing, set RMSD value as "NA"
                echo "${model} is missing in ${target}."
                csv_line="$csv_line,NA"
                nativefrac_csvline="${nativefrac_csvline},NA"
                nativefrac_cutline="${nativefrac_cutline},NA"
            fi
        done

    # Append the line to the CSV file
    echo "$csv_line" >> "${output_csv}"
    echo "${nativefrac_csvline}" >> "${nativefrac_output_csv}"
    echo "${nativefrac_cutline}" >> "${nativefrac_cut}"
done < "$input_file" 
