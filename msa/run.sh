#!/bin/bash

# Check if the input file exists
input_file=$1
if [ ! -e "$input_file" ]; then
    echo "Targetlist file '$input_file' not found."
    exit 1
fi

# This function calculates the number of effective sequences in the msa file using hmmbuild
calculate_neff() {
    msa=$1
    #echo $msa
    msa_output=$(hmmbuild out.hmm $msa | grep -A2 "eff_nseq" | awk 'NR==3 {print $7}')
    echo "${msa_output}"
}

# This function calculates the total number of sequences in the msa file using hmmbuild
calculate_nseq() {
    msa=$1
    #echo $msa
    nseq=$(hmmbuild out.hmm $msa | grep -A2 "eff_nseq" | awk 'NR==3 {print $3}')
    echo "${nseq}"
}

# This function replaces all U with T in the first sequence of the a3m msa file
correct_firstseq() {
    file=$1
    # Find the line numbers of the lines starting with '>'
    pos1=$(awk '/^>/{print NR; exit}' $file)
    pos2=$(awk -v p1=$pos1 'NR>p1 && /^>/{print NR; exit}' $file)
    # Use awk to replace 'U' with 'T' between pos1 and pos2
    awk -v p1=$pos1 -v p2=$pos2 'NR>p1 && NR<p2 {gsub(/U/,"T")};1' $file > $file.tmp
}

# This function replaces all U with T in the a3m msa file
correct_allseq() {
    file=$1
    #replace all U with T in the lines not starting with a '>'
    awk '/^>/{print;next}{gsub(/U/,"T")}1' $file > $file.tmp
}

# This function pads the sequences in the a3m msa file to the same length
pad_sequences() {
    file=$1
    dirpath=$(dirname $file)
    filename=$(basename $file .afa.tmp)
    outfile=$dirpath/${filename%.*}.afa
    max_length=0

    # Find the maximum length
    while IFS= read -r line; do
        if [[ $line != \>* ]] && [[ -n $line ]]; then
            len=${#line}
            if (( len > max_length )); then
                max_length=$len
            fi
        fi
    done < "$file"

    # Pad the sequences and write to the new file
    while IFS= read -r line; do
        if [[ $line == \>* ]]; then
            echo "$line" >> "$outfile"
        elif [[ -n $line ]]; then
            len=${#line}
            if (( len < max_length )); then
                padding=$(printf '%0.s-' $(seq $len $((max_length - 1))))
                echo "$line$padding" >> "$outfile"
            else
                echo "$line" >> "$outfile"
            fi
        fi
    done < "$file"
}

# This function formats the a3m file to afa format
reformat() {
    file=$1
    dirpath=$(dirname $file)
    filename=$(basename $file .a3m.tmp)
    reformat.pl fas a3m $file $dirpath/${filename%.*}.afa.tmp
}

# CSV header
csv_header="Target,df,rf,rhofold,trr"

# Output CSV file name
n_csv="n_comparison.csv"
neff_csv="neff_comparison.csv"

work_dir=$(readlink -f "$(dirname "$1")")

# Create or clear the output CSV csvfile
echo "$csv_header" > "$n_csv"

#
echo "$csv_header" > "${neff_csv}"


# Loop through each directory listed in the input file
while IFS= read -r target || [[ -n "$target" ]]; do

    # Check if the directory exists
    if [ ! -d "$target" ]; then
        echo "Target '$target' not found. Skipping."
        continue
    fi

    n_line="$target"
    neff_line="$target"

    #process rhofold.a3m file
    rhofold_a3m=$(find "$target" -type f -name "rhofold.a3m" | head -n 1)
    echo ${rhofold_a3m}
    if [ -f "$rhofold_a3m" ]; then
        correct_firstseq ${rhofold_a3m}
        reformat ${rhofold_a3m}.tmp
        python pad.py ${target}/rhofold.afa.tmp
        rm ${rhofold_a3m}.tmp

    fi

    #process trr.a3m file
    trr_a3m=$(find "$target" -type f -name "trr.a3m" | head -n 1)
    if [ -f "$trr_a3m" ]; then
        correct_allseq ${trr_a3m}
        reformat ${trr_a3m}.tmp
        python pad.py ${target}/trr.afa.tmp
        rm ${trr_a3m}.tmp
    fi

    # Array of msa files in each directory
    msafiles=("df.afa" "rf.afa" "rhofold.afa" "trr.afa")

    for msa in "${msafiles[@]}"
        do
            msa_file=$(find "$target" -type f -name "${msa}" | head -n 1)
            echo $msa_file
            # Check if model PDB files were found
            if [ -f "$msa_file" ]; then
                
                echo "Calculating Neff of the msa for ${msa_file} for ${target}"
                neff_value=$(calculate_neff ${msa_file})
                echo ${neff_value}
                nseq_value=$(calculate_nseq ${msa_file})
                neff_line="${neff_line},${neff_value}"
                n_line="${n_line},${nseq_value}"
            
            else
                # If either native or model PDB file is missing, set RMSD value as "NA"
                echo "${msa} is missing in ${target}."
                n_line="${n_line},NA"
                neff_line="${neff_line},NA"
            fi
        done

    # Append the line to the CSV file
    echo "${n_line}" >> "${n_csv}"
    echo "${neff_line}" >> "${neff_csv}"
done < "$input_file" 
