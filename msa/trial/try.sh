#!/usr/bin/env bash

pad_sequences() {
    file=$1
    outfile="${file%.*}_equal.${file##*.}"
    max_length=0

    # Find the maximum length
    while IFS= read -r line; do
        if [[ $line != \>* ]]; then
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
        else
            len=${#line}
            if (( len < max_length )); then
                padding=$(printf '%0.s-' $(seq $len $max_length))
                echo "$line$padding" >> "$outfile"
            else
                echo "$line" >> "$outfile"
            fi
        fi
    done < "$file"
}

pad_sequences test.a3m

