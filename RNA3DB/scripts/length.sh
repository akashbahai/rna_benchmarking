#!/bin/bash

# Read each line from targetlist.txt
while IFS= read -r targetdir; do
    # Check if the directory exists
    if [ -d "$targetdir" ]; then
        fastafile="$targetdir/${targetdir}.fasta"

        # Check if the FASTA file exists
        if [ -f "$fastafile" ]; then
            # Get the targetname from the directory name
            targetname=$(basename "$targetdir")

            # Get the length of the second line in the FASTA file
            length=$(awk 'NR==2 { print length }' "$fastafile")

            # Print targetname and sequence length
            echo "$targetdir $length"
        else
            echo "Error: $fastafile not found"
        fi
    else
        echo "Error: $targetdir is not a directory"
    fi
done < $1
