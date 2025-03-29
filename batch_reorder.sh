#!/bin/bash

# Loop through all files matching mol*.pdb
for file in mol*.pdb; do
    # Define output file by appending "_reordered" before ".pdb"
    output="${file%.pdb}_reordered.pdb"
    
    # Run reorder-atoms.py script
    ./reorder_atoms.py --input "$file" --output "$output" --template template.pdb

    echo "Processed: $file -> $output"
done

echo "All files processed!"

