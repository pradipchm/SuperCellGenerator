#!/bin/bash

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cif) cif_base="$2"; shift ;;
        --matrix) matrix="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$cif_base" || -z "$matrix" ]]; then
    echo "Usage: ./nontemplate_morethanone.sh --cif <cif_basename> --matrix <a,b,c>"
    exit 1
fi

# Define filenames based on the basename
cif_file="${cif_base}.cif"
pdb_file="${cif_base}.pdb"
sdf_file="${cif_base}.sdf"

# Convert CIF to PDB and clean up
codcif2sdf "$cif_file" > "$sdf_file"
obabel "$sdf_file" -O "$pdb_file"
sed -i '/^CONECT/d' "$pdb_file"

# Convert CIF to SDF and split into separate molecules
obabel "$sdf_file" -O "${cif_base}_.pdb" -m --separate

# Remove CONECT lines from all generated PDB files
for pdb in ${cif_base}_*.pdb; do
    sed -i '/^CONECT/d' "$pdb"
done

# Reorder atoms in the first molecule
./reorder_atoms.py --input "${cif_base}_1.pdb" --output "${cif_base}_1_reorder.pdb" --template template.pdb

# Generate supercell
./ASE_cif_to_pymatgen_supercell_cif.py --input "$cif_file" --output supercell.cif --matrix "$matrix"
codcif2sdf supercell.cif > supercell.sdf

# Convert supercell SDF to separate PDB files and clean up
obabel supercell.sdf -O mol.pdb -m --separate
for pdb in mol*.pdb; do
    sed -i '/^CONECT/d' "$pdb"
done

# Batch reorder atoms
./batch_reorder.sh

# Merge reordered PDB files into a single supercell PDB
obabel mol*_reordered.pdb -O supercell.pdb --join
rm mol*.pdb

# Map sequence using the reordered molecule as a template
./mapping_sequence.py --input supercell.pdb --output supercell_reorder.pdb --template "${cif_base}_1_reorder.pdb" --input2 "$cif_file" --matrix "$matrix"

echo "Processing completed for $cif_base with matrix $matrix."

