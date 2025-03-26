#!/bin/bash

# Default values (optional)
CIF_FILE=""
MATRIX=""

# Function to display usage
usage() {
    echo "Usage: $0 --cif <input_identifier> --matrix <a,b,c>"
    echo "Example: $0 --cif 1241883 --matrix 10,3,4"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cif)
            CIF_FILE="$2"
            shift 2
            ;;
        --matrix)
            MATRIX="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if both arguments are provided
if [[ -z "$CIF_FILE" || -z "$MATRIX" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Define filenames based on input CIF
PDB_FILE="${CIF_FILE}.pdb"
UNIQUE_PDB_FILE="${CIF_FILE}_unique.pdb"
SUPERCELL_CIF="supercell.cif"
SUPERCELL_SDF="supercell.sdf"
SUPERCELL_PDB="supercell.pdb"
REORDERED_SUPERCELL_PDB="supercell_reordered.pdb"

echo "Processing CIF: $CIF_FILE.cif with matrix: $MATRIX"

# Step 1: Convert CIF to PDB
obabel "${CIF_FILE}.cif" -O "$PDB_FILE"

# Step 2: Remove CONECT lines from PDB
sed -i '/^CONECT/d' "$PDB_FILE"

# Step 3: Rename unique atom names
./unique_atom_name.py --input "$PDB_FILE" --output "$UNIQUE_PDB_FILE"

# Step 4: Copy unique PDB as template
cp "$UNIQUE_PDB_FILE" template.pdb

# Step 5: Generate supercell CIF using ASE script with matrix
./ASE_cif_to_pymatgen_supercell_cif.py --input "${CIF_FILE}.cif" --output "$SUPERCELL_CIF" --matrix "$MATRIX"

# Step 6: Convert supercell CIF to SDF
codcif2sdf "$SUPERCELL_CIF" > "$SUPERCELL_SDF"

# Step 7: Convert SDF to multiple PDBs
obabel "$SUPERCELL_SDF" -O mol.pdb -m --separate

# Step 8: Remove CONECT lines from generated PDBs
sed -i '/^CONECT/d' mol*.pdb

# Step 9: Run batch reorder script
./batch_reorder.sh

# Step 10: Join reordered PDB files into a single supercell PDB
obabel mol*_reordered.pdb -O "$SUPERCELL_PDB" --join

# Step 11: Clean up intermediate PDB files
rm -f mol*.pdb

# Step 12: Perform final mapping sequence using template with matrix
./mapping_sequence.py --input "$SUPERCELL_PDB" --output "$REORDERED_SUPERCELL_PDB" --template template.pdb --matrix "$MATRIX"

echo "Processing completed. Final reordered supercell PDB: $REORDERED_SUPERCELL_PDB"

