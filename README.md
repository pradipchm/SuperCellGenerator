# Super Cell Constructor
Tool for setting up molecular crystal supercells from CIF files

## Authors:
- [Pradip Si](https://www.valsson.info/members/pradip-si), University of North Texas

## Requirments
- ASE
- Pymatgen
- Rdkit
- Open babel
- NumPy
- COD Tools [cod-tools](https://wiki.crystallography.net/cod-tools/)

## Instructions 

### Generate a PBC supercell from a CIF file collected from the CCDC database
1. Save the CIF files for all polymorphs as PDB (obabel file.cif -O file.pdb) and add missing Hs if required ( `obabel file.pdb -O file.pdb -h`) using openbabel and delete CONECT lines (sed -i '/^CONECT/d' input.pdb)
2. Select one as a template (preferably choose a unit cell with one molecule) and give a unique atom name using unique_atom_name.py (that will be used to match the atom sequence for other PDBs).
3. Use `ASE_cif_to_pymatgen_supercell_cif.py` to generate the supercell cif file for the template. (Change the space group to "P1" in the cif file if it shows an error due to the space group) 
4. Convert the supercell.cif to sdf (codcif2sdf supercell.cif > supercell.sdf)
5. Split all the molecules (obabel supercell.sdf -O mol.pdb -m --separate) and delete CONECT for all PDB files (sed -i '/^CONECT/d' mol*.pdb)
6. Reorder all the PDBs with a template.pdb (./batch_reorder.sh (There should be a file template.pdb in the path))
7. Merge all the reordered PDBs (obabel mol*_reordered.pdb -O supercell.pdb --join) and delete those files (rm mol*.pdb)  
8. Use `mapping_sequence.py` to match the atom sequence to the template PDB( ```./mapping_sequence_new.py --input supercell.pdb --output supercell_reordered.pdb --template template.pdb --matrix 6,5,4```. 

### If the molecule is not selected as the template file, before proceeding to the step3

- First, reorder the PDB file to match the atom sequence with the template PDB file using `reorder-atoms.py`. (Follow steps 4,5 and 6 for unit cells having more than one molecule to split and reorder.)
- It should pass the validation ( might visualize to ensure that the connectivity is similar to the template).
- Return to step 3.


**Check the new file after each step to ensure that the code is doing the correct task**

### Use the Example codes to run all the steps simultaneously for template and non-template (with one and more than one molecule) crystals to make the supercell 
```
   ./template.sh --cif input.cif --matrix a,b,c 
   ./nontemplate.sh --cif input.cif --matrix a,b,c 
   ./nontemplate_morethanone.sh --cif input.cif --matrix a,b,c
``` 

## Acknowledgements
The development of this tool was supported by a DOE Early Career Award (BES Condensed Phase and Interfacial Molecular Science (CPIMS) / DE-SC0024283)


