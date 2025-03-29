#!/usr/bin/env python

import numpy as np
import rdkit.Chem
import re

def reorder_atoms(mol_pdb_fname, template_pdb_fname, output_pdb_fname):
    from rdkit.Chem import rdmolfiles
    """
    function reorders the atom sequence of a PDB file to match the atom sequence of another PDB file.

    Parametars:
    ----------
    mol_pdb_fname : str
        pdb file which atom sequence need to be matched with template pdb
    template_pdb_fname : str
        pdb file which will be use as template to match the atom sequence of other pdb
    output_pdb_fname : str
        pdb file with the reordered atom sequence
        
    """
    mol_to_transform = rdkit.Chem.rdmolfiles.MolFromPDBFile(mol_pdb_fname, removeHs=False)
    transform_order = list(rdmolfiles.CanonicalRankAtoms(mol_to_transform, includeChirality=False))

    mol_template = rdkit.Chem.rdmolfiles.MolFromPDBFile(template_pdb_fname, removeHs=False)
    template_order = list(rdmolfiles.CanonicalRankAtoms(mol_template, includeChirality=False))

    if len(template_order) != len(transform_order):
        raise RuntimeError('Number of atoms differs between template and molecule to transform.')

    i_transform_order = [int(i) for i in np.argsort(transform_order)]
    i_template_order =  [int(i) for i in np.argsort(template_order)]

    N_atoms = len(template_order)

    pos_to_transform = mol_to_transform.GetConformers()[0].GetPositions()
    lines = [None]*N_atoms
    for _, (otr, ote) in enumerate(zip(i_transform_order, i_template_order)):
        #print(mol_to_transform.GetAtoms()[otr].GetPDBResidueInfo().GetName(),
        #      mol_template.GetAtoms()[ote].GetPDBResidueInfo().GetName())
        pdb_entry_template = mol_template.GetAtoms()[ote].GetPDBResidueInfo()
        lines[ote] = '{ATOM:<6}{serial_number:>5} {atom_name:<4}{alt_loc_indicator:<1}{res_name:<3} {chain_id:<1}{res_seq_number:>4}{insert_code:<1}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp_factor:6.2f}           {atom_type:<3}'.format(
                    ATOM='HETATM',
                    serial_number=pdb_entry_template.GetSerialNumber(),
                    atom_name=pdb_entry_template.GetName(),
                    alt_loc_indicator=' ',
                    res_name=pdb_entry_template.GetResidueName(),
                    chain_id=pdb_entry_template.GetChainId(),
                    res_seq_number=pdb_entry_template.GetResidueNumber(),
                    insert_code=' ',
                    x= pos_to_transform[otr, 0],
                    y= pos_to_transform[otr, 1],
                    z= pos_to_transform[otr, 2],
                    occupancy=1.0,
                    temp_factor=0.0,
                    atom_type = re.split(r'(\s*\d+)', pdb_entry_template.GetName().strip())[0]
            )
    with open(mol_pdb_fname, 'r') as infile, open(output_pdb_fname, 'w') as f:
        for line in infile:
            if line.startswith('HETATM'):
                break
            f.write(line)
        for line in lines:
            if line is not None:
                print(line, file=f)
        print("END", file=f)


def validate(template_pdb_fname, output_pdb_fname):
    import rdkit.Chem.rdPartialCharges
    """
    Function compares two pdb files and returns 'Validation OK' if the the PDBs are similar
    
    Parametars:
    ----------
    template_pdb_fname : str
        template pdb which is being compared with
    output_pdb_fname : str
        pdb file that needs to match with template
        
    """
    
    mol_template = rdkit.Chem.rdmolfiles.MolFromPDBFile(template_pdb_fname, removeHs=False)
    finished_molecule = rdkit.Chem.rdmolfiles.MolFromPDBFile(output_pdb_fname, removeHs=False)
    N_atoms = finished_molecule.GetNumAtoms()
    for i in range(N_atoms):
        assert finished_molecule.GetAtoms()[i].GetAtomicNum() == mol_template.GetAtoms()[i].GetAtomicNum()
        assert finished_molecule.GetAtoms()[i].GetDegree() == mol_template.GetAtoms()[i].GetDegree()
        assert finished_molecule.GetAtoms()[i].GetTotalDegree() == mol_template.GetAtoms()[i].GetTotalDegree(), i
        assert finished_molecule.GetAtoms()[i].GetHybridization() == mol_template.GetAtoms()[i].GetHybridization(), i
        assert finished_molecule.GetAtoms()[i].GetFormalCharge() == mol_template.GetAtoms()[i].GetFormalCharge(), i
        assert finished_molecule.GetAtoms()[i].GetTotalValence() == mol_template.GetAtoms()[i].GetTotalValence(), i    
    rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(finished_molecule, throwOnParamFailure=True)
    charges_finished = np.array([float(finished_molecule.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(N_atoms)])
    rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol_template, throwOnParamFailure=True)
    charges_template = np.array([float(mol_template.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(N_atoms)])
    assert np.allclose(charges_finished, charges_template)
    print('Validation OK')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='reorder ligand atoms according to template',
                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', metavar='pdb file', default='in.pdb',
                    help='(in) pdb file name of the molecule to be reoredered')
    parser.add_argument('--template', metavar='pdb file', default='template.pdb',
                    help='(in) pdb file name of the template molecule with the desired order')
    parser.add_argument('--output', metavar='pdb file', default='out.pdb',
                    help='(out) pdb file name for the reordered molecule')
    parser.add_argument('--novalidation', action='store_true',
                    help='Skip validation of output.')

    args = parser.parse_args()
    reorder_atoms(args.input, args.template, args.output)
    if not args.novalidation:
        validate(args.template, args.output)
