# ligand_tools.py
# This module handles loading ligands and generating 3D conformers

from rdkit import Chem
from rdkit.Chem import AllChem

def load_ligand(smiles):
    """
    Convert SMILES to an RDKit molecule and add hydrogens.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    return mol

def generate_conformers(mol, num_confs=10):
    """
    Generate 3D conformers for the ligand.
    """
    AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=42)
    return mol
