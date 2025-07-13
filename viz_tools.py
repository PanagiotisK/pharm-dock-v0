# viz_tools.py
# Tools for pharmacophore + ligand visualization in PyMOL

from rdkit import Chem
from rdkit.Chem import AllChem
import os

'''
You can expand this to:

    Draw lines from ligand atom → matched pharmacophore
    Label only matched features (using EmbedLib.EmbedPharmacophore)
    Show the protein pocket in background (protein.pdb from Fpocket)
'''

def save_ligand_conformer(mol, conf_id, filename="ligand_conf.pdb"):
    """
    Save one conformer of a molecule to PDB for PyMOL.
    """
    writer = Chem.rdmolfiles.PDBWriter(filename)
    writer.write(mol, confId=conf_id)
    writer.close()

def save_pharmacophore_as_pml(pharmacophore, filename="pharmacophore.pml"):
    """
    Save pharmacophore features as PyMOL commands (spheres + labels).
    """
    with open(filename, 'w') as f:
        for i, feat in enumerate(pharmacophore.GetFeatures()):
            pos = feat.GetPos()
            name = feat.GetFamily()
            color = "red" if name == "Donor" else "blue" if name == "Acceptor" else "yellow"
            f.write(f"pseudoatom pharm_{i}, pos=[{pos.x},{pos.y},{pos.z}], color={color}, label={name}_{i}\n")
        f.write("show spheres\n")
