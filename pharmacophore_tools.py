# pharmacophore_tools.py
# This module uses RDKit to generate pharmacophore models

from rdkit.Chem.Pharm3D import EmbedLib

def generate_pharmacophore(mol):
    """
    Generate a basic pharmacophore model using RDKit's default features.
    """
    factory = EmbedLib.PharmacophoreFactory()
    factory.AddFeaturesFromMol(mol, includeOnly="Donor,Acceptor,Hydrophobe")
    return factory.GetPharmacophore()
