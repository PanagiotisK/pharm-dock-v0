# docking_tools.py
# This module aligns ligand conformers to the pharmacophore and scores them

from rdkit.Chem.Pharm3D import EmbedLib

def score_ligand_conformers(mol, pharmacophore):
    """
    Score all conformers by how well they match the pharmacophore.
    Returns a list of (conf_id, score).
    """
    scores = []
    for conf_id in mol.GetConformerIds():
        mapping = EmbedLib.EmbedPharmacophore(mol, pharmacophore, confId=conf_id, maxMatches=5)
        score = len(mapping) if mapping else 0
        scores.append((conf_id, score))
    return scores
