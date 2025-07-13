# scoring_tools.py
# Implements distance-weighted pharmacophore scoring
import numpy as np
from rdkit.Chem.Pharm3D import EmbedLib

# Define weights for different pharmacophore feature types
FEATURE_WEIGHTS = {
    'Donor': 2.0,
    'Acceptor': 1.5,
    'Hydrophobe': 1.0
}

def compute_distance(p1, p2):
    """
    Compute Euclidean distance between two 3D points.
    """
    return np.linalg.norm(np.array(p1) - np.array(p2))

def score_pharmacophore_match(mol, pharmacophore, conf_id):
    """
    Score a single conformer based on distance-weighted feature matches.
    """
    total_score = 0.0
    mapping = EmbedLib.EmbedPharmacophore(mol, pharmacophore, confId=conf_id, maxMatches=5)
    
    if mapping is None:
        return 0.0  # no matches

    # mapping is list of (ligand_atom_id, pharm_feature)
    conf = mol.GetConformer(conf_id)

    for lig_idx, pharm_feat in mapping:
        lig_coord = conf.GetAtomPosition(lig_idx)
        d = compute_distance(lig_coord, pharm_feat.GetPos())
        w = FEATURE_WEIGHTS.get(pharm_feat.GetFamily(), 1.0)
        score = w / (1 + d**2)
        total_score += score

    return total_score

def score_all_conformers(mol, pharmacophore):
    """
    Score all conformers of a molecule and return (conf_id, score) list.
    """
    results = []
    for conf_id in mol.GetConformerIds():
        score = score_pharmacophore_match(mol, pharmacophore, conf_id)
        results.append((conf_id, score))
    return sorted(results, key=lambda x: -x[1])
