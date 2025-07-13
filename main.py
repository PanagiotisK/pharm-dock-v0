# main.py
# This is the main orchestrator script that puts all modules together

from ligand_tools import load_ligand, generate_conformers
from protein_tools import run_fpocket
from pharmacophore_tools import generate_pharmacophore
from docking_tools import score_ligand_conformers
from scoring_tools import score_all_conformers
from viz_tools import save_ligand_conformer, save_pharmacophore_as_pml

# STEP 1: Run Fpocket to detect protein pockets
pockets = run_fpocket("1AKE.pdb")
print("Detected pockets:", pockets)

# STEP 2: Load ligand and generate conformers
smiles = "CCOc1ccc2nc(S(N)(=O)=O)sc2c1"
ligand = load_ligand(smiles)
ligand = generate_conformers(ligand)

# STEP 3: Build pharmacophore from ligand
pharmacophore = generate_pharmacophore(ligand)

# STEP 4: Align ligand conformers to pharmacophore
scores = score_ligand_conformers(ligand, pharmacophore)

# STEP 5: Rank poses
for conf_id, score in sorted(scores, key=lambda x: -x[1]):
    print(f"Conformer {conf_id} matched {score} pharmacophore points")

# STEP 6: Extended score conformers using distance-weighted pharmacophore scoring
dw_scores = score_all_conformers(ligand, pharmacophore)
for dw_conf_id, dw_score in dw_scores:
    print(f"Conformer {dw_conf_id}: score = {dw_score:.3f}")

# STEP 7: Save top conformer and pharmacophore
top_conf = scores[0][0]
save_ligand_conformer(ligand, top_conf, "ligand_top.pdb")
save_pharmacophore_as_pml(pharmacophore, "pharmacophore.pml")

print("Open in PyMOL:")
print(" pymol ligand_top.pdb pharmacophore.pml")
