# Pharmacophore-Guided Docking Pipeline

## Overview
A modular Python pipeline to detect protein pockets, generate ligand conformers, build a pharmacophore model, align ligands, and score docking poses.

## Installation

```bash
git clone https://github.com/PanagiotisK/pharm-dock-v0
cd pharm-dock-v0
```

# Install Python dependencies
```bash
pip install -r requirements.txt
```

# Install fpocket
```bash
sudo apt-get install fpocket
```

# Install PyMol
```bash
sudo apt-get install pymol-open-source
```

## Usage
```bash
python main.py
```

This will:
- Detect pockets in examples/protein.pdb
- Load the ligand from examples/ligand.smiles
- Generate conformers
- Build a pharmacophore model
- Score each conformer

## Files

- protein_tools.py – run fpocket and extract pocket files
- ligand_tools.py – load SMILES and generate conformers
- pharmacophore_tools.py – build pharmacophore from ligand
- docking_tools.py – align conformers and score
- main.py – orchestrator of modules

## Next Steps

- Filter poses to a specific pocket region
- Integrate docking engines like AutoDock or Open Babel
- Add physics-based scoring (H-bonds, VDW, electrostatics)
- Validate using active/inactive ligands


