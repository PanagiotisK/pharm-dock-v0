# protein_tools.py
# This module handles protein loading and pocket detection using Fpocket

import subprocess
import os

def run_fpocket(pdb_path):
    """
    Run fpocket on a given PDB file and return the path to detected pocket PDB files.
    """
    subprocess.run(['fpocket', '-f', pdb_path])
    pocket_dir = pdb_path.replace('.pdb', '_out/pockets/')
    pocket_files = [os.path.join(pocket_dir, f) for f in os.listdir(pocket_dir) if f.endswith('_atm.pdb')]
    return sorted(pocket_files)
