import os
from typing import List

import MDAnalysis as mda
from MDAnalysis import Universe


def read_pdb_file(file: str) -> Universe:
    return mda.Universe(file)


def read_pdb_folder(folder: str) -> List[Universe]:
    pdb_files = [f for f in os.listdir(folder) if f.endswith(".pdb")]
    universes = []

    for pdb_file in pdb_files:
        pdb_path = os.path.join(folder, pdb_file)
        universes.append(mda.Universe(pdb_path))

    return universes
