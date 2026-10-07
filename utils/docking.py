from pathlib import Path
from vina import Vina

# AChE active site center (E20 donepezil centroid from 1EVE crystal structure)
BOX_CENTER = [2.8, 64.4, 68.0]
BOX_SIZE   = [25.0, 25.0, 25.0]
RECEPTOR   = Path("data/protein_ready.pdbqt")


def run_docking(ligand_pdbqt: str, out_pdbqt: str, exhaustiveness: int = 8, n_poses: int = 5) -> list[float]:
    if not RECEPTOR.exists():
        raise FileNotFoundError("Protein not prepared. Run setup_protein.py first.")

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(RECEPTOR))
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(out_pdbqt, n_poses=n_poses, overwrite=True)

    energies = v.energies(n_poses=n_poses)
    return [round(e[0], 2) for e in energies]


def read_pdbqt(path: str) -> str:
    with open(path) as f:
        return f.read()
