from pathlib import Path
from vina import Vina

TARGETS = {
    "AChE — Acetylcholinesterase (Alzheimer's)": {
        "receptor": "data/protein_ready.pdbqt",
        "protein_pdb": "data/protein_clean.pdb",
        "box_center": [2.8, 64.4, 68.0],
        "box_size":   [25.0, 25.0, 25.0],
        "pdb_id": "1EVE",
    },
    "EGFR — Epidermal Growth Factor Receptor (Cancer)": {
        "receptor": "data/egfr_ready.pdbqt",
        "protein_pdb": "data/egfr_clean.pdb",
        "box_center": [22.0, 0.3, 52.8],
        "box_size":   [25.0, 25.0, 25.0],
        "pdb_id": "1M17",
    },
    "BuChE — Butyrylcholinesterase (AChE homolog)": {
        "receptor": "data/buche_ready.pdbqt",
        "protein_pdb": "data/buche_clean.pdb",
        "box_center": [139.6, 116.1, 40.9],
        "box_size":   [25.0, 25.0, 25.0],
        "pdb_id": "4BDS",
    },
    "COX-2 — Cyclooxygenase-2 (off-target)": {
        "receptor": "data/cox2_ready.pdbqt",
        "protein_pdb": "data/cox2_clean.pdb",
        "box_center": [165.4, 203.3, 205.5],
        "box_size":   [25.0, 25.0, 25.0],
        "pdb_id": "5IKT",
    },
    "VEGFR2 — Vascular Endothelial Growth Factor Receptor (kinase off-target)": {
        "receptor": "data/vegfr2_ready.pdbqt",
        "protein_pdb": "data/vegfr2_clean.pdb",
        "box_center": [-24.6, -0.4, -10.9],
        "box_size":   [25.0, 25.0, 25.0],
        "pdb_id": "4ASD",
    },
}


def run_docking(ligand_pdbqt: str, out_pdbqt: str,
                target: str = "AChE — Acetylcholinesterase (Alzheimer's)",
                exhaustiveness: int = 8, n_poses: int = 5) -> list[float]:
    cfg = TARGETS[target]
    receptor = Path(cfg["receptor"])
    if not receptor.exists():
        raise FileNotFoundError(f"Receptor not found: {receptor}. Run setup_protein.py first.")

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(receptor))
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=cfg["box_center"], box_size=cfg["box_size"])
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(out_pdbqt, n_poses=n_poses, overwrite=True)

    scores = []
    with open(out_pdbqt) as f:
        for line in f:
            if "VINA RESULT" in line:
                scores.append(round(float(line.split()[3]), 2))
    return scores[:n_poses]


def read_pdbqt(path: str) -> str:
    with open(path) as f:
        return f.read()
