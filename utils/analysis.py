import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import io
import pandas as pd


def score_chart(scores: dict[str, float]) -> bytes:
    """Bar chart comparing docking scores. scores = {name: kcal/mol}"""
    names = list(scores.keys())
    vals  = list(scores.values())
    colors = ["#c0392b" if i == 0 else "#2563eb" for i in range(len(names))]

    fig, ax = plt.subplots(figsize=(max(5, len(names) * 1.4), 4))
    fig.patch.set_facecolor("#0f1320")
    ax.set_facecolor("#1a2035")

    bars = ax.bar(names, vals, color=colors, width=0.5, zorder=3)
    ax.set_ylabel("Binding energy (kcal/mol)", color="#e2e8f0", fontsize=11)
    ax.set_title("Docking score comparison", color="#e2e8f0", fontsize=13, pad=12)
    ax.tick_params(colors="#e2e8f0")
    ax.spines[:].set_color("#253060")
    ax.yaxis.grid(True, color="#253060", linestyle="--", alpha=0.5, zorder=0)

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, val - 0.15,
                f"{val:.1f}", ha="center", va="top", color="white",
                fontsize=10, fontweight="bold")

    ax.annotate("More negative = stronger binding →",
                xy=(0.99, 0.04), xycoords="axes fraction",
                ha="right", fontsize=8, color="#9ca3af", style="italic")

    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=140, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    return buf.read()


def get_interactions(protein_pdb: str, ligand_pdbqt: str, smiles: str = None) -> pd.DataFrame | None:
    try:
        import warnings, MDAnalysis as mda, prolif as plf
        from rdkit import Chem
        from rdkit.Chem import AllChem
        warnings.filterwarnings("ignore")

        # Load protein (needs hydrogens — use protein_ready.pdb)
        u = mda.Universe(protein_pdb, guess_bonds=False)
        u.guess_TopologyAttrs(to_guess=["elements", "bonds"])
        protein_mol = plf.Molecule.from_mda(u.select_atoms("protein"))

        # Build ligand from SMILES + docked coordinates from PDBQT
        if smiles:
            coords = []
            with open(ligand_pdbqt) as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                    elif line.startswith("ENDMDL"):
                        break
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)
            conf = mol.GetConformer()
            heavy_idx = [i for i, a in enumerate(mol.GetAtoms()) if a.GetAtomicNum() != 1]
            for i, xyz in zip(heavy_idx[:len(coords)], coords):
                conf.SetAtomPosition(i, xyz)
            ligand_mol = plf.Molecule.from_rdkit(mol)
        else:
            lig_u = mda.Universe(ligand_pdbqt, guess_bonds=False)
            lig_u.guess_TopologyAttrs(to_guess=["elements", "bonds"])
            ligand_mol = plf.Molecule.from_mda(lig_u.atoms, inferrer=None)

        fp = plf.Fingerprint()
        fp.run_from_iterable([ligand_mol], protein_mol)
        df = fp.to_dataframe()

        # Flatten multi-level columns to readable strings
        rows = []
        for col in df.columns:
            residue, interaction = col[1], col[2]
            if df[col].iloc[0]:
                rows.append({"Residue": residue, "Interaction": interaction})
        return pd.DataFrame(rows) if rows else pd.DataFrame()
    except Exception as e:
        raise e
