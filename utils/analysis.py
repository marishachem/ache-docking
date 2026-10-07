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
        protein_mol = plf.Molecule.from_mda(u.select_atoms("protein"), force=True)

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

        rows = []
        for col in df.columns:
            residue, interaction = col[1], col[2]
            if df[col].iloc[0]:
                rows.append({"Residue": residue, "Interaction": interaction})
        return pd.DataFrame(rows) if rows else pd.DataFrame()
    except Exception as e:
        raise e


INTERACTION_COLORS = {
    "HBDonor":      "#60a5fa",
    "HBAcceptor":   "#34d399",
    "Hydrophobic":  "#fbbf24",
    "PiStacking":   "#a78bfa",
    "PiCation":     "#f472b6",
    "CationPi":     "#f472b6",
    "Anionic":      "#f87171",
    "Cationic":     "#fb923c",
    "VdWContact":   "#94a3b8",
    "EdgeToFace":   "#c084fc",
    "FaceToFace":   "#818cf8",
}


def interaction_chart(idf: pd.DataFrame) -> bytes:
    """Dot-matrix chart: residues (y) × interaction types (x), colored by type."""
    if idf is None or idf.empty:
        return b""

    residues = idf["Residue"].unique().tolist()
    inter_types = idf["Interaction"].unique().tolist()

    fig, ax = plt.subplots(figsize=(max(5, len(inter_types) * 1.1),
                                    max(3, len(residues) * 0.45 + 1.5)))
    fig.patch.set_facecolor("#0f1320")
    ax.set_facecolor("#0f1320")

    for yi, res in enumerate(residues):
        res_ints = idf[idf["Residue"] == res]["Interaction"].tolist()
        for xi, itype in enumerate(inter_types):
            if itype in res_ints:
                color = INTERACTION_COLORS.get(itype, "#e2e8f0")
                ax.scatter(xi, yi, s=220, color=color, zorder=3,
                           linewidths=0.5, edgecolors="#1e293b")

    ax.set_xticks(range(len(inter_types)))
    ax.set_xticklabels(inter_types, rotation=35, ha="right",
                       color="#e2e8f0", fontsize=9)
    ax.set_yticks(range(len(residues)))
    ax.set_yticklabels(residues, color="#e2e8f0", fontsize=9,
                       fontfamily="monospace")
    ax.set_xlim(-0.6, len(inter_types) - 0.4)
    ax.set_ylim(-0.6, len(residues) - 0.4)
    ax.tick_params(length=0)
    ax.spines[:].set_color("#253060")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#1e293b", linewidth=0.8)
    ax.xaxis.grid(True, color="#1e293b", linewidth=0.8)
    ax.set_title("Protein–Ligand Interactions", color="#e2e8f0",
                 fontsize=12, pad=10)

    # Legend
    seen = sorted(idf["Interaction"].unique())
    patches = [mpatches.Patch(color=INTERACTION_COLORS.get(t, "#e2e8f0"), label=t)
               for t in seen]
    ax.legend(handles=patches, loc="upper left", bbox_to_anchor=(1.01, 1),
              framealpha=0.15, labelcolor="#e2e8f0", fontsize=8,
              facecolor="#1e293b", edgecolor="#253060")

    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=140, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    buf.seek(0)
    return buf.read()
