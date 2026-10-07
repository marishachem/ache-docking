import streamlit as st
import streamlit.components.v1 as components
import tempfile, os
from pathlib import Path

st.set_page_config(page_title="Protein–Ligand Docking", page_icon="🔬", layout="wide")

# ── Per-target molecules ───────────────────────────────────────────────────────
TARGET_MOLECULES = {
    "AChE — Acetylcholinesterase (Alzheimer's)": {
        "🎓 Diploma molecule (piperidinone)": "O=C1N(c2ccc(OC)cc2)C(c2ccccc2)(C(=O)OCC)C(CC(=O)OCC)CC1c1ccccc1",
        "💊 Donepezil (reference drug)":      "O=C(Cc1ccc(OC)c(OC)c1)N1CCC(Cc2ccccc2)CC1",
        "💊 Rivastigmine":                     "CCN(C)C(=O)Oc1cccc(C(C)N(C)CC)c1",
        "💊 Galantamine":                      "COc1ccc2c(c1)C[C@H]1[C@@H](O)CC[N@@+]3(C)CC=C[C@@H]1[C@@H]23",
    },
    "EGFR — Epidermal Growth Factor Receptor (Cancer)": {
        "💊 Gefitinib (1st gen)":   "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "💊 Erlotinib (1st gen)":   "C#Cc1cccc(Nc2ncnc3cc(OCCO)c(OC)cc23)c1",
        "💊 Lapatinib (2nd gen)":   "CS(=O)(=O)CCNCc1ccc(-c2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)o1",
        "💊 Osimertinib (3rd gen)": "C=CC(=O)Nc1cc2c(Nc3ccc(F)c(NC(=O)/C=C/CN(C)C)c3)ncnc2cc1N(C)CCN1CCOCC1",
    },
}

TARGET_PROTEIN_PDB = {
    "AChE — Acetylcholinesterase (Alzheimer's)": "data/protein_clean.pdb",
    "EGFR — Epidermal Growth Factor Receptor (Cancer)": "data/egfr_clean.pdb",
}

# Hydrogenated PDBs needed for ProLIF interaction detection
TARGET_PROTEIN_H_PDB = {
    "AChE — Acetylcholinesterase (Alzheimer's)": "data/protein_ready.pdb",
    "EGFR — Epidermal Growth Factor Receptor (Cancer)": "data/egfr_ready.pdb",
}

TARGET_DESCRIPTIONS = {
    "AChE — Acetylcholinesterase (Alzheimer's)": (
        "Docking into **Acetylcholinesterase** (PDB: 1EVE), the Alzheimer's drug target. "
        "The diploma molecule is a piperidine-2-one derivative synthesised during a chemistry "
        "degree — structurally related to donepezil, the approved AChE inhibitor."
    ),
    "EGFR — Epidermal Growth Factor Receptor (Cancer)": (
        "Docking into **EGFR kinase domain** (PDB: 1M17), a key cancer drug target. "
        "EGFR mutations drive lung, breast and colorectal cancers. The three generations of "
        "approved inhibitors (gefitinib → osimertinib) show how drug resistance shapes drug design."
    ),
}

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🔬 Docking Lab")
    st.divider()

    target = st.selectbox("Target protein", list(TARGET_MOLECULES.keys()))
    molecules = TARGET_MOLECULES[target]

    st.divider()
    mode = st.radio("Input mode", ["Preset molecules", "Custom SMILES"])
    if mode == "Preset molecules":
        selected = st.selectbox("Choose molecule", list(molecules.keys()))
        smiles = molecules[selected]
        mol_name = selected.split("(")[0].strip().lstrip("💊🎓").strip()
    else:
        smiles = st.text_area("SMILES", placeholder="Paste SMILES here…", height=100)
        mol_name = "Custom molecule"

    if smiles:
        st.code(smiles, language=None)

    st.divider()
    exhaustiveness = st.slider("Exhaustiveness", 4, 32, 8,
                               help="Higher = more accurate but slower")
    n_poses = st.slider("Number of poses", 1, 10, 5)

    run_btn = st.button("▶ Run Docking", type="primary", use_container_width=True)

    st.divider()
    st.caption("**First-time setup:**")
    st.code("python3 setup_protein.py\npip install -r requirements.txt", language="bash")

# ── Main area ──────────────────────────────────────────────────────────────────
st.title("Protein–Ligand Docking")
st.markdown(TARGET_DESCRIPTIONS[target])

if not smiles:
    st.info("Select a molecule in the sidebar and click **Run Docking**.")
    st.stop()

# ── 2D structure + properties ──────────────────────────────────────────────────
col1, col2 = st.columns([1, 2])
with col1:
    st.subheader("2D Structure")
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        import base64

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            drawer = rdMolDraw2D.MolDraw2DSVG(320, 240)
            drawer.drawOptions().padding = 0.15
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            b64 = base64.b64encode(svg.encode()).decode()
            st.markdown(
                f'<img src="data:image/svg+xml;base64,{b64}" style="width:100%;border-radius:10px;">',
                unsafe_allow_html=True
            )
        else:
            st.error("Invalid SMILES")
            st.stop()
    except Exception as e:
        st.warning(f"Could not render 2D structure: {e}")

with col2:
    st.subheader("Molecular Properties")
    try:
        from utils.ligand import get_mol_props
        props = get_mol_props(smiles)
        if props:
            lipinski_pass = (
                props["MW"] <= 500 and props["LogP"] <= 5 and
                props["HBD"] <= 5 and props["HBA"] <= 10
            )
            p1, p2, p3 = st.columns(3)
            p1.metric("MW (g/mol)", props["MW"], help="Lipinski: ≤500")
            p2.metric("LogP", props["LogP"], help="Lipinski: ≤5")
            p3.metric("TPSA (Å²)", props["TPSA"], help="Veber: ≤140")
            p4, p5, p6 = st.columns(3)
            p4.metric("HB Donors", props["HBD"], help="Lipinski: ≤5")
            p5.metric("HB Acceptors", props["HBA"], help="Lipinski: ≤10")
            p6.metric("Rot. Bonds", props["RotBonds"], help="Veber: ≤10")
            if lipinski_pass:
                st.success("✅ Passes Lipinski Rule of Five")
            else:
                st.warning("⚠️ Fails one or more Lipinski rules")
    except Exception as e:
        st.warning(f"Could not compute properties: {e}")

st.divider()

# ── Docking ────────────────────────────────────────────────────────────────────
if run_btn:
    from utils.docking import TARGETS
    receptor = TARGETS[target]["receptor"]
    if not Path(receptor).exists():
        st.error(f"Receptor not prepared: `{receptor}`. Run `python3 setup_protein.py` first.")
        st.stop()

    with st.spinner("Preparing ligand…"):
        try:
            from utils.ligand import smiles_to_pdbqt
            lig_pdbqt = str(Path(tempfile.mkdtemp()) / "ligand.pdbqt")
            smiles_to_pdbqt(smiles, lig_pdbqt)
        except Exception as e:
            st.error(f"Ligand preparation failed: {e}")
            st.stop()

    with st.spinner("Running AutoDock Vina… (this takes ~30–60s)"):
        try:
            from utils.docking import run_docking, read_pdbqt
            out_pdbqt = str(Path(tempfile.mkdtemp()) / "docked.pdbqt")
            energies = run_docking(lig_pdbqt, out_pdbqt, target=target,
                                   exhaustiveness=exhaustiveness, n_poses=n_poses)
            docked_str = read_pdbqt(out_pdbqt)
            st.session_state["docked"] = {
                "energies": energies,
                "pdbqt": docked_str,
                "lig_pdbqt": lig_pdbqt,
                "out_pdbqt": out_pdbqt,
                "name": mol_name,
                "smiles": smiles,
                "target": target,
            }
        except Exception as e:
            st.error(f"Docking failed: {e}")
            st.stop()

# ── Results ────────────────────────────────────────────────────────────────────
if "docked" not in st.session_state:
    st.stop()

d = st.session_state["docked"]
energies = d["energies"]
best = energies[0]

st.subheader(f"Results — {d['name']}")

col_a, col_b, col_c = st.columns(3)
col_a.metric("Best binding energy", f"{best} kcal/mol",
             help="More negative = stronger binding")
col_b.metric("Poses found", len(energies))
col_c.metric("Target", d["target"].split("—")[0].strip())

tab1, tab2, tab3 = st.tabs(["🧬 3D Viewer", "📊 Score Comparison", "🔗 Interactions"])

# ── Tab 1: 3D Viewer ───────────────────────────────────────────────────────────
with tab1:
    st.caption("Protein = light blue cartoon · Your molecule = green sticks")
    try:
        protein_pdb_path = TARGET_PROTEIN_PDB[d["target"]]
        with open(protein_pdb_path) as f:
            protein_str = f.read()

        viewer_html = f"""
        <html><head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
        <style>body{{margin:0;background:#0f1320;}}#viewer{{width:100%;height:520px;}}</style>
        </head><body>
        <div id="viewer"></div>
        <script>
        let viewer = $3Dmol.createViewer('viewer', {{backgroundColor:'#0f1320'}});
        viewer.addModel(`{protein_str.replace("`","'")}`, 'pdb');
        viewer.setStyle({{model:0}}, {{cartoon:{{color:'#90cdf4', opacity:0.85}}}});
        viewer.addModel(`{d["pdbqt"].replace("`","'")}`, 'pdbqt');
        viewer.setStyle({{model:1}}, {{stick:{{colorscheme:'greenCarbon', radius:0.25}}}});
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:0.12, color:'#a78bfa'}},
                          {{model:0}});
        viewer.zoomTo({{model:1}});
        viewer.render();
        </script>
        </body></html>
        """
        components.html(viewer_html, height=530)
    except Exception as e:
        st.error(f"3D viewer error: {e}")

# ── Tab 2: Score comparison ────────────────────────────────────────────────────
with tab2:
    st.caption("Add more molecules via sidebar to compare scores.")

    score_key = f"all_scores_{d['target']}"
    if score_key not in st.session_state:
        st.session_state[score_key] = {}
    st.session_state[score_key][d["name"]] = best

    try:
        from utils.analysis import score_chart
        img = score_chart(st.session_state[score_key])
        st.image(img, use_container_width=True)
    except Exception as e:
        st.error(f"Chart error: {e}")

    st.subheader("All poses")
    import pandas as pd
    df = pd.DataFrame({
        "Pose": [f"Pose {i+1}" for i in range(len(energies))],
        "Binding energy (kcal/mol)": energies
    })
    st.dataframe(df, use_container_width=True, hide_index=True)

# ── Tab 3: Interactions ────────────────────────────────────────────────────────
with tab3:
    st.caption("Hydrogen bonds, hydrophobic contacts, π-stacking and more")

    # 3D viewer
    try:
        protein_pdb_path = TARGET_PROTEIN_PDB[d["target"]]
        with open(protein_pdb_path) as f:
            protein_str = f.read()
        viewer_html = f"""
        <html><head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
        <style>body{{margin:0;background:#0f1320;}}#viewer{{width:100%;height:420px;}}</style>
        </head><body>
        <div id="viewer"></div>
        <script>
        let viewer = $3Dmol.createViewer('viewer', {{backgroundColor:'#0f1320'}});
        viewer.addModel(`{protein_str.replace("`","'")}`, 'pdb');
        viewer.setStyle({{model:0}}, {{cartoon:{{color:'#90cdf4', opacity:0.85}}}});
        viewer.addModel(`{d["pdbqt"].replace("`","'")}`, 'pdbqt');
        viewer.setStyle({{model:1}}, {{stick:{{colorscheme:'greenCarbon', radius:0.3}}}});
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:0.12, color:'#a78bfa'}}, {{model:0}});
        viewer.zoomTo({{model:1}});
        viewer.render();
        </script>
        </body></html>
        """
        components.html(viewer_html, height=430)
    except Exception as e:
        st.error(f"3D viewer error: {e}")

    st.divider()

    # Interaction diagram
    try:
        from utils.analysis import get_interactions, interaction_chart
        h_pdb = TARGET_PROTEIN_H_PDB[d["target"]]
        with st.spinner("Calculating interactions…"):
            idf = get_interactions(h_pdb, d["out_pdbqt"], smiles=d["smiles"])
        if idf is not None and not idf.empty:
            img = interaction_chart(idf)
            st.image(img, use_container_width=True)
            with st.expander(f"Raw data ({len(idf)} interactions)"):
                st.dataframe(idf, use_container_width=True, hide_index=True)
        else:
            st.info("No interactions detected in the docked pose.")
    except Exception as e:
        st.error(f"Interaction analysis failed: {e}")
