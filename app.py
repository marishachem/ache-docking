import streamlit as st
import streamlit.components.v1 as components
import tempfile, os, re
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

# ── Name → SMILES helpers ──────────────────────────────────────────────────────
import urllib.request, urllib.parse, urllib.error, json as _json

def _try_pubchem(name):
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           f"{urllib.parse.quote(name)}/property/"
           f"IsomericSMILES,IUPACName,MolecularFormula,MolecularWeight/JSON")
    with urllib.request.urlopen(url, timeout=8) as r:
        d = _json.loads(r.read())["PropertyTable"]["Properties"][0]
    return d.get("IsomericSMILES",""), d.get("IUPACName", name), d.get("MolecularFormula",""), d.get("MolecularWeight","")

def _try_cir(name):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{urllib.parse.quote(name)}/smiles"
    with urllib.request.urlopen(url, timeout=8) as r:
        return r.read().decode().strip(), name, "", ""

def _try_opsin(name):
    url = f"https://opsin.ch.cam.ac.uk/opsin/{urllib.parse.quote(name)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=10) as r:
        data = _json.loads(r.read())
    return data.get("smiles",""), name, "", ""

def name_to_smiles(name):
    from rdkit import Chem
    for fn, label in [(_try_pubchem, "PubChem"), (_try_cir, "NCI CIR"), (_try_opsin, "OPSIN")]:
        try:
            smi, iupac, formula, mw = fn(name)
            if smi and Chem.MolFromSmiles(smi):
                return smi, iupac, formula, mw, label
        except Exception:
            continue
    return None, None, None, None, None

# ── Sidebar ────────────────────────────────────────────────────────────────────
with st.sidebar:
    st.title("🔬 Docking Lab")
    st.divider()

    target = st.selectbox("Target protein", list(TARGET_MOLECULES.keys()))
    molecules = TARGET_MOLECULES[target]

    st.divider()
    mode = st.radio("Input mode", ["Preset molecules", "Search by name", "Custom SMILES"])

    if mode == "Preset molecules":
        selected = st.selectbox("Choose molecule", list(molecules.keys()))
        smiles = molecules[selected]
        mol_name = selected.split("(")[0].strip().lstrip("💊🎓").strip()

    elif mode == "Search by name":
        name_query = st.text_input("Compound name", placeholder="e.g. ibuprofen, erlotinib…")
        smiles = ""
        mol_name = name_query or "Custom molecule"
        if name_query:
            with st.spinner("Searching…"):
                smi, iupac, formula, mw, source = name_to_smiles(name_query)
            if smi:
                smiles = smi
                mol_name = name_query
                st.success(f"Found via **{source}**")
                if formula and mw:
                    st.caption(f"{formula} · {mw} g/mol")
            else:
                st.error("Not found. Try the Name→SMILES tool below for complex IUPAC names.")

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

# ── Name → SMILES converter ────────────────────────────────────────────────────
with st.expander("🔤 Name → SMILES converter", expanded=False):
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import Descriptors
    import base64

    lookup_name = st.text_input(
        "Compound name",
        placeholder="e.g. aspirin, erlotinib, caffeine, N-[(4'-methoxy)phenyl]-3-phenyl-…",
        key="lookup_input"
    )
    if lookup_name:
        with st.spinner("Searching PubChem → NCI CIR → OPSIN…"):
            found_smi, found_iupac, found_formula, found_mw, found_source = name_to_smiles(lookup_name)

        if found_smi:
            mol = Chem.MolFromSmiles(found_smi)
            mw_calc = round(Descriptors.MolWt(mol), 2)
            st.success(f"Found via **{found_source}**")

            rc1, rc2 = st.columns([1, 1])
            with rc1:
                drawer = rdMolDraw2D.MolDraw2DSVG(340, 260)
                drawer.drawOptions().padding = 0.15
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                b64 = base64.b64encode(svg.encode()).decode()
                st.markdown(
                    f'<img src="data:image/svg+xml;base64,{b64}" '
                    f'style="width:100%;border-radius:10px;border:1px solid #e2e8f0;">',
                    unsafe_allow_html=True
                )
            with rc2:
                st.markdown("**SMILES**")
                st.code(found_smi, language=None)
                if found_formula:
                    st.markdown(f"**Formula:** {found_formula}")
                st.markdown(f"**MW:** {found_mw or mw_calc} g/mol")
                if found_iupac and found_iupac != lookup_name:
                    st.markdown(f"**IUPAC:** {found_iupac}")
                st.caption("Copy the SMILES above, then use **Custom SMILES** mode in the sidebar to dock it.")
        else:
            st.error("Not found in PubChem, NCI CIR or OPSIN. Try a slightly different name or spelling.")
    else:
        st.caption("Searches PubChem → NCI CIR → OPSIN in order. Works with common names, drug names, brand names, and IUPAC names.")

st.divider()

if not smiles:
    st.info("Select a molecule in the sidebar and click **Run Docking**.")
    st.stop()

# ── 2D structure + properties ──────────────────────────────────────────────────
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import base64

col1, col2 = st.columns([1, 2])
with col1:
    st.subheader("2D Structure")
    try:
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
col_a.metric("Best binding energy", f"{best} kcal/mol", help="More negative = stronger binding")
col_b.metric("Poses found", len(energies))
col_c.metric("Target", d["target"].split("—")[0].strip())

# ── Score interpretation ───────────────────────────────────────────────────────
BENCHMARKS = {
    "AChE — Acetylcholinesterase (Alzheimer's)": [
        ("Donepezil", -10.49), ("Rivastigmine", -7.74), ("Galantamine", -7.70),
    ],
    "EGFR — Epidermal Growth Factor Receptor (Cancer)": [
        ("Lapatinib", -8.40), ("Osimertinib", -8.11),
        ("Gefitinib", -8.07), ("Erlotinib", -7.66),
    ],
}

def score_label(score):
    if score < -10:  return "🟢 Excellent", "Stronger than most approved drugs"
    elif score < -8: return "🟢 Very good",  "Comparable to approved drugs"
    elif score < -6: return "🟡 Moderate",   "Promising hit, worth optimizing"
    elif score < -4: return "🟠 Weak",        "Some affinity, unlikely to be active"
    else:            return "🔴 Poor",         "Essentially no binding"

label, desc = score_label(best)
benchmarks = BENCHMARKS.get(d["target"], [])
worst_ref = max(benchmarks, key=lambda x: x[1])[1] if benchmarks else None

with st.expander("📖 How to read this score", expanded=True):
    ic1, ic2 = st.columns([1, 2])
    with ic1:
        st.markdown(f"**Your score:** `{best} kcal/mol`")
        st.markdown(f"**Rating:** {label}")
        st.markdown(f"*{desc}*")
        if benchmarks and worst_ref:
            diff = best - worst_ref
            sign = "better" if diff < 0 else "weaker"
            st.markdown(f"**vs. weakest ref drug ({worst_ref}):** {abs(diff):.2f} kcal/mol {sign}")
    with ic2:
        st.markdown("**Reference scale:**")
        for rng, dot, meaning in [
            ("< −10",     "🟢", "Excellent — stronger than most approved drugs"),
            ("−8 to −10", "🟢", "Very good — comparable to approved drugs"),
            ("−6 to −8",  "🟡", "Moderate — promising, worth optimizing"),
            ("−4 to −6",  "🟠", "Weak — unlikely to be active"),
            ("> −4",      "🔴", "Poor — essentially no binding"),
        ]:
            here = (
                (rng == "< −10"     and best < -10) or
                (rng == "−8 to −10" and -10 <= best < -8) or
                (rng == "−6 to −8"  and -8  <= best < -6) or
                (rng == "−4 to −6"  and -6  <= best < -4) or
                (rng == "> −4"      and best >= -4)
            )
            st.markdown(f"{dot} `{rng}` — {meaning} {'◀ **your score**' if here else ''}")
        if benchmarks:
            st.markdown("**Approved drugs on this target:**")
            for name, score in benchmarks:
                st.markdown(f"- {name}: `{score}` kcal/mol")

tab1, tab2, tab3 = st.tabs(["🧬 3D Viewer", "📊 Score Comparison", "🔗 Interactions"])

# ── Tab 1: 3D Viewer ───────────────────────────────────────────────────────────
with tab1:
    st.caption("Protein colored by secondary structure · Ligand = green sticks")
    try:
        protein_pdb_path = TARGET_PROTEIN_PDB[d["target"]]
        with open(protein_pdb_path) as f:
            protein_str = f.read()
        viewer_html = f"""
        <html><head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
        <style>body{{margin:0;background:#e8edf5;}}#v1{{width:100%;height:520px;}}</style>
        </head><body><div id="v1"></div>
        <script>
        let viewer = $3Dmol.createViewer('v1', {{backgroundColor:'#e8edf5'}});
        viewer.addModel(`{protein_str.replace("`","'")}`, 'pdb');
        viewer.setStyle({{model:0}}, {{cartoon:{{color:'spectrum'}}}});
        viewer.addModel(`{d["pdbqt"].replace("`","'")}`, 'pdbqt');
        viewer.setStyle({{model:1}}, {{stick:{{colorscheme:'greenCarbon', radius:0.3}}}});
        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity:0.08, color:'white'}}, {{model:0}});
        viewer.zoomTo({{model:1}});
        viewer.render();
        </script></body></html>
        """
        components.html(viewer_html, height=530)
    except Exception as e:
        st.error(f"3D viewer error: {e}")

# ── Tab 2: Score comparison ────────────────────────────────────────────────────
with tab2:
    st.caption("Dock more molecules to add bars to the chart.")
    score_key = f"all_scores_{d['target']}"
    if score_key not in st.session_state:
        st.session_state[score_key] = {}
    st.session_state[score_key][d["name"]] = best

    try:
        from utils.analysis import score_chart
        st.image(score_chart(st.session_state[score_key]), use_container_width=True)
    except Exception as e:
        st.error(f"Chart error: {e}")

    st.subheader("All poses")
    import pandas as pd
    st.dataframe(pd.DataFrame({
        "Pose": [f"Pose {i+1}" for i in range(len(energies))],
        "Binding energy (kcal/mol)": energies
    }), use_container_width=True, hide_index=True)

# ── Tab 3: Interactions ────────────────────────────────────────────────────────
with tab3:
    st.caption("Hydrogen bonds, hydrophobic contacts, π-stacking and more")
    try:
        from utils.analysis import get_interactions, interaction_chart

        h_pdb = TARGET_PROTEIN_H_PDB[d["target"]]
        with st.spinner("Calculating interactions…"):
            idf = get_interactions(h_pdb, d["out_pdbqt"], smiles=d["smiles"])

        protein_pdb_path = TARGET_PROTEIN_PDB[d["target"]]
        with open(protein_pdb_path) as f:
            protein_str = f.read()

        itype_colors = {
            "HBDonor": "#60a5fa", "HBAcceptor": "#34d399",
            "Hydrophobic": "#fbbf24", "PiStacking": "#a78bfa",
            "PiCation": "#f472b6", "CationPi": "#f472b6",
            "VdWContact": "#94a3b8", "EdgeToFace": "#c084fc",
        }
        resi_js = "[]"
        if idf is not None and not idf.empty:
            resi_list = []
            for _, row in idf.iterrows():
                m = re.match(r'[A-Z]+(\d+)\.(\w+)', row["Residue"])
                if m:
                    color = itype_colors.get(row["Interaction"], "#fb923c")
                    resi_list.append(
                        f'{{resi:{m.group(1)},chain:"{m.group(2)}",'
                        f'color:"{color}",label:"{row["Residue"]}"}}'
                    )
            resi_js = "[" + ",".join(resi_list) + "]"

        viewer_html3 = f"""
        <html><head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
        <style>body{{margin:0;background:#f1f5f9;}}#v3{{width:100%;height:440px;}}</style>
        </head><body><div id="v3"></div>
        <script>
        setTimeout(function() {{
            let viewer = $3Dmol.createViewer('v3', {{backgroundColor:'#f1f5f9'}});
            viewer.addModel(`{protein_str.replace("`","'")}`, 'pdb');
            viewer.setStyle({{model:0}}, {{cartoon:{{color:'#c8d0dc', opacity:0.5}}}});
            let residues = {resi_js};
            residues.forEach(r => {{
                viewer.addStyle({{model:0, resi:r.resi, chain:r.chain}},
                    {{stick:{{color:r.color, radius:0.22}},
                     cartoon:{{color:r.color, opacity:0.9}}}});
                viewer.addLabel(r.label, {{
                    position:{{resi:r.resi, chain:r.chain}},
                    backgroundColor:'rgba(0,0,0,0.55)', fontColor:'white',
                    fontSize:10, borderRadius:3
                }});
            }});
            viewer.addModel(`{d["pdbqt"].replace("`","'")}`, 'pdbqt');
            viewer.setStyle({{model:1}}, {{stick:{{colorscheme:'greenCarbon', radius:0.32}}}});
            viewer.zoomTo({{model:1}});
            viewer.render();
        }}, 150);
        </script></body></html>
        """
        components.html(viewer_html3, height=450)

        st.divider()

        if idf is not None and not idf.empty:
            st.image(interaction_chart(idf), use_container_width=True)
            with st.expander(f"Raw data ({len(idf)} interactions)"):
                st.dataframe(idf, use_container_width=True, hide_index=True)
        else:
            st.info("No interactions detected in the docked pose.")

    except Exception as e:
        st.error(f"Interaction analysis failed: {e}")
