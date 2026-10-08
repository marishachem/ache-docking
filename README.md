# Protein–Ligand Docking Lab

Computational docking of small molecules into multiple drug targets using **AutoDock Vina** — with an interactive Streamlit UI for visualization, interaction analysis, and selectivity evaluation.

The main molecule is a **piperidine-2-one derivative** synthesised during a chemistry degree (2013). This project explores whether it could bind to the same active sites as approved drugs — computationally.

## Targets

| Protein | PDB | Disease / Role |
|---------|-----|---------------|
| **AChE** — Acetylcholinesterase | 1EVE | Alzheimer's — primary target |
| **BuChE** — Butyrylcholinesterase | 4BDS | AChE closest homolog — selectivity test |
| **EGFR** — Epidermal Growth Factor Receptor | 1M17 | Lung/breast cancer — primary target |
| **VEGFR2** — Vascular Endothelial Growth Factor Receptor | 4ASD | Kinase off-target for EGFR |
| **COX-2** — Cyclooxygenase-2 | 5IKT | Common off-target |

## Features

- **Two primary targets** — switch between AChE and EGFR in the sidebar
- **Name → SMILES converter** — type any compound name (common, IUPAC, brand) and get SMILES via PubChem → NCI CIR → OPSIN fallback chain
- **AutoDock Vina** docking with live progress bar
- **3D interactive viewer** — protein colored by secondary structure, ligand in green (3Dmol.js)
- **Score comparison chart** — your molecule vs. approved drugs for the selected target
- **Score interpretation** — color-coded rating (🟢🟡🟠🔴) with reference scale and benchmark comparison
- **Interaction viewer** — interacting residues highlighted in 3D by interaction type with residue labels
- **Interaction dot-matrix chart** — residues × interaction types, color-coded (ProLIF)
- **🎯 Selectivity Panel** — dock the same molecule against multiple proteins simultaneously, get a selectivity verdict (✅ Selective / ⚠️ Borderline / ❌ Not selective)
- **Molecular properties** — MW, LogP, TPSA, Lipinski rule check
- **Custom SMILES input** — dock any molecule

## Getting Started

```bash
git clone https://github.com/marishachem/ache-docking
cd ache-docking

# Install dependencies (conda recommended for AutoDock Vina)
conda create -n docking python=3.11
conda activate docking
conda install -c conda-forge vina
pip install -r requirements.txt

# One-time protein setup (downloads PDBs, prepares PDBQT)
python3 setup_protein.py

# Run the app
streamlit run app.py
```

Opens at `http://localhost:8501`

## Score Reference

| Score (kcal/mol) | Rating |
|-----------------|--------|
| < −10 | 🟢 Excellent |
| −8 to −10 | 🟢 Very good — comparable to approved drugs |
| −6 to −8 | 🟡 Moderate — promising hit |
| −4 to −6 | 🟠 Weak |
| > −4 | 🔴 Poor |

Scores are meaningful for **comparing molecules on the same target**. They are not experimental IC₅₀ or Kd values.

## Selectivity

A molecule is considered **selective** when it binds its primary target **2+ kcal/mol stronger** than any off-target. The selectivity panel docks against all relevant proteins and gives an automatic verdict.

## Preset Molecules

**AChE target:**
| Molecule | Role |
|----------|------|
| Piperidinone derivative | Diploma molecule (synthesised 2013) |
| Donepezil | Approved AChE inhibitor (reference, −10.49) |
| Rivastigmine | Approved AChE inhibitor |
| Galantamine | Approved AChE inhibitor |

**EGFR target:**
| Molecule | Generation | Score |
|----------|-----------|-------|
| Gefitinib | 1st gen | −8.07 |
| Erlotinib | 1st gen | −7.66 |
| Lapatinib | 2nd gen | −8.40 |
| Osimertinib | 3rd gen | −8.11 |

Any custom molecule can be docked by name or SMILES.

## Project Structure

```
ache-docking/
├── app.py                  # Streamlit UI
├── setup_protein.py        # One-time protein preparation
├── utils/
│   ├── ligand.py           # SMILES → PDBQT (meeko + RDKit)
│   ├── docking.py          # AutoDock Vina wrapper (multi-target)
│   └── analysis.py         # Score chart, interaction chart, selectivity chart
├── data/                   # PDB and PDBQT files
└── requirements.txt
```

## Tech

Python · Streamlit · RDKit · AutoDock Vina · meeko · 3Dmol.js · BioPython · ProLIF · MDAnalysis · matplotlib · pdbfixer · OpenMM
