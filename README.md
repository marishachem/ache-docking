# Protein–Ligand Docking Lab

Computational docking of small molecules into two drug targets using **AutoDock Vina** — with an interactive Streamlit UI for visualization and analysis.

The main molecule is a **piperidine-2-one derivative** synthesised during a chemistry degree (2013). This project explores whether it could bind to the same active sites as approved drugs — computationally.

## Targets

| Protein | PDB | Disease | Reference drugs |
|---------|-----|---------|-----------------|
| **AChE** — Acetylcholinesterase | 1EVE | Alzheimer's | Donepezil (−10.49), Rivastigmine (−7.74), Galantamine (−7.70) |
| **EGFR** — Epidermal Growth Factor Receptor | 1M17 | Cancer (lung, breast) | Gefitinib (−8.07), Erlotinib (−7.66), Lapatinib (−8.40), Osimertinib (−8.11) |

Scores in kcal/mol — more negative = stronger binding.

## Features

- **Two protein targets** — switch between AChE and EGFR in the sidebar
- **AutoDock Vina** docking via Python API
- **3D interactive viewer** — full protein cartoon + ligand sticks, colored by secondary structure (3Dmol.js)
- **Score comparison chart** — your molecule vs. approved drugs for the selected target
- **Score interpretation panel** — color-coded rating (🟢🟡🟠🔴) with reference scale and benchmark comparison
- **Interaction viewer** — interacting residues highlighted in 3D by interaction type (H-bond, hydrophobic, π-stacking…) with labels
- **Interaction dot-matrix chart** — residues × interaction types, color-coded (ProLIF)
- **Molecular properties** — MW, LogP, TPSA, Lipinski rule check
- **Streamlit UI** — input any SMILES, get results in one click

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

## Project Structure

```
ache-docking/
├── app.py                  # Streamlit UI
├── setup_protein.py        # One-time protein preparation
├── utils/
│   ├── ligand.py           # SMILES → PDBQT (meeko + RDKit)
│   ├── docking.py          # AutoDock Vina wrapper (multi-target)
│   └── analysis.py         # Score chart, ProLIF interactions, dot-matrix chart
├── data/                   # Generated files (PDB, PDBQT)
└── requirements.txt
```

## Preset Molecules

**AChE target:**
| Molecule | Role |
|----------|------|
| Piperidinone derivative | Diploma molecule (synthesised 2013) |
| Donepezil | Approved AChE inhibitor (reference) |
| Rivastigmine | Approved AChE inhibitor |
| Galantamine | Approved AChE inhibitor |

**EGFR target:**
| Molecule | Generation |
|----------|-----------|
| Gefitinib | 1st gen EGFR inhibitor |
| Erlotinib | 1st gen EGFR inhibitor |
| Lapatinib | 2nd gen EGFR inhibitor |
| Osimertinib | 3rd gen EGFR inhibitor |

Any custom molecule can be docked by pasting its SMILES into the sidebar.

## Tech

Python · Streamlit · RDKit · AutoDock Vina · meeko · 3Dmol.js · BioPython · ProLIF · MDAnalysis · matplotlib · pdbfixer · OpenMM
