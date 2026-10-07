# AChE Docking

Computational docking of small molecules into **Acetylcholinesterase** (AChE, PDB: 1EVE) — the Alzheimer's disease drug target.

The main molecule is a **piperidine-2-one derivative** synthesised during a chemistry degree, structurally related to donepezil (the approved AChE inhibitor). This project explores whether it could bind to the same active site — computationally.

## Features

- **AutoDock Vina** docking via Python API
- **3D interactive viewer** — protein cartoon + ligand sticks + surface (3Dmol.js)
- **Score comparison chart** — your molecule vs. donepezil, rivastigmine, galantamine
- **Protein–ligand interaction table** — hydrogen bonds, hydrophobic contacts (ProLIF)
- **Molecular properties** — MW, LogP, TPSA, Lipinski rule check
- **Streamlit UI** — input any SMILES, get results in one click

## Getting Started

```bash
git clone https://github.com/marishachem/ache-docking
cd ache-docking
pip install -r requirements.txt

# One-time protein setup (downloads 1EVE, prepares PDBQT)
python3 setup_protein.py

# Run the app
streamlit run app.py
```

Opens at `http://localhost:8501`

## Project Structure

```
ache-docking/
├── app.py                  # Streamlit UI
├── setup_protein.py        # One-time protein preparation
├── utils/
│   ├── ligand.py           # SMILES → PDBQT
│   ├── docking.py          # AutoDock Vina wrapper
│   └── analysis.py         # Score chart + ProLIF interactions
├── data/                   # Generated files (PDB, PDBQT)
└── requirements.txt
```

## Molecules Included

| Molecule | Role |
|----------|------|
| Piperidinone derivative | Diploma molecule (synthesised 2013) |
| Donepezil | Approved AChE inhibitor (reference) |
| Rivastigmine | Approved AChE inhibitor |
| Galantamine | Approved AChE inhibitor |

## Requirements

- Python 3.10+
- AutoDock Vina (installed via `pip install vina`)
- meeko (ligand/receptor PDBQT preparation)
- pdbfixer + OpenMM (protein hydrogens)
- ProLIF + MDAnalysis (interactions, optional)

## Tech

Python · Streamlit · RDKit · AutoDock Vina · py3Dmol · BioPython · ProLIF · matplotlib
