"""
Run once before using the app:
    python3 setup_protein.py

Downloads 1EVE, cleans it, adds hydrogens, converts to PDBQT.
"""

import urllib.request
import subprocess
import sys
from pathlib import Path

DATA = Path("data")
DATA.mkdir(exist_ok=True)


def download_pdb():
    out = DATA / "1EVE.pdb"
    if out.exists():
        print("1EVE.pdb already exists, skipping download.")
        return
    print("Downloading 1EVE from RCSB...")
    urllib.request.urlretrieve(
        "https://files.rcsb.org/download/1EVE.pdb", out
    )
    print(f"Saved to {out}")


def clean_protein():
    from Bio.PDB import PDBParser, PDBIO, Select

    class ProteinOnly(Select):
        def accept_residue(self, residue):
            return residue.id[0] == " "

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("AChE", DATA / "1EVE.pdb")
    io = PDBIO()
    io.set_structure(structure)
    out = DATA / "protein_clean.pdb"
    io.save(str(out), ProteinOnly())
    print(f"Cleaned protein saved to {out}")


def add_hydrogens():
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    out = DATA / "protein_ready.pdb"
    if out.exists():
        print("protein_ready.pdb already exists, skipping.")
        return
    print("Adding missing atoms and hydrogens...")
    fixer = PDBFixer(filename=str(DATA / "protein_clean.pdb"))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    with open(out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"Saved to {out}")


def convert_to_pdbqt():
    out = DATA / "protein_ready.pdbqt"
    if out.exists():
        print("protein_ready.pdbqt already exists, skipping.")
        return
    print("Converting protein to PDBQT...")
    import shutil
    mk = shutil.which("mk_prepare_receptor.py") or \
         "/opt/homebrew/Caskroom/miniforge/base/envs/docking/bin/mk_prepare_receptor.py"
    result = subprocess.run(
        [mk, "-i", str(DATA / "protein_clean.pdb"),
         "-a", "--compute_charges", "-p", str(out)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("mk_prepare_receptor.py failed, trying obabel fallback...")
        result2 = subprocess.run(
            ["obabel", str(DATA / "protein_ready.pdb"), "-O", str(out), "-xr"],
            capture_output=True, text=True
        )
        if result2.returncode != 0:
            print("ERROR: Could not convert to PDBQT.")
            print("Install openbabel: brew install open-babel")
            sys.exit(1)
    print(f"Saved to {out}")


if __name__ == "__main__":
    download_pdb()
    clean_protein()
    add_hydrogens()
    convert_to_pdbqt()
    print("\nProtein setup complete. You can now run: streamlit run app.py")
