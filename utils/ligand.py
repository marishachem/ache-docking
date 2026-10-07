import subprocess
import tempfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def smiles_to_pdbqt(smiles: str, out_path: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) == -1:
        raise RuntimeError("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol)

    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        sdf_path = f.name
    Chem.MolToMolFile(mol, sdf_path)

    result = subprocess.run(
        ["mk_prepare_ligand.py", "-i", sdf_path, "-o", out_path],
        capture_output=True, text=True
    )
    Path(sdf_path).unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"mk_prepare_ligand.py failed:\n{result.stderr}")
    return True


def get_mol_props(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "MW": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "RotBonds": Descriptors.NumRotatableBonds(mol),
    }
