from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


def smiles_to_pdbqt(smiles: str, out_path: str) -> bool:
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) == -1:
        raise RuntimeError("3D embedding failed")
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_setups[0])

    with open(out_path, "w") as f:
        f.write(pdbqt_string)
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
