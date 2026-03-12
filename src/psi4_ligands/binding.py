"""Binding energy utilities using Psi4."""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

from Bio.PDB import PDBParser
from rdkit import Chem


def _guess_element(name: str) -> str:
    name = name.strip()
    if not name:
        return "X"
    if len(name) >= 2 and name[:2].isalpha():
        candidate = name[:2].title()
        if candidate in {
            "Cl",
            "Br",
            "Na",
            "Mg",
            "Al",
            "Si",
            "Ca",
            "Li",
            "Zn",
            "Fe",
            "Cu",
            "Mn",
            "Co",
            "Ni",
        }:
            return candidate
    return name[0].upper()


@dataclass(frozen=True)
class AtomCoord:
    element: str
    x: float
    y: float
    z: float


def _atom_sort_key(atom) -> tuple:
    chain_id = atom.get_parent().get_parent().id
    residue = atom.get_parent()
    resseq = residue.get_id()[1]
    icode = residue.get_id()[2]
    serial = getattr(atom, "serial_number", 0)
    return (chain_id, resseq, icode, atom.get_name(), atom.get_altloc(), serial)


def load_protein_coords(pdb_path: str) -> List[AtomCoord]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    coords_with_keys = []
    for atom in sorted(structure.get_atoms(), key=_atom_sort_key):
        element = atom.element.strip() if atom.element else ""
        if not element:
            element = _guess_element(atom.get_name())
        x, y, z = atom.get_coord()
        coords_with_keys.append(AtomCoord(element, float(x), float(y), float(z)))
    return coords_with_keys


def load_ligand_coords(sdf_path: str) -> List[AtomCoord]:
    supplier = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError(f"No valid molecule found in {sdf_path}")
    if mol.GetNumConformers() == 0:
        raise ValueError(f"Ligand in {sdf_path} has no conformers")
    conformer = mol.GetConformer()
    coords: List[AtomCoord] = []
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        coords.append(AtomCoord(atom.GetSymbol(), pos.x, pos.y, pos.z))
    return coords


def _coords_to_psi4_geometry(coords: Iterable[AtomCoord], charge: int, multiplicity: int) -> str:
    lines = [f"{charge} {multiplicity}"]
    for atom in coords:
        lines.append(f"{atom.element} {atom.x:.8f} {atom.y:.8f} {atom.z:.8f}")
    return "\n".join(lines)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _candidate_psi4_paths() -> list[Path]:
    env_path = os.environ.get("PSI4_PYTHONPATH")
    candidates: list[Path] = []
    if env_path:
        candidates.extend(Path(part) for part in env_path.split(os.pathsep) if part)

    root = _repo_root()
    candidates.extend(
        [
            root / "local" / "psi4" / "lib",
            root / "build" / "psi4" / "stage" / "lib",
        ]
    )
    return candidates


def _ensure_psi4_importable() -> None:
    attempted: list[str] = []
    for candidate in _candidate_psi4_paths():
        if not candidate.exists():
            continue
        candidate_str = str(candidate)
        attempted.append(candidate_str)
        if candidate_str not in sys.path:
            sys.path.insert(0, candidate_str)
        try:
            __import__("psi4")
            return
        except ImportError:
            continue

    attempted_paths = ", ".join(attempted) if attempted else "no local Psi4 paths found"
    raise ImportError(
        "Psi4 is not importable. Build the vendored Psi4 with "
        "`./scripts/bootstrap.sh` or set PSI4_PYTHONPATH. "
        f"Searched: {attempted_paths}"
    )


def binding_energy(
    protein_coords: Iterable[AtomCoord],
    ligand_coords: Iterable[AtomCoord],
    method: str,
    basis: str,
    charge: int = 0,
    multiplicity: int = 1,
    threads: int = 1,
) -> Tuple[float, float, float, float]:
    _ensure_psi4_importable()
    import psi4

    protein_list = list(protein_coords)
    ligand_list = list(ligand_coords)

    protein_geom = _coords_to_psi4_geometry(protein_list, charge, multiplicity)
    ligand_geom = _coords_to_psi4_geometry(ligand_list, charge, multiplicity)
    complex_geom = _coords_to_psi4_geometry(
        protein_list + ligand_list, charge, multiplicity
    )

    level = method if "/" in method else f"{method}/{basis}"

    psi4.set_num_threads(threads)

    protein_mol = psi4.geometry(protein_geom)
    ligand_mol = psi4.geometry(ligand_geom)
    complex_mol = psi4.geometry(complex_geom)

    e_protein = psi4.energy(level, molecule=protein_mol)
    e_ligand = psi4.energy(level, molecule=ligand_mol)
    e_complex = psi4.energy(level, molecule=complex_mol)

    delta = e_complex - (e_protein + e_ligand)
    return e_complex, e_protein, e_ligand, delta
