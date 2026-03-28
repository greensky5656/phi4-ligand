"""Binding energy utilities using Psi4 and xTB."""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from Bio.PDB import PDBParser
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


HARTREE_TO_KCAL_MOL = 627.509474


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


@dataclass(frozen=True)
class ResidueRef:
    chain_id: str
    hetflag: str
    resseq: int
    icode: str
    resname: str


@dataclass(frozen=True)
class BoxSelection:
    residues: tuple[ResidueRef, ...]
    atoms: tuple[AtomCoord, ...]
    cap_atoms: tuple[AtomCoord, ...]


def _atom_sort_key(atom) -> tuple:
    chain_id = atom.get_parent().get_parent().id
    residue = atom.get_parent()
    hetflag = residue.get_id()[0]
    resseq = residue.get_id()[1]
    icode = residue.get_id()[2]
    serial = getattr(atom, "serial_number", 0)
    return (chain_id, hetflag, resseq, icode, atom.get_name(), atom.get_altloc(), serial)


def _residue_sort_key(residue) -> tuple:
    chain_id = residue.get_parent().id
    hetflag, resseq, icode = residue.get_id()
    return (chain_id, hetflag, resseq, icode, residue.resname)


def _atomcoord_from_biopython(atom) -> AtomCoord:
    element = atom.element.strip() if atom.element else ""
    if not element:
        element = _guess_element(atom.get_name())
    x, y, z = atom.get_coord()
    return AtomCoord(element, float(x), float(y), float(z))


def _find_atom(residue, atom_name: str):
    for atom in residue.get_atoms():
        if atom.get_name().strip() == atom_name:
            return atom
    return None


def _normalized(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    if norm < 1.0e-8:
        return None
    return vector / norm


def _orthonormal_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, reference))) > 0.8:
        reference = np.array([0.0, 1.0, 0.0])
    first = np.cross(axis, reference)
    first = first / np.linalg.norm(first)
    second = np.cross(axis, first)
    second = second / np.linalg.norm(second)
    return first, second


def _tetrahedral_methyl_hydrogens(
    carbon_position: np.ndarray,
    bonded_position: np.ndarray,
    bond_length: float = 1.09,
) -> list[AtomCoord]:
    axis = _normalized(bonded_position - carbon_position)
    if axis is None:
        return []
    first, second = _orthonormal_basis(axis)
    lateral = (2.0 * np.sqrt(2.0)) / 3.0
    hydrogen_positions: list[AtomCoord] = []
    for angle_deg in (0.0, 120.0, 240.0):
        angle = np.deg2rad(angle_deg)
        direction = (
            (-1.0 / 3.0) * axis
            + lateral * (np.cos(angle) * first + np.sin(angle) * second)
        )
        position = carbon_position + direction * bond_length
        hydrogen_positions.append(
            AtomCoord("H", float(position[0]), float(position[1]), float(position[2]))
        )
    return hydrogen_positions


def _amide_hydrogen(
    nitrogen_position: np.ndarray,
    carbonyl_carbon_position: np.ndarray,
    methyl_carbon_position: np.ndarray,
    bond_length: float = 1.01,
) -> list[AtomCoord]:
    direction = _normalized(
        2.0 * nitrogen_position - carbonyl_carbon_position - methyl_carbon_position
    )
    if direction is None:
        return []
    position = nitrogen_position + direction * bond_length
    return [AtomCoord("H", float(position[0]), float(position[1]), float(position[2]))]


def _ace_cap(residue, previous_residue) -> list[AtomCoord]:
    previous_carbon = _find_atom(previous_residue, "C")
    previous_oxygen = _find_atom(previous_residue, "O")
    previous_alpha = _find_atom(previous_residue, "CA")
    if previous_carbon is None or previous_oxygen is None or previous_alpha is None:
        return []

    carbon_position = previous_carbon.get_coord()
    oxygen_position = previous_oxygen.get_coord()
    methyl_position = previous_alpha.get_coord()
    cap_atoms = [
        AtomCoord(
            "C", float(carbon_position[0]), float(carbon_position[1]), float(carbon_position[2])
        ),
        AtomCoord(
            "O", float(oxygen_position[0]), float(oxygen_position[1]), float(oxygen_position[2])
        ),
        AtomCoord(
            "C", float(methyl_position[0]), float(methyl_position[1]), float(methyl_position[2])
        ),
    ]
    cap_atoms.extend(_tetrahedral_methyl_hydrogens(methyl_position, carbon_position))
    return cap_atoms


def _nme_cap(residue, next_residue) -> list[AtomCoord]:
    carbonyl_carbon = _find_atom(residue, "C")
    next_nitrogen = _find_atom(next_residue, "N")
    next_alpha = _find_atom(next_residue, "CA")
    if carbonyl_carbon is None or next_nitrogen is None or next_alpha is None:
        return []

    nitrogen_position = next_nitrogen.get_coord()
    methyl_position = next_alpha.get_coord()
    carbonyl_position = carbonyl_carbon.get_coord()

    cap_atoms = [
        AtomCoord(
            "N", float(nitrogen_position[0]), float(nitrogen_position[1]), float(nitrogen_position[2])
        ),
        AtomCoord(
            "C", float(methyl_position[0]), float(methyl_position[1]), float(methyl_position[2])
        ),
    ]
    cap_atoms.extend(_amide_hydrogen(nitrogen_position, carbonyl_position, methyl_position))
    cap_atoms.extend(_tetrahedral_methyl_hydrogens(methyl_position, nitrogen_position))
    return cap_atoms


def infer_protein_fragment_charge(
    residues: Sequence[ResidueRef],
    capped: bool = True,
) -> int:
    charged_residues = {
        "ARG": 1,
        "ASP": -1,
        "GLU": -1,
        "HIP": 1,
        "HIPH": 1,
        "HSP": 1,
        "LYS": 1,
    }
    charge = sum(charged_residues.get(residue.resname.upper(), 0) for residue in residues)

    if not capped and residues:
        chains = {residue.chain_id for residue in residues}
        charge += len(chains)
        charge -= len(chains)
    return charge


def infer_ligand_formal_charge(sdf_path: str) -> int:
    supplier = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError(f"No valid molecule found in {sdf_path}")
    return int(Chem.GetFormalCharge(mol))


def infer_ligand_rotatable_bonds(sdf_path: str) -> int:
    supplier = Chem.SDMolSupplier(sdf_path, sanitize=True, removeHs=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError(f"No valid molecule found in {sdf_path}")
    return int(rdMolDescriptors.CalcNumRotatableBonds(mol))


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


def select_protein_box(
    pdb_path: str,
    ligand_coords: Sequence[AtomCoord],
    padding_angstroms: float,
    cap_residues: bool = True,
) -> BoxSelection:
    if padding_angstroms < 0:
        raise ValueError("Box padding must be non-negative")
    if not ligand_coords:
        raise ValueError("Ligand coordinates are required for box selection")

    xs = [atom.x for atom in ligand_coords]
    ys = [atom.y for atom in ligand_coords]
    zs = [atom.z for atom in ligand_coords]

    x_min = min(xs) - padding_angstroms
    x_max = max(xs) + padding_angstroms
    y_min = min(ys) - padding_angstroms
    y_max = max(ys) + padding_angstroms
    z_min = min(zs) - padding_angstroms
    z_max = max(zs) + padding_angstroms

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    ordered_residues = [
        residue
        for residue in sorted(structure.get_residues(), key=_residue_sort_key)
        if list(residue.get_atoms())
    ]

    selected_indices: set[int] = set()
    for index, residue in enumerate(ordered_residues):
        for atom in residue.get_atoms():
            x, y, z = atom.get_coord()
            if (
                x_min <= float(x) <= x_max
                and y_min <= float(y) <= y_max
                and z_min <= float(z) <= z_max
            ):
                selected_indices.add(index)
                break

    bridged_indices = set(selected_indices)
    for index in sorted(selected_indices):
        if index + 2 in selected_indices:
            bridged_indices.add(index + 1)
    selected_indices = bridged_indices

    selected_residues = []
    selected_atoms: list[AtomCoord] = []
    cap_atoms: list[AtomCoord] = []
    for index, residue in enumerate(ordered_residues):
        if index not in selected_indices:
            continue

        previous_selected = index - 1 in selected_indices
        next_selected = index + 1 in selected_indices

        if cap_residues and index > 0 and not previous_selected:
            new_cap_atoms = _ace_cap(residue, ordered_residues[index - 1])
            cap_atoms.extend(new_cap_atoms)
            selected_atoms.extend(new_cap_atoms)

        atoms = sorted(residue.get_atoms(), key=_atom_sort_key)
        selected_atoms.extend(_atomcoord_from_biopython(atom) for atom in atoms)

        if cap_residues and index + 1 < len(ordered_residues) and not next_selected:
            new_cap_atoms = _nme_cap(residue, ordered_residues[index + 1])
            cap_atoms.extend(new_cap_atoms)
            selected_atoms.extend(new_cap_atoms)

        hetflag, resseq, icode = residue.get_id()
        selected_residues.append(
            ResidueRef(
                chain_id=residue.get_parent().id,
                hetflag=hetflag,
                resseq=resseq,
                icode=icode,
                resname=residue.resname,
            )
        )

    if not selected_atoms:
        raise ValueError(
            f"No protein residues found within a {padding_angstroms:.3f} A box around the ligand"
        )

    return BoxSelection(
        residues=tuple(selected_residues),
        atoms=tuple(selected_atoms),
        cap_atoms=tuple(cap_atoms),
    )


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


def _coords_to_tblite_inputs(
    coords: Sequence[AtomCoord],
) -> tuple[np.ndarray, np.ndarray]:
    angstrom_to_bohr = 1.8897259886
    numbers = np.asarray(
        [_element_to_atomic_number(atom.element) for atom in coords],
        dtype=np.int32,
    )
    positions = np.asarray(
        [[atom.x, atom.y, atom.z] for atom in coords],
        dtype=float,
    )
    return numbers, positions * angstrom_to_bohr


def _element_to_atomic_number(symbol: str) -> int:
    symbol = symbol.strip().title()
    periodic_table = {
        "H": 1,
        "He": 2,
        "Li": 3,
        "Be": 4,
        "B": 5,
        "C": 6,
        "N": 7,
        "O": 8,
        "F": 9,
        "Ne": 10,
        "Na": 11,
        "Mg": 12,
        "Al": 13,
        "Si": 14,
        "P": 15,
        "S": 16,
        "Cl": 17,
        "Ar": 18,
        "K": 19,
        "Ca": 20,
        "Sc": 21,
        "Ti": 22,
        "V": 23,
        "Cr": 24,
        "Mn": 25,
        "Fe": 26,
        "Co": 27,
        "Ni": 28,
        "Cu": 29,
        "Zn": 30,
        "Ga": 31,
        "Ge": 32,
        "As": 33,
        "Se": 34,
        "Br": 35,
        "Kr": 36,
        "Rb": 37,
        "Sr": 38,
        "Y": 39,
        "Zr": 40,
        "Nb": 41,
        "Mo": 42,
        "Tc": 43,
        "Ru": 44,
        "Rh": 45,
        "Pd": 46,
        "Ag": 47,
        "Cd": 48,
        "In": 49,
        "Sn": 50,
        "Sb": 51,
        "Te": 52,
        "I": 53,
        "Xe": 54,
        "Cs": 55,
        "Ba": 56,
        "La": 57,
        "Ce": 58,
        "Pr": 59,
        "Nd": 60,
        "Pm": 61,
        "Sm": 62,
        "Eu": 63,
        "Gd": 64,
        "Tb": 65,
        "Dy": 66,
        "Ho": 67,
        "Er": 68,
        "Tm": 69,
        "Yb": 70,
        "Lu": 71,
        "Hf": 72,
        "Ta": 73,
        "W": 74,
        "Re": 75,
        "Os": 76,
        "Ir": 77,
        "Pt": 78,
        "Au": 79,
        "Hg": 80,
        "Tl": 81,
        "Pb": 82,
        "Bi": 83,
        "Po": 84,
        "At": 85,
        "Rn": 86,
    }
    if symbol not in periodic_table:
        raise ValueError(f"Unsupported element for xTB: {symbol}")
    return periodic_table[symbol]


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


def _run_xtb_energy(
    coords: Sequence[AtomCoord],
    method: str,
    charge: int,
    multiplicity: int,
    threads: int,
    solvent: str | None = None,
    solvent_state: str | None = None,
) -> float:
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

    try:
        from tblite.interface import Calculator
    except ImportError as exc:
        raise ImportError(
            "GFN-xTB support requires the `tblite` package in the repo venv. "
            "Install prerequisites with `brew install libomp pkgconf`, then install "
            "`tblite` into `.venv311`."
        ) from exc

    method_lookup = {
        "xtb": "GFN2-xTB",
        "gfn1-xtb": "GFN1-xTB",
        "gfn2-xtb": "GFN2-xTB",
        "ipea1-xtb": "IPEA1-xTB",
    }
    method_key = method.strip().lower()
    if method_key not in method_lookup:
        raise ValueError(f"Unsupported xTB method: {method}")

    numbers, positions = _coords_to_tblite_inputs(coords)
    calc = Calculator(
        method_lookup[method_key],
        numbers,
        positions,
        charge=float(charge),
        uhf=max(multiplicity - 1, 0),
    )
    calc.set("verbosity", 0)
    if solvent is not None:
        if solvent_state is None:
            calc.add("alpb-solvation", solvent)
        else:
            calc.add("alpb-solvation", solvent, solvent_state)
    result = calc.singlepoint()
    return float(result.get("energy"))


def binding_energy(
    protein_coords: Iterable[AtomCoord],
    ligand_coords: Iterable[AtomCoord],
    method: str,
    basis: str,
    protein_charge: int = 0,
    ligand_charge: int = 0,
    complex_charge: int | None = None,
    multiplicity: int = 1,
    threads: int = 1,
    xtb_alpb_solvent: str | None = None,
    xtb_alpb_state: str | None = None,
) -> Tuple[float, float, float, float]:
    protein_list = list(protein_coords)
    ligand_list = list(ligand_coords)
    method_normalized = method.strip().lower()
    if complex_charge is None:
        complex_charge = protein_charge + ligand_charge

    if method_normalized.endswith("-xtb") or method_normalized == "xtb":
        e_protein = _run_xtb_energy(
            protein_list,
            method=method,
            charge=protein_charge,
            multiplicity=multiplicity,
            threads=threads,
            solvent=xtb_alpb_solvent,
            solvent_state=xtb_alpb_state,
        )
        e_ligand = _run_xtb_energy(
            ligand_list,
            method=method,
            charge=ligand_charge,
            multiplicity=multiplicity,
            threads=threads,
            solvent=xtb_alpb_solvent,
            solvent_state=xtb_alpb_state,
        )
        e_complex = _run_xtb_energy(
            protein_list + ligand_list,
            method=method,
            charge=complex_charge,
            multiplicity=multiplicity,
            threads=threads,
            solvent=xtb_alpb_solvent,
            solvent_state=xtb_alpb_state,
        )
    else:
        _ensure_psi4_importable()
        import psi4

        protein_geom = _coords_to_psi4_geometry(
            protein_list, protein_charge, multiplicity
        )
        ligand_geom = _coords_to_psi4_geometry(
            ligand_list, ligand_charge, multiplicity
        )
        complex_geom = _coords_to_psi4_geometry(
            protein_list + ligand_list, complex_charge, multiplicity
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


def estimate_binding_free_energy(
    interaction_energy_hartree: float,
    interaction_scale: float = 0.35,
    entropy_penalty_kcal_mol: float = 10.0,
) -> float:
    if interaction_scale < 0:
        raise ValueError("Interaction scale must be non-negative")
    return (
        interaction_energy_hartree * interaction_scale
        + entropy_penalty_kcal_mol / HARTREE_TO_KCAL_MOL
    )


def estimate_entropy_penalty_kcal_mol(
    ligand_rotatable_bonds: int,
    base_penalty_kcal_mol: float = 6.0,
    per_rotor_penalty_kcal_mol: float = 0.8,
) -> float:
    if ligand_rotatable_bonds < 0:
        raise ValueError("Rotatable bond count must be non-negative")
    return base_penalty_kcal_mol + per_rotor_penalty_kcal_mol * ligand_rotatable_bonds
