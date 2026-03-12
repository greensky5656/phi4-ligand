"""CLI for Psi4 ligand binding energy calculations."""

from __future__ import annotations

import argparse
import sys

from .binding import binding_energy, load_ligand_coords, load_protein_coords


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute binding energy from protein PDB and ligand SDF."
    )
    parser.add_argument("--protein-pdb", required=True, help="Path to protein PDB")
    parser.add_argument("--ligand-sdf", required=True, help="Path to ligand SDF")
    parser.add_argument("--method", default="b3lyp", help="QM method (default: b3lyp)")
    parser.add_argument(
        "--basis", default="6-31g*", help="Basis set (default: 6-31g*)"
    )
    parser.add_argument(
        "--charge", type=int, default=0, help="Total charge (default: 0)"
    )
    parser.add_argument(
        "--multiplicity", type=int, default=1, help="Spin multiplicity (default: 1)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Psi4 threads (default: 1 for reproducibility)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    protein_coords = load_protein_coords(args.protein_pdb)
    ligand_coords = load_ligand_coords(args.ligand_sdf)

    e_complex, e_protein, e_ligand, delta = binding_energy(
        protein_coords,
        ligand_coords,
        method=args.method,
        basis=args.basis,
        charge=args.charge,
        multiplicity=args.multiplicity,
        threads=args.threads,
    )

    print(f"E_complex: {e_complex:.8f} Eh")
    print(f"E_protein: {e_protein:.8f} Eh")
    print(f"E_ligand : {e_ligand:.8f} Eh")
    print(f"Delta    : {delta:.8f} Eh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
