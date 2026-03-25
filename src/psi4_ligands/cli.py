"""CLI for protein/ligand binding energy calculations."""

from __future__ import annotations

import argparse
import sys

from .binding import (
    binding_energy,
    infer_ligand_formal_charge,
    infer_protein_fragment_charge,
    load_ligand_coords,
    load_protein_coords,
    select_protein_box,
)

HARTREE_TO_KCAL_MOL = 627.509474


def _parse_charge_argument(value: str) -> str | int:
    lowered = value.strip().lower()
    if lowered == "auto":
        return "auto"
    return int(value)


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
        "--protein-charge",
        type=_parse_charge_argument,
        default="auto",
        help="Protein fragment charge or auto (default: auto)",
    )
    parser.add_argument(
        "--ligand-charge",
        type=_parse_charge_argument,
        default="auto",
        help="Ligand fragment charge or auto (default: auto)",
    )
    parser.add_argument(
        "--complex-charge",
        type=_parse_charge_argument,
        help="Complex charge or auto (default: protein-charge + ligand-charge)",
    )
    parser.add_argument(
        "--multiplicity", type=int, default=1, help="Spin multiplicity (default: 1)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Backend threads (default: 1 for reproducibility)",
    )
    parser.add_argument(
        "--box-padding-angstroms",
        type=float,
        help="Select complete protein residues with any atom inside a ligand-centered axis-aligned box padded by this many angstroms",
    )
    parser.add_argument(
        "--cap-box-residues",
        dest="cap_box_residues",
        action="store_true",
        default=True,
        help="Cap cut peptide bonds in boxed protein selections (default: enabled)",
    )
    parser.add_argument(
        "--no-cap-box-residues",
        dest="cap_box_residues",
        action="store_false",
        help="Do not cap cut peptide bonds in boxed protein selections",
    )
    parser.add_argument(
        "--xtb-alpb-solvent",
        help="Optional ALPB solvent for xTB methods, for example water",
    )
    parser.add_argument(
        "--xtb-alpb-state",
        help="Optional ALPB solution state for xTB methods, for example bar1mol or reference",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    ligand_coords = load_ligand_coords(args.ligand_sdf)
    if args.box_padding_angstroms is None:
        protein_coords = load_protein_coords(args.protein_pdb)
        selected_residues = None
        cap_atoms = ()
    else:
        box = select_protein_box(
            args.protein_pdb,
            ligand_coords,
            args.box_padding_angstroms,
            cap_residues=args.cap_box_residues,
        )
        protein_coords = list(box.atoms)
        selected_residues = box.residues
        cap_atoms = box.cap_atoms

    if args.protein_charge == "auto":
        if selected_residues is not None:
            protein_charge = infer_protein_fragment_charge(
                selected_residues,
                capped=args.cap_box_residues,
            )
        else:
            protein_charge = 0
    else:
        protein_charge = args.protein_charge

    if args.ligand_charge == "auto":
        ligand_charge = infer_ligand_formal_charge(args.ligand_sdf)
    else:
        ligand_charge = args.ligand_charge

    if args.complex_charge in (None, "auto"):
        complex_charge = protein_charge + ligand_charge
    else:
        complex_charge = args.complex_charge

    e_complex, e_protein, e_ligand, delta = binding_energy(
        protein_coords,
        ligand_coords,
        method=args.method,
        basis=args.basis,
        protein_charge=protein_charge,
        ligand_charge=ligand_charge,
        complex_charge=complex_charge,
        multiplicity=args.multiplicity,
        threads=args.threads,
        xtb_alpb_solvent=args.xtb_alpb_solvent,
        xtb_alpb_state=args.xtb_alpb_state,
    )

    if selected_residues is not None:
        print(f"Selected residues: {len(selected_residues)}")
        print(f"Selected atoms   : {len(protein_coords)}")
        print(f"Cap atoms        : {len(cap_atoms)}")
        print(f"Capping         : {'on' if args.cap_box_residues else 'off'}")
    print(
        f"Charges        : protein={protein_charge}, ligand={ligand_charge}, "
        f"complex={complex_charge}"
    )
    if args.xtb_alpb_solvent:
        print(
            f"xTB solvent    : ALPB {args.xtb_alpb_solvent}"
            + (f" ({args.xtb_alpb_state})" if args.xtb_alpb_state else "")
        )
    print("Model note     : relative local interaction score, not a binding free energy")
    print(f"E_complex: {e_complex * HARTREE_TO_KCAL_MOL:.2f} kcal/mol")
    print(f"E_protein: {e_protein * HARTREE_TO_KCAL_MOL:.2f} kcal/mol")
    print(f"E_ligand : {e_ligand * HARTREE_TO_KCAL_MOL:.2f} kcal/mol")
    print(f"Score    : {delta * HARTREE_TO_KCAL_MOL:.2f} kcal/mol")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
