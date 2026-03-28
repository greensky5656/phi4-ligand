#!/usr/bin/env python3
"""Run all benchmark protein/ligand pairs and write a report."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from psi4_ligands.binding import (
    HARTREE_TO_KCAL_MOL,
    binding_energy,
    estimate_binding_free_energy,
    estimate_entropy_penalty_kcal_mol,
    infer_ligand_formal_charge,
    infer_ligand_rotatable_bonds,
    infer_protein_fragment_charge,
    load_ligand_coords,
    select_protein_box,
)


@dataclass(frozen=True)
class BenchmarkResult:
    run: str
    selected_residues: int
    selected_atoms: int
    cap_atoms: int
    protein_charge: int
    ligand_charge: int
    complex_charge: int
    ligand_rotors: int
    docking_affinity_kcal_mol: float | None
    e_complex_kcal_mol: float
    e_protein_kcal_mol: float
    e_ligand_kcal_mol: float
    dE_bind_kcal_mol: float
    entropy_penalty_kcal_mol: float
    interaction_scale: float
    dG_bind_kcal_mol: float


def _parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Run all benchmark receptor/ligand pairs and write a report."
    )
    parser.add_argument(
        "--benchmarks-dir",
        default=str(root / "samples" / "benchmarks"),
        help="Directory containing benchmark subfolders and selection_manifest.tsv",
    )
    parser.add_argument("--method", default="gfn2-xtb", help="Method to use")
    parser.add_argument("--basis", default="6-31g*", help="Basis for non-xTB methods")
    parser.add_argument(
        "--box-padding-angstroms",
        type=float,
        default=4.0,
        help="Ligand-centered box padding in angstroms",
    )
    parser.add_argument(
        "--cap-box-residues",
        dest="cap_box_residues",
        action="store_true",
        default=True,
        help="Cap cut peptide bonds in boxed selections",
    )
    parser.add_argument(
        "--no-cap-box-residues",
        dest="cap_box_residues",
        action="store_false",
        help="Disable capping for boxed selections",
    )
    parser.add_argument("--threads", type=int, default=1, help="Backend threads")
    parser.add_argument(
        "--xtb-alpb-solvent",
        default="water",
        help="Optional ALPB solvent for xTB methods",
    )
    parser.add_argument(
        "--xtb-alpb-state",
        help="Optional ALPB solution state for xTB methods",
    )
    parser.add_argument(
        "--interaction-scale",
        type=float,
        default=0.35,
        help="Compression factor for heuristic dG_bind estimation",
    )
    parser.add_argument(
        "--entropy-penalty-kcal-mol",
        type=float,
        help="Explicit entropy penalty in kcal/mol",
    )
    parser.add_argument(
        "--entropy-base-kcal-mol",
        type=float,
        default=6.0,
        help="Base entropy penalty when auto mode is used",
    )
    parser.add_argument(
        "--entropy-per-rotor-kcal-mol",
        type=float,
        default=0.8,
        help="Additional entropy penalty per rotatable bond in auto mode",
    )
    parser.add_argument(
        "--report-prefix",
        default="benchmark_report",
        help="Output filename prefix inside the benchmarks directory",
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        help="Optional subset of benchmark folder names to run",
    )
    return parser.parse_args()


def _load_manifest(manifest_path: Path) -> dict[str, dict[str, str]]:
    if not manifest_path.exists():
        return {}
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["run"]: row for row in reader}


def _kcal_mol(energy_hartree: float) -> float:
    return energy_hartree * HARTREE_TO_KCAL_MOL


def _run_one(
    run_dir: Path,
    manifest_row: dict[str, str] | None,
    args: argparse.Namespace,
) -> BenchmarkResult:
    ligand_path = run_dir / "ligand.sdf"
    receptor_path = run_dir / "receptor.pdb"

    ligand_coords = load_ligand_coords(str(ligand_path))
    box = select_protein_box(
        str(receptor_path),
        ligand_coords,
        args.box_padding_angstroms,
        cap_residues=args.cap_box_residues,
    )
    protein_coords = list(box.atoms)

    protein_charge = infer_protein_fragment_charge(
        box.residues,
        capped=args.cap_box_residues,
    )
    ligand_charge = infer_ligand_formal_charge(str(ligand_path))
    ligand_rotors = infer_ligand_rotatable_bonds(str(ligand_path))
    complex_charge = protein_charge + ligand_charge

    e_complex, e_protein, e_ligand, delta = binding_energy(
        protein_coords,
        ligand_coords,
        method=args.method,
        basis=args.basis,
        protein_charge=protein_charge,
        ligand_charge=ligand_charge,
        complex_charge=complex_charge,
        threads=args.threads,
        xtb_alpb_solvent=args.xtb_alpb_solvent,
        xtb_alpb_state=args.xtb_alpb_state,
    )

    if args.entropy_penalty_kcal_mol is None:
        entropy_penalty = estimate_entropy_penalty_kcal_mol(
            ligand_rotors,
            base_penalty_kcal_mol=args.entropy_base_kcal_mol,
            per_rotor_penalty_kcal_mol=args.entropy_per_rotor_kcal_mol,
        )
    else:
        entropy_penalty = args.entropy_penalty_kcal_mol

    dG_bind = estimate_binding_free_energy(
        delta,
        interaction_scale=args.interaction_scale,
        entropy_penalty_kcal_mol=entropy_penalty,
    )

    docking_affinity = None
    if manifest_row is not None:
        raw_value = manifest_row.get("score_only_affinity_kcal_per_mol", "").strip()
        if raw_value:
            docking_affinity = float(raw_value)

    return BenchmarkResult(
        run=run_dir.name,
        selected_residues=len(box.residues),
        selected_atoms=len(protein_coords),
        cap_atoms=len(box.cap_atoms),
        protein_charge=protein_charge,
        ligand_charge=ligand_charge,
        complex_charge=complex_charge,
        ligand_rotors=ligand_rotors,
        docking_affinity_kcal_mol=docking_affinity,
        e_complex_kcal_mol=_kcal_mol(e_complex),
        e_protein_kcal_mol=_kcal_mol(e_protein),
        e_ligand_kcal_mol=_kcal_mol(e_ligand),
        dE_bind_kcal_mol=_kcal_mol(delta),
        entropy_penalty_kcal_mol=entropy_penalty,
        interaction_scale=args.interaction_scale,
        dG_bind_kcal_mol=_kcal_mol(dG_bind),
    )


def _write_tsv(results: list[BenchmarkResult], output_path: Path) -> None:
    fieldnames = [
        "run",
        "selected_residues",
        "selected_atoms",
        "cap_atoms",
        "protein_charge",
        "ligand_charge",
        "complex_charge",
        "ligand_rotors",
        "docking_affinity_kcal_mol",
        "e_complex_kcal_mol",
        "e_protein_kcal_mol",
        "e_ligand_kcal_mol",
        "dE_bind_kcal_mol",
        "entropy_penalty_kcal_mol",
        "interaction_scale",
        "dG_bind_kcal_mol",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row.__dict__)


def _fmt(value: float | None, digits: int = 2) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def _write_markdown(
    results: list[BenchmarkResult],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    lines = [
        "# Benchmark Report",
        "",
        f"Method: `{args.method}`",
        f"Box padding: `{args.box_padding_angstroms}` A",
        f"Capping: `{'on' if args.cap_box_residues else 'off'}`",
        f"xTB solvent: `{args.xtb_alpb_solvent or 'none'}`",
        f"dG model: `dG_bind = interaction_scale * dE_bind + entropy_penalty`",
        "",
        "| run | docking_affinity | dE_bind | entropy_penalty | dG_bind |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in results:
        lines.append(
            f"| {row.run} | {_fmt(row.docking_affinity_kcal_mol, 3)} | "
            f"{_fmt(row.dE_bind_kcal_mol)} | {_fmt(row.entropy_penalty_kcal_mol)} | "
            f"{_fmt(row.dG_bind_kcal_mol)} |"
        )
    output_path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = _parse_args()
    benchmarks_dir = Path(args.benchmarks_dir).resolve()
    manifest = _load_manifest(benchmarks_dir / "selection_manifest.tsv")
    run_dirs = sorted(path for path in benchmarks_dir.iterdir() if path.is_dir())
    if args.runs:
        requested = set(args.runs)
        run_dirs = [path for path in run_dirs if path.name in requested]

    results: list[BenchmarkResult] = []
    for run_dir in run_dirs:
        result = _run_one(run_dir, manifest.get(run_dir.name), args)
        results.append(result)
        report_tsv = benchmarks_dir / f"{args.report_prefix}.tsv"
        report_md = benchmarks_dir / f"{args.report_prefix}.md"
        _write_tsv(results, report_tsv)
        _write_markdown(results, report_md, args)
        print(
            f"{result.run}: docking={_fmt(result.docking_affinity_kcal_mol, 3)} "
            f"dE_bind={result.dE_bind_kcal_mol:.2f} dG_bind={result.dG_bind_kcal_mol:.2f}",
            flush=True,
        )

    print(f"Wrote {report_tsv}")
    print(f"Wrote {report_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
