#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$ROOT/.venv311/bin/python}"

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --protein-pdb FILE --ligand-sdf FILE [options]

Required:
  --protein-pdb FILE              Protein PDB path
  --ligand-sdf FILE               Ligand SDF path

Options:
  --cutoff FLOAT                  Ligand-centered box padding in angstroms (default: 4.0)
  --method NAME                   Method to use (default: gfn2-xtb)
  --threads INT                   Threads (default: 1)
  --protein-charge INT|auto       Protein fragment charge (default: auto)
  --ligand-charge INT|auto        Ligand charge (default: auto)
  --complex-charge INT|auto       Complex charge (default: auto)
  --multiplicity INT              Spin multiplicity (default: 1)
  --solvent NAME                  xTB ALPB solvent, e.g. water (default: water)
  --solvent-state NAME            Optional ALPB state, e.g. bar1mol
  --cap                           Cap boxed residue cuts (default: on)
  --no-cap                        Disable capping for boxed residue cuts
  --basis NAME                    Basis for non-xTB methods (default: 6-31g*)
  --python-bin FILE               Python interpreter to use (default: .venv311/bin/python)
  -h, --help                      Show this help

Example:
  $(basename "$0") \\
    --protein-pdb "$ROOT/samples/B49_receptor.pdb" \\
    --ligand-sdf "$ROOT/samples/B49_ligandl.sdf"
EOF
}

protein_pdb=""
ligand_sdf=""
cutoff="4.0"
method="gfn2-xtb"
threads="1"
protein_charge="auto"
ligand_charge="auto"
complex_charge="auto"
multiplicity="1"
solvent="water"
solvent_state=""
cap_mode="on"
basis="6-31g*"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --protein-pdb)
      protein_pdb="$2"
      shift 2
      ;;
    --ligand-sdf)
      ligand_sdf="$2"
      shift 2
      ;;
    --cutoff)
      cutoff="$2"
      shift 2
      ;;
    --method)
      method="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --protein-charge)
      protein_charge="$2"
      shift 2
      ;;
    --ligand-charge)
      ligand_charge="$2"
      shift 2
      ;;
    --complex-charge)
      complex_charge="$2"
      shift 2
      ;;
    --multiplicity)
      multiplicity="$2"
      shift 2
      ;;
    --solvent)
      solvent="$2"
      shift 2
      ;;
    --solvent-state)
      solvent_state="$2"
      shift 2
      ;;
    --cap)
      cap_mode="on"
      shift
      ;;
    --no-cap)
      cap_mode="off"
      shift
      ;;
    --basis)
      basis="$2"
      shift 2
      ;;
    --python-bin)
      PYTHON_BIN="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$protein_pdb" || -z "$ligand_sdf" ]]; then
  usage >&2
  exit 1
fi

cli_args=(
  -m psi4_ligands.cli
  --protein-pdb "$protein_pdb"
  --ligand-sdf "$ligand_sdf"
  --method "$method"
  --basis "$basis"
  --box-padding-angstroms "$cutoff"
  --threads "$threads"
  --protein-charge "$protein_charge"
  --ligand-charge "$ligand_charge"
  --complex-charge "$complex_charge"
  --multiplicity "$multiplicity"
)

if [[ "$cap_mode" == "on" ]]; then
  cli_args+=(--cap-box-residues)
else
  cli_args+=(--no-cap-box-residues)
fi

if [[ -n "$solvent" ]]; then
  cli_args+=(--xtb-alpb-solvent "$solvent")
fi

if [[ -n "$solvent_state" ]]; then
  cli_args+=(--xtb-alpb-state "$solvent_state")
fi

exec "$PYTHON_BIN" "${cli_args[@]}"
