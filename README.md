# psi4-ligands

Local setup for protein/ligand binding energy calculations using Psi4.

## Notes
- Psi4 is vendored in this repo under `vendor/psi4`.
- Installation is source-based and local to this repo; no conda is required.
- The supported bootstrap path builds Psi4 into `local/psi4` and uses a repo-local Python 3.11 venv at `.venv311`.

## Bootstrap

Run the bootstrap script from the repo root:

```bash
./scripts/bootstrap.sh
```

Requirements:
- `python3.11`
- `cmake`
- `gfortran`

Optional environment overrides:
- `PYTHON_BIN=/path/to/python3.11`
- `JOBS=8`
- `CC=clang`
- `CXX=clang++`
- `FC=gfortran`
- `VENV_DIR=/custom/venv`

The bootstrap script:
- creates `.venv311`
- installs the Python dependencies needed by this package and Psi4
- installs `psi4-ligands` in editable mode
- configures and builds the vendored Psi4 source into `local/psi4`

At runtime, `psi4_ligands` looks for Psi4 in:
- `PSI4_PYTHONPATH`
- `local/psi4/lib`
- `build/psi4/stage/lib`

Compatibility note:
- `rdkit-pypi` in this setup requires `numpy<2`.

## CLI

Example:

```bash
.venv311/bin/python -m psi4_ligands.cli \
  --protein-pdb path/to/protein.pdb \
  --ligand-sdf path/to/ligand.sdf
```

Reproducibility:
- Protein atoms are sorted deterministically by chain/residue/atom identifiers.
- Psi4 threads default to 1 (override with `--threads`).
