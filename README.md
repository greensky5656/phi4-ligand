# psi4-ligands

Local setup for protein/ligand binding energy calculations using Psi4.

Install instructions are in [INSTALL.md](/Users/jasondeckman/Code/affinra/psi4-ligands/INSTALL.md).

Quick install on macOS with Homebrew:

```bash
brew install python@3.11 cmake gcc libomp pkgconf
git clone <repo-url>
cd psi4-ligands
./scripts/bootstrap.sh
```

## CLI

Example:

```bash
.venv311/bin/python -m psi4_ligands.cli \
  --protein-pdb path/to/protein.pdb \
  --ligand-sdf path/to/ligand.sdf
```

GFN2-xTB example with a ligand-centered `4 A` protein box:

```bash
.venv311/bin/python -m psi4_ligands.cli \
  --protein-pdb samples/B49_receptor.pdb \
  --ligand-sdf samples/B49_ligandl.sdf \
  --method gfn2-xtb \
  --box-padding-angstroms 4 \
  --xtb-alpb-solvent water \
  --threads 1
```

Example with explicit `dG_bind` tuning:

```bash
.venv311/bin/python -m psi4_ligands.cli \
  --protein-pdb samples/B49_receptor.pdb \
  --ligand-sdf samples/B49_ligandl.sdf \
  --method gfn2-xtb \
  --box-padding-angstroms 4 \
  --xtb-alpb-solvent water \
  --interaction-scale 0.35 \
  --entropy-base-kcal-mol 6.0 \
  --entropy-per-rotor-kcal-mol 0.8
```

Equivalent wrapper script:

```bash
./scripts/run_score.sh \
  --protein-pdb samples/B49_receptor.pdb \
  --ligand-sdf samples/B49_ligandl.sdf
```

## Notes

- Psi4 is vendored in this repo under `vendor/psi4`.
- Installation is source-based and local to this repo.
- The supported bootstrap path builds Psi4 into `local/psi4` and uses a repo-local Python 3.11 venv at `.venv311`.
- `rdkit-pypi` in this setup requires `numpy<2`.

Reproducibility:
- Protein atoms are sorted deterministically by chain/residue/atom identifiers.
- Boxed protein selection includes whole residues in deterministic chain/residue order.
- Single-residue gaps are bridged before capping to avoid overlapping cap construction.
- Psi4 threads default to 1 (override with `--threads`).
- CLI energies are reported in `kcal/mol`.

GFN-xTB support:
- The repo uses `tblite` for `gfn1-xtb`, `gfn2-xtb`, and `ipea1-xtb`.
- On Homebrew-based macOS installs, `tblite` needs `libomp` and `pkgconf` available during installation.
- `./scripts/bootstrap.sh` installs `tblite` automatically when `libomp` is present.
- Optional ALPB solvent can be enabled with `--xtb-alpb-solvent`.
- Boxed protein fragments are capped by default with backbone-based ACE/NME-style caps.
- Capping can be disabled with `--no-cap-box-residues` or `./scripts/run_score.sh --no-cap`.
- Protein and ligand charges default to `auto` for CLI runs.
- `dE_bind` is the raw local interaction energy.
- `dG_bind` is a heuristic estimate computed from the interaction energy with a compression factor and an entropy penalty.
- The default entropy penalty is `6.0 + 0.8 * N_rotatable_bonds` kcal/mol, but it can be overridden from the CLI.
