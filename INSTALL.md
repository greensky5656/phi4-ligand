# Install

## macOS (Homebrew)

Install system prerequisites:

```bash
brew install python@3.11 cmake gcc libomp pkgconf
```

Clone and bootstrap:

```bash
git clone <repo-url>
cd psi4-ligands
./scripts/bootstrap.sh
```

What `./scripts/bootstrap.sh` does:
- creates the repo-local `.venv311`
- installs Python dependencies into the repo venv
- installs `tblite` into `.venv311` when `libomp` is available
- installs `psi4-ligands` in editable mode
- builds the vendored Psi4 into `local/psi4`

## Verify

Check that `tblite` is installed in the repo venv:

```bash
.venv311/bin/python -m pip show tblite
```

Run the sample scoring workflow:

```bash
./scripts/run_score.sh \
  --protein-pdb samples/B49_receptor.pdb \
  --ligand-sdf samples/B49_ligandl.sdf
```

## Requirements

Required system tools:
- `python3.11`
- `cmake`
- `gfortran` from `gcc`
- `libomp`
- `pkgconf`

## Notes

- The default wrapper uses `gfn2-xtb`, `4.0 Å` cutoff, `ALPB water`, and `1` thread.
- Boxed residue capping is enabled by default.
- Protein and ligand charges default to `auto`.
- Output is reported in `kcal/mol`.
- The reported `Score` is a relative local interaction score, not a binding free energy.

## Useful Commands

Disable capping:

```bash
./scripts/run_score.sh \
  --protein-pdb your_protein.pdb \
  --ligand-sdf your_ligand.sdf \
  --no-cap
```

Change cutoff:

```bash
./scripts/run_score.sh \
  --protein-pdb your_protein.pdb \
  --ligand-sdf your_ligand.sdf \
  --cutoff 6.0
```

Run the CLI directly:

```bash
.venv311/bin/python -m psi4_ligands.cli \
  --protein-pdb your_protein.pdb \
  --ligand-sdf your_ligand.sdf \
  --method gfn2-xtb \
  --box-padding-angstroms 4 \
  --xtb-alpb-solvent water \
  --threads 1
```
