#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${VENV_DIR:-$ROOT/.venv311}"
PYTHON_BIN="${PYTHON_BIN:-/opt/homebrew/bin/python3.11}"
BUILD_DIR="${BUILD_DIR:-$ROOT/build/psi4}"
INSTALL_DIR="${INSTALL_DIR:-$ROOT/local/psi4}"
JOBS="${JOBS:-4}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  if command -v python3.11 >/dev/null 2>&1; then
    PYTHON_BIN="$(command -v python3.11)"
  else
    echo "python3.11 is required. Set PYTHON_BIN or install Python 3.11." >&2
    exit 1
  fi
fi

if ! command -v cmake >/dev/null 2>&1; then
  echo "cmake is required." >&2
  exit 1
fi

if ! command -v gfortran >/dev/null 2>&1; then
  echo "gfortran is required." >&2
  exit 1
fi

"$PYTHON_BIN" -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install --upgrade pip setuptools wheel
"$VENV_DIR/bin/python" -m pip install \
  "numpy<2" \
  biopython \
  scipy \
  pydantic \
  msgpack \
  networkx \
  pint \
  psutil \
  py-cpuinfo \
  rdkit-pypi

if [[ -d /opt/homebrew/opt/libomp ]]; then
  export LDFLAGS="-L/opt/homebrew/opt/libomp/lib ${LDFLAGS:-}"
  export CPPFLAGS="-I/opt/homebrew/opt/libomp/include ${CPPFLAGS:-}"
  export PKG_CONFIG_PATH="/opt/homebrew/opt/libomp/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  "$VENV_DIR/bin/python" -m pip install tblite
else
  cat <<EOF

Skipping tblite install because libomp was not found.
For GFN-xTB support on Homebrew-based macOS installs:
  brew install libomp pkgconf
  export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
  export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
  export PKG_CONFIG_PATH="/opt/homebrew/opt/libomp/lib/pkgconfig"
  $VENV_DIR/bin/python -m pip install tblite
EOF
fi

"$VENV_DIR/bin/python" -m pip install -e "$ROOT" --no-deps

CXX_FLAGS="${CXX_FLAGS:-}"
if [[ -d /opt/homebrew/include ]]; then
  CXX_FLAGS="-I/opt/homebrew/include ${CXX_FLAGS}"
fi
CXX_FLAGS="${CXX_FLAGS#"${CXX_FLAGS%%[![:space:]]*}"}"
CXX_FLAGS="${CXX_FLAGS%"${CXX_FLAGS##*[![:space:]]}"}"

CONFIGURE_ARGS=(
  -S "$ROOT/vendor/psi4"
  -B "$BUILD_DIR"
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR"
  -DPython_EXECUTABLE="$VENV_DIR/bin/python"
  -DCMAKE_C_COMPILER="${CC:-clang}"
  -DCMAKE_CXX_COMPILER="${CXX:-clang++}"
  -DCMAKE_Fortran_COMPILER="${FC:-gfortran}"
)

if [[ -n "$CXX_FLAGS" ]]; then
  CONFIGURE_ARGS+=("-DCMAKE_CXX_FLAGS=$CXX_FLAGS")
fi

cmake "${CONFIGURE_ARGS[@]}"
cmake --build "$BUILD_DIR" --target install -j "$JOBS"

cat <<EOF

Bootstrap complete.

Use Psi4 from this repo with:
  $VENV_DIR/bin/python -m psi4_ligands.cli --protein-pdb path/to/protein.pdb --ligand-sdf path/to/ligand.sdf

If you need to override the local Psi4 location:
  export PSI4_PYTHONPATH="$INSTALL_DIR/lib"

For GFN-xTB support, ensure tblite is present in the repo venv:
  $VENV_DIR/bin/python -m pip show tblite
EOF
