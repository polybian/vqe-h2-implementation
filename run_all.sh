#!/usr/bin/env bash
# Runs every implementation in order.
set -euo pipefail
cd "$(dirname "$0")"

# Prefer the project's venv (where qiskit lives); fall back to python3.
if [ -x ".venv/bin/python" ]; then
  PY=".venv/bin/python"
else
  PY="python3"
fi

echo "=== 01 — H2 Hamiltonian and exact ground state ==="
$PY scripts/01_h2_hamiltonian.py
echo
echo "=== 02 — VQE ground state ==="
$PY scripts/02_vqe_ground_state.py
echo
echo "=== 03 — Dissociation curve ==="
$PY scripts/03_dissociation_curve.py
echo
echo "All three implementations completed successfully."
