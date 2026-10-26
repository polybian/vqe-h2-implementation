"""
01 - The H2 Hamiltonian: from Pauli coefficients to the exact answer.

Builds the 2-qubit H2 Hamiltonian (STO-3G basis, equilibrium bond length
R = 0.735 Å) from its literature Pauli coefficients - the parity-mapped
Hamiltonian with two-qubit reduction:

    H = -0.4804·II + 0.3935·ZI − 0.3935·IZ − 0.0113·ZZ + 0.1812·XX

The script prints the Hamiltonian, diagonalizes it exactly, and checks
the ground-state energy against the published FCI/STO-3G value
(-1.1373 Ha). This exact value is the reference every later
implementation is measured against.

Run:
    python3 scripts/01_h2_hamiltonian.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (
    H2_FCI_LITERATURE_HA,
    exact_two_qubit_energy,
    h2_two_qubit_matrix,
    h2_two_qubit_paulis,
)


def main() -> None:
    paulis = h2_two_qubit_paulis()
    print("H2 Hamiltonian (STO-3G, R = 0.735 Å, 2 qubits):")
    for term in ("II", "ZI", "IZ", "ZZ", "XX"):
        print(f"  {paulis[term]:+8.4f} · {term}")

    H = h2_two_qubit_matrix()
    print("\nDense matrix:")
    with np.printoptions(precision=4, suppress=True):
        print(H)

    eigvals = np.linalg.eigvalsh(H)
    print(f"\nEigenvalues: {np.round(eigvals, 6)}")

    e0 = exact_two_qubit_energy()
    print(f"\nExact ground-state energy : {e0:.6f} Ha")
    print(f"Literature FCI/STO-3G     : {H2_FCI_LITERATURE_HA:.4f} Ha")
    print(f"Difference                : {abs(e0 - H2_FCI_LITERATURE_HA) * 1000:.2f} mHa")

    assert abs(e0 - H2_FCI_LITERATURE_HA) < 0.005, \
        "Hamiltonian does not match the published H2 FCI energy!"
    print("\nOK - the Pauli Hamiltonian reproduces the published H2 "
          "ground-state energy.")


if __name__ == "__main__":
    main()
