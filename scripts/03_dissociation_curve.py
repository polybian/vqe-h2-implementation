"""
03 — The H2 dissociation curve: mapping the molecule's energy landscape.

For each bond length R from 0.4 to 2.5 Å this implementation:

  1. Runs a restricted Hartree-Fock STO-3G calculation (PySCF) to get
     the molecular-orbital integrals at that geometry.
  2. Builds the second-quantized Hamiltonian and parity-maps it to
     qubits (the shared core in common.py).
  3. Diagonalizes it EXACTLY (full CI in the STO-3G basis).

The energies below are therefore exact diagonalization results, not
VQE approximations — labeled honestly. The curve must show a minimum
near the equilibrium bond length (~0.74 Å): stretched or compressed,
the molecule is less bound.

Run:
    python3 scripts/03_dissociation_curve.py [--points 13]
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import exact_ground_state_energy

try:
    from pyscf import ao2mo, gto, scf
except ImportError as exc:
    raise SystemExit("PySCF is required for script 03: pip install pyscf") from exc


def energy_at_bond_length(r_ang: float) -> float:
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {r_ang}", basis="sto-3g",
                unit="Angstrom", verbose=0)
    mf = scf.RHF(mol)
    mf.max_cycle = 100
    mf.kernel()
    C = mf.mo_coeff
    h1 = C.T @ mf.get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), 2)
    return exact_ground_state_energy(h1, eri, mol.energy_nuc())


def ascii_plot(rs: np.ndarray, es: np.ndarray, width: int = 46) -> str:
    lo, hi = es.min(), es.max()
    lines = []
    for r, e in zip(rs, es):
        pos = int(width * (e - lo) / (hi - lo)) if hi > lo else 0
        marker = "*" if r == rs[np.argmin(es)] else "·"
        lines.append(f"  {r:4.2f} Å  {e:+.4f} Ha  "
                     + " " * pos + marker)
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=13)
    args = parser.parse_args()

    rs = np.linspace(0.4, 2.5, args.points)
    print(f"Computing exact FCI/STO-3G energies at {args.points} bond lengths...\n")
    es = np.array([energy_at_bond_length(r) for r in rs])

    print("  R        E(R)            dissociation curve (·), minimum (*)")
    print(ascii_plot(rs, es))

    i_min = int(np.argmin(es))
    r_min, e_min = rs[i_min], es[i_min]
    print(f"\nMinimum: R = {r_min:.3f} Å, E = {e_min:.6f} Ha")
    print(f"Dissociation limit (2.5 Å): E = {es[-1]:.6f} Ha  "
          f"(bound by {(es[-1] - e_min) * 1000:.1f} mHa)")

    assert 0.6 <= r_min <= 0.9, \
        f"Minimum at {r_min:.2f} Å is far from the known 0.74 Å!"
    assert e_min < es[0] and e_min < es[-1], \
        "Curve has no interior minimum — unphysical!"
    print("\nOK — the dissociation curve bottoms out at the equilibrium "
          "bond length, as quantum chemistry demands.")


if __name__ == "__main__":
    main()
