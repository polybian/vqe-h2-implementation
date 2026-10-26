"""
02 - VQE from scratch: the variational ground state of H2.

Implements the Variational Quantum Eigensolver without a black-box
VQE routine:

  * Ansatz (hardware-efficient): Ry(θ0)⊗Ry(θ1) → CNOT(0,1) →
    Ry(θ2)⊗Ry(θ3) on 2 qubits.
  * Expectation values ⟨ψ(θ)|H|ψ(θ)⟩ evaluated on the Qiskit Aer
    statevector simulator (exact, noise-free).
  * Classical optimization loop: scipy COBYLA, with the energy of
    every iteration printed.

The run converges to the exact ground-state energy found in script 01.
The final assertion requires agreement within 2% - the variational
principle guarantees we can never go *below* the exact value.

Run:
    python3 scripts/02_vqe_ground_state.py [--maxiter 200]
"""

import argparse
import os
import sys

import numpy as np
from scipy.optimize import minimize

from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit_aer import AerSimulator
from qiskit.quantum_info import SparsePauliOp, Statevector

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import exact_two_qubit_energy, h2_two_qubit_paulis


def build_ansatz() -> QuantumCircuit:
    theta = ParameterVector("θ", 4)
    qc = QuantumCircuit(2)
    qc.ry(theta[0], 0)
    qc.ry(theta[1], 1)
    qc.cx(0, 1)
    qc.ry(theta[2], 0)
    qc.ry(theta[3], 1)
    return qc


def make_cost_fn(ansatz: QuantumCircuit, hamiltonian: SparsePauliOp,
                 simulator: AerSimulator):
    def cost(params: np.ndarray) -> float:
        bound = ansatz.assign_parameters(params)
        bound.save_statevector()
        result = simulator.run(bound).result()
        psi = Statevector(result.data(0)["statevector"])
        return float(psi.expectation_value(hamiltonian).real)
    return cost


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maxiter", type=int, default=200)
    args = parser.parse_args()

    paulis = h2_two_qubit_paulis()
    hamiltonian = SparsePauliOp.from_list(
        [(p, c) for p, c in paulis.items()])
    e_exact = exact_two_qubit_energy()
    print(f"Target (exact) ground-state energy: {e_exact:.6f} Ha\n")

    ansatz = build_ansatz()
    simulator = AerSimulator(method="statevector", seed_simulator=123)
    cost = make_cost_fn(ansatz, hamiltonian, simulator)

    history = []

    def callback(params):
        e = cost(params)
        history.append(e)
        print(f"  iter {len(history):>3}: energy = {e:.6f} Ha")

    rng = np.random.default_rng(42)
    x0 = rng.uniform(-np.pi, np.pi, 4)
    print(f"Initial parameters: {np.round(x0, 4)}")
    print(f"Initial energy     : {cost(x0):.6f} Ha\n")
    print("Optimizing with COBYLA...")

    res = minimize(cost, x0, method="COBYLA",
                   options={"maxiter": args.maxiter, "disp": False},
                   callback=callback)

    e_vqe = float(res.fun)
    print(f"\nVQE ground-state energy : {e_vqe:.6f} Ha")
    print(f"Exact value             : {e_exact:.6f} Ha")
    print(f"Iterations (callback)     : {len(history)}")
    print(f"Function evaluations    : {res.nfev}")
    err_pct = 100 * abs(e_vqe - e_exact) / abs(e_exact)
    print(f"Relative error          : {err_pct:.3f}%")

    assert e_vqe >= e_exact - 1e-9, \
        "Variational principle violated: VQE below the exact energy!"
    assert err_pct < 2.0, f"VQE did not converge: {err_pct:.2f}% error"
    print("\nOK - VQE converged to the H2 ground state within 2% of exact.")


if __name__ == "__main__":
    main()
