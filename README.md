# vqe-h2-implementation

A from-scratch Variational Quantum Eigensolver (VQE) for the hydrogen
molecule (the "hello world" of quantum chemistry) built with Qiskit,
SciPy, and PySCF. No black-box VQE routine: the ansatz, the expectation
values, and the optimization loop are all explicit.

This one connects my quantum hobby with my interest in molecular engineering.

## The three implementations

| # | Script | What it does |
|---|--------|--------------|
| 1 | `01_h2_hamiltonian.py` | Builds the 2-qubit H2 Hamiltonian from its literature Pauli coefficients (hardcoded in `scripts/common.py`); exact diagonalization vs the published FCI value |
| 2 | `02_vqe_ground_state.py` | VQE from scratch: hardware-efficient Ry ansatz, Aer statevector expectations, scipy COBYLA loop |
| 3 | `03_dissociation_curve.py` | Exact FCI/STO-3G energies across bond lengths 0.4–2.5 Å (PySCF integrals, own parity mapping, exact diagonalization, **not** VQE) |

## Architecture

```
 PySCF (RHF/STO-3G)          shared core                 Qiskit Aer
┌──────────────────┐   ┌──────────────────────┐   ┌──────────────────────┐
│ 1e/2e integrals  │──▶│ parity-mapped qubit  │──▶│ statevector simulator│
│ at each geometry │   │ Hamiltonian          │   │ evaluates ⟨ψ(θ)|H|ψ⟩ │
└──────────────────┘   └──────────────────────┘   └──────────────────────┘
        │                        │                            │
        │                        ▼                            ▼
        │               exact diagonalization        scipy COBYLA updates θ
        │               (scripts 01, 03)             (script 02)
```

The 2-qubit Hamiltonian (R = 0.735 Å, parity mapping + two-qubit
reduction; II term includes nuclear repulsion):

    H = −0.332404·II + 0.397937·ZI − 0.397937·IZ − 0.011280·ZZ + 0.180931·XX

Coefficients are the published literature values for H2/STO-3G at
equilibrium (parity mapping + two-qubit reduction), hardcoded in
`scripts/common.py`, not recomputed from integrals in scripts 01/02.
See [docs/architecture.md](docs/architecture.md).

## Quickstart

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

bash run_all.sh
```

## Validation

Measured on Qiskit Aer 2.6.0 (statevector, seeded) + PySCF 2.14.0, 2026-10-21:

| Implementation | Prediction | Measured |
|------------|-----------|----------|
| H2 Hamiltonian (script 01) | exact energy matches published FCI/STO-3G (−1.1373 Ha) | **−1.137305 Ha**: 0.01 mHa from literature, assertion passed |
| VQE ground state (script 02) | converges within 2% of exact, never below it | **−1.137305 Ha** vs exact −1.137305 Ha: **0.000%** relative error, 49 iterations, 70 function evaluations, assertion passed. Note: expectations are evaluated on the **exact statevector simulator** (no shot noise), which is why the error is 0.000%. On real hardware it would not be. The 2% assert threshold is loose (≈22.7 mHa); chemical accuracy is ~1.6 mHa |
| Dissociation curve (script 03) | minimum near the equilibrium bond length (~0.74 Å) | minimum at **R = 0.750 Å, E = −1.137117 Ha**; 2.5 Å limit −0.936055 Ha (bound by 201.1 mHa): assertion passed |

Script 03 does **not** run VQE at any geometry: it takes PySCF
integrals, applies this repo's own parity mapping, and diagonalizes
exactly (FCI/STO-3G). VQE is the approximate method under test in
script 02; exact diagonalization is the yardstick.

## Why this matters (the one-paragraph version)

Finding molecular ground states is the core computational task of
chemistry: reaction rates, battery voltages, and material properties
all flow from it. VQE splits the work the way near-term quantum
computing demands: the quantum processor prepares and measures trial
states, a classical optimizer steers. H2 is where every quantum
chemistry method gets debugged first, because it's the smallest
molecule with electron correlation and we know its exact answer.

## Plain-language guide

[simplified-explanation.md](simplified-explanation.md): VQE and the
hydrogen molecule explained without scary math: the guitar tuner, the
translation from chemistry to qubits, and stretching the molecule
(Spanish).
