# Architecture

## What this implementation is

A from-scratch Variational Quantum Eigensolver (VQE) for the hydrogen
molecule — the "hello world" of quantum chemistry — built with Qiskit,
SciPy, and PySCF. No black-box VQE routine: the ansatz, the expectation
values, and the optimization loop are all explicit.

## Layout

```
vqe-h2-implementation/
├── README.md
├── requirements.txt           # qiskit + qiskit-aer + numpy + scipy + pyscf
├── run_all.sh
├── scripts/
│   ├── common.py              # Pauli algebra, parity mapper, fermionic
│   │                          # Hamiltonian builder, exact solver,
│   │                          # literature 2-qubit H2 Hamiltonian
│   ├── 01_h2_hamiltonian.py   # 2-qubit H2 from Pauli coefficients;
│   │                          # exact diagonalization vs literature
│   ├── 02_vqe_ground_state.py # VQE: Ry ansatz + Aer statevector +
│   │                          # scipy COBYLA
│   └── 03_dissociation_curve.py  # FCI/STO-3G energies vs bond length
│                              # (exact diagonalization, PySCF integrals)
├── docs/
│   ├── architecture.md        # this file
│   └── interview-answers.md
└── simplified-explanation.md  # plain-language guide (Spanish)
```

## The pipeline

```
 PySCF (RHF/STO-3G)          common.py                     Qiskit Aer
┌──────────────────┐   ┌──────────────────────┐   ┌──────────────────────┐
│ 1e/2e integrals  │──▶│ parity-mapped qubit  │──▶│ statevector simulator│
│ at each geometry │   │ Hamiltonian          │   │ evaluates ⟨ψ(θ)|H|ψ⟩ │
└──────────────────┘   └──────────────────────┘   └──────────────────────┘
        │                        │                            │
        │                        ▼                            ▼
        │               exact diagonalization        scipy COBYLA updates θ
        │               (scripts 01, 03)             (script 02)
```

## Script 01 — the Hamiltonian

The 2-qubit H2 Hamiltonian at R = 0.735 Å, parity mapping with
two-qubit reduction (particle-number and spin Z2 symmetries):

    H = −0.332404·II + 0.397937·ZI − 0.397937·IZ − 0.011280·ZZ + 0.180931·XX

The coefficients were derived from first principles (PySCF integrals →
parity mapping → numerical tapering), not copied: the II term carries
the nuclear repulsion energy. Diagonalizing gives −1.137305 Ha,
matching the published FCI/STO-3G value (−1.1373 Ha) to 0.01 mHa.

## Script 02 — the VQE

- **Ansatz** (hardware-efficient): Ry(θ₀)⊗Ry(θ₁) → CNOT(0,1) →
  Ry(θ₂)⊗Ry(θ₃). Four parameters, one entangling layer.
- **Quantum part**: for each θ, the circuit runs on the Aer
  statevector simulator and ⟨ψ(θ)|H|ψ(θ)⟩ is evaluated exactly.
- **Classical part**: scipy COBYLA minimizes the energy; every
  iteration is printed.
- The variational principle is asserted, not assumed: the script
  checks the VQE energy never drops below the exact value and lands
  within 2% of it.

## Script 03 — the dissociation curve

For each bond length, PySCF recomputes the RHF/STO-3G integrals, the
shared core rebuilds the parity-mapped Hamiltonian, and exact
diagonalization gives the FCI energy. These points are **exact
diagonalization, not VQE** — the honest label matters: VQE is the
approximate method being validated, exact diagonalization is the
yardstick. The curve's minimum must sit at the equilibrium bond
length (~0.74 Å).

## Design choices

- **Parity mapping written by hand** (verified against the fermionic
  anticommutation relations as matrices) instead of importing a
  mapper: the mapping is the physics, and hiding it would defeat the
  purpose of rebuilding it.
- **Seeded RNG + seeded simulator**: script 02 reproduces bit-for-bit.
- **Assertions as verdicts**: each script encodes its physics claim
  (literature match, variational principle + 2%, interior minimum),
  so "it works" is a passing check, not commentary.

## What this is not

Not a production quantum-chemistry package: STO-3G is a minimal basis,
there is no error mitigation, and the simulator is noise-free. It is
a faithful operational model of the VQE loop — the same loop that
runs on real hardware, minus the noise.
