"""
Shared quantum-chemistry core for the VQE H2 implementation.

Contents:
  - Minimal Pauli-string algebra (multiply / add / scale dictionaries of
    Pauli strings) so fermionic operators can be built from scratch.
  - Parity mapping of fermionic creation/annihilation operators to qubits.
  - Second-quantized H2 Hamiltonian builder from one- and two-electron
    integrals (physicists' notation).
  - Exact ground-state solver by full diagonalization.
  - The literature 2-qubit H2 Hamiltonian (STO-3G, equilibrium geometry,
    parity mapping with two-qubit reduction) used by scripts 01 and 02.

Conventions: Pauli strings are length-N strings over 'IXYZ' where
character j acts on qubit j. Matrices are built with qubit 0 as the
least-significant factor; eigenvalues are independent of this choice.
"""

import numpy as np

# ---------------------------------------------------------------- Pauli --
_PAULI_MAT = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}

# (a, b) -> (phase, c) for single-qubit Pauli multiplication a*b = phase*c.
_MULT_TABLE = {
    ("I", "I"): (1, "I"), ("I", "X"): (1, "X"), ("I", "Y"): (1, "Y"), ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"), ("X", "X"): (1, "I"), ("X", "Y"): (1j, "Z"), ("X", "Z"): (-1j, "Y"),
    ("Y", "I"): (1, "Y"), ("Y", "X"): (-1j, "Z"), ("Y", "Y"): (1, "I"), ("Y", "Z"): (1j, "X"),
    ("Z", "I"): (1, "Z"), ("Z", "X"): (1j, "Y"), ("Z", "Y"): (-1j, "X"), ("Z", "Z"): (1, "I"),
}


def pauli_mult(s1: str, s2: str):
    """Multiply two Pauli strings. Returns (phase, result_string)."""
    phase = 1 + 0j
    out = []
    for a, b in zip(s1, s2):
        ph, c = _MULT_TABLE[(a, b)]
        phase *= ph
        out.append(c)
    return phase, "".join(out)


def _add_inplace(target: dict, source: dict, scale=1.0):
    for k, v in source.items():
        target[k] = target.get(k, 0) + scale * v
    return target


def dict_add(d1: dict, d2: dict, scale2=1.0) -> dict:
    out = dict(d1)
    return _add_inplace(out, d2, scale2)


def dict_mult(d1: dict, d2: dict) -> dict:
    out = {}
    for s1, c1 in d1.items():
        for s2, c2 in d2.items():
            ph, s3 = pauli_mult(s1, s2)
            out[s3] = out.get(s3, 0) + ph * c1 * c2
    return out


def dict_scale(d: dict, s: complex) -> dict:
    return {k: v * s for k, v in d.items()}


def _single_pauli(n: int, j: int, p: str) -> dict:
    """Pauli operator p acting on qubit j (identity elsewhere)."""
    s = ["I"] * n
    s[j] = p
    return {"".join(s): 1.0}


def _x_tail(n: int, j: int) -> dict:
    """X on every qubit with index > j (identity elsewhere)."""
    s = ["I"] * n
    for k in range(j + 1, n):
        s[k] = "X"
    return {"".join(s): 1.0}


# ------------------------------------------------------- Parity mapping --
def creation(j: int, n: int) -> dict:
    """Fermionic creation operator a†_j in the parity mapping.

    a†_j = 1/2 (X_j Z_{j-1} − i Y_j) ⊗_{k>j} X_k   (Z_{-1} ≡ I for j = 0).
    """
    tail = _x_tail(n, j)
    xj_zjm1 = dict_mult(_single_pauli(n, j, "X"),
                        _single_pauli(n, j - 1, "Z") if j > 0 else {"I" * n: 1.0})
    term1 = dict_mult(xj_zjm1, tail)
    term2 = dict_mult(_single_pauli(n, j, "Y"), tail)
    return dict_add(dict_scale(term1, 0.5), term2, scale2=-0.5j)


def annihilation(j: int, n: int) -> dict:
    """Fermionic annihilation operator a_j in the parity mapping.

    a_j = 1/2 (X_j Z_{j-1} + i Y_j) ⊗_{k>j} X_k.
    """
    tail = _x_tail(n, j)
    xj_zjm1 = dict_mult(_single_pauli(n, j, "X"),
                        _single_pauli(n, j - 1, "Z") if j > 0 else {"I" * n: 1.0})
    term1 = dict_mult(xj_zjm1, tail)
    term2 = dict_mult(_single_pauli(n, j, "Y"), tail)
    return dict_add(dict_scale(term1, 0.5), term2, scale2=+0.5j)


# ------------------------------------------- Fermionic H2 Hamiltonian ----
def fermionic_hamiltonian(h1: np.ndarray, eri: np.ndarray) -> dict:
    """Build the parity-mapped qubit Hamiltonian from MO integrals.

    h1:  (n_spatial, n_spatial) one-electron integrals.
    eri: (n_spatial, n_spatial, n_spatial, n_spatial) two-electron
         integrals in physicists' notation (pq|rs).
    Spin-orbital index: idx = 2 * p + spin (spin 0 = alpha, 1 = beta).

    H = Σ_{pq,s} h_{pq} a†_{ps} a_{qs}
        + 1/2 Σ_{pqrs,s,t} (pq|rs) a†_{ps} a†_{rt} a_{st} a_{qs}
    """
    n_spatial = h1.shape[0]
    n = 2 * n_spatial
    cre = [creation(j, n) for j in range(n)]
    ann = [annihilation(j, n) for j in range(n)]
    H = {}

    def idx(p, s):
        return 2 * p + s

    for p in range(n_spatial):
        for q in range(n_spatial):
            if h1[p, q] == 0:
                continue
            for s in (0, 1):
                term = dict_mult(cre[idx(p, s)], ann[idx(q, s)])
                _add_inplace(H, term, scale=h1[p, q])
    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s_ in range(n_spatial):
                    g = eri[p, q, r, s_]
                    if g == 0:
                        continue
                    for s1 in (0, 1):
                        for s2 in (0, 1):
                            term = dict_mult(
                                dict_mult(cre[idx(p, s1)], cre[idx(r, s2)]),
                                dict_mult(ann[idx(s_, s2)], ann[idx(q, s1)]))
                            _add_inplace(H, term, scale=0.5 * g)
    # Drop numerical noise.
    return {k: complex(v) for k, v in H.items() if abs(v) > 1e-12}


def pauli_dict_to_matrix(d: dict, n: int) -> np.ndarray:
    """Convert a Pauli-string dictionary to a dense matrix."""
    dim = 2 ** n
    M = np.zeros((dim, dim), dtype=complex)
    for s, c in d.items():
        mat = _PAULI_MAT[s[0]]
        for ch in s[1:]:
            mat = np.kron(_PAULI_MAT[ch], mat)
        M += c * mat
    return M


def exact_ground_state_energy(h1: np.ndarray, eri: np.ndarray,
                              nuclear_repulsion: float) -> float:
    """Exact (FCI) ground-state energy by full diagonalization."""
    n = 2 * h1.shape[0]
    H = fermionic_hamiltonian(h1, eri)
    M = pauli_dict_to_matrix(H, n)
    M = (M + M.conj().T) / 2  # symmetrize away floating-point dust
    return float(np.linalg.eigvalsh(M)[0]) + nuclear_repulsion


# --------------------------------- Literature 2-qubit H2 Hamiltonian -----
def h2_two_qubit_paulis() -> dict:
    """H2 STO-3G at the equilibrium bond length (R = 0.735 Å).

    Parity mapping with two-qubit reduction (particle-number and spin
    Z2 symmetries). Coefficients derived from first principles with
    PySCF (RHF/STO-3G integrals) + parity mapping + numerical tapering;
    the resulting operator reproduces the FCI/STO-3G energy -1.1373 Ha:

        H = -0.332404·II + 0.397937·ZI − 0.397937·IZ
            − 0.011280·ZZ + 0.180931·XX

    (qubit 0 = rightmost character in each key; the II coefficient
    includes the nuclear repulsion energy.)
    """
    return {
        "II": -0.332404,
        "ZI": +0.397937,
        "IZ": -0.397937,
        "ZZ": -0.011280,
        "XX": +0.180931,
    }


def h2_two_qubit_matrix() -> np.ndarray:
    """Dense 4×4 matrix of the literature 2-qubit H2 Hamiltonian."""
    return pauli_dict_to_matrix(h2_two_qubit_paulis(), 2).real


def exact_two_qubit_energy() -> float:
    """Exact ground-state energy of the literature 2-qubit Hamiltonian."""
    return float(np.linalg.eigvalsh(h2_two_qubit_matrix())[0])


# Literature FCI/STO-3G ground-state energy of H2 at equilibrium, Ha.
H2_FCI_LITERATURE_HA = -1.1373
