"""
Z_37 conference attack on srg(333,166,82,83)  ==  symmetric conference matrix C(334) == H(668).

Assume a fixed-point-free automorphism sigma of order p=37 acting on the 333 vertices.
333 = 9*37  => 9 orbits of size 37; identify each orbit with Z_37 (sigma = +1).
The adjacency matrix A is then a 9x9 array of 37x37 CIRCULANTS:
  block (i,j) is the circulant with connection set D[i][j] subset of Z_37,
  A[(i,a),(j,b)] = 1  iff  (b-a) mod 37  in  D[i][j].
Symmetry of A (undirected, no loops):  D[j][i] = -D[i][j]  and  0 not in D[i][i],  D[i][i] = -D[i][i].

SRG(333,166,82,83) identity (derived):  A^2 + A - 83 I = 83 J.
Fourier over Z_37 block-diagonalizes: for character t in Z_37 define the 9x9 matrix
  M_t[i][j] = sum_{d in D[i][j]} omega^{t d},  omega = exp(2 pi i/37).
Then A is the SRG  <=>
  (t=0)   R := M_0 (integer orbit matrix, row/col sums 166) satisfies R^2 + R - 83 I_9 = 83*37 J_9 / ... (checked numerically)
  (t!=0)  M_t is Hermitian with  M_t^2 + M_t - 83 I_9 = 0   (all eigenvalues in {(-1 +- sqrt(333))/2}).
This module builds A, checks the SRG identity directly AND via characters, and validates on Paley graphs.
"""
import numpy as np

P = 37
NORB = 9
V = P * NORB   # 333

# ---------- core construction ----------
def build_adjacency(D):
    """D: NORB x NORB list of sets (subsets of Z_P). Returns V x V adjacency (numpy int8)."""
    A = np.zeros((V, V), dtype=np.int64)
    for i in range(NORB):
        for j in range(NORB):
            Dij = D[i][j]
            if not Dij:
                continue
            for a in range(P):
                base = i*P + a
                for d in Dij:
                    b = (a + d) % P
                    A[base, j*P + b] = 1
    return A

def check_symmetry(D):
    for i in range(NORB):
        for j in range(NORB):
            need = set((-d) % P for d in D[j][i])
            if set(D[i][j]) != need:
                return False, (i, j)
        if 0 in D[i][i]:
            return False, ("loop", i)
    return True, None

def orbit_matrix(D):
    return np.array([[len(D[i][j]) for j in range(NORB)] for i in range(NORB)], dtype=np.int64)

# ---------- SRG checks ----------
def srg_identity_direct(A, k, lam, mu):
    """check A^2 + (mu-lam) A - (k-mu) I == mu J  (general SRG identity)."""
    n = A.shape[0]
    lhs = A @ A + (mu - lam)*A - (k - mu)*np.eye(n, dtype=np.int64)
    rhs = mu*np.ones((n, n), dtype=np.int64)
    return np.array_equal(lhs, rhs)

def char_matrices(D, k, lam, mu):
    """build M_t for all t, check Hermitian SRG condition for t!=0 and orbit condition for t=0."""
    w = np.exp(2j*np.pi/P)
    ok = True
    details = []
    I9 = np.eye(NORB)
    for t in range(P):
        M = np.zeros((NORB, NORB), dtype=complex)
        for i in range(NORB):
            for j in range(NORB):
                M[i, j] = sum(w**((t*d) % P) for d in D[i][j])
        if t == 0:
            R = np.rint(M.real).astype(np.int64)
            # t=0 condition: R^2 +(mu-lam)R -(k-mu)I = mu*P*J_9  (J_333 quotient = P*J_9)
            lhs = R @ R + (mu-lam)*R - (k-mu)*np.eye(NORB, dtype=np.int64)
            rhs = mu*P*np.ones((NORB, NORB), dtype=np.int64)
            cond = np.array_equal(lhs, rhs)
        else:
            herm = np.allclose(M, M.conj().T, atol=1e-9)
            resid = M @ M + (mu-lam)*M - (k-mu)*I9
            cond = herm and np.allclose(resid, 0, atol=1e-7)
        details.append((t, cond))
        ok = ok and cond
    return ok, details

# ---------- Paley-graph validation builders ----------
def paley_single_orbit(p):
    """Paley graph P(p) as a single circulant over Z_p (1 orbit). Returns connection set = QR(p)."""
    qr = set((x*x) % p for x in range(1, p))
    return qr

if __name__ == "__main__":
    import sys
    print("=== VALIDATION 1: Paley P(37) single circulant is srg(37,18,8,9) ===")
    # one orbit, p=37: build directly
    qr = paley_single_orbit(37)
    A1 = np.zeros((37, 37), dtype=np.int64)
    for a in range(37):
        for d in qr:
            A1[a, (a+d) % 37] = 1
    sym = np.array_equal(A1, A1.T)
    print(f"  |QR(37)|={len(qr)} (k should be 18), symmetric={sym}")
    print(f"  srg(37,18,8,9) identity A^2+(mu-lam)A-(k-mu)I==muJ : "
          f"{srg_identity_direct(A1,18,8,9)}")

    print("\n=== VALIDATION 2: character/Fourier equivalence on a constructed example ===")
    # build a random symmetric D (NOT an SRG) and confirm direct-check == character-check (both should agree)
    rng = np.random.default_rng(0)
    D = [[set() for _ in range(NORB)] for _ in range(NORB)]
    for i in range(NORB):
        for j in range(i, NORB):
            if i == j:
                # symmetric circulant, no 0: pick pairs {d,-d}
                S = set()
                for d in range(1, P//2+1):
                    if rng.random() < 0.4:
                        S |= {d, (P-d) % P}
                D[i][j] = S
            else:
                S = set(int(x) for x in rng.choice(P, size=rng.integers(0, P), replace=False))
                D[i][j] = S
                D[j][i] = set((-d) % P for d in S)
    symok, where = check_symmetry(D)
    A = build_adjacency(D)
    direct = srg_identity_direct(A, 166, 82, 83)
    via, _ = char_matrices(D, 166, 82, 83)
    print(f"  symmetry ok={symok}; random example is SRG(direct)={direct}; via-characters={via}; "
          f"AGREE={direct==via}  (expect both False, agreeing)")
    # also: A symmetric?
    print(f"  A symmetric={np.array_equal(A,A.T)}; row sums all equal={len(set(A.sum(1)))==1} "
          f"(={A.sum(1)[0]})")
