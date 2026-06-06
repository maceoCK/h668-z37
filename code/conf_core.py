"""
Parametric block-circulant SRG / conference-graph machinery.
A vertex set = NORB orbits of size P (an order-P semiregular automorphism, sigma=+1 on each orbit).
Adjacency A is a NORB x NORB array of P x P circulants with connection sets D[i][j] subset Z_P.
A[(i,a),(j,b)] = 1 iff (b-a) mod P in D[i][j].
Undirected, no loops:  D[j][i] = -D[i][j];  0 not in D[i][i]; D[i][i] = -D[i][i].

SRG(v,k,lam,mu) identity:  A^2 + (mu-lam)A - (k-mu)I = mu J.
Fourier over Z_P:  M_t[i][j] = sum_{d in D[i][j]} w^{t d},  w=exp(2pi i/P).
 A is the SRG  <=>  (t=0) R=M_0 (integer) with R^2+(mu-lam)R-(k-mu)I_NORB = mu*P*J_NORB
                    (t!=0) M_t Hermitian with M_t^2+(mu-lam)M_t-(k-mu)I_NORB = 0.
"""
import numpy as np

def build_adjacency(D, p, norb):
    V = p*norb
    A = np.zeros((V, V), dtype=np.int64)
    for i in range(norb):
        for j in range(norb):
            for a in range(p):
                base = i*p + a
                row = A[base]
                jb = j*p
                for d in D[i][j]:
                    row[jb + (a+d) % p] = 1
    return A

def check_symmetry(D, p, norb):
    for i in range(norb):
        for j in range(norb):
            if set(D[i][j]) != set((-d) % p for d in D[j][i]):
                return False
        if 0 in D[i][i] or set(D[i][i]) != set((-d) % p for d in D[i][i]):
            return False
    return True

def orbit_matrix(D, norb):
    return np.array([[len(D[i][j]) for j in range(norb)] for i in range(norb)], dtype=np.int64)

def srg_identity_direct(A, k, lam, mu):
    n = A.shape[0]
    lhs = A @ A + (mu-lam)*A - (k-mu)*np.eye(n, dtype=np.int64)
    return np.array_equal(lhs, mu*np.ones((n, n), dtype=np.int64))

def srg_via_characters(D, p, norb, k, lam, mu, atol=1e-7):
    w = np.exp(2j*np.pi/p)
    I = np.eye(norb)
    for t in range(p):
        M = np.zeros((norb, norb), dtype=complex)
        for i in range(norb):
            for j in range(norb):
                if D[i][j]:
                    M[i, j] = sum(w**((t*d) % p) for d in D[i][j])
        if t == 0:
            R = np.rint(M.real).astype(np.int64)
            lhs = R @ R + (mu-lam)*R - (k-mu)*np.eye(norb, dtype=np.int64)
            if not np.array_equal(lhs, mu*p*np.ones((norb, norb), dtype=np.int64)):
                return False
        else:
            if not np.allclose(M, M.conj().T, atol=atol):
                return False
            if not np.allclose(M@M + (mu-lam)*M - (k-mu)*I, 0, atol=atol):
                return False
    return True

# ---------- GF(p^2) and Paley P(p^2) in block-circulant form (for validation) ----------
def gf2_setup(p):
    # smallest non-residue r0; element a+b*theta, theta^2=r0
    sq = set((x*x) % p for x in range(p))
    r0 = next(x for x in range(1, p) if x not in sq)
    def mul(u, v):
        (a, b), (c, d) = u, v
        return ((a*c + b*d*r0) % p, (a*d + b*c) % p)
    def is_square(u):
        if u == (0, 0):
            return False
        e = (p*p - 1)//2
        r = (1, 0)
        base = u
        while e:
            if e & 1: r = mul(r, base)
            base = mul(base, base); e >>= 1
        return r == (1, 0)
    return r0, mul, is_square

def paley_pp_blocks(p):
    """Paley graph on GF(p^2): vertex (i,a) ~ field elt theta*i + a, sigma = +1 (add 1 in 'a').
       Returns connection sets D[i][j] subset Z_p (norb=p)."""
    r0, mul, is_square = gf2_setup(p)
    # field elt of vertex (i,a) = a*(1,0) + i*(0,1) = (a, i)   [theta=(0,1)]
    D = [[set() for _ in range(p)] for _ in range(p)]
    for i in range(p):
        for j in range(p):
            for d in range(p):
                # difference of (j, ?) - (i, ?): elt diff = (d, j-i)  (a-component diff d, theta-comp j-i)
                diff = (d % p, (j - i) % p)
                if diff != (0, 0) and is_square(diff):
                    D[i][j].add(d)
    return D
