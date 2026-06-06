"""
Foundations for the LP(333) deep dive. Everything here is EXACT (integers / sympy).
Goal: rigorously verify the spectral framework and extract the cyclotomic constraints
for a Legendre pair of length ell = 333 = 9 * 37, the object equivalent to H(668).

Definitions:
  x in {+-1}^ell. PAF_x(k) = sum_i x_i x_{i+k}  (indices mod ell).
  Legendre pair (LP): u,v in {+-1}^ell with PAF_u(k)+PAF_v(k) = -2 for all k != 0.
  => row sums satisfy (sum u)^2 + (sum v)^2 = 2, i.e. each is +-1.
  PSD_x(s) = |sum_j x_j w^{js}|^2,  w = exp(2 pi i / ell).
  Flat-spectrum form of LP: PSD_u(s)+PSD_v(s) = 2 ell + 2 for all s != 0.
"""
import itertools, numpy as np

def paf(x):
    x = np.asarray(x, dtype=np.int64); L = len(x)
    return np.array([int(np.sum(x*np.roll(x, -k))) for k in range(L)], dtype=np.int64)

def is_lp(u, v):
    L = len(u)
    s = paf(u)+paf(v)
    return all(s[k] == -2 for k in range(1, L)) and s[0] == 2*L

def legendre_seq(p):
    # +-1 sequence: x_0 = 1 (or choice), x_i = legendre(i,p) for i!=0 mapped to +-1, x_0:=1
    qr = set((i*i) % p for i in range(1, p))
    return [1 if (i == 0 or i in qr) else -1 for i in range(p)]

def psd_exact(x, s, L):
    # exact |sum x_j zeta^{js}|^2 using sympy cyclotomics
    import sympy as sp
    w = sp.exp(2*sp.pi*sp.I*s/L)
    val = sum(x[j]*w**j for j in range(L))
    return sp.nsimplify(sp.expand(val*sp.conjugate(val)), rational=True)

def compress(x, m):
    """m-compression to length m: c_r = sum_{j == r (mod m)} x_j  (r in 0..m-1)."""
    c = [0]*m
    for j, xv in enumerate(x):
        c[j % m] += xv
    return c

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("=== sanity: prime length 37, u=v=Legendre is an LP ===")
    p = 37
    u = legendre_seq(p)
    print("  sum u =", sum(u), " PAF_u(1..) all -1? ",
          all(paf(u)[k] == -1 for k in range(1, p)))
    print("  is_lp(u,u) for ell=37:", is_lp(u, u))

    print("\n=== flat-spectrum constant: PSD_u(s)+PSD_v(s) should be 2*ell+2 ===")
    # use a small known LP to verify the spectral identity exactly
    import sympy as sp
    for s in [1, 5, 12]:
        tot = psd_exact(u, s, p) + psd_exact(u, s, p)  # v=u here
        print(f"  ell=37 s={s}: PSD_u+PSD_v = {tot}   (2*ell+2 = {2*p+2})")

    print("\n=== exhaustive LP(9) (9 = 3^2): do full LPs exist, and their 3-compressions ===")
    L = 9
    seqs = []
    for bits in itertools.product([1, -1], repeat=L):
        if sum(bits) in (1, -1):
            seqs.append(bits)
    pafs = {tuple(bits): tuple(int(v) for v in paf(bits)) for bits in seqs}
    target = tuple([2*L] + [-2]*(L-1))
    sols = []
    for a in seqs:
        pa = pafs[a]
        for b in seqs:
            pb = pafs[b]
            if all(pa[k]+pb[k] == target[k] for k in range(L)):
                sols.append((a, b))
    print(f"  #LP(9) ordered pairs (rowsum +-1): {len(sols)}")
    if sols:
        a, b = sols[0]
        print("  example u =", a, " v =", b)
        print("  3-compression of u:", compress(a, 3), " of v:", compress(b, 3))
