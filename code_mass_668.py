"""
Code-mass-formula machinery instantiated at H(668), 668 = 4*167 == 4 (mod 8).

THEOREMS BEING VERIFIED (proofs in analysis; elementary, hold for EVERY Hadamard
matrix H of order n = 4m, m odd):

 T1. (binary non-self-orthogonality) With H normalized (first row all +1),
     B = (H+J)/2 over F2 has pairwise row inner products  b_i . b_j = n/4 = m (mod 2)
     for i,j >= 2 distinct.  At 668: m = 167 == 1, so the binary row span is NOT
     self-orthogonal -> the Ozeki/Tonchev doubly-even pipeline (needs 8|n) cannot start.

 T2. (binary collapse) det((H+J)/2) = +-2 * m^(n/2)  [since det(H+J) = 2 det H for
     row-normalized H], so v_2(det B) = 1, so rank_F2(B) = n-1 exactly; the all-ones
     vector is in the right kernel; each row has even weight n/2.  Hence
     rowspan_F2((H+J)/2) = E_n, the even-weight code -- INDEPENDENT of H.

 T3. (design code trivial) The 2-(n-1, n/2-1, n/4-1) design (n-1 = 667) has
     N N^T = m I + (m-1) J == I (mod 2) for m odd, so its F2 code is ALL of F2^(n-1).

 T4. (Z4 collapse) The Z4 row span of H:
     - self-orthogonal since h_i.h_j = n delta_ij == 0 mod 4;
     - SNF of H is forced: d_1=1, d_n=n, d_i d_{n+1-i}=n, mod-2 rank of H is 1
       (H == J mod 2), v_2(det) = n  ==>  v_2 pattern (0, 1^(n-2), 2); for m prime
       the full SNF is (1, 2^(n/2-1), (2m)^(n/2-1), 4m);
     - hence |C| = 4 * 2^(n-2) = 2^n = 4^(n/2): C is SELF-DUAL over Z4;
     - every row lies in Klemm's code K_n = {x in Z4^n : all coords congruent mod 2,
       sum(x) == 0 mod 4} (entries odd; integer row sums are n or 0 == 0 mod 4);
       |K_n| = 2^n  ==>  C = K_n exactly.  INDEPENDENT of H.  Type I (rows have
       Euclidean weight n == 4 mod 8), so the 'Type II Z4 needs 8|n' gate never fires.

 CONTRAST (n == 0 mod 8, Sylvester H_16): there the binary span IS self-orthogonal
 doubly-even -- showing exactly which residue class the 8|n machinery lives in.
"""
import numpy as np
from math import gcd
from itertools import product
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form

# ---------------------------------------------------------------- constructions
def paley(q):
    """Paley type I Hadamard matrix of order n=q+1, q prime == 3 mod 4."""
    n = q + 1
    sq = {(i * i) % q for i in range(1, q)}
    def chi(x):
        x %= q
        return 0 if x == 0 else (1 if x in sq else -1)
    S = np.zeros((n, n), dtype=int)
    S[0, 1:] = 1
    S[1:, 0] = -1
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                S[i, j] = chi(j - i)
    H = S + np.eye(n, dtype=int)
    assert (H @ H.T == n * np.eye(n, dtype=int)).all(), "not Hadamard"
    return H

def sylvester(k):
    H = np.array([[1]])
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    return H

def normalize_first_row(H):
    """Negate columns so first row is all +1 (monomial equiv)."""
    return H @ np.diag(H[0])

def full_normalize(H):
    H = normalize_first_row(H)
    D = np.diag(H[:, 0])
    return D @ H  # first column now +1 too; first row untouched

# ---------------------------------------------------------------- F2 utilities
def gf2_rank(A):
    A = (A % 2).astype(np.int8).copy()
    nr, nc = A.shape
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nr) if A[i, c]), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(nr):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
    return r

# ---------------------------------------------------------------- per-order check
def check_order(H, label):
    n = H.shape[0]; m = n // 4
    assert n % 8 == 4 and m % 2 == 1
    H = normalize_first_row(H)
    B = (H + 1) // 2                      # (H+J)/2, 0/1 matrix
    # T1: Gram over F2
    G = (B @ B.T) % 2
    exp = np.full((n, n), m % 2, dtype=int)
    exp[0, :] = 0; exp[:, 0] = 0
    np.fill_diagonal(exp, 0)              # weights n (row0) and n/2 are even
    assert (G == exp).all(), "T1 Gram pattern failed"
    # T2: rank n-1, kernel = <1>, rows even weight  => code = E_n
    assert gf2_rank(B) == n - 1, "T2 rank failed"
    assert ((B.sum(axis=1)) % 2 == 0).all()
    assert ((B @ np.ones(n, dtype=int)) % 2 == 0).all()
    detH = int(Matrix(H.tolist()).det())
    detHpJ = int(Matrix((H + 1).tolist()).det())   # det(H+J) (J=all-ones)
    assert detHpJ == 2 * detH, "det(H+J) = 2 det H failed"
    # T3: design code full
    Hn = full_normalize(H)
    core = Hn[1:, 1:]
    N = (1 + core) // 2                  # +1 positions: incidence of 2-(n-1, 2m-1, m-1)
    NNt = N @ N.T
    assert (NNt == m * np.eye(n - 1, dtype=int) + (m - 1)).all(), "design params"
    assert gf2_rank(N) == n - 1, "T3 full 2-rank failed"
    # T4: SNF forced; |C_Z4| = 2^n; rows in K_n
    snf = smith_normal_form(Matrix(H.tolist()))
    d = [abs(int(snf[i, i])) for i in range(n)]
    expected = sorted([1] + [2] * (n // 2 - 1) + [2 * m] * (n // 2 - 1) + [4 * m])
    assert sorted(d) == expected, f"SNF {sorted(d)} != {expected}"
    size = 1
    for di in d:
        size *= 4 // gcd(di, 4)
    assert size == 2 ** n, "Z4 span size != 2^n"
    assert (np.abs(H) == 1).all()                          # all entries odd in Z4
    assert ((H.sum(axis=1)) % 4 == 0).all(), "row sums != 0 mod 4"
    print(f"  {label}: n={n}, m={m} | Gram=J-I pattern OK | rank2(B)={n-1} -> code=E_{n} | "
          f"design 2-rank={n-1} (full) | SNF=(1,2^{n//2-1},{2*m}^{n//2-1},{4*m}) | "
          f"|C_Z4|=2^{n}, rows in K_{n} -> C=K_{n}")

# ---------------------------------------------------------------- Z4 closure n=12
def z4_closure_check():
    H = normalize_first_row(paley(11))
    n = 12
    gens = [tuple(int(x) % 4 for x in row) for row in H]
    S = {tuple([0] * n)}
    for g in gens:
        S = {tuple((s[i] + k * g[i]) % 4 for i in range(n)) for s in S for k in range(4)}
    assert len(S) == 2 ** n, f"closure size {len(S)} != {2**n}"
    # K_n membership for every element
    for w in S:
        ps = {x % 2 for x in w}
        assert len(ps) == 1, "coords not parity-uniform"
        assert sum(w) % 4 == 0, "sum != 0 mod 4"
    odd = sum(1 for w in S if w[0] % 2 == 1)
    evens = [w for w in S if w[0] % 2 == 0]
    # swe check: even part = {2x : x even weight}: #(k coords equal 2) = C(12,k), k even
    from math import comb
    cnt = {}
    for w in evens:
        k = sum(1 for x in w if x == 2)
        assert all(x in (0, 2) for x in w)
        cnt[k] = cnt.get(k, 0) + 1
    assert all(cnt.get(k, 0) == (comb(12, k) if k % 2 == 0 else 0) for k in range(13))
    assert odd == 2 ** 11
    print(f"  Z4 closure of H_12: |C| = {len(S)} = 2^12; all words in K_12; "
          f"odd words = {odd} = 2^11 (each swe b^12); even part = 2(E_12) binomial -> "
          f"swe = 1/2[(a+c)^12 + (a-c)^12] + 2^11 b^12  VERIFIED == K_12")

# ---------------------------------------------------------------- contrast n=16
def contrast_16():
    H = normalize_first_row(sylvester(4))
    B = (H + 1) // 2
    G = (B @ B.T) % 2
    assert (G == 0).all(), "H_16 binary span should be self-orthogonal"
    r = gf2_rank(B)
    wts = sorted({int(w) for w in B.sum(axis=1)})
    assert r <= 8
    print(f"  CONTRAST Sylvester H_16 (n==0 mod 8): Gram == 0 (self-orthogonal!), "
          f"rank2 = {r} <= 8, row weights {wts} all == 0 mod 4 -> doubly-even pipeline "
          f"EXISTS at 0 mod 8 and provably NOT at 4 mod 8.")

# ---------------------------------------------------------------- 668 arithmetic
def at_668():
    n, m = 668, 167
    print("\n=== EXACT ARITHMETIC AT 668 ===")
    print(f"  668 = 4*167;  668 mod 8 = {n % 8};  m = 167, m mod 2 = {m % 2}")
    print(f"  T1: pairwise F2 inner product of rows of (H+J)/2 = n/4 mod 2 = {(-(-n // 4)) % 2}"
          f"  -> NOT self-orthogonal over F2")
    v2detB = 1 + (n // 2) * 2 - n
    print(f"  T2: det((H+J)/2) = +-2*167^334, v_2 = {v2detB} -> rank_F2 = {n-1} -> code = E_668 (forced)")
    print(f"  T3: N N^T = 167 I + 166 J == I (mod 2) [167 odd, 166 even] -> F2 design code = F2^667 (forced)")
    print(f"  T4: SNF(H) forced = (1, 2^333, 334^333, 668); |C_Z4| = 4*2^666 = 2^668 = 4^334")
    print(f"      rows: Euclidean weight = 668 == {n % 8} (mod 8) -> Type I, NOT Type II")
    print(f"      C_Z4 = K_668 (Klemm), swe = 1/2[(a+c)^668+(a-c)^668] + 2^667 b^668")
    print(f"  Gleason gates: doubly-even binary self-dual length 668 exist iff 8|668: {668 % 8 == 0}")
    print(f"                 Type II Z4 length 668 exist iff 8|668: {668 % 8 == 0}")
    print(f"      -> both families EMPTY at 668, but H(668) provably produces NEITHER (T1, T4).")
    # Type I binary self-dual at 668 exists: <(1,1)>^334. Mass:
    import sys
    sys.set_int_max_str_digits(100000)
    Nsd = 1
    for i in range(1, 334):
        Nsd *= (2 ** i + 1)
    print(f"  # binary self-dual codes length 668 = prod_(i=1..333)(2^i+1): {len(str(Nsd))} digits > 0 (family nonempty)")
    # F167 self-dual length 668: q=167==3 mod 4 -> exist iff 4|n
    print(f"  F167 self-dual [668,334]: q mod 4 = {167 % 4}, need 4|n: {n % 4 == 0} -> EXISTS;")
    a, b = next((a, b) for a in range(167) for b in range(167) if (a*a + b*b + 1) % 167 == 0)
    print(f"      explicit witness a^2+b^2=-1 mod 167: a={a}, b={b}; "
          f"[4,2] gen [I | A], A=[[{a},{b}],[{b},{-a % 167}]], direct-sum 167 copies -> length 668")
    print(f"  Construction A4(K_668) = D_668^+ (odd unimodular; even would need 8|n) -- fixed lattice, H-independent")

if __name__ == "__main__":
    print("=== SMALL-ORDER MACHINE VERIFICATION (n == 4 mod 8) ===")
    for q in (11, 19, 27):
        if q == 27:
            # q=27 not prime; use GF(27)? skip -> use q=43 instead handled below
            continue
        check_order(paley(q), f"Paley H_{q+1}")
    # order 28 via Paley q=27 needs GF(27); instead verify a third order with q=43 (n=44)
    check_order(paley(43), "Paley H_44")
    # scrambled equivalence test on H_20: random monomial transform, re-normalize
    rng = np.random.default_rng(668)
    H = paley(19)
    P = np.eye(20, dtype=int)[rng.permutation(20)]
    Q = np.eye(20, dtype=int)[rng.permutation(20)]
    D1 = np.diag(rng.choice([-1, 1], 20)); D2 = np.diag(rng.choice([-1, 1], 20))
    check_order(P @ D1 @ H @ D2 @ Q, "scrambled H_20 (random monomial equiv)")
    # transpose
    check_order(paley(19).T, "H_20 transpose")
    z4_closure_check()
    contrast_16()
    at_668()
    print("\nALL ASSERTIONS PASSED")
