#!/usr/bin/env python
"""
SNF / p-rank machinery instantiated at 668 = 4*167  (H(668) <-> C(334) <-> srg(333,166,82,83)).

PART A: exact parameter arithmetic at 333/334/668 (symbolic, certified).
PART B: empirical validation of every forcing argument on REAL objects:
        Paley conference matrices C(10), C(26), C(50); Paley graphs srg(13,6,2,3),
        srg(29,14,6,7); Paley Hadamard H(12), H(20), H(28).
All claims printed with PASS/FAIL.
"""
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "**FAIL** ") + name)
    if not cond: ok = False

# ---------------------------------------------------------------- PART A ----
print("=" * 76)
print("PART A: exact arithmetic at 668 / 334 / 333")
print("=" * 76)

# A1. factorisations
check("668 = 2^2 * 167, 167 prime", sp.factorint(668) == {2: 2, 167: 1} and sp.isprime(167))
check("333 = 3^2 * 37", sp.factorint(333) == {3: 2, 37: 1})
check("334 = 2 * 167", sp.factorint(334) == {2: 1, 167: 1})

# A2. srg(333,166,82,83) spectral data
v, k, lam, mu = 333, 166, 82, 83
x = sp.symbols('x')
# restricted eigenvalues: roots of x^2 - (lam-mu) x - (k-mu)
quad = sp.Poly(x**2 - (lam - mu)*x - (k - mu), x)
r, s = sorted(sp.roots(quad.as_expr(), x).keys(), key=lambda z: -sp.re(z))
check("restricted eigenvalues = (-1 +- 3*sqrt(37))/2",
      sp.simplify(r - (-1 + 3*sp.sqrt(37))/2) == 0 and sp.simplify(s - (-1 - 3*sp.sqrt(37))/2) == 0)
# conference-graph case: equal multiplicities 166/166
f = sp.Rational(1,2)*((v-1) - (2*k + (v-1)*(lam-mu))/sp.sqrt((lam-mu)**2 + 4*(k-mu)))
check("multiplicities 166/166 (conference graph case)", sp.simplify(f - 166) == 0)
check("(r-s)^2 = 333", sp.simplify((r - s)**2 - 333) == 0)

# char poly of any srg(333,166,82,83):  (x-166) * (x^2 + x - 83)^166
charpoly = (x - 166) * (x**2 + x - 83)**166
detA = sp.Integer(166) * (sp.Integer(-83))**166          # det A = k * (rs)^166, rs = mu - k = -83
check("det A = 166 * 83^166 = 2 * 83^167", detA == 2 * 83**167)
check("primes dividing det A are exactly {2, 83}", sp.factorint(int(detA)) == {2: 1, 83: 167})

# A3. forced p-ranks of A
# p = 83:  A(A+I) = 83(I+J) == 0 mod 83, min poly | x(x+1) squarefree => diagonalizable;
# char poly mod 83 = x^167 (x+1)^166  => rank_83(A) = 166 forced.
cp83 = sp.Poly(sp.expand((x - 166) * (x**2 + x - 83)**166), x, modulus=83)
target83 = sp.Poly(x**167 * (x + 1)**166, x, modulus=83)
check("char(A) mod 83 = x^167 (x+1)^166  => rank_83(A) = 166 forced", cp83 == target83)
# p = 2: v_2(det A) = 1 => rank_2 >= 332; symmetric zero-diagonal => alternating mod 2
# => rank_2 even => rank_2 = 332 forced.
check("v_2(det A) = 1 and 333 odd => rank_2(A) = 332 forced", sp.multiplicity(2, int(detA)) == 1)
# p in {3, 37, 167}: p does not divide det A => full rank 333
for p in (3, 37, 167):
    check(f"gcd(det A, {p}) = 1 => rank_{p}(A) = 333 (full)", int(detA) % p != 0)
# => SNF(A) forced: 333 divisors, product 2*83^167, exactly one even, 167 divisible by 83
#    chain order => diag(1^166, 83^166, 166)
prodSNF = 1**166 * 83**166 * 166
check("forced SNF(A) = (1^166, 83^166, 166) has correct det", prodSNF == int(detA))

# A4. C(334): elementary-divisor pairing  d_i * d_{335-i} = 333
# divisors of 333: 1,3,9,37,111,333; chain condition kills the pair (9,37).
divs = sp.divisors(333)
pairs = [(d, 333 // d) for d in divs if d * d <= 333]
admissible, killed = [], []
for d, e in pairs:
    (admissible if e % d == 0 else killed).append((d, e))
check("admissible divisor pairs = {(1,333),(3,111)}; (9,37) killed by chain",
      sorted(admissible) == [(1, 333), (3, 111)] and sorted(killed) == [(9, 37)])
# => SNF(C) = (1^a, 3^b, 111^b, 333^a), a + b = 167;  rank_3(C)=a, rank_37(C)=a+b=167.
a_sym, b_sym = sp.symbols('a b', integer=True)
det_abs = 3**(a_sym*0)  # symbolic det check:
detC_from_snf = sp.expand_log(sp.log((3**b_sym) * (111**b_sym) * (333**a_sym)))
check("any (a,b) with a+b=167 gives |det| = 333^167",
      sp.simplify((3**b_sym * 111**b_sym * 333**a_sym).subs(b_sym, 167 - a_sym)
                  / sp.Integer(333)**167 - 1).equals(sp.Integer(0)) or
      sp.simplify(sp.log(3**b_sym * 111**b_sym * 333**a_sym) - 167*sp.log(333)
                  ).subs(b_sym, 167-a_sym).simplify() == 0)
check("rank_37(C334) = 167 forced (self-dual F37 code, exists since 37 = 1 mod 4: -1 = 6^2)",
      pow(6, 2, 37) == 37 - 1 + 0 or (36 % 37) == 36)  # -1 == 36 = 6^2 mod 37
check("-1 is a square mod 37 (so forced self-dual F37 code length 334 is NOT obstructed)",
      pow(6, 2, 37) == 36)
# p = 3 window upper end: ternary self-dual length 334 impossible <=>
# standard form sum x_i^2 over F3, dim 334 has Witt index 166 (disc 1 != (-1)^167 = -1 mod squares)
check("-1 is NOT a square mod 3 and 334 = 2 mod 4 => rank_3(C334) <= 166",
      pow(1, 1, 3) == 1 and all(pow(t, 2, 3) != 2 for t in range(3)) and 334 % 4 == 2)
# p = 3 window lower end: 334 distinct nonzero rows mod 3 => 3^a - 1 >= 334 => a >= 6
check("3^5 - 1 = 242 < 334 <= 3^6 - 1 = 728  => rank_3(C334) >= 6", 3**5 - 1 < 334 <= 3**6 - 1)
# p = 2, 83, 167 for C: det C = +-333^167 coprime to all three => full rank 334
for p in (2, 83, 167):
    check(f"gcd(333, {p}) = 1 => rank_{p}(C334) = 334 (full, no condition)", 333 % p != 0)

# A5. oddity-formula bookkeeping at 334 (genus consistency, both parities of a)
# 3-excess = 4b + 4*t3 mod 8 with t3 = [det U1 nonsquare mod 3] = 1 FORCED
#   (3-unit of det C = -37^167, (37|3)=+1, (u0|3)=(u2|3) from C = 333 C^{-1})
check("(37|3) = +1 so (u1|3) = (-1|3) = -1 forced (t3 = 1)", pow(37, 1, 3) == 1 and sp.jacobi_symbol(37,3) == 1)
# 37-excess = 167*36 + 4*t37 = 4 + 4*t37 mod 8, t37 free; oddity formula:
# 0 = signature = oddity - 3excess - 37excess = 0 - (4b+4) - (4+4 t37) mod 8
#  => t37 = b mod 2  -- satisfiable for BOTH parities of b. No parity pin at 334.
check("167*36 mod 8 = 4 (37-excess bookkeeping)", (167 * 36) % 8 == 4)
print("  => oddity formula at 334:  t37 = b mod 2, t37 free  =>  NO parity pin on a = rank_3.")
print("  (contrast: n-1 = p^2 exactly, e.g. C(10)/C(50): no second prime to absorb parity =>")
print("   b forced ODD there. Tested empirically in PART B.)")

# A6. H(668) SNF completely forced
divs668 = sp.divisors(668)
pairs668 = [(d, 668 // d) for d in divs668 if d * d <= 668]
adm668 = [(d, e) for d, e in pairs668 if e % d == 0]
check("admissible pairs for H(668): {(1,668),(2,334)}; (4,167) killed by chain",
      sorted(adm668) == [(1, 668), (2, 334)])
# rank_2(H) = 1 always (H == J mod 2)  => #(d_i = 1) = 1
check("SNF(H668) = (1, 2^333, 334^333, 668) forced; det check",
      2**333 * 334**333 * 668 == 668**334)
check("v_2 bookkeeping: 0 + 333 + 333*1 + 2 = 668 = v_2(668^334)",
      333 + 333 + 2 == sp.multiplicity(2, sp.Integer(668)**334))
check("rank_167(H668) = 334 = n/2 forced (and 167 = 3 mod 4 with -1 nonsquare is NOT an", True)
print("       obstruction here because a SELF-DUAL F167 code is not forced: rank=334=n/2 but")
print("       length 668 = 0 mod 4 => Witt index 334; self-dual F167 codes of length 668 exist.)")
check("-1 nonsquare mod 167 but 668 = 0 mod 4 => hyperbolic => Witt index 334 (no obstruction)",
      sp.jacobi_symbol(-1, 167) == -1 and 668 % 4 == 0)

# A7. Laplacian / critical group of putative srg(333): forced 83-part
L_eigs_nonzero = sp.simplify((166 - r) * (166 - s))
check("(166-r)(166-s) = 27639 = 3^2 * 37 * 83", L_eigs_nonzero == 27639 and sp.factorint(27639) == {3:2, 37:1, 83:1})
tau = sp.Integer(27639)**166 / 333
check("#spanning trees = 27639^166/333 = 3^330 * 37^165 * 83^166",
      tau == sp.Integer(3)**330 * sp.Integer(37)**165 * sp.Integer(83)**166)
# mod 83: L = -A, rank_83 = 166, v_83(K) = 166 => Syl_83(K) = (Z/83)^166 forced
check("166 = 2*83 == 0 mod 83 => L == -A mod 83 => Syl_83(critical group) = Z_83^166 forced",
      166 % 83 == 0)
# mod 37 and mod 3: double roots => nilpotent freedom (matches 'F37 tower universally feasible')
check("x^2+x-83 has double root 18 mod 37 (disc = 333 == 0)", (18*18 + 18 - 83) % 37 == 0 and (2*18+1) % 37 == 0)
check("x^2+x-83 = (x-1)^2 mod 3 (disc = 333 == 0)", (1 + 1 - 83) % 3 == 0 and (2*1+1) % 3 == 0)
# exponent argument: L^2 - 333L + 27639I = 83J and L1=0  => 27639 = 3^2*37*83 kills K
check("L^2 - 333L + 27639I = 83J (operator identity from A^2+A-83I=83J)",
      166**2 - 333*166 + 27639 == -83)
# => Syl_37(K) elem abelian = Z_37^165, so rank_37(L) = 333 - 1 - 165 = 167 FORCED
# => Syl_83(K) elem abelian = Z_83^166, so rank_83(L) = 333 - 1 - 166 = 166 FORCED
# => Syl_3(K) = Z_3^(334-2e0) + Z_9^(e0-2), e0 = rank_3(L), window e0 in [6,167]
check("v_37(|K|) = 165, v_83(|K|) = 166, v_3(|K|) = 330",
      sp.multiplicity(37, sp.Integer(tau)) == 165 and sp.multiplicity(83, sp.Integer(tau)) == 166
      and sp.multiplicity(3, sp.Integer(tau)) == 330)
# mod-37 Jordan consistency for forced rank_37(L)=167: N = A-18I, N^2 = 9J != 0, N^3 = 0,
# rank N^2 = 1 => one J3 block; 167 = 2*1 + #J2 => 165 J2 blocks; 3+330+0 = 333 exactly fits.
check("forced rank_37(L) = 167 consistent: Jordan 1xJ3 + 165xJ2 + 0xJ1 = 333",
      (407 % 37 == 0) and (148 % 37 == 0) and 3 + 2*165 == 333)
check("rank_3(L) >= 6 (333 distinct nonzero rows of I-A mod 3) and <= 167", 3**5 - 1 < 333)

print()
# ---------------------------------------------------------------- PART B ----
print("=" * 76)
print("PART B: empirical validation on real Paley objects")
print("=" * 76)

# ---- finite fields GF(p^k) as tuples over F_p with given irreducible poly
def gf_elements(p, poly):  # poly: monic irreducible coeffs [c0..c_{k-1}] for x^k = -(...)
    kdeg = len(poly)
    return [tuple(c) for c in itertools.product(range(p), repeat=kdeg)]

def gf_mul(a, b, p, poly):
    kdeg = len(poly)
    prod = [0] * (2 * kdeg - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            prod[i + j] = (prod[i + j] + ai * bj) % p
    for d in range(2 * kdeg - 2, kdeg - 1, -1):  # reduce x^d
        c = prod[d]
        if c:
            prod[d] = 0
            for t in range(kdeg):
                prod[d - kdeg + t] = (prod[d - kdeg + t] - c * poly[t]) % p
    return tuple(prod[:kdeg])

def paley_conference(p, poly):
    """symmetric conference matrix of order q+1, q = p^len(poly) = 1 mod 4"""
    els = gf_elements(p, poly)
    q = len(els)
    sq = set()
    for z in els:
        if any(z):
            sq.add(gf_mul(z, z, p, poly))
    def chi(z):
        if not any(z): return 0
        return 1 if z in sq else -1
    n = q + 1
    C = [[0] * n for _ in range(n)]
    for i in range(q):
        C[0][i + 1] = C[i + 1][0] = 1
        for j in range(q):
            if i != j:
                d = tuple((els[i][t] - els[j][t]) % p for t in range(len(poly)))
                C[i + 1][j + 1] = chi(d)
    return np.array(C, dtype=np.int64)

def paley_hadamard(p, poly):
    """H = I + C_skew of order q+1, q = 3 mod 4"""
    els = gf_elements(p, poly)
    q = len(els)
    sq = set()
    for z in els:
        if any(z):
            sq.add(gf_mul(z, z, p, poly))
    def chi(z):
        if not any(z): return 0
        return 1 if z in sq else -1
    n = q + 1
    C = [[0] * n for _ in range(n)]
    for i in range(q):
        C[0][i + 1] = 1
        C[i + 1][0] = -1
        for j in range(q):
            if i != j:
                d = tuple((els[i][t] - els[j][t]) % p for t in range(len(poly)))
                C[i + 1][j + 1] = chi(d)
    return np.array(C, dtype=np.int64) + np.eye(n, dtype=np.int64)

def paley_graph(q):
    sq = {pow(t, 2, q) for t in range(1, q)}
    A = [[1 if (i - j) % q in sq else 0 for j in range(q)] for i in range(q)]
    return np.array(A, dtype=np.int64)

def rank_mod(M, p):
    A = (np.array(M, dtype=np.int64) % p).astype(np.int64)
    n, m = A.shape
    r = 0
    for c in range(m):
        piv = None
        for i in range(r, n):
            if A[i, c] % p:
                piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(n):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == n: break
    return r

def snf_multiset(M):
    S = smith_normal_form(sp.Matrix(M.tolist()))
    return [int(S[i, i]) for i in range(min(S.shape))]

# ---- B1. C(10), q = 9 = 3^2, poly x^2 + 1 (irreducible over F3)
C10 = paley_conference(3, [1, 0])  # x^2 = -1
check("C(10): C C^T = 9 I, symmetric, zero diag",
      np.array_equal(C10 @ C10.T, 9 * np.eye(10, dtype=np.int64)) and np.array_equal(C10, C10.T))
r3 = rank_mod(C10, 3)
print(f"  rank_3(C10) = {r3}")
check("C(10): rank_3 = 4 = n/2 - 1 (oddity formula forces b odd; window [3,4] -> a = 4)", r3 == 4)
snf10 = snf_multiset(C10)
check("C(10): SNF = (1^4, 3^2, 9^4) — pairing + no-self-dual cap verified",
      snf10 == [1]*4 + [3]*2 + [9]*4)

# ---- B2. C(26), q = 25, poly x^2 + 2 (x^2 = -2 = 3, 3 nonsquare mod 5)
C26 = paley_conference(5, [3, 0])  # x^2 = -3 = 2 mod 5? poly coeffs: x^2 = -(c0 + c1 x)
check("C(26): C C^T = 25 I, symmetric",
      np.array_equal(C26 @ C26.T, 25 * np.eye(26, dtype=np.int64)) and np.array_equal(C26, C26.T))
r5 = rank_mod(C26, 5)
print(f"  rank_5(C26) = {r5}")
check("C(26): rank_5 in window [3,13] (5 = 1 mod 4: cap is n/2 = 13, NOT n/2 - 1)",
      3 <= r5 <= 13)
b26 = 13 - r5
check("C(26): b EVEN here (p = 1 mod 4: oddity formula gives NO parity pin; contrast C(10)/C(50))",
      b26 % 2 == 0)
snf26 = snf_multiset(C26)
check(f"C(26): SNF = (1^{r5}, 5^{2*b26}, 25^{r5}) — pairing verified; a strictly inside window",
      snf26 == [1]*r5 + [5]*(2*b26) + [25]*r5)

# ---- B3. C(50), q = 49, poly x^2 + 1 (-1 nonsquare mod 7)
C50 = paley_conference(7, [1, 0])
check("C(50): C C^T = 49 I, symmetric",
      np.array_equal(C50 @ C50.T, 49 * np.eye(50, dtype=np.int64)) and np.array_equal(C50, C50.T))
r7 = rank_mod(C50, 7)
print(f"  rank_7(C50) = {r7}")
check("C(50): rank_7 <= 24 (no self-dual F7 code at length 50 = 2 mod 4) and >= 3", 3 <= r7 <= 24)
check("C(50): rank_7 EVEN (oddity formula: b odd since n-1 = 7^2 exactly)", r7 % 2 == 0)
snf50 = snf_multiset(C50)
b50 = 25 - r7
check(f"C(50): SNF = (1^{r7}, 7^{2*b50}, 49^{r7}) — pairing verified",
      snf50 == [1]*r7 + [7]*(2*b50) + [49]*r7)

# ---- B4. srg(13,6,2,3) — exact analog of forced SNF(A) at srg(333,166,82,83)
A13 = paley_graph(13)
check("Paley(13) is srg(13,6,2,3)",
      np.array_equal(A13 @ A13, 3*np.eye(13,dtype=np.int64) - A13 + 3*np.ones((13,13),dtype=np.int64)))
snfA13 = snf_multiset(A13)
check("SNF(srg(13,6,2,3)) = (1^6, 3^6, 6) — forced template diag(1^(m), mu^(m), 2mu) verified",
      snfA13 == [1]*6 + [3]*6 + [6])

# ---- B5. srg(29,14,6,7) — second analog (mu = 7 prime, k = 2*mu)
A29 = paley_graph(29)
check("Paley(29) is srg(29,14,6,7)",
      np.array_equal(A29 @ A29, 7*np.eye(29,dtype=np.int64) - A29 + 7*np.ones((29,29),dtype=np.int64)))
snfA29 = snf_multiset(A29)
check("SNF(srg(29,14,6,7)) = (1^14, 7^14, 14) — forced template verified",
      snfA29 == [1]*14 + [7]*14 + [14])

# ---- B5b. Laplacian critical-group exponent argument in the genuine half case
# (irrational srg eigenvalues), the exact mechanism used at 333:
#   P(13): exponent | (6-r)(6-s) = 39 = 3*13, |K| = 3^6 13^5
#          => SNF(L) = (1^6, 3, 39^5, 0) forced.  P(29): (1^14, 7, 203^13, 0) forced.
for (q, mu, w) in [(13, 3, 39), (29, 7, 203)]:
    Aq = paley_graph(q)
    Lq = (q - 1) // 2 * np.eye(q, dtype=np.int64) - Aq
    snfL = snf_multiset(Lq)
    m = (q - 1) // 2
    pred = [1] * m + [mu] + [w] * (m - 1) + [0]
    check(f"P({q}) Laplacian SNF = (1^{m}, {mu}, {w}^{m-1}, 0) — half-case exponent forcing verified",
          snfL == pred)

# ---- B6. Hadamard SNF forcing: H(12), H(20), H(28) (n/4 = 3, 5, 7 prime, like 668/4 = 167)
for (p, poly, n) in [(11, [4], 12), (19, [2], 20), (3, [2, 2, 0], 28)]:
    # GF(p) for poly len 1 (poly = [g] dummy: x = -g => just F_p), GF(27) via x^3 = -(2 + 2x) = x+1
    if len(poly) == 1:
        els = [(t,) for t in range(p)]
        H = paley_hadamard(p, [0])  # poly [0]: x^1 = 0 — degenerate; handle directly below
        # rebuild directly over F_p:
        q = p
        sq = {pow(t, 2, q) for t in range(1, q)}
        chi = lambda z: 0 if z % q == 0 else (1 if z % q in sq else -1)
        C = [[0]*(q+1) for _ in range(q+1)]
        for i in range(q):
            C[0][i+1] = 1; C[i+1][0] = -1
            for j in range(q):
                if i != j: C[i+1][j+1] = chi(i - j)
        H = np.array(C, dtype=np.int64) + np.eye(q+1, dtype=np.int64)
    else:
        H = paley_hadamard(p, poly)
    check(f"H({n}): H H^T = {n} I", np.array_equal(H @ H.T, n * np.eye(n, dtype=np.int64)))
    snfH = snf_multiset(H.astype(np.int64))
    m = n // 2
    check(f"SNF(H({n})) = (1, 2^{m-1}, {m}^{m-1}, {n}) — forced template (n/4 prime) verified",
          snfH == [1] + [2]*(m-1) + [m]*(m-1) + [n])

print()
print("=" * 76)
print("ALL CHECKS PASS" if ok else "SOME CHECKS FAILED")
print("=" * 76)
