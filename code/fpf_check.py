#!/usr/bin/env python
"""
fpf_check.py -- numerical/exact sanity checks for the theorem:

    Every automorphism of order 37 of a strongly regular graph with
    parameters (333,166,82,83) is fixed-point-free.

Proof being verified (see fpf_verification.md): orbit-partition quotient
matrix Q, char(Q) = (x-166)(x^2+x-83)^{166-18s}, the two trace identities
tr(Q), tr(Q^2), plus Cauchy-Schwarz on the orbit-orbit block B, give a
contradiction for every fixed-point count f = 333-37s with 1 <= s <= 8.

The script verifies, with EXACT arithmetic:
  PART 0: spectral arithmetic of srg(333,166,82,83);
  PART 1: every algebraic identity in the proof, symbolically in s;
  PART 2: the pincer infeasibility for s=1..8 exhaustively over integer G,
          and its (necessary!) feasibility at s=9, cross-checked against the
          actual orbit matrix R found in conf_orbit3.out;
  PART 3: end-to-end validation of the whole framework (all-or-nothing lemma,
          equitability, char(Q) factorization, both trace identities, the
          Cauchy-Schwarz direction) on REAL conference graphs with REAL
          odd-prime-order automorphisms with fixed points:
            Paley(13) with x->9x  (p=3, f=1, s=4)
            Paley(29) with x->16x (p=7, f=1, s=4)
            Paley(37) with x->x+1 (p=37, f=0, s=1; degenerate FPF case)
            Paley(125) with Frobenius x->x^5 (p=3, f=5, s=40; fixed subgraph C5)

Runs in seconds. Exit code 0 iff all checks pass.
"""

from fractions import Fraction
import sys

import numpy as np
import sympy as sp

FAILURES = []


def exact_rank(M):
    """Rank of an integer sympy/list matrix, exact Fraction Gaussian elimination."""
    rows = [[Fraction(int(M[i, j])) for j in range(M.cols)] for i in range(M.rows)]
    rank, ncols, r = 0, M.cols, 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pr = rows[r]
        for i in range(r + 1, len(rows)):
            if rows[i][c]:
                fac = rows[i][c] / pr[c]
                rows[i] = [a - fac * b for a, b in zip(rows[i], pr)]
        r += 1
        rank += 1
        if r == len(rows):
            break
    return rank


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not cond:
        FAILURES.append(name)


# ----------------------------------------------------------------------
# PART 0: spectral arithmetic for srg(333,166,82,83)
# ----------------------------------------------------------------------
print("PART 0: spectral facts for srg(333,166,82,83)")
n, k, lam, mu = 333, 166, 82, 83
x = sp.Symbol("x")
quad = x**2 - (lam - mu) * x - (k - mu)          # x^2 + x - 83
theta = sp.Rational(-1, 2) + 3 * sp.sqrt(37) / 2
tau = sp.Rational(-1, 2) - 3 * sp.sqrt(37) / 2

check("333 = 9*37", n == 9 * 37)
check("quad = x^2+x-83", sp.expand(quad - (x**2 + x - 83)) == 0)
check("theta root of quad", sp.simplify(quad.subs(x, theta)) == 0)
check("tau root of quad", sp.simplify(quad.subs(x, tau)) == 0)
check("theta+tau = -1, theta*tau = -83",
      sp.simplify(theta + tau + 1) == 0 and sp.simplify(theta * tau + 83) == 0)
check("theta^2+tau^2 = 167", sp.simplify(theta**2 + tau**2 - 167) == 0)
check("disc 333 = 9*37 not a perfect square (quad irreducible /Q)",
      sp.Poly(quad, x).is_irreducible and sp.sqrt(333).is_rational is not True)
# conference condition: 2k + (n-1)(lam-mu) = 0  ==> equal multiplicities (n-1)/2
check("conference condition 2k+(n-1)(lam-mu)=0", 2 * k + (n - 1) * (lam - mu) == 0)
check("multiplicities (n-1)/2 = 166", (n - 1) // 2 == 166)
check("166 mod 37 = 18 (fixed-degree congruence)", k % 37 == 18)
check("82 mod 37 = 8, 83 mod 37 = 9", lam % 37 == 8 and mu % 37 == 9)

# ----------------------------------------------------------------------
# PART 1: symbolic verification of every identity in the proof
# ----------------------------------------------------------------------
print("\nPART 1: symbolic identities (in the number of 37-orbits s)")
s_, G_, T_ = sp.symbols("s G T")
f_ = 333 - 37 * s_                       # fixed points     (Lemma 1)
b_ = (f_ + s_ - 1) / 2                   # pair multiplicity in char(Q)
check("b = (f+s-1)/2 = 166-18s", sp.simplify(b_ - (166 - 18 * s_)) == 0)

trQ = 166 + b_ * (theta + tau)           # = 166 - b
check("tr(Q) = 18s  ==> tr(B) = 18s (since tr(A_Delta)=0)",
      sp.simplify(trQ - 18 * s_) == 0)

trQ2 = 166**2 + b_ * (theta**2 + tau**2)
check("tr(Q^2) = 55278 - 3006 s", sp.simplify(trQ2 - (55278 - 3006 * s_)) == 0)

two_m = 166 * f_ - 37 * G_               # Lemma 2: sum of Delta-degrees
SB2_from_traces = sp.simplify(trQ2 - two_m - 74 * G_)
check("Sum B_ij^2 = tr(Q^2) - 2m - 74G = 3136 s - 37 G",
      sp.simplify(SB2_from_traces - (3136 * s_ - 37 * G_)) == 0)

# row sums: sum_ij B_ij = 166 s - G ; off-diagonal sum T = 166s - G - 18s
G_of_T = 148 * s_ - T_
SB2_in_T = sp.simplify(SB2_from_traces.subs(G_, G_of_T))
check("Sum B_ij^2 = 37 T - 2340 s", sp.simplify(SB2_in_T - (37 * T_ - 2340 * s_)) == 0)

# Cauchy-Schwarz pincer: 37T - 2340s >= 324 s + T^2/(s(s-1))  (s>=2)
# <=>  T^2 - 37 s(s-1) T + 2664 s^2 (s-1) <= 0
pincer_quad = sp.expand((324 * s_ + T_**2 / (s_ * (s_ - 1)) - SB2_in_T) * s_ * (s_ - 1))
target_quad = T_**2 - 37 * s_ * (s_ - 1) * T_ + 2664 * s_**2 * (s_ - 1)
check("pincer quadratic T^2 - 37s(s-1)T + 2664 s^2(s-1)",
      sp.simplify(pincer_quad - target_quad) == 0)
disc = sp.expand((37 * s_ * (s_ - 1))**2 - 4 * 2664 * s_**2 * (s_ - 1))
check("discriminant = s^2 (s-1) (1369(s-1) - 10656)",
      sp.simplify(disc - s_**2 * (s_ - 1) * (1369 * (s_ - 1) - 10656)) == 0)

print("  discriminant sign across the full range:")
for sv in range(1, 10):
    d = sv * sv * (sv - 1) * (1369 * (sv - 1) - 10656)
    tag = "  <-- infeasible (proof case)" if 1 <= sv <= 8 else "  <-- feasible (FPF case, as required)"
    if sv == 1:
        tag = "  <-- s=1 handled separately (T=0 forced)"
    print(f"    s={sv} (f={333-37*sv:3d}): disc = {d:>10d}{tag}")
    if 2 <= sv <= 8:
        check(f"disc<0 at s={sv}", d < 0)
check("disc>0 at s=9 (pincer correctly does NOT exclude FPF)",
      9 * 9 * 8 * (1369 * 8 - 10656) > 0)

# ----------------------------------------------------------------------
# PART 2: exhaustive integer infeasibility for s=1..8; feasibility at s=9
# ----------------------------------------------------------------------
print("\nPART 2: exhaustive pincer check over all integer G, s=1..9")
for sv in range(1, 9):
    f = 333 - 37 * sv
    if sv == 1:
        # B is 1x1: tr forces B11=18, row sum forces gamma1 = G = 148
        B11 = 18
        G = 148
        required = 3136 * sv - 37 * G          # = Sum B_ij^2 demanded by tr(Q^2)
        check(f"s=1 (f={f}): required Sum B^2 = {required} < 0 <= {B11**2} = B11^2",
              required < 0 and required != B11**2,
              f"identity forces SumB^2={required}, but SumB^2=324")
        continue
    # s >= 2: for every integer G, required SumB^2 must be >= CS lower bound; show it never is
    worst = None
    for G in range(0, 4 * f + 1):              # a_u <= 4 so G <= 4f (generous)
        T = 148 * sv - G
        required = Fraction(3136 * sv - 37 * G)
        lower = Fraction(324 * sv) + Fraction(T * T, sv * (sv - 1))
        margin = lower - required               # must be > 0 for ALL G (infeasibility)
        if worst is None or margin < worst[0]:
            worst = (margin, G)
    check(f"s={sv} (f={f}): infeasible for every integer G in [0,{4*f}]",
          worst[0] > 0, f"min (CS_lower - required) = {worst[0]} at G={worst[1]}")

# s = 9 (f = 0): G = 0 forced, T = 1332 forced; pincer must be satisfiable
sv, G = 9, 0
T = 148 * sv - G
required = 3136 * sv - 37 * G
lower = Fraction(324 * sv) + Fraction(T * T, sv * (sv - 1))
check("s=9 (f=0): required Sum B^2 = 28224 >= CS lower bound 27558 (margin 666)",
      Fraction(required) >= lower and required == 28224 and lower == 27558)

# cross-check against the actual orbit matrix found in conf_orbit3.out
R = sp.Matrix([
    [10, 19, 19, 19, 19, 20, 20, 20, 20],
    [19, 14, 22, 22, 20, 15, 15, 18, 21],
    [19, 22, 14, 22, 20, 15, 18, 21, 15],
    [19, 22, 22, 14, 20, 15, 21, 15, 18],
    [19, 20, 20, 20, 18, 24, 15, 15, 15],
    [20, 15, 15, 15, 24, 20, 19, 19, 19],
    [20, 15, 18, 21, 15, 19, 24, 17, 17],
    [20, 18, 21, 15, 15, 19, 17, 24, 17],
    [20, 21, 15, 18, 15, 19, 17, 17, 24]])
J9, I9 = sp.ones(9, 9), sp.eye(9)
check("R symmetric, row sums 166", R.T == R and all(sum(R.row(i)) == 166 for i in range(9)))
check("R^2 + R - 83 I = 3071 J (= 83*37 J)", R * R + R - 83 * I9 == 3071 * J9)
check("tr(R) = 162 = 18 s", R.trace() == 162)
check("Sum R_ij^2 = 28224 (matches tr(Q^2) identity at s=9, G=0)",
      sum(e * e for e in R) == 28224)
M1 = R - 166 * I9
M2 = R * R + R - 83 * I9
check("(R-166I)(R^2+R-83I) = 0 and rank(R-166I)=8  ==> char(R)=(x-166)(x^2+x-83)^4",
      M1 * M2 == sp.zeros(9, 9) and exact_rank(M1) == 8)

# ----------------------------------------------------------------------
# PART 3: end-to-end framework validation on real graphs
# ----------------------------------------------------------------------
print("\nPART 3: framework validation on real conference graphs")


def paley_prime(q):
    QR = {(i * i) % q for i in range(1, q)}
    A = np.zeros((q, q), dtype=np.int64)
    for i in range(q):
        for j in range(q):
            if i != j and (i - j) % q in QR:
                A[i, j] = 1
    return A


def paley125():
    """Paley graph on GF(125) = GF(5)[t]/(t^3+3t+3); vertices are field elements."""
    p = 5

    def mul(a, b):
        # multiply two cubics' residues: coefficient tuples (c0,c1,c2)
        c = [0] * 5
        for i in range(3):
            for j in range(3):
                c[i + j] = (c[i + j] + a[i] * b[j]) % p
        # reduce: t^3 = 2t + 2, t^4 = 2t^2 + 2t   (mod t^3+3t+3 over GF(5))
        for d in (4, 3):
            if c[d]:
                c[d - 3] = (c[d - 3] + 2 * c[d]) % p
                c[d - 2] = (c[d - 2] + 2 * c[d]) % p
                c[d] = 0
        return (c[0], c[1], c[2])

    elems = [(a, b, c) for a in range(p) for b in range(p) for c in range(p)]
    idx = {e: i for i, e in enumerate(elems)}
    squares = {mul(e, e) for e in elems if e != (0, 0, 0)}
    assert len(squares) == 62
    A = np.zeros((125, 125), dtype=np.int64)
    for i, u in enumerate(elems):
        for j, v in enumerate(elems):
            if i != j:
                d = tuple((u[t] - v[t]) % p for t in range(3))
                if d in squares:
                    A[i, j] = 1
    # Frobenius x -> x^5
    def pw5(e):
        e2 = mul(e, e)
        e4 = mul(e2, e2)
        return mul(e4, e)
    perm = [idx[pw5(e)] for e in elems]
    return A, perm


def validate(name, A, perm, p, srg):
    nn, kk, ll, mm = srg
    print(f"  --- {name}: srg{srg}, automorphism of order {p} ---")
    A = np.asarray(A, dtype=np.int64)
    J = np.ones((nn, nn), dtype=np.int64)
    I = np.eye(nn, dtype=np.int64)
    check(f"{name}: srg identity A^2 = kI + lam A + mu(J-I-A)",
          np.array_equal(A @ A, kk * I + ll * A + mm * (J - I - A)))
    P = np.zeros((nn, nn), dtype=np.int64)
    for v in range(nn):
        P[perm[v], v] = 1
    check(f"{name}: perm is an automorphism", np.array_equal(P @ A @ P.T, A))
    # order p
    Pi = np.linalg.matrix_power(P, p)
    check(f"{name}: perm has order {p}", np.array_equal(Pi, I) and not np.array_equal(P, I))

    # orbits
    seen, orbits, fixed = set(), [], []
    for v in range(nn):
        if v in seen:
            continue
        orb, w = [v], perm[v]
        while w != v:
            orb.append(w)
            w = perm[w]
        seen.update(orb)
        (fixed if len(orb) == 1 else orbits).append(orb)
    fixed = [o[0] for o in fixed]
    f, s = len(fixed), len(orbits)
    check(f"{name}: all non-fixed orbits have size {p}; f={f} = n mod-{p}-congruent",
          all(len(o) == p for o in orbits) and (nn - f) % p == 0)

    # all-or-nothing lemma
    aon = all(int(A[u, o].sum()) in (0, p) for u in fixed for o in [np.array(o) for o in orbits])
    check(f"{name}: all-or-nothing (fixed vertex adjacent to all or none of each orbit)", aon)

    # constant adjacency orbit->orbit and orbit->fixed
    constant = True
    for oi in orbits:
        for oj in orbits:
            vals = {int(A[w, oj].sum()) for w in oi}
            if len(vals) != 1:
                constant = False
        for u in fixed:
            vals = {int(A[w, u]) for w in oi}
            if len(vals) != 1:
                constant = False
    check(f"{name}: orbit partition is equitable (constant block counts)", constant)

    # quotient matrix Q (classes: fixed singletons then orbits)
    classes = [[u] for u in fixed] + orbits
    t = len(classes)
    Q = np.zeros((t, t), dtype=np.int64)
    for a, Ca in enumerate(classes):
        v0 = Ca[0]
        for bb, Cb in enumerate(classes):
            Q[a, bb] = int(A[v0, list(Cb)].sum())
    S = np.zeros((nn, t), dtype=np.int64)
    for a, Ca in enumerate(classes):
        for v in Ca:
            S[v, a] = 1
    check(f"{name}: A S = S Q (equitability, matrix form)", np.array_equal(A @ S, S @ Q))
    check(f"{name}: Q row sums = k", all(Q[a].sum() == kk for a in range(t)))

    # char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^b with b=(f+s-1)/2, via minimal-poly + rank
    Qs = sp.Matrix(Q.tolist())
    It = sp.eye(t)
    M1 = Qs - kk * It
    M2 = Qs * Qs - (ll - mm) * Qs - (kk - mm) * It
    check(f"{name}: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0", M1 * M2 == sp.zeros(t, t))
    check(f"{name}: rank(Q-kI) = f+s-1 (k simple in Q)", exact_rank(M1) == t - 1)
    check(f"{name}: f+s odd", (f + s) % 2 == 1)
    b = (f + s - 1) // 2
    disc_q = (ll - mm) ** 2 + 4 * (kk - mm)
    check(f"{name}: eigenvalue quadratic irreducible (disc {disc_q} non-square)",
          sp.sqrt(disc_q).is_rational is not True)
    # direct exact charpoly comparison when t is small
    if t <= 12:
        cp = Qs.charpoly(x).as_expr()
        target = sp.expand((x - kk) * (x**2 - (ll - mm) * x - (kk - mm)) ** b)
        check(f"{name}: char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^{b} (direct)",
              sp.expand(cp - target) == 0)

    # trace identities
    Adel = A[np.ix_(fixed, fixed)]
    two_m = int(Adel.sum())
    Nblock = Q[:f, f:]                      # = p * N (0/1 N)
    G = int((Nblock // p).sum())
    B = Q[f:, f:]
    SB2 = int((B * B).sum())
    trB = int(np.trace(B))
    trQ2 = int(np.trace(Q @ Q))
    check(f"{name}: tr(B) = k + b(lam-mu) = {kk + b*(ll-mm)}", trB == kk + b * (ll - mm))
    check(f"{name}: tr(Q^2) = k^2 + b[(lam-mu)^2+2(k-mu)] = {kk*kk + b*((ll-mm)**2+2*(kk-mm))}",
          trQ2 == kk * kk + b * ((ll - mm) ** 2 + 2 * (kk - mm)))
    check(f"{name}: tr(Q^2) = 2m + 2pG + Sum B^2  ({trQ2} = {two_m}+{2*p*G}+{SB2})",
          trQ2 == two_m + 2 * p * G + SB2)
    check(f"{name}: 2m = kf - pG  ({two_m} = {kk*f} - {p*G})", two_m == kk * f - p * G)
    check(f"{name}: B row sums k - gamma_i",
          all(int(B[i].sum()) == kk - int(Nblock[:, i].sum() // p) for i in range(s)))

    # Cauchy-Schwarz pincer direction (must HOLD on a real example -- no false kill)
    D = trB
    T = int(B.sum()) - D
    lower = Fraction(D * D, s) + (Fraction(T * T, s * (s - 1)) if s >= 2 else 0)
    check(f"{name}: pincer C-S inequality satisfied (SumB^2={SB2} >= {lower})",
          Fraction(SB2) >= lower)
    print(f"      f={f}, s={s}, b={b}, G={G}, 2m={two_m}, trB={D}, T={T}, "
          f"SumB^2={SB2}, CS lower={lower}")


validate("Paley(13)/x->9x", paley_prime(13),
         [(9 * v) % 13 for v in range(13)], 3, (13, 6, 2, 3))
validate("Paley(29)/x->16x", paley_prime(29),
         [(16 * v) % 29 for v in range(29)], 7, (29, 14, 6, 7))
validate("Paley(37)/x->x+1", paley_prime(37),
         [(v + 1) % 37 for v in range(37)], 37, (37, 18, 8, 9))
A125, perm125 = paley125()
validate("Paley(125)/Frobenius", A125, perm125, 3, (125, 62, 30, 31))

# ----------------------------------------------------------------------
print()
if FAILURES:
    print(f"RESULT: {len(FAILURES)} FAILURE(S):")
    for nm in FAILURES:
        print("  -", nm)
    sys.exit(1)
print("RESULT: ALL CHECKS PASS — the pincer kills f = 37c for all 1<=c<=8; "
      "only f=0 (and the identity, f=333) survive. Theorem verified.")
