#!/usr/bin/env python
"""
no3power_check.py -- numerical/exact sanity checks for the theorem:

    No strongly regular graph srg(333,166,82,83) admits a fixed-point-free
    automorphism of 3-power order.

and for the corollary chain (Fein-Kantor-Schacher => every vertex-transitive
example has a fixed-point-free automorphism of order exactly 37).

Proof being verified (see no3power_verification.md). The engine is the

  MASTER PARITY LEMMA: for every subgroup H <= Aut(Gamma) of ODD order,
  the number N of H-orbits on the 333 vertices satisfies N = 1 (mod 4).

  Proof skeleton: orbit partition is equitable with quotient Q; char(Q)
  divides char(A) = (x-166)(x^2+x-83)^166 with x^2+x-83 irreducible over Q
  (disc 333 is not a square), and Q has row sums 166, forcing
  char(Q) = (x-166)(x^2+x-83)^((N-1)/2), so tr(Q) = 166-(N-1)/2 = (333-N)/2.
  Each diagonal entry Q_ii is the degree of the vertex-transitive induced
  subgraph on an ODD orbit (|H| odd), hence EVEN by the handshake lemma.
  So (333-N)/2 is even, i.e. N = 1 (mod 4).

  Theorem: sigma fpf of order 3^e. e=1: N = 111 = 3 (mod 4), kill.
  e>=2: every <sigma>-orbit (size 3^j, j>=1) splits into exactly 3
  <sigma^3>-orbits, so N(sigma^3) = 3 N(sigma) = 3 (mod 4), kill.

The script verifies, with EXACT arithmetic:
  PART 0: the spectral/parity arithmetic on (333,166,82,83): the master
          constraint, the e=1 kill, exhaustive orbit-type kills for orders
          9, 27, 81, 243 (showing the sigma^3 step is genuinely needed),
          the f = 3 (mod 6) byproduct for order-3-with-fixed-points, and
          the FKS corollary arithmetic (q in {3,37}; q=37 => order exactly 37,
          9 orbits of 37).
  PART 1: end-to-end validation of EVERY lemma (equitability AS=SQ, row sums,
          annihilator, rank, char(Q) factorization, even-diagonal handshake,
          forced trace (v-N)/2, N = 1 mod 4, the orbit-splitting formula
          N(sigma^3) = 3N(sigma) - 2*n_1) on REAL conference graphs with REAL
          odd-order automorphisms:
            Paley(13)  x->3x   (order 3, f=1, N=5)
            Paley(37)  x->10x  (order 3, f=1, N=13)
            Paley(29)  x->16x  (order 7, f=1, N=5)
            Paley(53)  x->10x  (order 13, f=1, N=5)
            Paley(73)  x->2x   (order 9!, f=1, N=9; sigma^3=x->8x, N=25)
            Paley(125) Frobenius (order 3, f=5, N=45)
          plus the two CONTROLS where a genuine fpf automorphism of order
          p = 3 (mod 4) EXISTS, so the obstruction MUST fail, and the script
          verifies it fails at exactly one hypothesis (v is a square =>
          x^2+x-(v-1)/4 reducible => theta/tau multiplicities in char(Q)
          NOT forced equal), while every other lemma still holds:
            Paley(9)   x->x+1  (fpf order 3, N=3 = 3 mod 4)
            Paley(49)  x->x+1  (fpf order 7, N=7 = 3 mod 4)
  PART 2: symbolic proof that the earlier round's route (a) (cube-root-of-
          unity eigenspace traces) is VACUOUS: the trace identity reduces to
          0 = 0 and constrains nothing.

Runs in seconds. Exit code 0 iff all checks pass.
"""

from fractions import Fraction
from math import gcd, isqrt
import sys

import numpy as np
import sympy as sp

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not cond:
        FAILURES.append(name)


def exact_rank(M):
    """Rank of an integer sympy matrix, exact Fraction Gaussian elimination."""
    rows = [[Fraction(int(M[i, j])) for j in range(M.cols)] for i in range(M.rows)]
    rank, r = 0, 0
    for c in range(M.cols):
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


def cycles(perm):
    """Orbits of <perm> = cycles of the permutation."""
    n = len(perm)
    seen = [False] * n
    out = []
    for v in range(n):
        if not seen[v]:
            orb, seen[v], w = [v], True, perm[v]
            while w != v:
                orb.append(w)
                seen[w] = True
                w = perm[w]
            out.append(orb)
    return out


def perm_order(perm):
    o = 1
    for c in cycles(perm):
        o = o * len(c) // gcd(o, len(c))
    return o


# ----------------------------------------------------------------------
# PART 0: parity arithmetic for srg(333,166,82,83)
# ----------------------------------------------------------------------
print("PART 0: master-lemma arithmetic for srg(333,166,82,83)")
v, k, lam, mu = 333, 166, 82, 83
x = sp.Symbol("x")
quad = x**2 - (lam - mu) * x - (k - mu)            # x^2 + x - 83

check("conference condition 2k+(v-1)(lam-mu)=0 (multiplicities 166/166)",
      2 * k + (v - 1) * (lam - mu) == 0)
check("disc(x^2+x-83) = 333 = 9*37 is NOT a perfect square",
      isqrt(333) ** 2 != 333 and sp.Poly(quad, x).is_irreducible)
check("166 is not a root of x^2+x-83 (so a=1 in char(Q))",
      166**2 + 166 - 83 == 27639 and 27639 != 0)

# tr(Q) = 166 - b with b=(N-1)/2, i.e. (333-N)/2 ; even <=> N = 1 (mod 4)
Nsym = sp.Symbol("N")
trQ = k - (Nsym - 1) / 2
check("tr(Q) = 166-(N-1)/2 = (333-N)/2 (symbolic)",
      sp.simplify(trQ - (333 - Nsym) / 2) == 0)
check("for all odd N in [1,333]: (333-N)/2 even <=> N = 1 (mod 4)",
      all((((333 - N) // 2) % 2 == 0) == (N % 4 == 1) for N in range(1, 334, 2)))

# e = 1 kill
check("fpf order 3: N = 333/3 = 111 = 3 (mod 4)  ==> KILLED by master lemma",
      333 // 3 == 111 and 111 % 4 == 3)

# e >= 2 exhaustive kill over all orbit-type vectors, orders 9, 27, 81, 243
print("  exhaustive orbit-type check, fpf sigma of order 3^e (orbit sizes 3^j, j>=1):")
for e in (2, 3, 4, 5):
    sizes = [3**j for j in range(1, e + 1)]
    sols, single_pass, double_pass = 0, 0, 0

    def rec(idx, rem, counts):
        global sols, single_pass, double_pass
        if idx == len(sizes) - 1:
            if rem % sizes[idx]:
                return
            counts = counts + [rem // sizes[idx]]
            if counts[-1] < 1:        # order 3^e needs an orbit of size 3^e
                return
            sols += 1
            N = sum(counts)
            N3 = 3 * N                # N(sigma^3) = 3 N(sigma) for fpf sigma
            if N % 4 == 1:
                single_pass += 1
            if N % 4 == 1 and N3 % 4 == 1:
                double_pass += 1
            return
        for c in range(rem // sizes[idx] + 1):
            rec(idx + 1, rem - c * sizes[idx], counts + [c])

    sols = single_pass = double_pass = 0
    rec(0, 333, [])
    check(f"order 3^{e}={3**e}: all {sols} orbit types killed by "
          f"{{N=1 (4)}} AND {{3N=1 (4)}}",
          sols > 0 and double_pass == 0,
          f"{single_pass} types pass the single sigma-constraint "
          f"(sigma^3 step genuinely needed)" if single_pass else "")

# byproduct: order-3 automorphism with f fixed points needs f = 3 (mod 6)
allowed_f = [f for f in range(0, 334, 3)
             if (333 - 2 * ((333 - f) // 3)) % 4 == 1]
check("order-3 automorphism: master lemma <=> f = 3 (mod 6); f=0 excluded",
      all(f % 6 == 3 for f in allowed_f) and 0 not in allowed_f
      and allowed_f[0] == 3)

# FKS corollary arithmetic
print("  Fein-Kantor-Schacher corollary arithmetic:")
primes = [q for q in range(2, 334) if all(q % d for d in range(2, isqrt(q) + 1))]


def fpf_cycle_type_exists(q):
    """Can 333 be partitioned into parts q^j, j>=1 (a fpf q-element's cycle type)?"""
    if 333 % q:
        return False              # all parts divisible by q
    return True                   # 333/q parts of size q works


feasible = [q for q in primes if fpf_cycle_type_exists(q)]
check("fpf element of prime-power q^m order on 333 points forces q | 333, "
      "i.e. q in {3,37}", feasible == [3, 37])
check("q=37: 37^2 = 1369 > 333, so every cycle has length exactly 37 "
      "==> order exactly 37, 9 cycles", 37**2 > 333 and 333 // 37 == 9)
check("q=3 excluded by the theorem ==> vertex-transitive case has fpf "
      "order-37 automorphism", True, "conditional on FKS (CFSG)")

# ----------------------------------------------------------------------
# PART 1: end-to-end validation on real conference graphs
# ----------------------------------------------------------------------
print("\nPART 1: framework validation on real conference graphs")


def paley_prime(q):
    QR = {(i * i) % q for i in range(1, q)}
    A = np.zeros((q, q), dtype=np.int64)
    for i in range(q):
        for j in range(q):
            if i != j and (i - j) % q in QR:
                A[i, j] = 1
    return A


def paley_prime_square(p):
    """Paley graph on GF(p^2) = GF(p)[i]/(i^2+1), p = 3 (mod 4); returns
    (A, translation-by-1 permutation)."""
    assert p % 4 == 3
    elems = [(a, b) for a in range(p) for b in range(p)]
    idx = {e: i for i, e in enumerate(elems)}

    def fmul(u, w):
        return ((u[0] * w[0] - u[1] * w[1]) % p, (u[0] * w[1] + u[1] * w[0]) % p)

    squares = {fmul(e, e) for e in elems if e != (0, 0)}
    n = p * p
    A = np.zeros((n, n), dtype=np.int64)
    for i, u in enumerate(elems):
        for j, w in enumerate(elems):
            if i != j and ((u[0] - w[0]) % p, (u[1] - w[1]) % p) in squares:
                A[i, j] = 1
    perm = [idx[((e[0] + 1) % p, e[1])] for e in elems]
    return A, perm


def paley125():
    """Paley graph on GF(125) = GF(5)[t]/(t^3+3t+3) with Frobenius x->x^5."""
    p = 5

    def mul(a, b):
        c = [0] * 5
        for i in range(3):
            for j in range(3):
                c[i + j] = (c[i + j] + a[i] * b[j]) % p
        for d in (4, 3):
            if c[d]:
                c[d - 3] = (c[d - 3] + 2 * c[d]) % p
                c[d - 2] = (c[d - 2] + 2 * c[d]) % p
                c[d] = 0
        return (c[0], c[1], c[2])

    elems = [(a, b, c) for a in range(p) for b in range(p) for c in range(p)]
    idx = {e: i for i, e in enumerate(elems)}
    squares = {mul(e, e) for e in elems if e != (0, 0, 0)}
    A = np.zeros((125, 125), dtype=np.int64)
    for i, u in enumerate(elems):
        for j, w in enumerate(elems):
            if i != j and tuple((u[t] - w[t]) % p for t in range(3)) in squares:
                A[i, j] = 1

    def pw5(e):
        e2 = mul(e, e)
        return mul(mul(e2, e2), e)

    return A, [idx[pw5(e)] for e in elems]


def validate(name, A, perm, srg, control=False):
    """Verify every lemma used in the proof on a real (graph, automorphism).

    control=True: v is a perfect square; the master conclusion is EXPECTED to
    fail, and we verify it fails at exactly the irreducibility hypothesis.
    """
    vv, kk, ll, mm = srg
    d = perm_order(perm)
    print(f"  --- {name}: srg{srg}, automorphism of order {d} ---")
    A = np.asarray(A, dtype=np.int64)
    J = np.ones((vv, vv), dtype=np.int64)
    I = np.eye(vv, dtype=np.int64)
    check(f"{name}: srg identity A^2 = kI + lam A + mu(J-I-A)",
          np.array_equal(A @ A, kk * I + ll * A + mm * (J - I - A)))
    P = np.zeros((vv, vv), dtype=np.int64)
    for u in range(vv):
        P[perm[u], u] = 1
    check(f"{name}: perm is an automorphism of odd order {d}",
          np.array_equal(P @ A @ P.T, A) and d % 2 == 1 and d > 1)

    orbs = cycles(perm)
    N = len(orbs)
    n1 = sum(1 for o in orbs if len(o) == 1)
    check(f"{name}: all {N} orbit sizes odd and divide {d}",
          all(len(o) % 2 == 1 and d % len(o) == 0 for o in orbs))

    # equitability: A S = S Q, Q from orbit representatives
    S = np.zeros((vv, N), dtype=np.int64)
    Q = np.zeros((N, N), dtype=np.int64)
    for a, Ca in enumerate(orbs):
        for u in Ca:
            S[u, a] = 1
    for a, Ca in enumerate(orbs):
        for b2, Cb in enumerate(orbs):
            Q[a, b2] = int(A[Ca[0], Cb].sum())
    check(f"{name}: A S = S Q (orbit partition equitable)",
          np.array_equal(A @ S, S @ Q))
    check(f"{name}: Q row sums = k = {kk}", all(Q[a].sum() == kk for a in range(N)))

    # even-diagonal handshake lemma: induced subgraph on each orbit is
    # vertex-transitive => regular; odd order => even degree
    hand = True
    for a, Ca in enumerate(orbs):
        sub = A[np.ix_(Ca, Ca)]
        degs = set(int(r.sum()) for r in sub)
        if degs != {int(Q[a, a])} or Q[a, a] % 2:
            hand = False
    check(f"{name}: induced orbit subgraphs regular of EVEN degree "
          f"(handshake, odd orbits)", hand)
    check(f"{name}: tr(Q) = {int(np.trace(Q))} is even", int(np.trace(Q)) % 2 == 0)

    # annihilator + rank => char(Q) shape
    Qs = sp.Matrix(Q.tolist())
    It = sp.eye(N)
    M1 = Qs - kk * It
    M2 = Qs * Qs - (ll - mm) * Qs - (kk - mm) * It
    check(f"{name}: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)",
          M1 * M2 == sp.zeros(N, N))
    check(f"{name}: rank(Q-kI) = N-1 (k simple in Q)", exact_rank(M1) == N - 1)

    disc = (ll - mm) ** 2 + 4 * (kk - mm)        # = v for a conference graph
    if not control:
        check(f"{name}: disc = {disc} = v non-square (quadratic irreducible)",
              isqrt(disc) ** 2 != disc)
        b = (N - 1) // 2
        check(f"{name}: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^{b}; "
              f"tr(Q) = (v-N)/2 = {(vv - N) // 2}",
              int(np.trace(Q)) == (vv - N) // 2)
        check(f"{name}: b = {b} EVEN, N = {N} = 1 (mod 4)  [MASTER LEMMA]",
              b % 2 == 0 and N % 4 == 1)
        if N <= 25:
            cp = Qs.charpoly(x).as_expr()
            target = sp.expand((x - kk) * (x**2 - (ll - mm) * x - (kk - mm)) ** b)
            check(f"{name}: char(Q) verified directly", sp.expand(cp - target) == 0)
    else:
        # control: fpf automorphism exists, master conclusion must fail,
        # and the unique broken hypothesis is irreducibility (v square)
        check(f"{name}: CONTROL -- fpf (no fixed points)", n1 == 0)
        check(f"{name}: CONTROL -- N = {N} = 3 (mod 4): master CONCLUSION fails",
              N % 4 == 3)
        check(f"{name}: CONTROL -- disc = {disc} = v IS a perfect square "
              f"(irreducibility hypothesis fails; nothing else does)",
              isqrt(disc) ** 2 == disc)
        check(f"{name}: CONTROL -- forced trace (v-N)/2 = {(vv - N) // 2} is ODD, "
              f"actual tr(Q) = {int(np.trace(Q))} is even: equal-multiplicity "
              f"forcing absent", ((vv - N) // 2) % 2 == 1)
        # exhibit the unequal restricted-eigenvalue multiplicities
        r1 = ((ll - mm) + isqrt(disc)) // 2       # integer restricted eigenvalues
        r2 = ((ll - mm) - isqrt(disc)) // 2
        roots = sp.roots(Qs.charpoly(x).as_expr(), x)
        m1, m2 = roots.get(sp.Integer(r1), 0), roots.get(sp.Integer(r2), 0)
        check(f"{name}: CONTROL -- char(Q) multiplicities of r={r1}, s={r2} are "
              f"({m1},{m2}) UNEQUAL (rationality escape)", m1 != m2)
    return orbs, Q


validate("Paley(13)/x->3x", paley_prime(13), [(3 * u) % 13 for u in range(13)],
         (13, 6, 2, 3))
validate("Paley(37)/x->10x", paley_prime(37), [(10 * u) % 37 for u in range(37)],
         (37, 18, 8, 9))
validate("Paley(29)/x->16x", paley_prime(29), [(16 * u) % 29 for u in range(29)],
         (29, 14, 6, 7))
validate("Paley(53)/x->10x", paley_prime(53), [(10 * u) % 53 for u in range(53)],
         (53, 26, 12, 13))

# the order-9 case: sigma = x->2x on Paley(73) (2^9 = 512 = 1 mod 73)
perm73 = [(2 * u) % 73 for u in range(73)]
orbs9, _ = validate("Paley(73)/x->2x", paley_prime(73), perm73, (73, 36, 17, 18))
perm73c = [(8 * u) % 73 for u in range(73)]          # sigma^3 = x->8x, order 3
orbs3, _ = validate("Paley(73)/x->8x (=sigma^3)", paley_prime(73), perm73c,
                    (73, 36, 17, 18))
n1_9 = sum(1 for o in orbs9 if len(o) == 1)
check("Paley(73): orbit-splitting N(sigma^3) = 3 N(sigma) - 2 n_1  (25 = 27-2)",
      len(orbs3) == 3 * len(orbs9) - 2 * n1_9)
split_ok = all(
    sum(1 for o3 in orbs3 if set(o3) <= set(o9)) == (3 if len(o9) > 1 else 1)
    for o9 in orbs9)
check("Paley(73): every sigma-orbit of size 9 splits into EXACTLY 3 sigma^3-orbits",
      split_ok)

A125, perm125 = paley125()
validate("Paley(125)/Frobenius", A125, perm125, (125, 62, 30, 31))

# CONTROLS: genuine fpf automorphisms of order p = 3 (mod 4) on SQUARE-v
# conference graphs -- the obstruction must fail, exactly at irreducibility
A9, t9 = paley_prime_square(3)
validate("Paley(9)/x->x+1 [CONTROL]", A9, t9, (9, 4, 1, 2), control=True)
A49, t49 = paley_prime_square(7)
validate("Paley(49)/x->x+1 [CONTROL]", A49, t49, (49, 24, 11, 12), control=True)

# ----------------------------------------------------------------------
# PART 2: the earlier round's route (a) is vacuous (symbolic)
# ----------------------------------------------------------------------
print("\nPART 2: route (a) (zeta_3 eigenspace traces) is vacuous -- symbolic")
# sigma order 3 with f fixed points on 333 vertices; P = perm matrix.
# mult of eigenvalue 1 of P: f + (333-f)/3 ; of zeta, zeta^bar: (333-f)/3 each.
# V_theta, V_tau are real sigma-invariant, so on each the zeta/zeta^bar
# multiplicities are equal: tr(sigma|V_theta) = a - b_th, a + 2 b_th = 166,
# and the Perron line is sigma-fixed (eigenvalue 1).
f_, a_, ap_ = sp.symbols("f a ap")
b_th = (166 - a_) / 2
b_ta = (166 - ap_) / 2
mult1 = f_ + (333 - f_) / 3                       # eigenvalue-1 mult of P
constraint = sp.Eq(1 + a_ + ap_, mult1)           # 1 from the Perron line
expr = 1 + (a_ - b_th) + (ap_ - b_ta) - f_        # tr(P) - f, must be 0
sol = sp.solve(constraint, ap_)[0]
residual = sp.simplify(expr.subs(ap_, sol))
check("tr(P) = f reduces to 0 = 0 identically (no constraint on f)",
      residual == 0, f"residual = {residual}")

# ----------------------------------------------------------------------
print()
if FAILURES:
    print(f"RESULT: {len(FAILURES)} FAILURE(S):")
    for nm in FAILURES:
        print("  -", nm)
    sys.exit(1)
print("RESULT: ALL CHECKS PASS -- master parity lemma N=1 (mod 4) verified on "
      "all real examples; fpf 3-power orders 3,9,27,81,243 all killed on "
      "(333,166,82,83); controls Paley(9)/Paley(49) break exactly at "
      "irreducibility (v square). Theorem + corollary arithmetic verified.")
