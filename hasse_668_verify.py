#!/usr/bin/env python3
"""
Hasse-Minkowski / quadratic-form / lattice-genus obstruction check for H(668).

Question: H H^T = 668 I_668 with H in {+-1}^{668x668}. The quadratic-form functor
sees only: "the form 668*I_668 is represented by / equivalent to I_668 over Q (resp. Z)".
We verify, at 668, with explicit local computations:

  (A) I_668 ~ 668*I_668 over Q (Hasse-Minkowski: disc + Hasse invariant at every p + signature)
  (B) the p=167 subtlety: (668,668)_167 = -1 (167 == 3 mod 4), but it enters with EVEN
      exponent C(668,2), so it cancels;  same at p=2.
  (C) BRC for the associated 2-(667,333,166) design: x^2 = 167 y^2 - 166 z^2 has (1,1,1).
  (D) an EXPLICIT integral witness: M = I_167 (x) A, A a 4x4 Lipschitz-quaternion matrix
      from 668 = 25^2+5^2+3^2+3^2, M M^T = 668*I_668, M integral.
      => the lattice statement "exists sublattice L of Z^668 isometric to sqrt(668)Z^668,
         with 668 Z^668 < L < Z^668, index 668^334" is TRUE constructively.
      => genus theory cannot obstruct: any lattice with Gram 668*I in an orthogonal basis
         is literally isometric to sqrt(668)Z^668 (one isometry class, one genus class).
All arithmetic is exact (int). Hilbert symbol implementation is self-tested against
the product formula and bimultiplicativity on random inputs.
"""
import random
from math import gcd

# ---------- Hilbert symbol (a,b)_p over Q_p, p prime or p=-1 for R ----------
def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    t = pow(a, (p - 1) // 2, p)
    return 1 if t == 1 else -1

def val_unit(a, p):
    v = 0
    while a % p == 0:
        a //= p; v += 1
    return v, a

def hilbert(a, b, p):
    """Hilbert symbol (a,b)_p, a,b nonzero ints; p prime, or p == -1 for the real place."""
    assert a != 0 and b != 0
    if p == -1:
        return -1 if (a < 0 and b < 0) else 1
    if p == 2:
        al, u = val_unit(a, 2); be, v = val_unit(b, 2)
        eps = lambda x: ((x - 1) // 2) % 2
        om  = lambda x: ((x * x - 1) // 8) % 2
        e = eps(u) * eps(v) + al * om(v) + be * om(u)
        return -1 if e % 2 else 1
    al, u = val_unit(a, p); be, v = val_unit(b, p)
    e = (al * be * ((p - 1) // 2)) % 2
    s = (-1) ** e
    if be % 2:
        s *= legendre(u, p)
    if al % 2:
        s *= legendre(v, p)
    return s

def prime_factors(n):
    n = abs(n); fs = set(); d = 2
    while d * d <= n:
        while n % d == 0:
            fs.add(d); n //= d
        d += 1
    if n > 1:
        fs.add(n)
    return fs

# ---------- self-tests ----------
def selftest():
    rng = random.Random(668)
    # product formula: prod over all places of (a,b)_p = 1
    for _ in range(4000):
        a = rng.randint(-400, 400); b = rng.randint(-400, 400)
        if a == 0 or b == 0:
            continue
        places = {2} | prime_factors(a) | prime_factors(b)
        prod = hilbert(a, b, -1)
        for p in places:
            prod *= hilbert(a, b, p)
        assert prod == 1, (a, b)
    # bimultiplicativity and symmetry, random small primes
    primes = [2, 3, 5, 7, 11, 13, 37, 167]
    for _ in range(4000):
        a = rng.randint(1, 300); b = rng.randint(1, 300); c = rng.randint(1, 300)
        p = rng.choice(primes)
        assert hilbert(a, b, p) == hilbert(b, a, p)
        assert hilbert(a * b, c, p) == hilbert(a, c, p) * hilbert(b, c, p)
    # known values
    assert hilbert(-1, -1, 2) == -1
    assert hilbert(-1, -1, -1) == -1
    assert hilbert(2, 3, 3) == legendre(2, 3) == -1
    assert hilbert(5, 5, 5) == legendre(-1, 5) == 1
    print("[selftest] Hilbert symbol: product formula (4000 random), bimultiplicativity"
          " (4000 random), known values -- ALL PASS")

# ---------- (A)+(B): Hasse-Minkowski for I_668 vs 668*I_668 ----------
def hasse_check(n):
    print(f"\n=== (A) Rational equivalence I_{n} ~ {n}*I_{n} (Hasse-Minkowski) ===")
    # discriminant of n*I_n is n^n mod squares
    disc_sq = (n % 2 == 0) or False  # n^n square iff n even (exponent n even)
    exp_disc = n  # n^n
    print(f"disc(n*I_n) = {n}^{n}; exponent {n} even => square => disc matches disc(I_n)=1: "
          f"{exp_disc % 2 == 0}")
    # Hasse invariant c_p(n*I_n) = prod_{i<j} (n,n)_p = (n,n)_p^{C(n,2)}
    C2 = n * (n - 1) // 2
    print(f"C({n},2) = {C2}, parity: {'EVEN' if C2 % 2 == 0 else 'ODD'}")
    relevant = sorted(prime_factors(n) | {2})
    ok = True
    for p in relevant + [-1]:
        s = hilbert(n, n, p)
        cp = s ** (C2 % 2) if C2 % 2 else 1
        tag = "R" if p == -1 else f"p={p}"
        print(f"  ({n},{n})_{tag} = {s:+d}   -> c_{tag}(n*I_n) = ({s:+d})^{C2} = {cp:+d}"
              f"   (c of I_n is +1) {'MATCH' if cp == 1 else 'MISMATCH'}")
        ok &= (cp == 1)
    # all other p: (n,n)_p = 1 automatically (units, odd p not dividing n)
    print(f"  all other p: ({n},{n})_p = +1 (both units, p odd) -> c_p = +1, MATCH")
    # signature: both positive definite rank n
    print(f"  signature: both forms positive definite of rank {n}: MATCH")
    verdict = ok and (exp_disc % 2 == 0)
    print(f"==> Hasse-Minkowski: I_{n} and {n}*I_{n} are Q-equivalent: {verdict}")
    return verdict

# the p=167 subtlety, explicitly
def p167_subtlety():
    print("\n=== (B) The p=167 'subtlety' (167 == 3 mod 4), traced ===")
    print(f"167 mod 4 = {167 % 4}")
    print(f"(-1|167) = {legendre(-1, 167)}  (== -1 since 167 == 3 mod 4)")
    s2, s167 = hilbert(668, 668, 2), hilbert(668, 668, 167)
    print(f"(668,668)_2   = {s2:+d}   [= (167,-1)_2,  167 == 3 mod 4]")
    print(f"(668,668)_167 = {s167:+d}   [= (-1|167) = -1: NONTRIVIAL local symbol]")
    C2 = 668 * 667 // 2
    print(f"BUT exponent C(668,2) = {C2} = 334*667 is EVEN because 668 == 0 mod 4,")
    print(f"so c_2 = c_167 = (-1)^{C2} = +1. The nontrivial symbols cancel in pairs.")
    print("CONTRAST: for n == 2 mod 4 the exponent C(n,2) is ODD and the condition")
    print("(n,n)_p = +1 for all p <=> n = sum of two squares becomes a real gate")
    print("(this is exactly the conference-matrix/BRC gate; at 334: 333 = 18^2+3^2 passes).")
    C2_334 = 334 * 333 // 2
    print(f"  check: C(334,2) = {C2_334} is {'ODD' if C2_334 % 2 else 'EVEN'}; "
          f"(333,333)_3 = {hilbert(333,333,3):+d}, (333,333)_37 = {hilbert(333,333,37):+d} "
          f"(both +1 since 333 = 18^2+3^2)")
    return s167 == -1 and s2 == -1 and C2 % 2 == 0

# ---------- (C) BRC at the Hadamard 2-design ----------
def brc_check():
    print("\n=== (C) Bruck-Ryser-Chowla at the 2-(667,333,166) design ===")
    v, k, lam = 667, 333, 166
    nn = k - lam
    print(f"H(668) <=> symmetric 2-({v},{k},{lam}) design; order n = k - lambda = {nn} = 167")
    print(f"v = {v} odd, (v-1)/2 = {(v-1)//2} (odd) => BRC equation: x^2 = {nn} y^2 - {lam} z^2")
    x, y, z = 1, 1, 1
    lhs, rhs = x * x, nn * y * y - lam * z * z
    print(f"  witness (x,y,z)=(1,1,1): {lhs} = {nn}*1 - {lam}*1 = {rhs}: {lhs == rhs}")
    print("==> BRC passes IDENTICALLY for every Hadamard design (n - lambda = m - (m-1) = 1).")
    return lhs == rhs

# ---------- (D) explicit integral witness: quaternion blocks ----------
def quaternion_witness():
    print("\n=== (D) Explicit INTEGRAL witness M with M M^T = 668*I_668 ===")
    # four squares for 668
    sol = None
    for a in range(26, 0, -1):
        for b in range(a, -1, -1):
            r = 668 - a * a - b * b
            if r < 0:
                continue
            for c in range(b, -1, -1):
                d2 = r - c * c
                if d2 < 0:
                    continue
                d = int(d2 ** 0.5)
                for dd in (d - 1, d, d + 1):
                    if 0 <= dd <= c and dd * dd == d2:
                        sol = (a, b, c, dd)
                        break
                if sol:
                    break
            if sol:
                break
        if sol:
            break
    a, b, c, d = sol
    assert a*a + b*b + c*c + d*d == 668
    print(f"Lagrange: 668 = {a}^2 + {b}^2 + {c}^2 + {d}^2")
    A = [[ a,  b,  c,  d],
         [-b,  a, -d,  c],
         [-c,  d,  a, -b],
         [-d, -c,  b,  a]]
    # verify A A^T = 668 I_4 exactly
    AAT = [[sum(A[i][k] * A[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    ok4 = all(AAT[i][j] == (668 if i == j else 0) for i in range(4) for j in range(4))
    print(f"4x4 quaternion matrix A (rows of the left-regular rep of q=a+bi+cj+dk):")
    for row in A:
        print("   ", row)
    print(f"A A^T = 668*I_4 verified exactly: {ok4}")
    # M = I_167 (x) A : 668x668 block diagonal; verify M M^T = 668 I_668 exactly
    n = 668
    ok_full = True
    # block-diagonal product check without materializing 668x668 floats: exact, by blocks
    # (i-block, j-block) of M M^T is A A^T if i==j else 0 -- structural, but verify by
    # building the full integer matrix anyway (cheap at this size).
    M = [[0] * n for _ in range(n)]
    for blk in range(167):
        for i in range(4):
            for j in range(4):
                M[4 * blk + i][4 * blk + j] = A[i][j]
    # exact full verification of M M^T (integer arithmetic)
    for i in range(n):
        Mi = M[i]
        bi = i // 4
        for j in range(4 * bi, 4 * bi + 4):  # only same block can be nonzero
            s = sum(Mi[k] * M[j][k] for k in range(4 * bi, 4 * bi + 4))
            if s != (668 if i == j else 0):
                ok_full = False
        # off-block rows: supports disjoint => inner product 0 (verified structurally)
    # spot-verify 2000 random off-block pairs fully
    rng = random.Random(37)
    for _ in range(2000):
        i = rng.randrange(n); j = rng.randrange(n)
        if i // 4 == j // 4:
            continue
        s = sum(M[i][k] * M[j][k] for k in range(n))
        if s != 0:
            ok_full = False
    print(f"M = I_167 (x) A  (668x668, integer entries in {{0,+-{d},+-{c},+-{b},+-{a}}}):")
    print(f"M M^T = 668*I_668 verified exactly (all same-block entries + 2000 random "
          f"off-block pairs): {ok_full}")
    # lattice consequences
    print("\nLattice consequences (genus theory):")
    print(f"  L = M Z^668 is a sublattice of Z^668 with Gram 668*I_668,")
    print(f"  i.e. L is ISOMETRIC to sqrt(668)*Z^668 (orthogonal basis, equal norms 668;")
    print(f"  a lattice with Gram n*I in some basis IS sqrt(n)Z^n -- one isometry class,")
    print(f"  hence one genus class: genus theory has nothing to separate).")
    print(f"  det(Gram) = 668^668, so [Z^668 : L] = 668^334 (= sqrt(668^668)).")
    print(f"  668*Z^668 = M M^T Z^668 subset M Z^668 = L: chain 668 Z^668 < L < Z^668 holds.")
    print(f"  This is EXACTLY the lattice-level statement a Hadamard H(668) would produce;")
    print(f"  it is TRUE with witness M. What M lacks is +-1 ENTRIES -- a property the")
    print(f"  quadratic-form/lattice functor does not see (M's entries are 0,+-3,+-5,+-25).")
    return ok4 and ok_full

if __name__ == "__main__":
    selftest()
    r_a = hasse_check(668)
    r_b = p167_subtlety()
    r_c = brc_check()
    r_d = quaternion_witness()
    print("\n" + "=" * 70)
    print(f"(A) Q-equivalence I_668 ~ 668 I_668 (Hasse-Minkowski): {'PASSES' in ('PASSES',) and r_a}")
    print(f"(B) p=167 local symbol nontrivial but cancels (even exponent): {r_b}")
    print(f"(C) BRC at 2-(667,333,166) passes with (1,1,1): {r_c}")
    print(f"(D) integral/lattice witness M M^T = 668 I_668 exists: {r_d}")
    print("VERDICT: the quaternion/Hasse machinery PROVABLY CANNOT FIRE at 668.")
