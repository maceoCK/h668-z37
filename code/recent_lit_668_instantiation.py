#!/usr/bin/env python3
"""
Instantiation AT 668 of the 2018-2026 'recent literature' nonexistence machinery sweep.

1. Kotsireas-Koutschan (JCD 2021, arXiv:2101.03116) Algorithm 1 at ell=333:
   complete spectrum of the (ell/3)=111-th PSD value of a putative LP(333).
   If the spectrum were EMPTY -> no LP(333) -> no two-circulant-core H(668).
2. Barrera Acevedo - O Cathain - Dietrich (JCD 2019, arXiv:1904.11460) hypotheses at p=167:
   p mod 4, Williamson rowsum decompositions w^2+x^2+y^2+z^2=668 (odd), Paley arithmetic.
3. Field-descent self-conjugacy (Turyn/Schmidt, used by all recent group-ring tech):
   is 2 (resp 167) self-conjugate mod 333? mod 167?
4. BRC / sum-of-two-squares gates at 334/333; Krein + absolute + Delsarte clique
   bounds at srg(333,166,82,83) (the parameter conditions every recent SRG
   nonexistence paper must pass through).
"""
from math import isqrt, gcd

def factorint(n):
    f, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def n_order(a, n):
    assert gcd(a, n) == 1
    t, k = a % n, 1
    while t != 1:
        t, k = (t * a) % n, k + 1
    return k

ELL = 333
M = ELL // 3            # 111
TARGET = 2 * ELL + 2    # 668  = PSD_A(111) + PSD_B(111)
QSUM = (2 * TARGET + 2) // 3 if False else 4 * M + 2  # 446 = sum of the two odd-square sums

print("=" * 72)
print("1. KOTSIREAS-KOUTSCHAN mod-3 LP CONSTRAINT, Algorithm 1 at ell = 333")
print("=" * 72)
print(f"ell = {ELL} = 3*{M}; constraint: Ahat+Bhat = {TARGET}, Ahat == 4 (mod 12),")
print(f"A1+A2+A3 = 1 (all odd, |Ai| <= {M}), A1^2+A2^2+A3^2 = (2*Ahat+1)/3,")
print(f"and (A-sum-of-squares) + (B-sum-of-squares) = 4m+2 = {4*M+2}")

def odd_triples_with_sum_squares(Q, bound):
    """all multisets of odd integers (x1,x2,x3), |xi|<=bound, sum=1, sum of squares=Q"""
    sols = []
    # enumerate |x1|>=|x2|>=|x3| over odd absolute values, then sign patterns
    for a in range(1, bound + 1, 2):
        if 3 * a * a < Q:
            continue
        if a * a > Q - 2:
            break
        for b in range(1, a + 1, 2):
            c2 = Q - a * a - b * b
            if c2 < 1 or c2 > b * b:
                continue
            c = isqrt(c2)
            if c * c != c2 or c % 2 == 0:
                continue
            # sign patterns giving sum exactly +1  (KK normalize rowsum to +1)
            for sa in (a, -a):
                for sb in (b, -b):
                    for sc in (c, -c):
                        if sa + sb + sc == 1:
                            t = tuple(sorted((sa, sb, sc), reverse=True))
                            if t not in sols:
                                sols.append(t)
    return sols

surviving = []
total_pairs_checked = 0
for Ahat in range(4, TARGET + 1, 12):
    Bhat = TARGET - Ahat
    if Bhat < 0:
        break
    if Ahat > Bhat:
        break  # unordered pairs only
    assert Bhat % 12 == 4, "consistency: Bhat must also be == 4 mod 12"
    total_pairs_checked += 1
    QA = (2 * Ahat + 1)
    QB = (2 * Bhat + 1)
    if QA % 3 or QB % 3:
        continue
    QA //= 3
    QB //= 3
    assert QA + QB == 4 * M + 2
    solsA = odd_triples_with_sum_squares(QA, M)
    solsB = odd_triples_with_sum_squares(QB, M)
    if solsA and solsB:
        surviving.append((Ahat, Bhat, len(solsA), len(solsB), solsA[0], solsB[0]))

print(f"\ncandidate [Ahat,Bhat] pairs with Ahat<=Bhat, Ahat==4 mod 12 : {total_pairs_checked}")
print(f"SURVIVING pairs (both Diophantine systems solvable, all-odd, sum=1): {len(surviving)}")
for Ahat, Bhat, nA, nB, exA, exB in surviving:
    print(f"  [Ahat,Bhat]=[{Ahat:3d},{Bhat:3d}]  QA={(2*Ahat+1)//3:3d} ({nA} odd sols, e.g. {exA})"
          f"  QB={(2*Bhat+1)//3:3d} ({nB} odd sols, e.g. {exB})")

verdict1 = "NON-EMPTY -> constraint PASSES (prunes, does not obstruct)" if surviving else \
           "EMPTY -> LP(333) DOES NOT EXIST (obstruction fires!)"
print(f"\n>>> KK mod-3 spectrum at ell=333: {verdict1}")

print()
print("=" * 72)
print("2. COCYCLIC H(4p) THEOREM (Barrera Acevedo-O Cathain-Dietrich) at p=167")
print("=" * 72)
p = 167
print(f"p = {p}, p mod 4 = {p % 4}  -> Thm 4.5b applies: every cocyclic H(668) is")
print("   Williamson-type or transposed-Ito over back-circulant 167x167 blocks.")
print(f"4p-1 = {4*p-1} = {factorint(4*p-1)}  prime power? {len(factorint(4*p-1))==1}")
print(f"2p-1 = {2*p-1} = {factorint(2*p-1)}  prime power? {len(factorint(2*p-1))==1}")
print("   -> Paley I (core 667) and Paley II (core 333) constructions both UNAVAILABLE.")
# Williamson rowsum decompositions: 4p = w^2+x^2+y^2+z^2 over POSITIVE ODD integers
decomps = []
for w in range(1, isqrt(4*p) + 1, 2):
    for x in range(w, isqrt(4*p) + 1, 2):
        for y in range(x, isqrt(4*p) + 1, 2):
            z2 = 4*p - w*w - x*x - y*y
            if z2 < y*y:
                continue
            z = isqrt(z2)
            if z*z == z2 and z % 2 == 1:
                decomps.append((w, x, y, z))
print(f"odd rowsum decompositions w^2+x^2+y^2+z^2 = 668 (w<=x<=y<=z): {len(decomps)}")
print(f"   {decomps}")
print(">>> hypothesis layer non-empty -> theorem CONSTRAINS but does not kill cocyclic H(668);")
print(">>> Williamson(167) and Ito(167) remain open finite searches.")

print()
print("=" * 72)
print("3. FIELD DESCENT / SELF-CONJUGACY (Turyn; Schmidt 1999+; Leung-Schmidt updates)")
print("=" * 72)
def self_conjugate(q, n):
    """q self-conjugate mod n iff some power q^t == -1 mod n' (n' = n / q-part)"""
    n1 = n
    while n1 % q == 0:
        n1 //= q
    if n1 == 1:
        return True
    t = q % n1
    seen = set()
    while t not in seen:
        if t == n1 - 1:
            return True
        seen.add(t)
        t = (t * q) % n1
    return False

for (q, n) in [(2, 333), (167, 333), (2, 167), (3, 167), (37, 167)]:
    sc = self_conjugate(q, n)
    print(f"  {q} self-conjugate mod {n}?  {sc}")
print("  ord_9(2)=%d (2^3=8=-1 mod 9: t odd); ord_37(2)=%d, 2^18 mod 37 = %d (t even)"
      % (n_order(2, 9), n_order(2, 37), pow(2, 18, 37)))
print("  ord_9(167 mod 9 = 5)=%d, 5^3 mod 9=%d; ord_37(167 mod 37 = 19)=%d, 19^18 mod 37=%d"
      % (n_order(5, 9), pow(5, 3, 9), n_order(19, 37), pow(19, 18, 37)))
print("  ord_167(2) = %d (odd -> 2 not self-conjugate mod 167; confirms project fact)"
      % n_order(2, 167))
print(">>> parity clash (t odd mod 9 vs t even mod 37): NEITHER 2 NOR 167 is")
print(">>> self-conjugate mod 333 -> all self-conjugacy-based descent provably cannot fire.")

print()
print("=" * 72)
print("4. PARAMETER GATES AT srg(333,166,82,83) / C(334)  (what any SRG-tech must pass)")
print("=" * 72)
v, k, lam, mu = 333, 166, 82, 83
import math
disc = (lam - mu) ** 2 + 4 * (k - mu)
s_ = (lam - mu - math.sqrt(disc)) / 2
r_ = (lam - mu + math.sqrt(disc)) / 2
f = g = (v - 1) // 2
print(f"eigenvalues r,s = (-1 +- 3*sqrt(37))/2 = {r_:.4f}, {s_:.4f}; mult {f},{g} (irrational, conference)")
# BRC for conference matrix C(334): n-1 = 333 must be sum of two squares
reps = [(a, isqrt(333 - a * a)) for a in range(1, 19) if isqrt(333 - a*a)**2 == 333 - a*a]
print(f"BRC gate C(334): 333 sum of two squares? {reps}  -> PASSES")
# Krein conditions (pass iff (r+1)(k+r+2rs) <= (k+r)(s+1)^2 and the r<->s swap)
k1 = (1 + r_) * (k + r_ + 2 * r_ * s_) - (k + r_) * (s_ + 1) ** 2
k2 = (1 + s_) * (k + s_ + 2 * r_ * s_) - (k + s_) * (r_ + 1) ** 2
print(f"Krein 1: LHS-RHS = {k1:.2f} <= 0? {k1 <= 1e-9} (PASS);  "
      f"Krein 2: LHS-RHS = {k2:.2f} <= 0? {k2 <= 1e-9} (PASS)")
print(f"absolute bound: v={v} <= f(f+3)/2 = {f*(f+3)//2}? {v <= f*(f+3)//2}")
print(f"Delsarte clique bound: 1 + k/(-s) = {1 + k / (-s_):.3f} -> max clique <= {int(1 + k / (-s_))}")
print(f"vertex parity (Bondarenko-style local counting needs small lambda): lambda={lam} (dense local graph)")
print(">>> all classical+modern parameter gates PASS; Paley comparison: same-shape conference")
print(">>> parameters realized at every prime power q==1 mod 4, so no parameter-shape-only")
print(">>> method can distinguish 333; nonexistence would need 333-specific structure.")

print()
print("=" * 72)
print("5. ARITHMETIC HYPOTHESES OF THE REMAINING RECENT TECH")
print("=" * 72)
print(f"balancedly (multi-)splittable Hadamard (Kharaghani-Suda 2018/2022): negative results")
print(f"   live at orders n=4k^2 (36=4*3^2, 100=4*5^2, tied to PP(6),PP(10)); 668/4 = 167 square? "
      f"{isqrt(167)**2 == 167}")
print(f"Miller 2025 multiplier nonexistence needs p^2 | v: 667 = {factorint(667)} squarefree -> no")
print(f"Wang 2025 NNS(42),(44) nonexistence: Yang route needs composite 167? 167 prime -> N/A")
print(f"668 = 4*167: 167 == 3 mod 4 -> 668 != 4u^2, circulant-Hadamard tech (Leung-Schmidt) N/A")
print(f"ell=333 mod 5 = {333 % 5} -> Kotsireas et al. mod-5 paper (arXiv:2111.02105) N/A at 333")
