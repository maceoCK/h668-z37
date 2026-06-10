#!/usr/bin/env python3
"""
Regular two-graph / equiangular-line GLOBAL machinery instantiated at order 668.

Layers:
  A. Regular two-graphs DIRECTLY on 668 points (and on 1336 points): feasibility.
  B. The meaningful instantiation: conference two-graph on 334 points
     <-> C(334) <-> srg(333,166,82,83) <-> 334 equiangular lines at 1/sqrt(333) in R^167.
     All classical global conditions, computed exactly.
  C. Equiangular-line bounds at the EXACT parameters (d=167, alpha=1/sqrt(333), N=334):
     relative bound (exact), Neumann threshold, absolute bound, full Delsarte LP over
     Gegenbauer cone (exact rational arithmetic), Perron-conjugate obstruction showing
     the Jiang-Tidor-Yao-Zhang-Zhao k(lambda) = infinity, bracketing Paley witnesses.
  D. Sufficiency demo in miniature: srg(13,6,2,3) -> C(14) -> H(28) verified, the exact
     map that would carry srg(333,166,82,83) -> C(334) -> H(668).
All arithmetic exact (Fraction / Q(sqrt(333))).
"""
from fractions import Fraction as F
from math import comb, isqrt
from itertools import product

print("=" * 78)
print("LAYER A: regular two-graphs DIRECTLY on 668 points")
print("=" * 78)

def two_graph_integer_feasible(n):
    """Nontrivial regular two-graph on n points, integer Seidel eigenvalues
    rho1 > 0 > rho2, rho1*rho2 = 1-n, multiplicities m1=n(-rho2)/(rho1-rho2),
    m2=n*rho1/(rho1-rho2) must be positive integers."""
    out = []
    m = n - 1
    for r1 in range(1, m + 1):
        if m % r1:
            continue
        r2 = -(m // r1)
        num1, den = n * (-r2), r1 - r2
        num2 = n * r1
        feas = (num1 % den == 0) and (num2 % den == 0)
        m1 = F(num1, den); m2 = F(num2, den)
        trivial = (r1 == 1) or (-r2 == 1)
        out.append((r1, r2, m1, m2, feas, trivial))
    return out

for n in (668, 1336):
    print(f"\nn = {n}:  rho1*rho2 = {1-n}")
    any_nontrivial = False
    for r1, r2, m1, m2, feas, triv in two_graph_integer_feasible(n):
        tag = "TRIVIAL (complete/void two-graph)" if triv else (
              "FEASIBLE" if feas else "INFEASIBLE (multiplicities not integral)")
        print(f"  (rho1,rho2)=({r1},{r2}):  m1={m1}, m2={m2}  -> {tag}")
        if feas and not triv:
            any_nontrivial = True
    # irrational case: rho = +-sqrt(n-1) <=> symmetric conference matrix of order n
    print(f"  irrational case rho=+-sqrt({n-1}): requires symmetric conference matrix "
          f"of order {n}; n mod 4 = {n%4} (symmetric C(n) forces n=2 mod 4) -> "
          + ("allowed" if n % 4 == 2 else "IMPOSSIBLE"))
    print(f"  ==> nontrivial regular two-graph on {n} points: "
          + ("possible a priori" if any_nontrivial else "INFEASIBLE — machinery vacuous at this order"))

# brute-force sanity demo of the n%4==0 conference obstruction at smallest case n=4
def conf4_exists():
    # symmetric, zero diag, +-1 off-diag, S^2 = 3I; 6 free upper-triangle entries
    for bits in product((-1, 1), repeat=6):
        S = [[0]*4 for _ in range(4)]
        k = 0
        for i in range(4):
            for j in range(i+1, 4):
                S[i][j] = S[j][i] = bits[k]; k += 1
        ok = True
        for i in range(4):
            for j in range(4):
                v = sum(S[i][t]*S[j][t] for t in range(4))
                if v != (3 if i == j else 0):
                    ok = False; break
            if not ok: break
        if ok:
            return True
    return False
print(f"\nBrute force: symmetric conference matrix of order 4 exists? {conf4_exists()} "
      "(classical: symmetric C(n) needs n=2 mod 4 — Belevitch 1950 / van Lint-Seidel 1966)")

print()
print("=" * 78)
print("LAYER B: conference two-graph on 334 points <-> C(334) <-> srg(333,166,82,83)")
print("=" * 78)

n = 334
print(f"\nn = 334 = 2 mod 4? {n % 4 == 2}   (needed for symmetric conference matrix)")

# van Lint–Seidel necessary condition: n-1 must be a sum of two squares
m = n - 1
reps = [(a, b) for a in range(isqrt(m)+1) for b in range(a, isqrt(m)+1) if a*a + b*b == m]
print(f"van Lint–Seidel: 333 = sum of two squares? representations: {reps}  -> "
      + ("PASS" if reps else "FAIL (would be a nonexistence proof)"))

# eigenvalue / multiplicity structure of the conference two-graph
print(f"Seidel eigenvalues +-sqrt(333); sqrt(333) irrational (18^2=324 < 333 < 361=19^2)")
print(f"multiplicities forced equal by trace 0: 167/167; n=334 even ✓")
a_coh = (n - 2) // 2
print(f"two-graph regularity: every pair in a = (n-2)/2 = {a_coh} coherent triples (integer ✓)")
tot = n*(n-1)*a_coh
print(f"total coherent triples = n(n-1)a/6 = {tot}/6 = {tot//6} (divisible ✓: {tot % 6 == 0})")

# ---- exact field Q(sqrt(333)) ----
class Q333:
    """a + b*sqrt(333), a,b rational"""
    __slots__ = ("a", "b")
    def __init__(s, a=0, b=0): s.a, s.b = F(a), F(b)
    def __add__(s, o): o = coerce(o); return Q333(s.a+o.a, s.b+o.b)
    __radd__ = __add__
    def __neg__(s): return Q333(-s.a, -s.b)
    def __sub__(s, o): return s + (-coerce(o))
    def __rsub__(s, o): return coerce(o) + (-s)
    def __mul__(s, o):
        o = coerce(o)
        return Q333(s.a*o.a + 333*s.b*o.b, s.a*o.b + s.b*o.a)
    __rmul__ = __mul__
    def __le__(s, o):  # exact comparison via sign of a + b*sqrt(333)
        d = s - coerce(o)
        if d.b == 0: return d.a <= 0
        if d.a == 0: return d.b <= 0
        # sign(a + b sqrt333): compare a^2 vs 333 b^2 with signs
        if d.a > 0 and d.b > 0: return False
        if d.a < 0 and d.b < 0: return True
        if d.a > 0:  # b<0: a + b s <=0 iff a <= -b s iff a^2 <= 333 b^2
            return d.a*d.a <= 333*d.b*d.b
        else:        # a<0, b>0: a + b s <= 0 iff 333 b^2 <= a^2
            return 333*d.b*d.b <= d.a*d.a
    def __float__(s): return float(s.a) + float(s.b)*333**0.5
    def __repr__(s): return f"({s.a} + {s.b}*sqrt333 ~ {float(s):.4f})"
def coerce(x): return x if isinstance(x, Q333) else Q333(x)

sqrt333 = Q333(0, 1)
v, k, lam, mu = 333, 166, 82, 83
print(f"\nsrg({v},{k},{lam},{mu}) (descendant of the two-graph):")
print(f"  identity k(k-lam-1) = (v-k-1)mu: {k*(k-lam-1)} = {(v-k-1)*mu} ✓")
r = (Q333(-1) + sqrt333) * F(1, 2)
s = (Q333(-1) - sqrt333) * F(1, 2)
print(f"  eigenvalues r = {r}, s = {s}  (= (-1+-3sqrt37)/2, conference graph)")
print(f"  multiplicities equal 166/166 (conference case; v=333=1 mod 4 ✓)")

# Krein conditions (exact in Q(sqrt333)):
K1_l = (r+1) * (Q333(k) + r + 2*r*s)
K1_r = (Q333(k) + r) * (s+1) * (s+1)
K2_l = (s+1) * (Q333(k) + s + 2*r*s)
K2_r = (Q333(k) + s) * (r+1) * (r+1)
print(f"  Krein 1: (r+1)(k+r+2rs) <= (k+r)(s+1)^2 :  {float(K1_l):.4f} <= {float(K1_r):.4f}"
      f"  -> {'PASS' if K1_l <= K1_r else 'FAIL'}   slack = {float(K1_r)-float(K1_l):.4f}")
print(f"  Krein 2: (s+1)(k+s+2rs) <= (k+s)(r+1)^2 :  {float(K2_l):.4f} <= {float(K2_r):.4f}"
      f"  -> {'PASS' if K2_l <= K2_r else 'FAIL'}   slack = {float(K2_r)-float(K2_l):.4f}")
f_mult = 166
print(f"  absolute bound: v <= f(f+3)/2 : {v} <= {f_mult*(f_mult+3)//2}  -> PASS"
      f"  (slack factor {f_mult*(f_mult+3)/2/v:.1f}x)")

print()
print("=" * 78)
print("LAYER C: 334 equiangular lines at angle 1/sqrt(333) in R^167 — all known bounds")
print("=" * 78)

d = 167
alpha2 = F(1, 333)   # alpha^2, exact

# 1. relative bound (Lemmens–Seidel 1973), exact
rel = F(d) * (1 - alpha2) / (1 - d*alpha2)
print(f"\n1. RELATIVE BOUND: N <= d(1-a^2)/(1-d a^2) = 167*(332/333)/(166/333) = {rel}")
print(f"   N = 334 attains it with EQUALITY (gap = {rel - 334}).")
print(f"   Equality <=> lines form a tight frame <=> regular two-graph on 334 pts <=> C(334).")

# 2. Neumann threshold (in Lemmens–Seidel 1973): N > 2d  =>  1/alpha odd integer
print(f"\n2. NEUMANN: N = 334 = 2d exactly; hypothesis 'N > 2d' FAILS -> theorem cannot fire.")
print(f"   Conversely a 335th line: 335 > 2*167 would force 1/alpha = sqrt(333) to be an")
print(f"   odd integer; isqrt(333) = {isqrt(333)}, {isqrt(333)}^2 = {isqrt(333)**2} != 333 -> N <= 334.")
print(f"   So Neumann INDEPENDENTLY caps N at 334: the putative system is exactly maximal.")
print(f"   Threshold 2d is permanently sharp for irrational angles: icosahedron = 6 = 2*3")
print(f"   lines at 1/sqrt(5) in R^3 (= Paley C(6) system) exists.")

# verify icosahedron meets relative bound too
rel_ico = F(3)*(1-F(1,5))/(1-3*F(1,5))
print(f"   (check: relative bound at (d=3, a=1/sqrt5) = {rel_ico} = 6, attained by C(6))")

# 3. absolute bound
absb = d*(d+1)//2
print(f"\n3. ABSOLUTE BOUND: N <= d(d+1)/2 = {absb}; 334 <= {absb}, slack factor {absb/334:.1f}x.")

# 4. Full Delsarte LP over the Gegenbauer cone, EXACT rational arithmetic.
#    Lines <-> antipodal codes: WLOG even polynomials; LP value = 1 + 1/max_k(-G_k(alpha))
#    over even k >= 2 (single constraint point => LP optimum at best single k).
#    G_k = C_k^lambda(alpha)/C_k^lambda(1), lambda = (d-2)/2 = 165/2, alpha = 1/sqrt(333).
lamb = F(165, 2)
# w_k = C_k(alpha) * (sqrt333 if k odd else 1)  -> rational for all k
w = [F(1), 2*lamb]            # w_0 = 1, w_1 = 2*lambda*alpha*sqrt333 = 2*lambda
Gvals = {}                     # even k -> G_k(alpha) exact
for kk in range(2, 121):
    if kk % 2 == 0:
        wk = (2*(kk + lamb - 1) * w[-1] / 333 - (kk + 2*lamb - 2) * w[-2]) / kk
    else:
        wk = (2*(kk + lamb - 1) * w[-1]       - (kk + 2*lamb - 2) * w[-2]) / kk
    w.append(wk)
    if kk % 2 == 0:
        Ck1 = comb(kk + 164, kk)     # C_k^lambda(1) = binom(k+2lambda-1, k), 2lambda-1=164
        Gvals[kk] = wk / Ck1

assert Gvals[2] == F(-1, 333), f"G_2 sanity failed: {Gvals[2]}"
print(f"\n4. DELSARTE LP (exact, even Gegenbauer degrees up to 120):")
print(f"   G_2(alpha) = {Gvals[2]} (exactly -1/333 ✓)")
ranked = sorted(Gvals.items(), key=lambda kv: kv[1])[:5]
for kk, g in ranked:
    if g < 0:
        print(f"   k={kk:3d}: G_k(alpha) = {float(g):+.3e}  -> single-k LP bound = 1 + 1/(-G_k) = {float(1 - 1/g):,.2f}")
best_k, best_g = ranked[0]
lp_val = 1 - 1/best_g
print(f"   LP OPTIMUM: attained at k = {best_k}, value = {lp_val} "
      f"({'EXACTLY 334' if lp_val == 334 else 'NOT 334'})")
runner = ranked[1]
print(f"   runner-up k = {runner[0]} gives {float(1 - 1/runner[1]):,.1f} — "
      f"gap factor {float((1 - 1/runner[1]) / 334):,.1f}x; LP can NEVER certify < 334 here")
print(f"   (and equality at k=2 is tight: complementary slackness with the putative")
print(f"    tight frame; any LP improvement would also kill existing Paley systems).")

# 5. Perron-conjugate obstruction: lambda = (1-alpha)/(2alpha) = (sqrt333 - 1)/2
#    minimal polynomial x^2 + x - 83; conjugates lam ~ 8.624 and lam' ~ -9.624
print(f"\n5. FIXED-ANGLE ASYMPTOTICS (Jiang–Tidor–Yao–Zhang–Zhao, Annals 2021):")
print(f"   lambda := (1-a)/(2a) = (sqrt(333)-1)/2, minimal poly x^2 + x - 83 = 0")
import math as _m
laml = (_m.sqrt(333)-1)/2; lamc = (-_m.sqrt(333)-1)/2
print(f"   conjugates: {laml:.6f} and {lamc:.6f}; |conjugate| - lambda = {abs(lamc)-laml:.6f} = 1 > 0")
print(f"   => lambda is NOT the largest conjugate in absolute value => NO finite graph has")
print(f"   adjacency spectral radius lambda => k(lambda) = INFINITY => JTYZZ regime is")
print(f"   N_a(d) = d + o(d) for FIXED a as d -> inf: non-effective at d = 167, and our")
print(f"   alpha = 1/sqrt(2d-1) varies with d (the 'conference diagonal'), outside scope.")

# 6. bracketing witnesses on the conference diagonal
def is_prime(p):
    if p < 2: return False
    for q in range(2, isqrt(p)+1):
        if p % q == 0: return False
    return True
print(f"\n6. BRACKETING PALEY WITNESSES on the diagonal alpha = 1/sqrt(2d-1):")
for q in (313, 337):
    print(f"   q = {q}: prime? {is_prime(q)}, q mod 4 = {q%4} -> Paley C({q+1}) EXISTS:")
    print(f"        {q+1} equiangular lines at 1/sqrt({q}) in R^{(q+1)//2} — KNOWN system")
print(f"   333 = 9*37 is not a prime power -> Paley does not construct C(334); but any")
print(f"   UNIFORM analytic bound that kills (334, R^167, 1/sqrt333) must separate it from")
print(f"   the existing (314, R^157, 1/sqrt313) and (338, R^169, 1/sqrt337) systems.")

print()
print("=" * 78)
print("LAYER D: sufficiency map verified in miniature: srg(13,6,2,3)->C(14)->H(28)")
print("=" * 78)
# Paley graph on 13
QR13 = {(x*x) % 13 for x in range(1, 13)}
A = [[1 if (i - j) % 13 in QR13 else 0 for j in range(13)] for i in range(13)]
S = [[0 if i == j else 1 - 2*A[i][j] for j in range(13)] for i in range(13)]  # Seidel
# conference matrix C(14) = [[0, 1^T],[1, S]]
C = [[0]*14 for _ in range(14)]
for i in range(13):
    C[0][i+1] = C[i+1][0] = 1
    for j in range(13):
        C[i+1][j+1] = S[i][j]
def matmulT(M):
    nn = len(M)
    return [[sum(M[i][t]*M[j][t] for t in range(nn)) for j in range(nn)] for i in range(nn)]
CCt = matmulT(C)
ok_conf = all(CCt[i][j] == (13 if i == j else 0) for i in range(14) for j in range(14))
print(f"C(14) from srg(13,6,2,3): C C^T = 13 I ? {ok_conf}")
# H(28) = [[C+I, C-I],[C-I, -C-I]]
H = [[0]*28 for _ in range(28)]
for i in range(14):
    for j in range(14):
        e = 1 if i == j else 0
        H[i][j]       = C[i][j] + e
        H[i][j+14]    = C[i][j] - e
        H[i+14][j]    = C[i][j] - e
        H[i+14][j+14] = -C[i][j] - e
HHt = matmulT(H)
ok_h = all(HHt[i][j] == (28 if i == j else 0) for i in range(28) for j in range(28))
print(f"H(28) by symmetric-conference doubling: H H^T = 28 I ? {ok_h}")
print("The IDENTICAL doubling maps C(334) -> H(668). Direction is one-way:")
print("C(334) => H(668); nonexistence of C(334) would NOT touch H(668).")
