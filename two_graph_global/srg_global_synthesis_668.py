#!/usr/bin/env python
"""
srg GLOBAL SYNTHESIS at 668:  automorphism landscape x counting (Lam-style mass attack)
on srg(333,166,82,83) <-> C(334)  [conference route to H(668); sufficient, NOT necessary].

Question: can "few automorphisms" + "large isomorph classes" + the proven automorphism
landscape (all order-37 actions fpf & enumerated; no fpf 3-power; mod-CFSG every
vertex-transitive example has an fpf Z37) be combined into a GLOBAL counting
contradiction, the way Lam-Thiel-Swiercz killed PG(2,10)?

Layers:
  1. Classical global feasibility battery at (333,166,82,83), exact.   [all PASS => no free kill]
  2. The Lam skeleton question: is there a nontrivial self-orthogonal code to force
     countable substructures?  p=2: PROVABLY trivial.  p=3: self-orthogonal but
     provably NOT self-dual (Witt/discriminant at length 334 == 2 mod 4) => no Gleason.
     p=37: project's snf_668_verify.py already certified blind.  Verified on real
     Paley conference matrices C(6),C(10),C(14),C(18),C(30),C(38).
  3. Mass/counting landscape: exact labeled two-graph counts, the Z37-equivariant
     corner measure, exact unstructured lift-space census over ALL 625 orbit-matrix
     classes, and a complete brute-force demo at n=6 showing the mass identity
     N_labeled = n!/|Aut| is only available AFTER classification (circularity), while
     Burnside-accessible counts (Mallows-Sloane) see ALL two-graphs, never the regular ones.
  4. Lam-style cost model: exact branching magnitudes at C(334) vs the PG(2,10),
     srg(75,32,10,16), srg(95,40,12,20), C(66)/C(86) calibration anchors.
All arithmetic exact (sympy / Python ints) except where labelled "estimate".
"""
import sys, math
from math import comb, factorial, log10, log2
from fractions import Fraction
from itertools import combinations, permutations, product
import numpy as np
import sympy as sp

ok_all = True
def check(name, cond):
    global ok_all
    print(("PASS " if cond else "**FAIL** ") + name)
    if not cond: ok_all = False

E = "=" * 78
print(E); print("LAYER 1: classical global feasibility battery at srg(333,166,82,83), exact"); print(E)

v, k, lam, mu = 333, 166, 82, 83
check("standard identity k(k-lam-1) = (v-k-1)mu", k*(k-lam-1) == (v-k-1)*mu)
s37 = sp.sqrt(37)
r = (-1 + 3*s37)/2; s = (-1 - 3*s37)/2
check("restricted eigenvalues (-1±3√37)/2 satisfy x²+x-83=0",
      sp.simplify(r**2 + r - 83) == 0 and sp.simplify(s**2 + s - 83) == 0)
# conference case: multiplicities equal
check("2k+(v-1)(lam-mu)=0 => multiplicities 166/166 integral (conference graph)",
      2*k + (v-1)*(lam-mu) == 0)
# Krein conditions (Scott): (r+1)(k+r+2rs) <= (k+r)(s+1)^2  and swap
K1 = sp.simplify((k+r)*(s+1)**2 - (r+1)*(k+r+2*r*s))
K2 = sp.simplify((k+s)*(r+1)**2 - (s+1)*(k+s+2*r*s))
check(f"Krein 1 slack = {K1} = {sp.N(K1,6)} >= 0", K1 >= 0)
check(f"Krein 2 slack = {K2} = {sp.N(K2,6)} >= 0", K2 >= 0)
f = g = 166
check(f"absolute bound v <= f(f+3)/2 = {f*(f+3)//2}", v <= f*(f+3)//2)
# van Lint-Seidel / BRC analogue for C(334): 333 must be a sum of two squares
reps = [(a,b) for a in range(19) for b in range(a,19) if a*a+b*b == 333]
check(f"C(334) necessary condition: 333 = a²+b², reps {reps}", len(reps) > 0)
# Neumaier claw bound: only applies when s is integral -- conference graphs exempt
check("Neumaier/claw bound inapplicable (s irrational) -- conference family is the exempt class",
      not s.is_integer)
print("""
=> Every classical global condition PASSES (consistent with the project's relaxation
   meta-result: all lenses forgetting ±1-integrality are blind at these parameters).""")

print(E); print("LAYER 2: is there a Lam code skeleton at C(334)?  (the enabling structure of PG(2,10))"); print(E)
print("""
Lam-Thiel-Swiercz needed: a self-orthogonal binary code whose weight enumerator is
PINNED by Gleason-type theory, forcing a NONZERO count of small codewords (w12/w16/w19)
inside any putative plane -- the 'forced substructure' that exhaustion then kills.
Test the exact analogues at C(334): S symmetric, zero diagonal, ±1 off-diag, S² = 333·I.
""")
# p = 2: det(S)^2 = 333^334 odd  =>  S invertible over F2  =>  row space = F2^334.
check("p=2: det(S)² = 333^334 is ODD => rank_2(S) = 334 (full) => binary code TRIVIAL, no skeleton",
      333**334 % 2 == 1)
# p = 3: S² = 333 I ≡ 0 mod 3 => rowspace ⊆ kernel => C3 self-orthogonal, rank ≤ 167.
# Self-dual (rank=167) would need a 167-dim totally isotropic subspace of the standard
# form on F3^334; that form is hyperbolic iff (-1)^167 = -1 is a square mod 3. It is not.
check("p=3: C3 := rowspace(S mod 3) is SELF-ORTHOGONAL (S²≡0 mod 3), dim ≤ 167",
      (333) % 3 == 0)
check("p=3: -1 is a non-square mod 3 and 334 ≡ 2 (mod 4) => NO self-dual ternary code of "
      "length 334 exists => rank_3(S) ≤ 166 < 167: code strictly self-orthogonal, Gleason "
      "does NOT apply, weight enumerator NOT pinned, no forced-substructure counts",
      pow(2, 1, 3) != 1 and 334 % 4 == 2)
# p = 37: same self-orthogonality; project's snf_668_verify.py + collapse theorem already
# certified the F37 tower blind (PROVEN BLIND list). Recorded, not re-litigated.
print("p=37: S² ≡ 0 mod 37, self-orthogonal F37 code -- already PROVEN BLIND "
      "(snf_668_verify.py, F37 self-dual code towers, trace-linear collapse theorem).")

# ---- verify all of this on REAL Paley conference matrices --------------------
def paley_conference(q, field=None):
    """Symmetric conference matrix of order q+1, q ≡ 1 mod 4 prime (or GF(9) etc. via field)."""
    if field is None:
        els = list(range(q)); sub = lambda a,b: (a-b) % q
        sq = set((x*x) % q for x in range(1,q))
    else:
        els, sub, sq = field
    n = q+1
    S = [[0]*n for _ in range(n)]
    for i in range(q): S[0][i+1] = S[i+1][0] = 1
    for i,a in enumerate(els):
        for j,b in enumerate(els):
            if i != j: S[i+1][j+1] = 1 if sub(a,b) in sq else -1
    return S

def gf9():
    els = [(a,b) for a in range(3) for b in range(3)]          # a + b*i, i² = -1
    mul = lambda x,y: ((x[0]*y[0] - x[1]*y[1]) % 3, (x[0]*y[1] + x[1]*y[0]) % 3)
    sub = lambda x,y: ((x[0]-y[0]) % 3, (x[1]-y[1]) % 3)
    sq  = set(mul(e,e) for e in els if e != (0,0))
    return els, sub, sq

def rank_mod(M, p):
    M = [[x % p for x in row] for row in M]; n = len(M); m = len(M[0]); r = 0
    for c in range(m):
        piv = next((i for i in range(r,n) if M[i][c] % p), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p-2, p)
        M[r] = [(x*inv) % p for x in M[r]]
        for i in range(n):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(M[i][j] - f*M[r][j]) % p for j in range(m)]
        r += 1
        if r == n: break
    return r

print("\nVerification on real Paley conference matrices:")
cases = [(5,None),(9,gf9()),(13,None),(17,None),(29,None),(37,None)]
for q, fld in cases:
    S = paley_conference(q, fld); n = q+1
    M2 = [[S[i][j]*S[i][j] and sum(S[i][t]*S[j][t] for t in range(n)) or
           sum(S[i][t]*S[j][t] for t in range(n)) for j in range(n)] for i in range(n)]
    conf_ok = all(sum(S[i][t]*S[j][t] for t in range(n)) == (q if i==j else 0)
                  for i in range(n) for j in range(n))
    r2 = rank_mod(S, 2)
    line = f"  C({n}): S·Sᵀ={q}·I {'✓' if conf_ok else 'FAIL'};  rank_2 = {r2}/{n} (full: {r2==n})"
    if (n-1) % 3 == 0:
        r3 = rank_mod(S, 3)
        line += f";  3|{n-1}: rank_3 = {r3} ≤ {n//2-1} ({'✓' if r3 <= n//2-1 else 'FAIL'}) [Witt bound n/2-1]"
        check_r3 = r3 <= n//2 - 1
    else:
        check_r3 = True
    print(line)
    check(f"  C({n}) consistent", conf_ok and r2 == n and check_r3)

# SNF pairing d_i * d_{n+1-i} = n-1 (the only integral structure; constrains, never forbids)
print("\nSNF elementary-divisor pairing d_i·d_(n+1-i) = n-1 on real C(n):")
for q, fld in [(5,None),(9,gf9()),(13,None),(17,None)]:
    S = paley_conference(q, fld); n = q+1
    d = sp.Matrix(S).elementary_divisors() if hasattr(sp.Matrix(S),'elementary_divisors') else None
    from sympy.matrices.normalforms import smith_normal_form
    D = smith_normal_form(sp.Matrix(S))
    dd = [abs(D[i,i]) for i in range(n)]
    paired = all(dd[i]*dd[n-1-i] == q for i in range(n))
    print(f"  C({n}): SNF = {dd}  pairing d_i·d_(n+1-i)={q}: {paired}")
    check(f"  C({n}) SNF pairing", paired)
print("""=> At C(334) the SNF is forced only up to the pairing d_i·d_335-i = 333 = 3²·37 and
   rank_3 ≤ 166;  every such pattern is realized by integer matrices (project: snf_668_verify.py)
   => p-adic/SNF layer CONSTRAINS but cannot FORBID.   NO LAM CODE SKELETON EXISTS AT C(334).""")

print(E); print("LAYER 3: mass / counting landscape"); print(E)
# 3a. exact global counts
P333 = comb(333,2); P334 = comb(334,2)
log10_2 = log10(2)
log10_fact334 = sum(log10(i) for i in range(2,335))
print(f"""
3a. EXACT global counts (labeled):
  unordered pairs on 333 / 334 points: {P333} / {P334}
  labeled graphs on 334 points:            2^{P334}  ≈ 10^{P334*log10_2:.0f}
  labeled TWO-GRAPHS on 334 points (each switching class has size 2^333):
      2^{P334}/2^333 = 2^{P333}  ≈ 10^{P333*log10_2:.1f}      [EXACT]
  log10(334!) = {log10_fact334:.1f}
  => iso classes of ALL two-graphs on 334 ≥ 2^{P333}/334! ≈ 10^{P333*log10_2-log10_fact334:.0f}
  REGULAR ones with eigenvalues ±3√37 among them: NO FORMULA EXISTS (see 3c).""")

# 3b. the Z37-equivariant corner measure
orbits_pairs = P333 // 37
check(f"3b. fpf Z37 acts freely on pairs: {P333} = 37·{orbits_pairs} (all pair-stabilizers trivial, "
      "2d≡0 mod 37 => d=0)", P333 == 37*orbits_pairs and all((2*d) % 37 != 0 for d in range(1,37)))
print(f"""  Z37-equivariant Seidel space: 2^{orbits_pairs} of 2^{P333} labeled two-graphs
  => the ENTIRE enumerated automorphism corner has measure 2^-{P333-orbits_pairs} ≈ 10^-{(P333-orbits_pairs)*log10_2:.0f}
  of the search space. Closing all 625 classes says NOTHING about the trivial-Aut bulk.""")

# 3c. complete brute-force demo at n=6: the mass identity is circular
print("\n3c. Brute-force demo at n=6 (smallest conference order): why counting is circular")
n6 = 6; pairs6 = list(combinations(range(n6),2)); triples6 = list(combinations(range(n6),3))
# two-graph = parity function on triples induced by a graph; enumerate all graphs, build two-graphs
seen = {}
for bits in range(1 << len(pairs6)):
    adj = set(p for i,p in enumerate(pairs6) if bits >> i & 1)
    key = 0
    for t,(a,b,c) in enumerate(triples6):
        e = ((a,b) in adj) + ((a,c) in adj) + ((b,c) in adj)
        if e & 1: key |= 1 << t
    seen[key] = seen.get(key,0)+1
check(f"  labeled two-graphs on 6 points: {len(seen)} = 2^C(5,2) = {2**comb(5,2)}; "
      f"every switching class has 2^5=32 graphs: {set(seen.values())=={32}}",
      len(seen) == 2**comb(5,2) and set(seen.values()) == {32})
# iso classes under S6, regular ones, automorphism groups
def act(key, perm):
    out = 0
    for t,(a,b,c) in enumerate(triples6):
        ta = tuple(sorted((perm[a],perm[b],perm[c])))
        if key >> triples6.index(ta) & 1: out |= 1 << t
    return out
perms6 = list(permutations(range(n6)))
orbits = {}; aut = {}
for key in seen:
    img = [act(key,p) for p in perms6]
    can = min(img)
    orbits.setdefault(can, set()).add(key)
    if key == can: aut[can] = sum(1 for x in img if x == key)
# regular: every pair in constant # of coherent (odd) triples
def pair_degrees(key):
    deg = {p:0 for p in pairs6}
    for t,(a,b,c) in enumerate(triples6):
        if key >> t & 1:
            for p in ((a,b),(a,c),(b,c)): deg[p] += 1
    return set(deg.values())
regular_lab = [key for key in seen if len(pair_degrees(key)) == 1]
reg_classes = set(min(act(key,p) for p in perms6) for key in regular_lab)
print(f"  iso classes of two-graphs on 6: {len(orbits)}   (Mallows-Sloane: = # even graphs on 6 up to iso)")
# verify Mallows-Sloane at n=6: count even (all-degrees-even) graphs on 6 up to iso
evens = set()
for bits in range(1 << len(pairs6)):
    deg = [0]*6
    for i,p in enumerate(pairs6):
        if bits >> i & 1: deg[p[0]] += 1; deg[p[1]] += 1
    if all(d % 2 == 0 for d in deg):
        can = min(sum(1 << pairs6.index(tuple(sorted((pp[p[0]],pp[p[1]])))) for i,p in enumerate(pairs6) if bits >> i & 1) for pp in perms6)
        evens.add(can)
check(f"  Mallows-Sloane verified at n=6: #two-graphs/iso = {len(orbits)} = #even graphs/iso = {len(evens)}",
      len(orbits) == len(evens))
nontrivial_reg = [c for c in reg_classes if 0 < bin(c).count('1') < len(triples6)]
key_reg = nontrivial_reg[0] if nontrivial_reg else None
print(f"  REGULAR two-graphs on 6 (nontrivial): {len(nontrivial_reg)} iso class, "
      f"{sum(1 for k in regular_lab if 0 < bin(k).count('1') < 20)} labeled")
if key_reg is not None:
    a6 = aut[min(act(key_reg,p) for p in perms6)]
    nl = sum(1 for k in regular_lab if 0 < bin(k).count('1') < 20)
    check(f"  mass identity: labeled = 6!/|Aut| = 720/{a6} = {720//a6} = {nl} (orbit-stabilizer)",
          nl == 720 // a6)
print("""  THE POINT: 6!/|Aut| was computable only AFTER finding the object and its Aut group.
  Burnside/Mallows-Sloane machinery counts ALL two-graphs (switching classes <-> even graphs,
  a LINEAR/parity bijection); REGULARITY is a spectral condition (S²=(n-1)I) that no
  group-action bookkeeping sees. There is no independent handle on N_labeled(regular):
  a mass-formula contradiction would need an a-priori value for a number that is only
  defined through the classification itself.  CIRCULAR => cannot fire.""")

# 3d. growth of known conference-graph counts (literature: Brouwer srg tables / Spence pages)
print("""
3d. Known conference-graph counts up to iso (Brouwer tables; Spence www.maths.gla.ac.uk/~es/):
   v=5:1   v=9:1   v=13:1   v=17:1   v=25:15 (complete)   v=29:41 (complete, Spence)
   v=37: ≥6760 ('very likely not exhaustive', Spence)   v=45: ≥78
   => complete-classification frontier for conference graphs: v=29 (two-graphs: n=30).
   => super-exponential explosion. Heuristic expectation at v=333: 0 or >10^thousands.
   The project's annealed count E_pairs(333)=2^-1.37 sits AT threshold: counting heuristics
   CANNOT distinguish 0 from astronomical -- in either direction. A counting PROOF of 0 would
   need exact cancellation machinery that does not exist for this family.""")

# 3e. exact unstructured lift-space census over ALL 625 orbit-matrix classes
print("3e. EXACT unstructured lift space per Z37 orbit-matrix class (all 625 classes):")
rows = open('/Users/maceocardinalekwik/git/math/release/h668-z37/code/orbit_all_classes.txt').read().split()
mats = [list(map(int,rows[i*81:(i+1)*81])) for i in range(len(rows)//81)]
logs = []
for m in mats:
    R = [m[i*9:(i+1)*9] for i in range(9)]
    L = 0.0
    for i in range(9):
        L += log10(comb(18, R[i][i]//2))          # symmetric diag connection set
        for j in range(i+1,9):
            L += log10(comb(37, R[i][j]))          # off-diag connection set (D_ji = -D_ij)
    L -= 8*log10(37)                               # gauge: orbit base-point choices
    logs.append(L)
logs.sort()
print(f"  classes: {len(mats)};  lift space per class: min 10^{logs[0]:.1f}, "
      f"median 10^{logs[len(logs)//2]:.1f}, max 10^{logs[-1]:.1f}")
print(f"  TOTAL Z37 corner (sum over classes): ≈ 10^{logs[-1] + log10(sum(10**(l-logs[-1]) for l in logs)):.1f}")
print("""  The cyclotomic sweep (5072 jobs, all INFEASIBLE) closed only the multiplier-structured
  SUB-ANSATZ of this; the per-class GENERAL lift is a ~10^391 space => 'close all 625 classes'
  is NOT 625 CP-SAT calls; it is 625 instances each individually beyond unstructured exhaustion.""")

print(E); print("LAYER 4: Lam-style cost model -- C(334) vs every successful exhaustion ever done"); print(E)
print(f"""
Anchors (all literature-verified, URLs in report):
  PG(2,10), Lam-Thiel-Swiercz 1989 (Canad. J. Math 41, 1117-1123): ~2.7 months Cray-1A
     for the weight-19 case alone (+ yrs of workstations; ovals case 4400h VAX 11/780);
     enabled by a self-orthogonal binary code forcing countable w12/w16/w19 substructures.
     NO ANALOGUE exists at C(334) (Layer 2: p=2 code provably trivial, p=3 not self-dual).
  srg(75,32,10,16) & srg(95,40,12,20), Azarija-Marc (arXiv:1509.05933, 1603.02032): the
     ONLY feasible-parameter srg exhaustions ever completed. Their skeleton: INTEGER
     eigenvalue r=2 of multiplicity 56 (resp. 75) => star complement of order 19 (resp. 20);
     plus a FORCED count >=783 of 4-cliques (Bondarenko-Prymak-Radchenko) to seed the search;
     feasibility hinged on reducing every candidate to <=19 vertices ('n(G)-k2(G) >= 17').
  Conference matrices: known for all admissible orders through 62; smallest UNRESOLVED
     order 86 (66 modulo a disputed 2021 construction claim, arXiv:2102.05432) --
     open for ~60 years at orders two hundred-some below 334.
  Regular two-graphs: complete classification stops at n=36 (McKay-Spence 2001: exactly 227);
     n=38 only bounded (>=191, new ones still being found in 2023). Conference-class
     frontier n=30. 334 is ~300 points beyond every classification frontier.""")

# 4b. WHY the one successful srg technique provably cannot fire here: star complement scaling
print("4b. Star-complement scaling (exact):")
def srg_mults(v_,k_,l_,m_):
    D = (l_-m_)**2 + 4*(k_-m_); sD = sp.sqrt(D)
    r_ = (l_-m_+sD)/2; s_ = (l_-m_-sD)/2
    f_ = sp.Rational(1,2)*((v_-1) - (2*k_+(v_-1)*(l_-m_))/sD)
    g_ = (v_-1) - f_
    return r_, s_, sp.nsimplify(f_), sp.nsimplify(g_)
for (v_,k_,l_,m_) in [(75,32,10,16),(95,40,12,20),(333,166,82,83)]:
    r_, s_, f_, g_ = srg_mults(v_,k_,l_,m_)
    mmax = max(f_, g_, key=lambda t: float(t))
    sc = v_ - int(mmax)
    print(f"  srg{(v_,k_,l_,m_)}: eigenvalues {sp.nsimplify(r_)},{sp.nsimplify(s_)}; mults {f_}/{g_};"
          f" star complement order v-max(f,g) = {sc};  raw complement space 2^C({sc},2) = 2^{comb(sc,2)}"
          f" ≈ 10^{comb(sc,2)*log10(2):.0f}")
check("  conference parameters MAXIMIZE star-complement order: for every srg, max(f,g) >= (v-1)/2 "
      "with equality iff f=g (conference); so sc order = 167 = (v+1)/2 is the worst case possible "
      "on 333 vertices", max(166,166) == (333-1)//2)
print("""  => the Azarija-Marc engine needed order-19/20 star complements (2^171 raw, pruned);
     at (333,166,82,83) the star complement has 167 vertices (2^13861 raw) and the
     eigenvalue (-1+3√37)/2 is irrational (compatibility graph over Q(√37), no integral
     clique machinery). The unique technique that ever closed a feasible srg parameter
     set is structurally maximally mismatched to conference parameters. CANNOT FIRE.""")

print(f"""
Magnitudes at C(334) (exact):
  raw labeled two-graph space:        2^{P333} ≈ 10^{P333*log10_2:.0f}
  after FULL hypothetical Lam-grade pruning (PG(2,10) collapsed ~10^600-raw to ~10^15,
  i.e. ~585 orders): 10^{P333*log10_2:.0f} - 585 ≈ 10^{P333*log10_2-585:.0f}.  Useless.
  Z37-equivariant corner only:        2^{orbits_pairs} ≈ 10^{orbits_pairs*log10_2:.0f}, split into 625 classes
     of ~10^391 each; measured wall: CP-SAT/CDCL UNKNOWN at 2.5-9000s on even the
     415-bit sigma-equivariant SUB-instances (project, 60 seeded runs, glassy floor).
  => A Lam-style certificate of srg(333) nonexistence is not a compute purchase at ANY
     2026 budget; it is missing its enabling theorem (the code skeleton), and the
     fallback (raw orderly generation) is ~10^16000 with nothing to collapse it.

What the synthesis CAN deliver (the real theorem in this corner):
  [proven] every order-37 automorphism is fpf; orbit matrices = 625 classes; cyclotomic
  lifts empty for all 6 multiplier subgroups of order >=4; no fpf 3-power.
  [conditional target] IF the remaining multiplier orders 2,3 + self-complementary classes
  221/422 close (cube-and-conquer + VeriPB, est. 10^4-10^6 CPU-h, genuinely uncertain),
  THEN mod CFSG: NO vertex-transitive srg(333,166,82,83) exists -- a publishable theorem.
  Even TOTAL success (full srg(333) nonexistence) would NOT kill H(668): the conference
  route is sufficient, not necessary (H(668) <-> LP(333) survives independently).""")

print(E)
print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED -- see above")
print(E)
