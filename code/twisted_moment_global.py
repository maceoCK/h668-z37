"""
Twisted second moments -> GLOBAL system for srg(333,166,82,83) with fpf Z37 (H(668) route).

Mechanism under test ("twisted moment global"):
  For each nonzero character w of Z37 the 9x9 Hermitian block M_w (M_w^2 + M_w = 83 I)
  has r-eigenprojector P_w of rank a_w (a_w = 4 for w in QR, 5 for w in QNR; a=4 WLOG).
  Twisted overlaps  tau_{w,w'} = Tr(P_w P_{w'})  satisfy the PROVEN subspace box
      max(0, a_w + a_{w'} - 9) <= tau_{w,w'} <= min(a_w, a_{w'}).
  Galois-trace (Fourier) identity, derived and machine-validated below:
      sum_{w in Z37*} Tr(M_w M_{uw}) = 37*F_0(u) - S2,
  where F_0(u) = sum_{ij} |D_ij ∩ u D_ij| (near-multiplier count) and S2 = sum R_ij^2 = 28224
  (universal for every admissible orbit matrix).  Substituting the eigen-decomposition:
      X_u := sum_{w in Z37*} tau_{w,uw} = (F_0(u) - 36)/9          [MASTER IDENTITY]
  QR-restricted version (Gauss-sum split, also validated):
      Y_u := sum_{w in QR} tau_{w,uw} = (F_0(u) - 198)/18  (u in QR)
                                       = (F_0(u) -  36)/18  (u in QNR)
  Boxes on X_u give NEW necessary conditions:
      u in QR\{1}:  198 <= F_0(u) <= 1494
      u in QNR   :   36 <= F_0(u) <= 1332   (so a QNR multiplier is impossible, defect >= 162)
  plus  F_0(u) even,  F_0(u) = F_0(u^{-1}),  and (order-3 classes, H = {1,10,26})
  F_0 constant on uH with F_0 ≡ E (mod 3), E = #{(i,j): R_ij ≡ 1 mod 3} forced.
  Character-twisted global sums (Bochner cuts): for every character psi of Z37* ≅ Z36,
      sum_u psi(u) F_0(u) = sum_{ij} |Psi_ij(psi)|^2 >= 0.
  Exact chi (=quadratic) component:
      2*sum_{u in QR} F_0(u) = sum_ij (R_ij - e_ij)^2 + sum_ij x_ij^2 + 36E,
      sum_{u in Z37*} F_0(u) = sum_ij (R_ij - e_ij)^2 + 36E,
  with e_ij = [0 in D_ij] (forced by R mod 3 on order-3 classes) and
  x_ij = sum_{g in D_ij} chi(g) (integer; ≡ R-e mod 2; ≡ 0 mod 3 on order-3 classes;
  even on the diagonal).

  THE RUN: per skeleton class (4 order-3 classes + 3 self-complementary classes 109/221/422)
  assemble the integer system over the 11 (order-3 cosets) / 17 (inversion classes) F_0
  variables with all of the above + per-class combinatorial bounds, and decide feasibility
  with CP-SAT.  Also run the sharper per-(class, u-coset) PERIOD-LEVEL system: the 6 real
  degree-6 Gauss-period components f_m of the shifted counts F(k) reproduce tau_{1,u} and
  ALL its Galois conjugates exactly; every conjugate must land in its box.
"""
import itertools, math, cmath, json, sys
import numpy as np
from ortools.sat.python import cp_model

HERE = "/Users/maceocardinalekwik/git/math/release/h668-z37/code"
P = 37
QR = sorted({(x * x) % P for x in range(1, P)})           # quadratic residues mod 37
QRs = set(QR)
H3 = [1, 10, 26]                                          # order-3 multiplier subgroup
G = 2                                                     # primitive root mod 37
IND = {pow(G, k, P): k for k in range(36)}                # discrete log base 2
INV = {u: pow(u, P - 2, P) for u in range(1, P)}
SQ37 = math.sqrt(37.0)

def chi(u):  # Legendre symbol mod 37
    return 1 if u % P in QRs else (-1 if u % P else 0)

# ----------------------------------------------------------------------------------
# 0. exact constants in Z[sqrt(37)]/2  (pairs (a,b) meaning (a + b*sqrt37)/2)
# ----------------------------------------------------------------------------------
def half_mul(x, y):  # ((a+b s)/2)*((c+d s)/2) = ((ac+37bd) + (ad+bc) s)/4 ; keep /2 rep
    a, b = x; c, d = y
    num_r, num_i = a * c + 37 * b * d, a * d + b * c
    assert num_r % 2 == 0 and num_i % 2 == 0
    return (num_r // 2, num_i // 2)

r_ = (-1, 3)         # r = (-1+3sqrt37)/2
s_ = (-1, -3)        # s = (-1-3sqrt37)/2
rms = (0, 6)         # r - s = 3 sqrt37  -> (0+6s)/2
def hadd(x, y): return (x[0] + y[0], x[1] + y[1])
def hscale(x, k): return (x[0] * k, x[1] * k)

s_rms = half_mul(s_, rms)            # s(r-s)
s2 = half_mul(s_, s_)                # s^2
print("== exact constants ==")
print("  s(r-s) =", s_rms, "/2 each comp -> value (-333-3sqrt37)/2:", s_rms == (-333, -3))
print("  s^2    =", s2, "-> (167+3sqrt37)/2:", s2 == (167, 3))
# sum_w [s(r-s)(a_w+a_uw) + 9 s^2] over 36 w with sum(a_w+a_uw)=324:
tot = hadd(hscale(s_rms, 324), hscale(s2, 36 * 9))
print("  324*s(r-s) + 324*s^2 =", tot, " => rational -26892:", tot == (-26892 * 2, 0))
assert tot == (-2 * 26892, 0)
# beta_u = s(r-s)(a_1+a_u) + 9 s^2  (per-conjugate constant in the period system)
beta_qr = hadd(hscale(s_rms, 8), hscale(s2, 9))    # u in QR  : (-1161+3sqrt37)/2
beta_qnr = hadd(hscale(s_rms, 9), hscale(s2, 9))   # u in QNR : -747
print("  beta(u in QR)  =", beta_qr, " == (-1161+3sqrt37)/2:", beta_qr == (-1161, 3))
print("  beta(u in QNR) =", beta_qnr, " == -747:", beta_qnr == (-1494, 0))
assert beta_qr == (-1161, 3) and beta_qnr == (-1494, 0)
# Master identity constant: 333 X_u = 37 F0(u) - 28224 + 26892 = 37 F0(u) - 1332
#                       =>  X_u = (F0(u)-36)/9.   (28224 = S2 universal, checked per class.)

# ----------------------------------------------------------------------------------
# V1. machine validation of the Galois-trace identity at p=37 on RANDOM block systems
#     (identity is pure Fourier -- needs no srg property)
# ----------------------------------------------------------------------------------
print("\n== V1: random-D Galois-trace identity at p=37, m=9 ==")
rng = np.random.default_rng(668)
zeta = np.exp(2j * np.pi * np.arange(P) / P)
ok = True
for trial in range(12):
    m = 9
    D = [[None] * m for _ in range(m)]
    for i in range(m):
        full = list(range(1, P))
        rng.shuffle(full)
        # negation-closed diagonal without 0
        half = [g for g in range(1, 19)]
        rng.shuffle(half)
        take = half[: rng.integers(2, 9)]
        D[i][i] = sorted(set([g for g in take] + [(P - g) % P for g in take]))
    for i in range(m):
        for j in range(i + 1, m):
            size = int(rng.integers(5, 30))
            sel = list(rng.choice(P, size=size, replace=False))
            D[i][j] = sorted(sel)
            D[j][i] = sorted({(P - g) % P for g in sel})
    M = np.zeros((P, m, m), dtype=complex)
    for i in range(m):
        for j in range(m):
            for g in D[i][j]:
                M[:, i, j] += zeta[(np.arange(P) * g) % P]
    S2r = sum(len(D[i][j]) ** 2 for i in range(m) for j in range(m))
    for u in [1, 2, 5, 36]:
        lhs = sum(np.trace(M[w] @ M[(u * w) % P]) for w in range(1, P))
        F0 = sum(len(set(D[i][j]) & {(u * g) % P for g in D[i][j]})
                 for i in range(m) for j in range(m))
        rhs = P * F0 - S2r
        if abs(lhs - rhs) > 1e-6 * max(1, abs(rhs)):
            ok = False
            print(f"  FAIL trial {trial} u={u}: lhs={lhs}, rhs={rhs}")
print("  12 random systems x 4 twists:", "ALL PASS" if ok else "FAILED")
assert ok

# ----------------------------------------------------------------------------------
# V2. end-to-end machinery validation on Paley P(25) (true positive, exact projectors)
# ----------------------------------------------------------------------------------
print("\n== V2: Paley P(25) end-to-end (m=5 orbits of Z5, srg(25,12,5,6)) ==")
q = 5
F25 = [(a, b) for b in range(q) for a in range(q)]   # a + b*x, x^2 = 2 (2 QNR mod 5)
def f25_mul(u, v):
    a, b = u; c, d = v
    return ((a * c + 2 * b * d) % q, (a * d + b * c) % q)
squares = set()
for el in F25:
    if el != (0, 0):
        squares.add(f25_mul(el, el))
# orbits: index b; D[i][j] = {d in Z5 : (0,(j-i)) + (d,0) in squares}
Dp = [[sorted(d for d in range(q) if ((d) % q, (j - i) % q) in squares)
       for j in range(q)] for i in range(q)]
z5 = np.exp(2j * np.pi * np.arange(q) / q)
Mp = np.zeros((q, q, q), dtype=complex)
for i in range(q):
    for j in range(q):
        for g in Dp[i][j]:
            Mp[:, i, j] += z5[(np.arange(q) * g) % q]
ok = True
for w in range(1, q):
    err = np.abs(Mp[w] @ Mp[w] + Mp[w] - 6 * np.eye(q)).max()
    ok &= err < 1e-9
print("  M_w^2+M_w=6I for w=1..4:", "PASS" if ok else "FAIL")
assert ok
rP, sP = 2.0, -3.0
Pproj = [(Mp[w] - sP * np.eye(q)) / (rP - sP) for w in range(q)]   # r-eigenprojector
aP = [int(round(np.trace(Pproj[w]).real)) for w in range(q)]
for w in range(1, q):
    assert np.abs(Pproj[w] @ Pproj[w] - Pproj[w]).max() < 1e-9
print("  a_w (w=1..4):", aP[1:])
S2p = sum(len(Dp[i][j]) ** 2 for i in range(q) for j in range(q))
ok = True; boxok = True
for u in range(1, q):
    lhs = sum(np.trace(Mp[w] @ Mp[(u * w) % q]) for w in range(1, q)).real
    F0 = sum(len(set(Dp[i][j]) & {(u * g) % q for g in Dp[i][j]})
             for i in range(q) for j in range(q))
    rhs1 = q * F0 - S2p
    Xu = sum(np.trace(Pproj[w] @ Pproj[(u * w) % q]) for w in range(1, q)).real
    rhs2 = (rP - sP) ** 2 * Xu + sP * (rP - sP) * sum(aP[w] + aP[(u * w) % q]
            for w in range(1, q)) + (q - 1) * q * sP ** 2
    ok &= abs(lhs - rhs1) < 1e-8 and abs(lhs - rhs2) < 1e-8
    for w in range(1, q):
        t = np.trace(Pproj[w] @ Pproj[(u * w) % q]).real
        lo = max(0, aP[w] + aP[(u * w) % q] - q)
        hi = min(aP[w], aP[(u * w) % q])
        boxok &= (lo - 1e-9 <= t <= hi + 1e-9)
print("  per-u master identity (both sides):", "PASS" if ok else "FAIL")
print("  all overlap boxes hold:", "PASS" if boxok else "FAIL")
assert ok and boxok

# ----------------------------------------------------------------------------------
# V3. swap-regime sanity at p=13 (Paley(13): scalar blocks, a-pattern flips QR/QNR)
# ----------------------------------------------------------------------------------
print("\n== V3: Paley(13) swap-regime sanity ==")
p13 = 13
QR13 = sorted({(x * x) % p13 for x in range(1, p13)})
z13 = np.exp(2j * np.pi * np.arange(p13) / p13)
M13 = np.array([sum(z13[(w * g) % p13] for g in QR13) for w in range(p13)])
r13, s13 = (-1 + math.sqrt(13)) / 2, (-1 - math.sqrt(13)) / 2
a13 = [None] + [1 if abs(M13[w] - r13) < 1e-9 else 0 for w in range(1, p13)]
swap = all((a13[w] == 1) == (w in QR13) for w in range(1, p13))
print("  a_w = [w in QR] (Galois swap):", "PASS" if swap else "FAIL")
ok = True
for u in range(2, p13):
    lhs = sum(M13[w] * M13[(u * w) % p13] for w in range(1, p13)).real
    F0 = len(set(QR13) & {(u * g) % p13 for g in QR13})
    ok &= abs(lhs - (p13 * F0 - 36)) < 1e-9
print("  Galois-trace identity:", "PASS" if ok else "FAIL")
assert swap and ok

# ----------------------------------------------------------------------------------
# load skeleton classes
# ----------------------------------------------------------------------------------
def load_lines(path):
    return [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
            for ln in open(path) if ln.strip()]

order3 = load_lines(f"{HERE}/order3_classes.txt")
allcls = load_lines(f"{HERE}/orbit_all_classes.txt")
selfc = {k: allcls[k] for k in (109, 221, 422)}
J9 = np.ones((9, 9), dtype=np.int64)

def check_R(R, name):
    assert np.array_equal(R, R.T), name
    assert all(R.sum(1) == 166), name
    assert np.array_equal(R @ R + R - 83 * np.eye(9, dtype=np.int64), 3071 * J9), name
    assert int((R * R).sum()) == 28224, name   # S2 universal
    return True

print("\n== skeleton classes ==")
CLASSES = []
for idx, R in enumerate(order3):
    check_R(R, f"order3-{idx}")
    assert all(int(R[i, j]) % 3 in (0, 1) for i in range(9) for j in range(9))
    CLASSES.append((f"order3-{idx}", R, True))
    print(f"  order3-{idx}: diag {[int(R[i,i]) for i in range(9)]}, "
          f"E(forced)={sum(int(R[i,j])%3==1 for i in range(9) for j in range(9))}")
for k, R in selfc.items():
    check_R(R, f"selfcomp-{k}")
    CLASSES.append((f"selfcomp-{k}", R, False))
    print(f"  selfcomp-{k}: diag {[int(R[i,i]) for i in range(9)]}")

# ----------------------------------------------------------------------------------
# C1. global per-class MIP over the F_0(u) system
# ----------------------------------------------------------------------------------
def comb_lower(R, u, order3_flag):
    """provable lower bound for F_0(u) from |A∩B| >= |A|+|B|-|ambient|"""
    tot = 0
    minus_coset = ({(P - h) % P for h in H3} if order3_flag else {P - 1})
    for i in range(9):
        for j in range(9):
            Rij = int(R[i, j])
            if i == j:
                if u in minus_coset:
                    tot += Rij                      # uD_ii = -D_ii = D_ii (h-invariance)
                else:
                    tot += max(0, 2 * Rij - 36)     # D_ii ⊆ Z37*, negation-closed
            else:
                tot += max(0, 2 * Rij - 37)
    return tot

def canon_u(u, order3_flag):
    """canonical representative of the F_0-symmetry class of u"""
    orb = {u % P, INV[u % P]}
    if order3_flag:
        orb |= {(h * u) % P for h in H3} | {(h * INV[u % P]) % P for h in H3}
    return min(orb)

def build_and_solve(name, R, order3_flag, want_opt=True):
    md = cp_model.CpModel()
    reps = sorted({canon_u(u, order3_flag) for u in range(2, P)} - {1})  # nontrivial classes
    nvar = len(reps)
    # F variables on classes; F(trivial)=1494
    F = {}
    for u0 in reps:
        qr = u0 in QRs
        lo = max(198 if qr else 36, comb_lower(R, u0, order3_flag))
        hi = 1494 if qr else 1332
        F[u0] = md.NewIntVar(lo, hi, f"F_{u0}")
    Ffull = {1: 1494}
    for u in range(2, P):
        if order3_flag and u in H3:
            Ffull[u] = 1494
        else:
            Ffull[u] = F[canon_u(u, order3_flag)]
    # parity
    for u0 in reps:
        t = md.NewIntVar(0, 800, f"half_{u0}")
        md.Add(F[u0] == 2 * t)
    # e matrix
    e = {}
    for i in range(9):
        for j in range(i + 1, 9):
            Rij = int(R[i, j])
            if order3_flag:
                assert Rij % 3 in (0, 1)
                e[(i, j)] = Rij % 3                 # forced
            else:
                e[(i, j)] = md.NewBoolVar(f"e_{i}_{j}")
    Evar = md.NewIntVar(0, 72, "E")
    md.Add(Evar == 2 * sum(e.values()))
    if order3_flag:
        Eval = 2 * sum(e.values())                  # forced integer
        assert Eval % 3 == 1494 % 3 == 0            # consistency: F(1)=1494 ≡ E (mod 3)
        for u0 in reps:                             # F ≡ E (mod 3), all u
            z = md.NewIntVar(-200, 600, f"z3_{u0}")
            md.Add(F[u0] - Eval == 3 * z)
    # x (chi-sum) variables and squares
    x = {}; w2 = {}
    for i in range(9):
        for j in range(i, 9):
            Rij = int(R[i, j])
            eij = 0 if i == j else e[(i, j)]
            xv = md.NewIntVar(-Rij, Rij, f"x_{i}_{j}")
            if not order3_flag and i != j:
                md.Add(xv <= Rij - eij); md.Add(-xv <= Rij - eij)
            # parity: x ≡ R - e (mod 2)
            tp = md.NewIntVar(-Rij, Rij, f"tp_{i}_{j}")
            if order3_flag:
                eijc = 0 if i == j else e[(i, j)]
                md.Add(xv - (Rij - eijc) == 2 * tp)
                z3 = md.NewIntVar(-Rij, Rij, f"x3_{i}_{j}")
                md.Add(xv == 3 * z3)                       # x ≡ 0 mod 3 (H-coset structure)
            else:
                md.Add(xv - Rij + (0 if i == j else e[(i, j)]) == 2 * tp)
            if i == j:
                th = md.NewIntVar(-Rij, Rij, f"xh_{i}_{j}")
                md.Add(xv == 2 * th)                       # diagonal chi-sum even
            wv = md.NewIntVar(0, Rij * Rij, f"w_{i}_{j}")
            md.AddMultiplicationEquality(wv, [xv, xv])
            x[(i, j)] = xv; w2[(i, j)] = wv
    # sums of (R-e)^2 and weights: full matrix = diag + 2*(i<j)
    def lin_Rp2():   # sum_ij (R-e)^2  (linear in e since e binary: (R-e)^2 = R^2-(2R-1)e)
        expr = sum(int(R[i, i]) ** 2 for i in range(9))
        for i in range(9):
            for j in range(i + 1, 9):
                Rij = int(R[i, j])
                if order3_flag:
                    expr += 2 * (Rij - e[(i, j)]) ** 2
                else:
                    expr += 2 * (Rij * Rij) - 2 * (2 * Rij - 1) * e[(i, j)]
        return expr
    def lin_W():
        return sum(w2[(i, i)] for i in range(9)) + 2 * sum(
            w2[(i, j)] for i in range(9) for j in range(i + 1, 9))
    # global sum identity:  sum_{u in Z37*} F(u) = sum (R-e)^2 + 36 E
    md.Add(sum(Ffull[u] for u in range(1, P)) == lin_Rp2() + 36 * Evar)
    # QR-sum identity: 2*sum_{u in QR} F(u) = sum(R-e)^2 + sum x^2 + 36E
    md.Add(2 * sum(Ffull[u] for u in QR) == lin_Rp2() + lin_W() + 36 * Evar)
    # Bochner cuts (tightened inward => any SAT remains honest):
    # sum_u psi_j(u) F(u) >= 0 for characters psi_j of Z36; scaled 1e4, margin 3.0
    SCALE = 10 ** 4
    jlist = ([3, 6, 9, 12, 15] if order3_flag else list(range(1, 18)))
    for jj in jlist:
        coeffs = {}
        for u in range(1, P):
            cval = math.cos(2 * math.pi * jj * IND[u] / 36.0)
            coeffs[u] = int(round(SCALE * cval))
        md.Add(sum(coeffs[u] * Ffull[u] for u in range(1, P)) >= 3 * SCALE)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 60
    status = solver.Solve(md)
    out = {"name": name, "status": solver.StatusName(status)}
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = {u0: solver.Value(F[u0]) for u0 in reps}
        out["F"] = sol
        out["E"] = solver.Value(Evar)
        # post-check Bochner in floats on the full profile
        prof = [1494 if (u == 1 or (order3_flag and u in H3)) else sol[canon_u(u, order3_flag)]
                for u in range(1, P)]
        seq = [prof[pow(G, k, P) - 1] for k in range(36)]
        dft = np.fft.fft(np.array(seq, dtype=float))
        out["bochner_min"] = float(np.real(dft).min())
    return out, md, F, reps

print("\n== C1: per-class global MIP over the F_0(u) twisted-moment system ==")
c1_results = []
for name, R, o3 in CLASSES:
    out, md, F, reps = build_and_solve(name, R, o3)
    c1_results.append(out)
    print(f"  {name}: {out['status']}", end="")
    if "F" in out:
        vals = sorted(out["F"].values())
        print(f"  E={out['E']}  F-range=[{vals[0]},{vals[-1]}]  "
              f"Bochner min DFT={out['bochner_min']:.1f}", end="")
    print()

# slack quantification: min/max of F at a generic QR class and at the -1 class
print("\n  slack intervals (min..max attainable under ALL constraints):")
for name, R, o3 in CLASSES:
    row = {}
    for label, target in (("F(7) [QR]", 7), ("F(36) [-1]", 36),
                          ("F(2) [QNR]", 2)):
        rng_lohi = []
        for sense in (False, True):
            out, md, F, reps = build_and_solve(name, R, o3, want_opt=False)
            tu = canon_u(target, o3)
            if sense: md.Maximize(F[tu])
            else: md.Minimize(F[tu])
            sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = 60
            st = sv.Solve(md)
            rng_lohi.append(sv.Value(F[tu]) if st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
                            else None)
        row[label] = rng_lohi
    print(f"    {name}: " + "  ".join(f"{k}={v[0]}..{v[1]}" for k, v in row.items()))

# ----------------------------------------------------------------------------------
# C2. period-level per-(class, u-coset) system for the 4 order-3 classes
#     unknowns: F0, f_1..f_6 (degree-6 Gauss-period components of shifted counts F(k))
#     six real-embedding constraints: 333*sigma_j(tau_{1,u}) = F0 + sum_m f_m theta_{m+j} - beta_j
#     must lie in [0, 1332] for every j (all conjugate boxes are [0,4]).
# ----------------------------------------------------------------------------------
print("\n== C2: period-level system (order-3 classes), exact conjugate boxes ==")
# cosets of +-H = {1,10,26,36,27,11} (order 6): reps 2^m, m=0..5
pmH = sorted({h % P for h in H3} | {(P - h) % P for h in H3})
cosets6 = []
seen = set()
for mexp in range(36):
    rep = pow(G, mexp, P)
    cs = frozenset((rep * h) % P for h in pmH)
    if cs not in seen:
        seen.add(cs); cosets6.append(cs)
assert len(cosets6) == 6
theta = [sum(math.cos(2 * math.pi * g / P) for g in cs) for cs in cosets6]
print("  degree-6 Gauss periods theta_m:", [f"{t:.6f}" for t in theta])
assert abs(sum(theta) + 1.0) < 1e-9
# Galois: sigma_{2}: coset_m -> coset containing 2*rep
perm = []
for mi, cs in enumerate(cosets6):
    g0 = next(iter(cs))
    tgt = next(k for k, c2 in enumerate(cosets6) if (2 * g0) % P in c2)
    perm.append(tgt)
# chi on coset class: chi of any element
chi6 = [1 if next(iter(cs)) in QRs else -1 for cs in cosets6]

SC = 10 ** 6
results_c2 = []
for name, R, o3 in CLASSES:
    if not o3: continue
    E = sum(int(R[i, j]) % 3 == 1 for i in range(9) for j in range(9))
    ucos = sorted({canon_u(u, True) for u in range(2, P)} - {1})
    line = []
    for u0 in ucos:
        qr = u0 in QRs
        md = cp_model.CpModel()
        lo = max(198 if qr else 36, comb_lower(R, u0, True))
        F0 = md.NewIntVar(lo, 1494 if qr else 1332, "F0")
        t = md.NewIntVar(0, 800, "t"); md.Add(F0 == 2 * t)
        z = md.NewIntVar(-200, 600, "z"); md.Add(F0 - E == 3 * z)
        f = [md.NewIntVar(0, 1494, f"f{m}") for m in range(6)]
        md.Add(F0 + 6 * sum(f) == 28224)
        # six embeddings: apply sigma_{2^j}: theta_m -> theta_{perm^j(m)}, sqrt37 -> chi(2^j)*
        for j in range(6):
            # embedding j corresponds to v = 2^j; chi(v) = (-1)^j
            sgn = 1 if j % 2 == 0 else -1
            if qr:
                beta = (-1161 + sgn * 3 * SQ37) / 2.0
            else:
                beta = -747.0
            # compute permuted theta indices: m -> perm applied j times
            idx = list(range(6))
            for _ in range(j):
                idx = [perm[k] for k in idx]
            coeffs = [int(round(SC * theta[idx[m]])) for m in range(6)]
            expr = SC * F0 + sum(coeffs[m] * f[m] for m in range(6)) - int(round(SC * beta))
            # inward margin 0.05 real units (rounding error << this)
            md.Add(expr >= int(0.05 * SC))
            md.Add(expr <= int((1332 - 0.05) * SC))
        sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = 30
        st = sv.Solve(md)
        line.append((u0, sv.StatusName(st)))
    bad = [x for x in line if x[1] not in ("OPTIMAL", "FEASIBLE")]
    results_c2.append((name, line))
    print(f"  {name}: {len(line)} u-cosets -> " +
          ("ALL FEASIBLE" if not bad else f"INFEASIBLE at {bad}"))

# ----------------------------------------------------------------------------------
# informativeness: where does the overlap bound beat pure combinatorics?
# ----------------------------------------------------------------------------------
print("\n== informativeness of the overlap cuts per class ==")
for name, R, o3 in CLASSES:
    lb_gen = comb_lower(R, 7, o3)            # generic u (not in -coset)
    lb_m1 = comb_lower(R, 36, o3)            # u = -1 coset
    print(f"  {name}: comb LB generic u = {lb_gen} (overlap forces 198 on QR: "
          f"{'BEATS comb' if 198 > lb_gen else 'comb stronger'}), "
          f"comb LB at -1 = {lb_m1}; QNR cap 1332 < 1494 (defect >= 162: always new)")

print("\nDONE")
