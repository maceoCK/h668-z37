"""
General SOUND coset-level CP-SAT model for H-cyclotomic lifts of an orbit matrix.
(No gauge-fix anywhere; class-level equations: conv counts of H-invariant sets are constant on
H-classes, so each pair (i<=j) contributes nc+1 equations instead of p.)

Used two ways:
  * validate(): true-positive anchors (P(9), P(25), Petersen, Paley(13)) must come back FEASIBLE
    and verify; brute-force-confirmed empty cases must come back INFEASIBLE.  This is what makes a
    production INFEASIBLE verdict trustworthy.
  * solve() on srg(333,166,82,83), p=37, mult=10 (order 3), with the VALIDATED pinning a=4
    (WLOG: a QNR decimation maps a=5 solutions to a=4 on the same skeleton).

A FEASIBLE on production = candidate Z37 order-3-cyclotomic srg(333) -> verify -> H(668).
INFEASIBLE on all 28 compatible skeletons = order-3 cyclotomic case CLOSED.
"""
import numpy as np
from ortools.sat.python import cp_model

def cosets_of(mult, p):
    H = []; x = 1
    while True:
        H.append(x); x = x * mult % p
        if x == 1: break
    seen, cs = set(), []
    for g in range(1, p):
        if g in seen: continue
        c = tuple(sorted((g * h) % p for h in H))
        seen |= set(c); cs.append(c)
    return H, cs

def solve(p, norb, k, lam, mu, mult, R, a_pin=None, time_s=1800, workers=8, log=False):
    H, CS = cosets_of(mult, p)
    nc, e = len(CS), len(H)
    cidx = {g: ci for ci, c in enumerate(CS) for g in c}
    negc = [cidx[(-c[0]) % p] for c in CS]
    QR = {(x * x) % p for x in range(1, p)}
    chi = [1 if c[0] in QR else -1 for c in CS]
    cnt0 = [[e if negc[ci] == cj else 0 for cj in range(nc)] for ci in range(nc)]
    cntC = [[[0] * nc for _ in range(nc)] for _ in range(nc)]
    for ci in range(nc):
        for cj in range(nc):
            t = [0] * p
            for x in CS[ci]:
                for y in CS[cj]:
                    t[(x + y) % p] += 1
            for d in range(nc):
                g = CS[d][0]
                assert all(t[gg] == t[g] for gg in CS[d]), "class constancy"
                cntC[ci][cj][d] = t[g]
    m = cp_model.CpModel()
    Y, Z = {}, {}
    for i in range(norb):
        for j in range(i, norb):
            r = int(R[i][j])
            if i == j:
                need = 2 * e if e % 2 == 1 else e             # diagonal: +-paired cosets, no 0
                if r % need: return "INFEASIBLE", None
                Z[i, j] = 0
            elif e >= 2:
                if r % e > 1: return "INFEASIBLE", None       # incompatible size
                Z[i, j] = r % e
            else:
                Z[i, j] = m.NewBoolVar(f"z{i}_{j}")
            for c in range(nc):
                Y[i, j, c] = m.NewBoolVar(f"y{i}_{j}_{c}")
            zr = Z[i, j]
            if isinstance(zr, int):
                m.Add(sum(Y[i, j, c] for c in range(nc)) * e == r - zr)
            else:
                m.Add(sum(Y[i, j, c] for c in range(nc)) * e + zr == r)
    def y(i, j, c):
        return Y[i, j, c] if i <= j else Y[j, i, negc[c]]
    def z(i, j):
        return Z[i, j] if i <= j else Z[j, i]
    for i in range(norb):
        for c in range(nc):
            m.Add(Y[i, i, c] == Y[i, i, negc[c]])
    if a_pin is not None:
        assert p == 37 and norb == 9 and (k, lam, mu) == (166, 82, 83) and all(h in QR for h in H)
        for c in range(nc):
            t = 3 * a_pin - 9 if chi[c] == 1 else 18 - 3 * a_pin
            m.Add(sum(Y[i, i, c] for i in range(norb)) == t)
    prod = {}
    def AND(va, vb):
        if isinstance(va, int): return vb if va else 0
        if isinstance(vb, int): return va if vb else 0
        key = (va.Index(), vb.Index()) if va.Index() <= vb.Index() else (vb.Index(), va.Index())
        if key not in prod:
            v = m.NewBoolVar(f"p{len(prod)}")
            m.AddMultiplicationEquality(v, [va, vb]); prod[key] = v
        return prod[key]
    for i in range(norb):
        for j in range(i, norb):
            for d in list(range(nc)) + ["zero"]:
                terms = []
                for l in range(norb):
                    zil, zlj = z(i, l), z(l, j)
                    for c in range(nc):
                        for cq in range(nc):
                            w = cnt0[c][cq] if d == "zero" else cntC[c][cq][d]
                            if w:
                                terms.append(w * AND(y(i, l, c), y(l, j, cq)))
                    if d == "zero":
                        t = AND(zil, zlj) if not (isinstance(zil, int) and isinstance(zlj, int)) else zil * zlj
                        if not isinstance(t, int) or t: terms.append(t)
                    else:
                        if isinstance(zil, int):
                            if zil: terms.append(y(l, j, d))
                        else: terms.append(AND(zil, y(l, j, d)))
                        if isinstance(zlj, int):
                            if zlj: terms.append(y(i, l, d))
                        else: terms.append(AND(zlj, y(i, l, d)))
                rhs = mu + ((k - mu) if (i == j and d == "zero") else 0)
                ind = z(i, j) if d == "zero" else y(i, j, d)
                m.Add(sum(terms) == rhs + (lam - mu) * ind)
    sv = cp_model.CpSolver()
    sv.parameters.max_time_in_seconds = float(time_s)
    sv.parameters.num_search_workers = workers
    if log: sv.parameters.log_search_progress = True
    st = sv.Solve(m)
    name = sv.StatusName(st)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [[set() for _ in range(norb)] for _ in range(norb)]
        for i in range(norb):
            for j in range(norb):
                s = set()
                zv = z(i, j)
                if (zv if isinstance(zv, int) else sv.Value(zv)): s.add(0)
                for c in range(nc):
                    if sv.Value(y(i, j, c)): s |= set(CS[c])
                sol[i][j] = s
    return name, sol

def validate():
    import conf_core as cc
    ok = True
    def chk(name, cond, info=""):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {info}")
        if not cond: ok = False
    for p, mult in ((3, 2), (5, 4)):
        v = p * p; k = (v - 1) // 2; lam = (v - 5) // 4; mu = (v - 1) // 4
        D = cc.paley_pp_blocks(p)
        R = [[len(D[i][j]) for j in range(p)] for i in range(p)]
        st, sol = solve(p, p, k, lam, mu, mult, R, time_s=120)
        good = sol is not None and cc.srg_identity_direct(cc.build_adjacency(sol, p, p), k, lam, mu)
        chk(f"P({v}) mult={mult} FEASIBLE+verifies", st in ("OPTIMAL", "FEASIBLE") and good, st)
    st, sol = solve(5, 2, 3, 0, 1, 4, [[2, 1], [1, 2]], time_s=60)
    good = sol is not None and cc.srg_identity_direct(cc.build_adjacency(sol, 5, 2), 3, 0, 1)
    chk("Petersen mult=4 FEASIBLE+verifies", st in ("OPTIMAL", "FEASIBLE") and good, st)
    st, sol = solve(13, 1, 6, 2, 3, 3, [[6]], time_s=60)
    chk("Paley(13) mult=3 FEASIBLE", st in ("OPTIMAL", "FEASIBLE"), st)
    # brute-force-confirmed empty cases must be INFEASIBLE
    st, _ = solve(7, 3, 10, 4, 5, 6, [[2, 4, 4], [4, 2, 4], [4, 4, 2]], time_s=120)
    chk("srg(21) case A INFEASIBLE", st == "INFEASIBLE", st)
    st, _ = solve(7, 3, 10, 4, 5, 6, [[4, 3, 3], [3, 4, 3], [3, 3, 4]], time_s=120)
    chk("srg(21) case B INFEASIBLE", st == "INFEASIBLE", st)
    st, _ = solve(13, 2, 10, 3, 4, 3, [[6, 3], [3, 6]], time_s=120)
    chk("srg(26) z=0 case INFEASIBLE (brute=0)", st == "INFEASIBLE", st)
    print("MODEL VALIDATION:", "ALL PASS" if ok else "FAILURES")
    return ok

if __name__ == "__main__":
    validate()
