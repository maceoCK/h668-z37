"""
Modal: SOUND coset-level cyclotomic lift, long-timeout decidability test (then full sweep).
No gauge-fix (the gauge-fix is unsound with the multiplier constraint). Decides whether the per-R
order-e (mult) lift is SAT (-> H(668) candidate) or UNSAT (-> closes that subgroup over R), or remains
UNKNOWN. Validated locally: reconstructs P(9), P(25).

Run (decidability test):  modal run modal_cyclo_coset.py --mode test --time-s 1800
"""
import modal
app = modal.App("h668-cyclo-coset")
image = modal.Image.debian_slim(python_version="3.11").pip_install("ortools", "numpy")
P = 37

@app.function(image=image, cpu=8.0, timeout=43200, retries=1)
def coset_lift(R, mult, time_s):
    from ortools.sat.python import cp_model
    p = P; norb = 9; k, lam, mu = 166, 82, 83
    H = set(); x = 1
    while True:
        H.add(x); x = x*mult % p
        if x == 1: break
    seen = set(); cs = []
    for a in range(1, p):
        if a in seen: continue
        c = tuple(sorted({(a*h) % p for h in H}))
        for v in c: seen.add(v)
        cs.append(c)
    nc = len(cs); e = len(cs[0]); cidx = {}
    for ci, c in enumerate(cs):
        for g in c: cidx[g] = ci
    negc = [cidx[(-cs[ci][0]) % p] for ci in range(nc)]
    cnt = [[[0]*p for _ in range(nc)] for _ in range(nc)]
    for ci in range(nc):
        for cj in range(nc):
            t = cnt[ci][cj]
            for a in cs[ci]:
                for b in cs[cj]:
                    t[(a+b) % p] += 1
    m = cp_model.CpModel()
    yc = {(i, j, c): m.NewBoolVar(f"y{i}_{j}_{c}") for i in range(norb) for j in range(norb) for c in range(nc)}
    z0 = {(i, j): m.NewBoolVar(f"z{i}_{j}") for i in range(norb) for j in range(norb)}
    for i in range(norb):
        for j in range(norb):
            m.Add(z0[j, i] == z0[i, j])
            for c in range(nc):
                m.Add(yc[j, i, c] == yc[i, j, negc[c]])
    for i in range(norb):
        m.Add(z0[i, i] == 0)
        for c in range(nc):
            m.Add(yc[i, i, c] == yc[i, i, negc[c]])
    for i in range(norb):
        for j in range(norb):
            m.Add(e*sum(yc[i, j, c] for c in range(nc)) + z0[i, j] == int(R[i][j]))
    prod = {}
    def AND(a, b):
        key = (a, b) if id(a) <= id(b) else (b, a)
        if key not in prod:
            v = m.NewBoolVar(f"p{len(prod)}"); m.AddMultiplicationEquality(v, [a, b]); prod[key] = v
        return prod[key]
    def inD(pp, qq, g):
        return z0[pp, qq] if g == 0 else yc[pp, qq, cidx[g]]
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                terms = []
                for ip in range(norb):
                    for c in range(nc):
                        for cpp in range(nc):
                            w = cnt[c][cpp][g]
                            if w: terms.append(w * AND(yc[i, ip, c], yc[ip, j, cpp]))
                    terms.append(AND(z0[i, ip], inD(ip, j, g)))
                    terms.append(AND(z0[ip, j], inD(i, ip, g)))
                    if g == 0:
                        terms.append(-AND(z0[i, ip], z0[ip, j]))
                m.Add(sum(terms) == mu + ((k-mu) if (i == j and g == 0) else 0) + (lam-mu)*inD(i, j, g))
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = float(time_s); sv.parameters.num_search_workers = 8
    st = sv.Solve(m); name = sv.StatusName(st)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [[sorted(([0] if sv.Value(z0[i, j]) else []) + [g for c in range(nc) if sv.Value(yc[i, j, c]) for g in cs[c]]) for j in range(norb)] for i in range(norb)]
    return {"status": name, "sat": sol is not None, "sol": sol, "R": R}

def _load(path):
    out = []
    for line in open(path):
        v = [int(t) for t in line.split()]
        if len(v) == 81: out.append([v[9*i:9*i+9] for i in range(9)])
    return out

@app.local_entrypoint()
def main(mode: str = "test", time_s: int = 1800):
    import os
    here = os.path.dirname(__file__)
    r3 = _load(os.path.join(here, "orbit3_compatible.txt"))
    r2 = _load(os.path.join(here, "orbit_all_sym.txt"))
    if mode == "test":
        jobs = [(r3[0], 10, time_s), (r3[len(r3)//2], 10, time_s), (r3[-1], 10, time_s),
                (r2[0], 36, time_s), (r2[len(r2)//2], 36, time_s), (r2[-1], 36, time_s)]
        print(f"DECIDABILITY TEST: 3 order-3 + 3 order-2 reps, {time_s}s each, SOUND coset solver")
        for r in coset_lift.starmap(jobs):
            print(f"  mult? status={r['status']}{'  *** SAT -> H(668) candidate! ***' if r['sat'] else ''}", flush=True)
