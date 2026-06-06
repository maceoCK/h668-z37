"""
Modal CPU app: ORBIT-MATRIX-FREE cyclotomic-nonexistence proof for srg(333,166,82,83).
For each multiplier subgroup <mult> of Z_37*, search for ANY 166-regular block-circulant srg whose
every connection set is <mult>-invariant, with the orbit matrix LEFT FREE. UNSAT for all 8 subgroups
=> no cyclotomically-structured Z_37-symmetric srg(333) exists (complete; no skeleton sampling).
Encoding validated locally (cyclo_complete.py reconstructs Paley P(13)/P(17)).

Run:  modal run modal_cyclo.py --time-s 7200
"""
import modal

app = modal.App("h668-cyclo-complete")
image = modal.Image.debian_slim(python_version="3.11").pip_install("ortools")
P, NORB, K, LAM, MU = 37, 9, 166, 82, 83

@app.function(image=image, cpu=16.0, timeout=28800, retries=0)
def cyclo_check(mult, time_s):
    from ortools.sat.python import cp_model
    p, norb, k, lam, mu = P, NORB, K, LAM, MU
    m = cp_model.CpModel()
    x = {(i, j, g): m.NewBoolVar(f"x_{i}_{j}_{g}")
         for i in range(norb) for j in range(norb) for g in range(p)}
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                m.Add(x[j, i, g] == x[i, j, (-g) % p])
    for i in range(norb):
        m.Add(x[i, i, 0] == 0)
        m.Add(sum(x[i, j, g] for j in range(norb) for g in range(p)) == k)   # 166-regular, R free
    for i in range(norb):
        for j in range(norb):
            for g in range(1, p):
                m.Add(x[i, j, g] == x[i, j, (mult * g) % p])                  # <mult>-invariant blocks
    prod = {}
    def AND(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in prod:
            c = m.NewBoolVar(f"p{len(prod)}"); m.AddMultiplicationEquality(c, [x[a], x[b]]); prod[key] = c
        return prod[key]
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                terms = [AND((i, ip, d), (ip, j, (g - d) % p)) for ip in range(norb) for d in range(p)]
                if i == j and g == 0:
                    m.Add(sum(terms) == k)
                else:
                    m.Add(sum(terms) == mu + (lam - mu) * x[i, j, g])
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = float(time_s)
    sv.parameters.num_search_workers = 16
    st = sv.Solve(m); name = sv.StatusName(st)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [[sorted(g for g in range(p) if sv.Value(x[i, j, g])) for j in range(norb)] for i in range(norb)]
    return {"mult": mult, "status": name, "solution": sol}

@app.local_entrypoint()
def main(time_s: int = 7200):
    def primroot(p):
        for g in range(2, p):
            x, seen = 1, set()
            for _ in range(p - 1):
                x = x * g % p; seen.add(x)
            if len(seen) == p - 1: return g
    g = primroot(P)
    mults = {e: pow(g, (P - 1) // e, P) for e in (2, 3, 4, 6, 9, 12, 18, 36)}
    print(f"orbit-matrix-free cyclotomic check, {time_s}s each, multipliers (order->mult): {mults}")
    jobs = [(mlt, time_s) for mlt in mults.values()]
    res = {}
    for r in cyclo_check.starmap(jobs):
        res[r["mult"]] = r["status"]
        tag = "  *** SAT -> H(668)! ***" if r["solution"] else ""
        print(f"  mult={r['mult']:2d}: {r['status']}{tag}")
    allinf = all(v == "INFEASIBLE" for v in res.values())
    print("=== ALL INFEASIBLE: NO cyclotomically-structured Z_37 srg(333) exists (COMPLETE). ===" if allinf
          else f"=== not all resolved: {res} ===")
