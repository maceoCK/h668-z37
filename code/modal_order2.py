"""
Modal: long-running, spectrally-augmented order-2 (negation-invariant) closure attempt for
srg(333,166,82,83) via Z_37. One container per theta-multiplicity a in {3,4,5,6} (the values pinned by
the validated Gauss-sum/Galois reduction), each with the diagonal pair-counts N_g hard-constrained:
   N_g = 3a-9 (g a QR mod 37),  18-3a (g a QNR),   g = 1..18.
Orbit-matrix FREE, every block <36>=<-1>-invariant. UNSAT for all four a => no order-2-cyclotomic
Z_37 srg(333) => Theorem 1 -> 8/8. Any SAT is verified in-container (it would be an H(668)).

Run:  modal run modal_order2.py --time-s 21600     (6h per a, 4 a-values in parallel)
"""
import modal

app = modal.App("h668-order2-spectral")
image = modal.Image.debian_slim(python_version="3.11").pip_install("ortools")
P, NORB, K, LAM, MU = 37, 9, 166, 82, 83
NEG = P - 1

@app.function(image=image, cpu=16.0, timeout=86400, retries=0)
def solve_a(a, time_s, seed=0):
    from ortools.sat.python import cp_model
    qr = set((x*x) % P for x in range(1, P))
    isqr = {g: (g in qr) for g in range(1, (P-1)//2 + 1)}
    m = cp_model.CpModel()
    x = {(i, j, g): m.NewBoolVar(f"x{i}_{j}_{g}")
         for i in range(NORB) for j in range(NORB) for g in range(P)}
    for i in range(NORB):
        for j in range(NORB):
            for g in range(P):
                m.Add(x[j, i, g] == x[i, j, (-g) % P])
    for i in range(NORB):
        m.Add(x[i, i, 0] == 0)
        m.Add(sum(x[i, j, g] for j in range(NORB) for g in range(P)) == K)
    for i in range(NORB):
        for j in range(NORB):
            for g in range(1, P):
                m.Add(x[i, j, g] == x[i, j, (NEG*g) % P])      # negation-invariant blocks
    for g in range(1, (P-1)//2 + 1):                            # spectral N_g constraint
        target = (3*a - 9) if isqr[g] else (18 - 3*a)
        m.Add(sum(x[i, i, g] for i in range(NORB)) == target)
    prod = {}
    def AND(p, q):
        key = (p, q) if p <= q else (q, p)
        if key not in prod:
            c = m.NewBoolVar(f"p{len(prod)}"); m.AddMultiplicationEquality(c, [x[p], x[q]]); prod[key] = c
        return prod[key]
    for i in range(NORB):
        for j in range(NORB):
            for g in range(P):
                terms = [AND((i, ip, d), (ip, j, (g-d) % P)) for ip in range(NORB) for d in range(P)]
                if i == j and g == 0:
                    m.Add(sum(terms) == K)
                else:
                    m.Add(sum(terms) == MU + (LAM - MU) * x[i, j, g])
    sv = cp_model.CpSolver()
    sv.parameters.max_time_in_seconds = float(time_s)
    sv.parameters.num_search_workers = 16
    sv.parameters.random_seed = int(seed)
    st = sv.Solve(m); name = sv.StatusName(st)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [[sorted(g for g in range(P) if sv.Value(x[i, j, g])) for j in range(NORB)] for i in range(NORB)]
        try:
            modal.Dict.from_name("h668-order2", create_if_missing=True)[f"a{a}_s{seed}"] = {"a": a, "D": sol}
        except Exception:
            pass
    return {"a": a, "seed": seed, "status": name, "sat": sol is not None}

@app.local_entrypoint()
def main(time_s: int = 36000, seeds: int = 4):
    # a in {3,6} are ELIMINATED analytically (force an all-18 orbit-matrix diagonal, proven nonexistent);
    # only a in {4,5} remain. Portfolio: several random seeds per a (diversifies the SAT search; a SAT hit
    # at any seed would be an H(668), and the longest runs push the UNSAT frontier).
    AS = (4, 5)
    jobs = [(a, time_s, s) for a in AS for s in range(seeds)]
    print(f"order-2 frontier on Modal: a in {AS}, {seeds} seeds each, {time_s}s, 16 cores -> {len(jobs)} containers")
    res = {}
    for r in solve_a.starmap(jobs):
        res[(r["a"], r["seed"])] = r["status"]
        tag = "  *** SAT -> H(668)! ***" if r["sat"] else ""
        print(f"  a={r['a']} seed={r['seed']}: {r['status']}{tag}", flush=True)
    byA = {a: set(res[(aa, s)] for (aa, s) in res if aa == a) for a in AS}
    if any(st in ("OPTIMAL", "FEASIBLE") for st in res.values()):
        print(f"=== SAT FOUND -> verify immediately (possible H(668)). {res} ===")
    elif all("INFEASIBLE" in byA[a] for a in AS):
        print("=== both a in {4,5} INFEASIBLE (some seed): order 2 CLOSED -> Theorem 1 = 8/8. ===")
    else:
        print(f"=== unresolved: {res}. a in {{3,6}} eliminated; a in {{4,5}} still open. ===")
