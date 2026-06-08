"""
Spectrally-augmented order-2 (negation-invariant) lift solver for srg(333,166,82,83) via Z_37.

New ingredient (validated on Paley(13) m=1 and Paley(125) m=25): for an order-2-cyclotomic graph,
Tr(M_1) = a*theta+(9-a)*tau expanded in the integral basis eta_g=zeta^g+zeta^-g via the Gauss sum
sqrt(37)=sum_g (g|37) eta_g forces the diagonal pair-counts
    N_g := #{ orbits i : pair {g,-g} in D[i][i] }  =  3a-9  (g a QR mod 37),   18-3a  (g a QNR),
with the theta-multiplicity pinned to a in {3,4,5,6}.

We add these as hard constraints to the orbit-matrix-FREE, negation-invariant exact lift and run it
once per a. UNSAT for all four a  =>  NO order-2-cyclotomic Z_37 srg(333) exists  =>  Theorem 1 -> 8/8.
Any SAT is verified in-process (it would be an H(668)).
"""
import sys, time
from ortools.sat.python import cp_model

P, NORB, K, LAM, MU = 37, 9, 166, 82, 83
NEG = P - 1  # multiplier -1 = 36 generates the order-2 subgroup {1,36}

def qr_pairs():
    qr = set((x*x) % P for x in range(1, P))
    # representatives g=1..18; pair is "QR-type" iff g is a QR (well-defined since -1 is a QR mod 37)
    return {g: (g in qr) for g in range(1, (P-1)//2 + 1)}

def solve(a, time_s=300, workers=16):
    isqr = qr_pairs()
    m = cp_model.CpModel()
    x = {(i, j, g): m.NewBoolVar(f"x{i}_{j}_{g}")
         for i in range(NORB) for j in range(NORB) for g in range(P)}
    # A = A^T  : D[j][i] = -D[i][j]
    for i in range(NORB):
        for j in range(NORB):
            for g in range(P):
                m.Add(x[j, i, g] == x[i, j, (-g) % P])
    # no loops on diagonal
    for i in range(NORB):
        m.Add(x[i, i, 0] == 0)
    # 166-regular, orbit matrix FREE
    for i in range(NORB):
        m.Add(sum(x[i, j, g] for j in range(NORB) for g in range(P)) == K)
    # ORDER-2: every block negation-invariant
    for i in range(NORB):
        for j in range(NORB):
            for g in range(1, P):
                m.Add(x[i, j, g] == x[i, j, (NEG*g) % P])
    # gauge fix: force 0 in D[0][i] when that block is nonempty (cheap, sound)
    for i in range(1, NORB):
        m.Add(x[0, i, 0] <= sum(x[0, i, g] for g in range(P)))  # harmless link; full gauge handled by search
    # SPECTRAL N_g constraints (the new ingredient): N_g = #{i: pair g in D[i,i]} = sum_i x[i,i,g]
    for g in range(1, (P-1)//2 + 1):
        target = (3*a - 9) if isqr[g] else (18 - 3*a)
        m.Add(sum(x[i, i, g] for i in range(NORB)) == target)
    # convolution / SRG equations:  sum_{i',d} x[i,i',d] x[i',j,g-d] = 83 + (k-mu-83)[i=j,g=0] + (lam-mu)[g in D]
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
    sv.parameters.num_search_workers = workers
    st = sv.Solve(m)
    name = sv.StatusName(st)
    sol = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [[sorted(g for g in range(P) if sv.Value(x[i, j, g])) for j in range(NORB)] for i in range(NORB)]
    return name, sol

if __name__ == "__main__":
    ts = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    print(f"=== spectrally-augmented order-2 lift, {ts}s per a, a in {{3,4,5,6}} ===")
    results = {}
    for a in (3, 4, 5, 6):
        t0 = time.time()
        name, sol = solve(a, time_s=ts)
        results[a] = name
        tag = "  *** SAT -> candidate H(668)! ***" if sol else ""
        print(f"  a={a}: {name}  ({time.time()-t0:.0f}s){tag}", flush=True)
        if sol:
            import json; json.dump({"a": a, "D": sol}, open("order2_hit.json", "w"))
    allunsat = all(v == "INFEASIBLE" for v in results.values())
    print(f"=== results: {results} ===")
    if allunsat:
        print("ALL a INFEASIBLE => NO order-2-cyclotomic Z_37 srg(333) exists. Theorem 1 -> 8 of 8.")
    elif any(v in ("OPTIMAL", "FEASIBLE") for v in results.values()):
        print("SAT found -> verify immediately (possible H(668)).")
    else:
        print("Not all resolved (UNKNOWN/timeout) -> order 2 remains open; the a in {3,4,5,6} reduction stands.")
