"""
COMPLETE cyclotomic-nonexistence check (upgrades Theorem 1 from a skeleton-sample to ALL skeletons).
For a multiplier m in Z_p*, search for ANY k-regular block-circulant srg(v,k,lam,mu) whose every
connection set D[i][j] is <m>-invariant (a union of cyclotomic cosets) -- with the orbit matrix R
LEFT FREE (only k-regularity imposed). If UNSAT for every nontrivial subgroup of Z_p*, then NO
cyclotomically-structured Z_p-symmetric srg exists, over all orbit matrices simultaneously.

Sound (no gauge over-constraint): INFEASIBLE is a genuine nonexistence within the cyclotomic class.
Validated: on a single circulant (norb=1) it FINDS the Paley srg(p,...) under the QR multiplier.
"""
from ortools.sat.python import cp_model
import conf_core as cc

def primroot(p):
    for g in range(2, p):
        x, seen = 1, set()
        for _ in range(p - 1):
            x = x * g % p; seen.add(x)
        if len(seen) == p - 1: return g

def subgroup_gens(p):
    g = primroot(p); out = {}
    for e in range(2, p):
        if (p - 1) % e == 0:
            out[e] = pow(g, (p - 1) // e, p)
    return out  # order e -> generator of the order-e subgroup

def solve_cyclo(p, norb, k, lam, mu, mult, time_s=300, workers=8):
    m = cp_model.CpModel()
    x = {(i, j, g): m.NewBoolVar(f"x_{i}_{j}_{g}")
         for i in range(norb) for j in range(norb) for g in range(p)}
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                m.Add(x[j, i, g] == x[i, j, (-g) % p])
    for i in range(norb):
        m.Add(x[i, i, 0] == 0)
    # k-regularity (R left free; this is the only cardinality constraint)
    for i in range(norb):
        m.Add(sum(x[i, j, g] for j in range(norb) for g in range(p)) == k)
    # cyclotomic ansatz: blocks invariant under <mult>
    for i in range(norb):
        for j in range(norb):
            for g in range(1, p):
                m.Add(x[i, j, g] == x[i, j, (mult * g) % p])
    # convolution / SRG conditions (NO gauge-fix -> INFEASIBLE stays sound)
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
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = time_s; sv.parameters.num_search_workers = workers
    st = sv.Solve(m); name = sv.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        D = [[set(g for g in range(p) if sv.Value(x[i, j, g])) for j in range(norb)] for i in range(norb)]
        return name, D
    return name, None

def validate():
    # single circulant (norb=1): Paley srg(p,(p-1)/2,(p-5)/4,(p-1)/4) under the QR multiplier -> SAT
    for p in (13, 17):
        k = (p - 1) // 2; lam = (p - 5) // 4; mu = (p - 1) // 4
        gens = subgroup_gens(p); qr_mult = gens[(p - 1) // 2]   # order-(p-1)/2 subgroup = QRs
        st, D = solve_cyclo(p, 1, k, lam, mu, qr_mult, time_s=30)
        ok = D is not None and cc.srg_identity_direct(cc.build_adjacency(D, p, 1), k, lam, mu)
        print(f"   validate Paley P({p})=srg({p},{k},{lam},{mu}) under QR multiplier: status={st} found_srg={ok}")

if __name__ == "__main__":
    import sys
    if "validate" in sys.argv:
        print("=== validate cyclo-complete model finds Paley single circulants ===")
        validate()
    else:
        p, norb, k, lam, mu = 37, 9, 166, 82, 83
        gens = subgroup_gens(p)
        ts = int(sys.argv[1]) if len(sys.argv) > 1 else 300
        print(f"=== COMPLETE cyclotomic check srg(333): {len(gens)} subgroups, {ts}s each ===")
        allinf = True
        for e, mlt in sorted(gens.items()):
            st, D = solve_cyclo(p, norb, k, lam, mu, mlt, time_s=ts)
            print(f"   subgroup order {e:2d} (mult={mlt:2d}): {st}" + ("  *** SAT -> H(668)! ***" if D else ""))
            if st != "INFEASIBLE": allinf = False
            if D:
                import numpy as np; np.save("cyclo_solution.npy", np.array([[sorted(D[i][j]) for j in range(norb)] for i in range(norb)], dtype=object))
        print("ALL INFEASIBLE -> NO cyclotomically-structured Z_37 srg(333) exists (complete)." if allinf
              else "Not all INFEASIBLE (see above).")
