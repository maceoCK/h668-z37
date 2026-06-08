"""
SOUND + fast cyclotomic (H-invariant) lift solver, coset-level.

Background: the element-level exact_lift with a gauge-fix (force 0 in D[0][i]) is UNSOUND for the
H-cyclotomic ansatz -- an H-invariant set shifted by c!=0 is no longer H-invariant ((mult-1)c=0 => c=0),
so the cyclotomic ansatz has NO base-point/gauge freedom, and the gauge-fix falsely excludes valid graphs
(it reports P(25) -- a real negation-invariant srg -- as INFEASIBLE). Dropping the gauge-fix is sound but
the element-level model is then too slow (UNKNOWN).

This solver is sound (no gauge) and small: variables are which H-cosets lie in each block, and the
convolution uses precomputed coset-addition counts (no per-element AND products). Validated to reconstruct
the negation-invariant Paley graphs P(9), P(25).
"""
from ortools.sat.python import cp_model

def Hcosets(p, mult):
    H = set(); x = 1
    while True:
        H.add(x); x = x * mult % p
        if x == 1: break
    seen = set(); cs = []
    for a in range(1, p):
        if a in seen: continue
        c = tuple(sorted({(a*h) % p for h in H}))
        for v in c: seen.add(v)
        cs.append(c)
    return cs                     # size-e cosets partitioning Z_p^*

def solve_cyclo_coset(p, norb, k, lam, mu, R, mult, time_s=120, workers=8):
    cs = Hcosets(p, mult); nc = len(cs); e = len(cs[0])
    cidx = {}                     # element -> coset index
    for ci, c in enumerate(cs):
        for g in c: cidx[g] = ci
    negc = [cidx[(-cs[ci][0]) % p] for ci in range(nc)]      # coset index of -c
    cnt = [[[0]*p for _ in range(nc)] for _ in range(nc)]    # cnt[c][c'][g] = #{(a,b) in c x c': a+b=g}
    for ci in range(nc):
        for cj in range(nc):
            t = cnt[ci][cj]
            for a in cs[ci]:
                for b in cs[cj]:
                    t[(a+b) % p] += 1
    m = cp_model.CpModel()
    yc = {(i, j, c): m.NewBoolVar(f"y{i}_{j}_{c}") for i in range(norb) for j in range(norb) for c in range(nc)}
    z0 = {(i, j): m.NewBoolVar(f"z{i}_{j}") for i in range(norb) for j in range(norb)}
    # symmetry D[j][i] = -D[i][j]
    for i in range(norb):
        for j in range(norb):
            m.Add(z0[j, i] == z0[i, j])
            for c in range(nc):
                m.Add(yc[j, i, c] == yc[i, j, negc[c]])
    # diagonal: no loops, negation-invariant
    for i in range(norb):
        m.Add(z0[i, i] == 0)
        for c in range(nc):
            m.Add(yc[i, i, c] == yc[i, i, negc[c]])
    # cardinality |D[i][j]| = e*sum_c yc + z0 = R[i][j]
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
    # convolution per (i,j,g):  sum_{i'} [ sum_{c,c'} yc[i,i',c] yc[i',j,c'] cnt[c][c'][g]
    #   + z0[i,i'] inD(i',j,g) + z0[i',j] inD(i,i',g) - (g==0) z0[i,i'] z0[i',j] ] = mu + (k-mu)[i=j,g=0] + (lam-mu)[g in D]
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
                rhs_const = mu + ((k - mu) if (i == j and g == 0) else 0)
                m.Add(sum(terms) == rhs_const + (lam - mu) * inD(i, j, g))
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = float(time_s); sv.parameters.num_search_workers = workers
    st = sv.Solve(m); name = sv.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        D = [[set() for _ in range(norb)] for _ in range(norb)]
        for i in range(norb):
            for j in range(norb):
                s = set()
                if sv.Value(z0[i, j]): s.add(0)
                for c in range(nc):
                    if sv.Value(yc[i, j, c]): s |= set(cs[c])
                D[i][j] = s
        return name, D
    return name, None

def validate():
    import conf_core as cc
    for p, (k, lam, mu), neg in [(3, (4, 1, 2), 2), (5, (12, 5, 6), 4)]:
        v = p*p; k = (v-1)//2; lam = (v-5)//4; mu = (v-1)//4
        D = cc.paley_pp_blocks(p); R = cc.orbit_matrix(D, p)
        # confirm Paley(p^2) blocks are negation-invariant (order-2-cyclotomic)
        ninv = all({(-g) % p for g in D[i][j]} == set(D[i][j]) for i in range(p) for j in range(p))
        st, Ds = solve_cyclo_coset(p, p, k, lam, mu, R, neg, time_s=60)
        ok = Ds is not None and cc.srg_identity_direct(cc.build_adjacency(Ds, p, p), k, lam, mu)
        print(f"  P({v})=srg({v},{k},{lam},{mu}) order-2 (mult={neg}): neg-invariant={ninv}  status={st}  reconstructs={ok}")

if __name__ == "__main__":
    validate()
