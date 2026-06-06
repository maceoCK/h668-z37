"""
EXACT lift solver (SAT/CP) for the Z_p block-circulant srg(v,k,lam,mu) problem.
Given orbit matrix R (norb x norb), find connection sets D[i][j] subset Z_p with |D[i][j]|=R[i][j],
D[j][i]=-D[i][j], D[i][i] symmetric (0 not in), such that the block-circulant A is srg(v,k,lam,mu).

Exact condition (A^2 = kI + lam A + mu(J-I-A)) in block-circulant form, for all i,j in Z_norb, g in Z_p:
  conv(i,j,g) := sum_{i'} sum_d [d in D[i][i']]*[(g-d) in D[i'][j]]  ==  k*[i=j,g=0] + lam*x[i][j][g] + mu*(1 - x[i][j][g] - [i=j,g=0])
Encoded in CP-SAT: x[i,j,g] Bool (membership), AND-products, cardinality, symmetry.
Validated on Paley P(9)=srg(9,4,1,2) and P(25)=srg(25,12,5,6) before use on srg(333).
"""
from ortools.sat.python import cp_model
import numpy as np
import conf_core as cc

def solve_lift(p, norb, k, lam, mu, R, time_s=120, workers=8, seed=0, gauge=True, mult=None, log=False):
    m = cp_model.CpModel()
    x = {(i, j, g): m.NewBoolVar(f"x_{i}_{j}_{g}")
         for i in range(norb) for j in range(norb) for g in range(p)}
    # symmetry A=A^T : block(j,i) is transpose of block(i,j) -> D[j][i] = -D[i][j]
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                m.Add(x[j, i, g] == x[i, j, (-g) % p])
    # no loops on diagonal blocks
    for i in range(norb):
        m.Add(x[i, i, 0] == 0)
    # cardinality |D[i][j]| = R[i][j]
    for i in range(norb):
        for j in range(norb):
            m.Add(sum(x[i, j, g] for g in range(p)) == int(R[i][j]))
    # GAUGE FIXING: each orbit i has a free base point; re-basing shifts block (i,j) by (c_i-c_j).
    # This is a p^(norb-1) symmetry. Fix it WLOG by forcing 0 in D[0][i] for every i with an edge.
    if gauge:
        for i in range(1, norb):
            if R[0][i] > 0:
                m.Add(x[0, i, 0] == 1)
    # CYCLOTOMIC ANSATZ (optional, incomplete but sound): force every block invariant under the
    # multiplier 'mult' (so each D[i][j] is a union of <mult>-cosets). Tiny structured search.
    if mult is not None:
        for i in range(norb):
            for j in range(norb):
                for g in range(1, p):
                    m.Add(x[i, j, g] == x[i, j, (mult * g) % p])
    # convolution / SRG equations
    prod = {}
    def AND(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in prod:
            c = m.NewBoolVar(f"p{len(prod)}")
            m.AddMultiplicationEquality(c, [x[a], x[b]])
            prod[key] = c
        return prod[key]
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                terms = [AND((i, ip, d), (ip, j, (g - d) % p)) for ip in range(norb) for d in range(p)]
                if i == j and g == 0:
                    m.Add(sum(terms) == k)
                else:
                    m.Add(sum(terms) == mu + (lam - mu) * x[i, j, g])
    sv = cp_model.CpSolver()
    sv.parameters.max_time_in_seconds = time_s
    sv.parameters.num_search_workers = workers
    sv.parameters.random_seed = seed
    if log:
        sv.parameters.log_search_progress = True
    st = sv.Solve(m)
    name = sv.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        D = [[set(g for g in range(p) if sv.Value(x[i, j, g])) for j in range(norb)] for i in range(norb)]
        return name, D
    return name, None

def validate():
    for p, (k, lam, mu) in [(3, (4, 1, 2)), (5, (12, 5, 6))]:
        v = p * p; k = (v - 1) // 2; lam = (v - 5) // 4; mu = (v - 1) // 4
        Dtrue = cc.paley_pp_blocks(p)
        R = cc.orbit_matrix(Dtrue, p)
        st, D = solve_lift(p, p, k, lam, mu, R, time_s=60)
        ok = False
        if D is not None:
            A = cc.build_adjacency(D, p, p)
            ok = cc.srg_identity_direct(A, k, lam, mu)
        print(f"  P({v})=srg({v},{k},{lam},{mu}) via Z_{p} blocks: status={st}  found_valid_srg={ok}")

if __name__ == "__main__":
    import sys
    if "validate" in sys.argv or len(sys.argv) == 1:
        print("=== VALIDATION on Paley graphs ===")
        validate()
