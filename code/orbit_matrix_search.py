"""
Enumerate the 9x9 orbit matrices R for a fixed-point-free order-37 automorphism of
srg(333,166,82,83).  R[i][j] = |D[i][j]| (neighbours of an O_i-vertex in O_j).
Constraints (derived + validated framework):
  - R symmetric, 9x9, entries: diagonal in {0,2,...,36} (even, no loops, symmetric circ),
    off-diagonal in [0,37].
  - row sums = k = 166.
  - quotient SRG identity  R^2 + R - 83 I = 3071 J   (mu-lam=1, k-mu=83, mu*p=83*37=3071):
        i=j : sum_k R[i][k]^2          = 3154 - R[i][i]
        i!=j: sum_k R[i][k]*R[j][k]     = 3071 - R[i][j]
  (Equivalently R has eigenvalues 166 (x1) and (-1+-sqrt333)/2 (x4 each); trace R = 162.)
Symmetry breaking: orbits are interchangeable -> impose R[0][0]<=R[1][1]<=...<=R[8][8]
and lexical tie-break is skipped (we just want existence + a manageable list).
"""
from ortools.sat.python import cp_model
import numpy as np, sys, time

N = 9
K = 166
DIAG_RHS = 3154   # sum_k R[i][k]^2 = 3154 - R[i][i]
OFF_RHS = 3071    # sum_k R[i][k]*R[j][k] = 3071 - R[i][j]  (i!=j)

def build_model(break_sym=True):
    m = cp_model.CpModel()
    R = {}
    for i in range(N):
        for j in range(i, N):
            if i == j:
                # Cauchy-Schwarz on the 8 off-diagonal entries of row i (fixed sum 166-d_i,
                # fixed sum-of-squares 3154-d_i-d_i^2, each in [0,37]) forces d_i in {10,..,26} even.
                v = m.NewIntVar(10, 26, f"R{i}{j}")
                h = m.NewIntVar(5, 13, f"h{i}")
                m.Add(v == 2*h)
            else:
                v = m.NewIntVar(0, 37, f"R{i}{j}")
            R[i, j] = v
            R[j, i] = v
    # row sums
    for i in range(N):
        m.Add(sum(R[i, k] for k in range(N)) == K)
    # global derived constraints
    m.Add(sum(R[i, i] for i in range(N)) == 162)                       # trace
    m.Add(sum(R[i, j] for i in range(N) for j in range(i+1, N)) == 666) # off-diag sum
    # quadratic identities: product var for each needed pair (i,k),(j,k)
    prodcache = {}
    def mul(a, b):
        ka = (min(a, b), max(a, b))
        if ka not in prodcache:
            pv = m.NewIntVar(0, 37*37, f"m_{ka[0]}_{ka[1]}")
            m.AddMultiplicationEquality(pv, [R[ka[0][0], ka[0][1]], R[ka[1][0], ka[1][1]]])
            prodcache[ka] = pv
        return prodcache[ka]
    for i in range(N):
        for j in range(i, N):
            terms = []
            for k in range(N):
                terms.append(mul((i, k), (j, k)))
            if i == j:
                m.Add(sum(terms) == DIAG_RHS - R[i, i])
            else:
                m.Add(sum(terms) == OFF_RHS - R[i, j])
    if break_sym:
        for i in range(N-1):
            m.Add(R[i, i] <= R[i+1, i+1])
    return m, R

def enumerate_R(max_sol=200, time_s=300):
    m, R = build_model()
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_s
    solver.parameters.num_search_workers = 8
    class Cb(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__(); self.sols=[]; self.n=0
        def on_solution_callback(self):
            self.n += 1
            mat = np.array([[self.Value(R[i, j]) for j in range(N)] for i in range(N)], dtype=np.int64)
            self.sols.append(mat)
            if self.n >= max_sol:
                self.StopSearch()
    cb = Cb()
    solver.parameters.enumerate_all_solutions = True
    status = solver.Solve(m, cb)
    return cb.sols, solver.StatusName(status), solver.WallTime()

def verify_R(R):
    R = np.asarray(R)
    cond = np.array_equal(R@R + R - 83*np.eye(N, dtype=np.int64), 3071*np.ones((N, N), dtype=np.int64))
    eig = np.linalg.eigvalsh(R.astype(float))
    return cond, sorted(round(e, 4) for e in eig)

if __name__ == "__main__":
    t0 = time.time()
    sols, status, wt = enumerate_R(max_sol=int(sys.argv[1]) if len(sys.argv) > 1 else 200,
                                   time_s=int(sys.argv[2]) if len(sys.argv) > 2 else 300)
    print(f"status={status}  #orbit-matrices found={len(sols)}  walltime={wt:.1f}s")
    r = ((-1+333**0.5)/2, (-1-333**0.5)/2)
    print(f"  (SRG eigenvalues r,s = {r[0]:.4f}, {r[1]:.4f}; expect R-spectrum = 166, r x4, s x4)")
    for idx, R in enumerate(sols[:5]):
        cond, eig = verify_R(R)
        print(f"  R#{idx}: identity-holds={cond}  trace={int(np.trace(R))}  diag={list(np.diag(R))}")
        print(f"        spectrum={eig}")
    import os
    OUT = os.path.dirname(os.path.abspath(__file__))
    if sols:
        np.save(os.path.join(OUT, "orbit_matrices.npy"), np.array(sols))
        print(f"  saved {len(sols)} orbit matrices")
    if status == "INFEASIBLE":
        print("  *** PROVEN: NO orbit matrix exists => srg(333,166,82,83) has NO fixed-point-free")
        print("      order-37 automorphism. (rigorous nonexistence) ***")
    elif status in ("OPTIMAL",) and not sols:
        print("  search complete, none found (treat as infeasible only if status INFEASIBLE).")
    elif not sols:
        print(f"  status={status}: NO solution found but NOT proven infeasible (need longer / better model).")
