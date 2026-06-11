"""
CP-SAT with the Galois trace-condition pinning (sc_trace_condition.py) added to the
hostile-validated selfcomp_search.Model for the open self-complementary classes 221/422.

THEOREM-backed added constraints (necessary for any exact solution; verified separately):
  branch QR : D_77 = QR(37) exactly, and per <+-6>-orbit (QR-rep d):
              [d in S] - [6d in S] + [d in T] - [6d in T] = -1   (u,v)=(3,6)
  branch QNR: D_77 = QNR(37) exactly, and per-orbit sum = +1     (u,v)=(6,3)
A solution to the class must satisfy exactly one branch; INFEASIBLE on both = class
instance CLOSED (modulo CP-SAT trust, same standard as the validated class-109 closure).

Usage: sc_pinned_cpsat.py <cls> <perm_csv> <branch QR|QNR> <time_s>
"""
import sys, json
import numpy as np
sys.path.insert(0, "/Users/maceocardinalekwik/git/math/release/h668-z37/code")
from selfcomp_search import Model, verify_srg, verify_selfcomp
from ortools.sat.python import cp_model

P = 37
QR = sorted({(k*k) % P for k in range(1, P)})
QNR = [d for d in range(1, P) if d not in QR]

def main():
    cls_idx = int(sys.argv[1])
    perm = tuple(int(t) for t in sys.argv[2].split(","))
    branch = sys.argv[3]
    time_s = float(sys.argv[4])
    mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
            for ln in open("/Users/maceocardinalekwik/git/math/release/h668-z37/code/orbit_all_classes.txt")]
    R = [[int(x) for x in row] for row in mats[cls_idx]]
    M = Model(P, 9, 166, 82, 83, R, perm, 6)
    assert M.ok, "complement-size consistency failed"
    # identify diagonal reps
    diag = [rep for rep in M.orbits if rep[0] == rep[1]]
    fixed = [rep for rep in diag if perm[rep[0]] == rep[0]]
    cyc   = [rep for rep in diag if perm[rep[0]] != rep[0]]
    assert len(fixed) == 1 and len(cyc) == 2, (diag, fixed, cyc)
    (uorb,) = fixed
    Sorb, Torb = cyc
    print(f"cls={cls_idx} perm={perm} branch={branch}: U-block={uorb} "
          f"S={Sorb}(card {R[Sorb[0]][Sorb[0]]}) T={Torb}(card {R[Torb[0]][Torb[0]]})", flush=True)
    m = M.m
    Uset = QR if branch == "QR" else QNR
    eps = -1 if branch == "QR" else +1
    for d in range(1, P):
        m.Add(M.X[(uorb, d)] == (1 if d in Uset else 0))
    # per <+-6>-orbit constraints, QR-rep d
    seen = [False]*P
    norb_added = 0
    for x in range(1, P):
        if seen[x]: continue
        ob = [x, (6*x) % P, P-x, (P - (6*x) % P) % P]
        for z in ob: seen[z] = True
        d = x if x in QR else (6*x) % P          # QR half rep
        d6 = (6*d) % P                            # QNR half rep
        m.Add(M.X[(Sorb, d)] - M.X[(Sorb, d6)] + M.X[(Torb, d)] - M.X[(Torb, d6)] == eps)
        norb_added += 1
    assert norb_added == 9
    sv = cp_model.CpSolver()
    sv.parameters.max_time_in_seconds = time_s
    sv.parameters.num_search_workers = 8
    st = sv.Solve(m)
    name = sv.StatusName(st)
    print(f"RESULT cls={cls_idx} perm={','.join(map(str,perm))} branch={branch}: {name} "
          f"(wall {sv.WallTime():.1f}s)", flush=True)
    if name in ("OPTIMAL", "FEASIBLE"):
        D = M.extract(sv.Value)
        ok1 = verify_srg(D, P, 9, 166, 82, 83)
        ok2 = verify_selfcomp(D, P, 9, perm, 6)
        print(f"*** SOLUTION FOUND *** verify_srg={ok1} verify_selfcomp={ok2}")
        json.dump([[sorted(D[i][j]) for j in range(9)] for i in range(9)],
                  open(f"/tmp/SOLUTION_cls{cls_idx}_{branch}.json", "w"))

if __name__ == "__main__":
    main()
