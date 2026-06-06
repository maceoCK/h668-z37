"""
Cheap COMPLETE closure for finer multipliers: an order-e-invariant srg(333) forces its 9x9 orbit
matrix R to have coset-compatible entries -- off-diagonal R[i][j] mod e in {0,1} (0 = whether the
identity element 0 is in the block), diagonal R[i][i] = 0 mod e (e even) or mod 2e (e odd, since
cosets pair up under negation). If NO orbit matrix satisfies the SRG identity R^2+R-83I=3071J AND
these mod-e constraints, then no order-e-invariant Z_37 srg(333,166,82,83) exists (complete; sound).
"""
import orbit_matrix_search as oms
from ortools.sat.python import cp_model
import sys

def check(e, time_s=300):
    m, R = oms.build_model()
    offvals = [v for v in range(0, 38) if v % e in (0, 1)]
    diag_mod = e if e % 2 == 0 else 2 * e
    diagvals = [v for v in range(10, 27, 2) if v % diag_mod == 0]
    if not diagvals:
        return "INFEASIBLE", "(no diagonal value in {10..26} is 0 mod %d)" % diag_mod
    for i in range(9):
        for j in range(9):
            if i == j:
                m.AddAllowedAssignments([R[i, j]], [[v] for v in diagvals])
            elif i < j:
                m.AddAllowedAssignments([R[i, j]], [[v] for v in offvals])
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = time_s; sv.parameters.num_search_workers = 8
    st = sv.Solve(m)
    return sv.StatusName(st), f"diagvals={diagvals}"

if __name__ == "__main__":
    ts = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    # finer multiplier orders (2 is unrestricted -> can't close this way); also re-confirm coarse
    for e in (9, 6, 4, 3, 12, 18, 36, 2):
        st, info = check(e, ts)
        verdict = "CLOSED (no order-%d-invariant srg)" % e if st == "INFEASIBLE" else \
                  ("compatible orbit matrix EXISTS" if st in ("OPTIMAL", "FEASIBLE") else "undecided")
        print(f"  order {e:2d}: {st:11s} -> {verdict}   {info}")
