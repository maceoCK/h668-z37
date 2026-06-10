"""
Modal: order-3 cyclotomic closure sweep for H(668) / srg(333,166,82,83).

Runs the VALIDATED class-level pinned CP-SAT model (cyclo_cpsat.solve, validated on P(9)/P(25)/
Petersen true-positives and brute-force-confirmed INFEASIBLE anchors) over the 28 order-3
coset-compatible orbit-matrix representatives (complete enumeration orbit_all_sym.txt filtered by
diag = 0 mod 6, offdiag mod 3 in {0,1}), with the Galois trace pinning a=4 (WLOG: a QNR decimation
maps any a=5 solution to a=4 on the same skeleton).

Outcomes per skeleton: INFEASIBLE (closes it), FEASIBLE (verified in-container -> H(668)!), UNKNOWN.
All 28 INFEASIBLE  =>  no Z37 order-3-cyclotomic srg(333) exists at all (order 3 CLOSED).

Run:  modal run --detach modal_order3_cpsat.py --time-s 9000
"""
import modal

app = modal.App("h668-order3-cpsat")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("ortools", "numpy")
         .add_local_file("cyclo_cpsat.py", "/root/cyclo_cpsat.py")
         .add_local_file("order3_skeletons.txt", "/root/order3_skeletons.txt"))

@app.function(image=image, cpu=8.0, timeout=14400, retries=2)
def run_one(idx: int, time_s: int):
    import sys, time
    sys.path.insert(0, "/root")
    import numpy as np
    from cyclo_cpsat import solve
    Rs = [[int(t) for t in ln.split()] for ln in open("/root/order3_skeletons.txt")]
    R = [Rs[idx][9 * i: 9 * i + 9] for i in range(9)]
    t0 = time.time()
    st, sol = solve(37, 9, 166, 82, 83, 10, R, a_pin=4, time_s=time_s, workers=8)
    el = time.time() - t0
    out = {"idx": idx, "status": st, "sec": round(el, 1), "verified": None, "sol": None}
    if sol is not None:
        # independent in-container verification on the full 333-vertex adjacency matrix
        p, norb, k, lam, mu = 37, 9, 166, 82, 83
        V = p * norb
        A = np.zeros((V, V), dtype=np.int64)
        for i in range(norb):
            for j in range(norb):
                for a in range(p):
                    for d in sol[i][j]:
                        A[i * p + a][j * p + (a + d) % p] = 1
        lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(V, dtype=np.int64)
        ok = bool(np.array_equal(A, A.T) and np.array_equal(lhs, mu * np.ones((V, V), dtype=np.int64)))
        out["verified"] = ok
        out["sol"] = [[sorted(sol[i][j]) for j in range(9)] for i in range(9)]
    return out

@app.local_entrypoint()
def main(time_s: int = 9000):
    import json
    n = sum(1 for _ in open("order3_skeletons.txt"))
    print(f"order-3 closure sweep: {n} skeletons x (pinned a=4, class-level CP-SAT, {time_s}s, 8 workers)")
    results = []
    for r in run_one.starmap([(i, time_s) for i in range(n)]):
        results.append(r)
        tag = ""
        if r["status"] in ("OPTIMAL", "FEASIBLE"):
            tag = "  *** SAT -> H(668) CANDIDATE, verified=%s ***" % r["verified"]
        print(f"  [{len(results)}/{n}] skeleton {r['idx']}: {r['status']} ({r['sec']}s){tag}", flush=True)
    with open("order3_cpsat_results.json", "w") as f:
        json.dump(results, f, indent=1)
    inf = sum(1 for r in results if r["status"] == "INFEASIBLE")
    unk = sum(1 for r in results if r["status"] not in ("INFEASIBLE", "OPTIMAL", "FEASIBLE"))
    sat = [r for r in results if r["status"] in ("OPTIMAL", "FEASIBLE")]
    print(f"\n=== SWEEP DONE === INFEASIBLE {inf}/{n}, UNKNOWN {unk}, SAT {len(sat)}")
    if sat:
        print("*** SAT skeletons:", [(r["idx"], r["verified"]) for r in sat])
    elif inf == n:
        print("ORDER 3 CLOSED: no Z37 order-3-cyclotomic srg(333,166,82,83) exists (a=4 WLOG).")
