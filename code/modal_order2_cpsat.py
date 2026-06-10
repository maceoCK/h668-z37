"""
Modal: order-2 (negation-invariant) cyclotomic closure attempt for H(668) / srg(333,166,82,83).

Same validated class-level pinned CP-SAT model as the order-3 sweep (cyclo_cpsat.solve), with
mult=36 (H={+-1}), pinning a=4 (WLOG by QNR decimation; each QR pair {g,-g} in exactly 3 diagonal
blocks, each QNR pair in 6).  Sweep target: the 625 true isomorphism classes of admissible orbit
matrices (orbit_all_classes.txt, exact lexmin canonicalization of the complete 2901 enumeration).

Calibrate first:   modal run --detach modal_order2_cpsat.py --mode probe --time-s 9000   (6 classes)
Full sweep (only if probes return INFEASIBLE):  --mode full
"""
import modal

app = modal.App("h668-order2-cpsat")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("ortools", "numpy")
         .add_local_file("cyclo_cpsat.py", "/root/cyclo_cpsat.py")
         .add_local_file("orbit_all_classes.txt", "/root/orbit_all_classes.txt"))

@app.function(image=image, cpu=8.0, timeout=14400, retries=2)
def run_one(idx: int, time_s: int):
    import sys, time
    sys.path.insert(0, "/root")
    import numpy as np
    from cyclo_cpsat import solve
    Rs = [[int(t) for t in ln.split()] for ln in open("/root/orbit_all_classes.txt")]
    R = [Rs[idx][9 * i: 9 * i + 9] for i in range(9)]
    t0 = time.time()
    st, sol = solve(37, 9, 166, 82, 83, 36, R, a_pin=4, time_s=time_s, workers=8)
    el = time.time() - t0
    out = {"idx": idx, "status": st, "sec": round(el, 1), "verified": None, "sol": None}
    if sol is not None:
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
def main(mode: str = "probe", time_s: int = 9000):
    import json
    n = sum(1 for _ in open("orbit_all_classes.txt"))
    idxs = [0, n // 5, 2 * n // 5, 3 * n // 5, 4 * n // 5, n - 1] if mode == "probe" else list(range(n))
    print(f"order-2 CP-SAT {mode}: {len(idxs)} of {n} classes (pinned a=4, {time_s}s, 8 workers)")
    results = []
    for r in run_one.starmap([(i, time_s) for i in idxs]):
        results.append(r)
        tag = ""
        if r["status"] in ("OPTIMAL", "FEASIBLE"):
            tag = "  *** SAT -> H(668) CANDIDATE, verified=%s ***" % r["verified"]
        print(f"  [{len(results)}/{len(idxs)}] class {r['idx']}: {r['status']} ({r['sec']}s){tag}", flush=True)
    with open(f"order2_cpsat_{mode}.json", "w") as f:
        json.dump(results, f, indent=1)
    inf = sum(1 for r in results if r["status"] == "INFEASIBLE")
    print(f"\n=== DONE === INFEASIBLE {inf}/{len(idxs)}")
