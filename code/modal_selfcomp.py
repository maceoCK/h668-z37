"""
Modal: self-complementary Z37-symmetric srg(333,166,82,83) search (fleet Action 1).

For each of the 3 self-complementary orbit-matrix classes (109, 221, 422 of the 625) x each
complementing permutation (cycle type (4,4,1)), search for blocks D[i][j] with
sigma(i,g) = (pi(i), 6g) a complementing isomorphism (b=0 WLOG, twist analysis verified).
Model: selfcomp_search.Model (validated: Paley(13)/(17) anchors + exact brute-force cross-checks).

A verified SAT = a self-complementary conference graph srg(333,166,82,83) => H(668) EXISTS.

Run:  modal run --detach modal_selfcomp.py --time-s 7200
"""
import modal

app = modal.App("h668-selfcomp")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("ortools", "numpy")
         .add_local_file("selfcomp_search.py", "/root/selfcomp_search.py")
         .add_local_file("conf_core.py", "/root/conf_core.py")
         .add_local_file("orbit_all_classes.txt", "/root/orbit_all_classes.txt"))

JOBS = []  # filled by entrypoint

@app.function(image=image, cpu=8.0, timeout=10800, retries=3)
def run_one(cls_idx: int, perm: list, time_s: int):
    import sys, time
    sys.path.insert(0, "/root")
    import numpy as np
    from selfcomp_search import Model, verify_srg, verify_selfcomp
    mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
            for ln in open("/root/orbit_all_classes.txt")]
    R = [[int(x) for x in row] for row in mats[cls_idx]]
    t0 = time.time()
    M = Model(37, 9, 166, 82, 83, R, tuple(perm), 6)
    out = {"cls": cls_idx, "perm": perm, "status": "INCOMPATIBLE", "sec": 0,
           "verified": None, "sol": None}
    if M.ok:
        st, sols = M.solve(time_s=time_s, workers=8)
        out["status"] = st
        if sols:
            D = sols[0]
            ok = verify_srg(D, 37, 9, 166, 82, 83) and verify_selfcomp(D, 37, 9, tuple(perm), 6)
            out["verified"] = bool(ok)
            out["sol"] = [[sorted(D[i][j]) for j in range(9)] for i in range(9)]
    out["sec"] = round(time.time() - t0, 1)
    return out

@app.local_entrypoint()
def main(time_s: int = 7200):
    import json, itertools
    import numpy as np
    mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
            for ln in open("orbit_all_classes.txt")]
    jobs = []
    for k in (109, 221, 422):
        R = mats[k]
        Rc = np.where(np.eye(9, dtype=bool), 36 - R, 37 - R)
        for q in itertools.permutations(range(9)):
            if np.array_equal(Rc, R[np.ix_(q, q)]):
                jobs.append((k, list(q), time_s))
    print(f"self-complementary search: {len(jobs)} jobs (class, perm), {time_s}s x 8 workers each")
    results = []
    for r in run_one.starmap(jobs):
        results.append(r)
        tag = ""
        if r["status"] in ("OPTIMAL", "FEASIBLE"):
            tag = f"  *** SAT verified={r['verified']} -> H(668)!! ***"
        print(f"  [{len(results)}/{len(jobs)}] class {r['cls']} perm {r['perm']}: {r['status']} ({r['sec']}s){tag}", flush=True)
    with open("selfcomp_results.json", "w") as f:
        json.dump(results, f, indent=1)
    sat = [r for r in results if r["status"] in ("OPTIMAL", "FEASIBLE")]
    inf = sum(1 for r in results if r["status"] == "INFEASIBLE")
    print(f"\n=== DONE === SAT {len(sat)}, INFEASIBLE {inf}/{len(jobs)}, rest UNKNOWN")
    if sat:
        print("*** SELF-COMPLEMENTARY srg(333) FOUND -> CONSTRUCT H(668) ***")
    elif inf == len(jobs):
        print("THEOREM: no srg(333,166,82,83) admits a fixed-point-free Z37 normalized by a complementing map.")
