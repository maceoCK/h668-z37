"""
Modal: self-complementary 221/422 search in the (C1)/(C2)-pinned corner.

The verified theorem-grade necessary condition (c1c2_verification.md, SOUND):
  (C1) the fixed-orbit diagonal D[7][7] is exactly QR(37) or exactly QNR(37);
  (C2) per <±6>-orbit with QR-rep d (d6 = 6d):
       x[s0][s0][d] - x[s0][s0][d6] + x[t0][t0][d] - x[t0][t0][d6] = -1 (QR branch) / +1 (QNR),
       s0 = 0, t0 = 2 (class 221) / 3 (class 422).
Any solution MUST satisfy one branch; INFEASIBLE on both branches for all perms of a class
closes the class. 12 jobs: (221 x 2 perms + 422 x 4 perms) x 2 branches.

Run:  modal run --detach modal_selfcomp_pinned.py --time-s 9000
"""
import modal

app = modal.App("h668-selfcomp-pinned")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("ortools", "numpy")
         .add_local_file("selfcomp_search.py", "/root/selfcomp_search.py")
         .add_local_file("conf_core.py", "/root/conf_core.py")
         .add_local_file("orbit_all_classes.txt", "/root/orbit_all_classes.txt"))

@app.function(image=image, cpu=8.0, timeout=14400, retries=2)
def run_one(cls_idx: int, perm: list, branch: str, time_s: int):
    import sys, time
    sys.path.insert(0, "/root")
    import numpy as np
    from selfcomp_search import Model, verify_srg, verify_selfcomp
    mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
            for ln in open("/root/orbit_all_classes.txt")]
    R = [[int(x) for x in row] for row in mats[cls_idx]]
    M = Model(37, 9, 166, 82, 83, R, tuple(perm), 6)
    out = {"cls": cls_idx, "perm": perm, "branch": branch, "status": "INCOMPATIBLE",
           "sec": 0, "verified": None, "sol": None}
    if not M.ok:
        return out
    p = 37
    QR = {(x * x) % p for x in range(1, p)}
    s0, t0 = 0, (2 if cls_idx == 221 else 3)
    # (C1): pin fixed-orbit diagonal
    for d in range(1, p):
        want = 1 if ((d in QR) == (branch == "QR")) else 0
        M.m.Add(M.X[((7, 7), d)] == want)
    # (C2): per <±6>-orbit coupling
    seen = set()
    eps = -1 if branch == "QR" else 1
    for x in range(1, p):
        if x in seen: continue
        orb = {x, (6 * x) % p, (-x) % p, (-6 * x) % p}
        seen |= orb
        d = x if x in QR else (6 * x) % p
        assert d in QR
        d6 = (6 * d) % p
        M.m.Add(M.X[((s0, s0), d)] - M.X[((s0, s0), d6)]
                + M.X[((t0, t0), d)] - M.X[((t0, t0), d6)] == eps)
    t0_ = time.time()
    st, sols = M.solve(time_s=time_s, workers=8)
    out["status"] = st
    out["sec"] = round(time.time() - t0_, 1)
    if sols:
        D = sols[0]
        ok = verify_srg(D, 37, 9, 166, 82, 83) and verify_selfcomp(D, 37, 9, tuple(perm), 6)
        out["verified"] = bool(ok)
        out["sol"] = [[sorted(D[i][j]) for j in range(9)] for i in range(9)]
    return out

@app.local_entrypoint()
def main(time_s: int = 9000):
    import json, itertools
    mats = []
    for ln in open("orbit_all_classes.txt"):
        v = [int(t) for t in ln.split()]
        if len(v) == 81: mats.append([v[9 * i: 9 * i + 9] for i in range(9)])
    jobs = []
    for k in (221, 422):
        R = mats[k]
        Rc = [[(36 - R[i][j]) if i == j else (37 - R[i][j]) for j in range(9)] for i in range(9)]
        for q in itertools.permutations(range(9)):
            if all(Rc[i][j] == R[q[i]][q[j]] for i in range(9) for j in range(9)):
                for br in ("QR", "QNR"):
                    jobs.append((k, list(q), br, time_s))
    print(f"pinned-corner search: {len(jobs)} jobs ((class,perm) x branch), {time_s}s x 8w")
    results = []
    for r in run_one.starmap(jobs):
        results.append(r)
        tag = f"  *** SAT verified={r['verified']} -> H(668)!! ***" if r["status"] in ("OPTIMAL", "FEASIBLE") else ""
        print(f"  [{len(results)}/{len(jobs)}] cls {r['cls']} {r['branch']} perm {r['perm']}: {r['status']} ({r['sec']}s){tag}", flush=True)
    json.dump(results, open("selfcomp_pinned_results.json", "w"), indent=1)
    inf = sum(1 for r in results if r["status"] == "INFEASIBLE")
    sat = [r for r in results if r["status"] in ("OPTIMAL", "FEASIBLE")]
    print(f"\n=== DONE === INFEASIBLE {inf}/{len(jobs)}, SAT {len(sat)}")
    if inf == len(jobs):
        print("THEOREM: classes 221/422 closed => NO self-complementary Z37-symmetric srg(333) exists at all.")
