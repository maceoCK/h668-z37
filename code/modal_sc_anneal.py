"""
Modal: simulated-annealing SAT-side lane for the self-complementary Z37 srg(333) search.
6 templates (classes 221, 422 x complementing perms; emitted by the hostile-validated
sc_emit_template.py) x seeds.  Annealer validated end-to-end on Paley(13) (finds QR instantly)
and drift-checked at production shape.  A SOLUTION (objective 0) is verified in-container on the
full 333x333 adjacency matrix => H(668).

Run:  modal run --detach modal_sc_anneal.py --time-s 7000 --seeds 10
"""
import modal

app = modal.App("h668-sc-anneal")
image = (modal.Image.debian_slim(python_version="3.11")
         .apt_install("g++")
         .pip_install("numpy")
         .add_local_file("sc_anneal.cpp", "/root/sc_anneal.cpp")
         .add_local_file("sc_tpl_221_0.txt", "/root/sc_tpl_221_0.txt")
         .add_local_file("sc_tpl_221_1.txt", "/root/sc_tpl_221_1.txt")
         .add_local_file("sc_tpl_422_2.txt", "/root/sc_tpl_422_2.txt")
         .add_local_file("sc_tpl_422_3.txt", "/root/sc_tpl_422_3.txt")
         .add_local_file("sc_tpl_422_4.txt", "/root/sc_tpl_422_4.txt")
         .add_local_file("sc_tpl_422_5.txt", "/root/sc_tpl_422_5.txt"))

@app.function(image=image, cpu=1.0, timeout=8400, retries=2)
def run_one(tpl: str, seed: int, time_s: int):
    import subprocess, numpy as np
    subprocess.run(["g++", "-O3", "-o", "/root/sc", "/root/sc_anneal.cpp"], check=True)
    pr = subprocess.run(["/root/sc", f"/root/{tpl}", str(seed), str(time_s)],
                        capture_output=True, text=True, timeout=time_s + 600)
    out = {"tpl": tpl, "seed": seed, "best": None, "sol": None, "verified": None}
    best = None
    sol_lines = []
    in_sol = False
    for ln in pr.stdout.splitlines():
        if ln.startswith("BEST") or ln.startswith("RESULT"):
            for t in ln.split():
                if t.startswith("obj=") or t.startswith("best="):
                    best = int(t.split("=")[1])
        if ln.startswith("SOLUTION"): in_sol = True
        elif in_sol and ln.startswith("D "): sol_lines.append(ln)
    out["best"] = best if not sol_lines else 0
    if sol_lines:
        p, norb, k, lam, mu = 37, 9, 166, 82, 83
        D = [[set() for _ in range(norb)] for _ in range(norb)]
        for ln in sol_lines:
            head, _, elems = ln.partition(":")
            _, i, j = head.split()
            D[int(i)][int(j)] = set(int(x) for x in elems.split())
        V = p * norb
        A = np.zeros((V, V), dtype=np.int64)
        for i in range(norb):
            for j in range(norb):
                for a in range(p):
                    for d in D[i][j]:
                        A[i * p + a][j * p + (a + d) % p] = 1
        lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(V, dtype=np.int64)
        ok = bool(np.array_equal(A, A.T) and np.array_equal(lhs, mu * np.ones((V, V), dtype=np.int64)))
        out["verified"] = ok
        out["sol"] = [[sorted(D[i][j]) for j in range(9)] for i in range(9)]
    return out

@app.local_entrypoint()
def main(time_s: int = 7000, seeds: int = 10):
    import json
    tpls = ["sc_tpl_221_0.txt", "sc_tpl_221_1.txt", "sc_tpl_422_2.txt",
            "sc_tpl_422_3.txt", "sc_tpl_422_4.txt", "sc_tpl_422_5.txt"]
    jobs = [(t, s, time_s) for t in tpls for s in range(1, seeds + 1)]
    print(f"SA lane: {len(jobs)} jobs ({len(tpls)} templates x {seeds} seeds, {time_s}s, 1 cpu)")
    results = []
    for r in run_one.starmap(jobs):
        results.append(r)
        tag = f"  *** SOLUTION verified={r['verified']} -> H(668)!! ***" if r["sol"] else ""
        print(f"  [{len(results)}/{len(jobs)}] {r['tpl']} seed {r['seed']}: best={r['best']}{tag}", flush=True)
    with open("sc_anneal_results.json", "w") as f:
        json.dump(results, f, indent=1)
    hits = [r for r in results if r["sol"]]
    bests = sorted(r["best"] for r in results if r["best"] is not None)
    print(f"\n=== DONE === solutions {len(hits)}; best objectives {bests[:10]}")
