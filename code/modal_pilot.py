"""
Modal: Experiment 1 pilot v2 — incremental cube-and-conquer with conflict-budgeted warming.

Per shard: load the order-3 skeleton-0 CNF once (94MB, 404k vars, anchor-validated emitter),
ROUND A: skim `count` cubes with a conflict budget (build the shared learned-clause DB cheaply,
recording UNSAT-within-budget rate and times), ROUND B: re-solve the first `measure` cubes
unbudgeted with a wall cap to measure the WARM full-refutation distribution (the go/no-go number).

Run:  modal run --detach modal_pilot.py
"""
import modal

app = modal.App("h668-pilot-v2")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("python-sat")
         .add_local_file("order3_sk0.cnf", "/root/order3_sk0.cnf")
         .add_local_file("order3_sk0.map", "/root/order3_sk0.map")
         .add_local_file("cubes_sk0.txt", "/root/cubes_sk0.txt"))

@app.function(image=image, cpu=2.0, timeout=10800, retries=2)
def shard(start: int, count: int, measure: int, budget: int, cap_s: int):
    import json, time, threading
    from pysat.formula import CNF
    from pysat.solvers import Cadical195
    mp = json.load(open("/root/order3_sk0.map"))
    nc = len(mp["cosets"]); norb = 9
    dv = [[mp["y"][f"{i},{i},{c}"] for c in range(nc)] for i in range(norb)]
    cubes = []
    with open("/root/cubes_sk0.txt") as f:
        for ln, line in enumerate(f):
            if ln < start: continue
            if ln >= start + count: break
            cubes.append([int(t) for t in line.split()])
    cnf = CNF(from_file="/root/order3_sk0.cnf")
    s = Cadical195(bootstrap_with=cnf.clauses)
    def assume(masks):
        out = []
        for i in range(norb):
            for c in range(nc):
                v = dv[i][c]
                out.append(v if (masks[i] >> c) & 1 else -v)
        return out
    # ---- round A: budgeted skim
    a_stats = {"unsat": 0, "undet": 0, "sat": 0, "times": []}
    for masks in cubes:
        s.conf_budget(budget)
        t0 = time.time()
        r = s.solve_limited(assumptions=assume(masks))
        a_stats["times"].append(round(time.time() - t0, 3))
        if r is True: a_stats["sat"] += 1
        elif r is False: a_stats["unsat"] += 1
        else: a_stats["undet"] += 1
    # ---- round B: warm unbudgeted with wall cap
    b_times, b_status = [], []
    for masks in cubes[:measure]:
        timer = threading.Timer(cap_s, s.interrupt)
        timer.start()
        t0 = time.time()
        r = s.solve_limited(assumptions=assume(masks), expect_interrupt=True)
        timer.cancel()
        s.clear_interrupt()
        b_times.append(round(time.time() - t0, 3))
        b_status.append("SAT" if r is True else ("UNSAT" if r is False else "TIMEOUT"))
    w = sorted(b_times)
    return {"start": start,
            "roundA": {"unsat": a_stats["unsat"], "undet": a_stats["undet"], "sat": a_stats["sat"],
                       "mean_s": round(sum(a_stats["times"]) / max(1, len(a_stats["times"])), 3)},
            "roundB": {"status": b_status, "times": b_times,
                       "median": w[len(w)//2] if w else None,
                       "p90": w[int(len(w)*0.9)] if w else None}}

@app.local_entrypoint()
def main(shards: int = 8, count: int = 2000, measure: int = 60, budget: int = 50000, cap_s: int = 300):
    import json
    jobs = [(k * 18000, count, measure, budget, cap_s) for k in range(shards)]
    print(f"pilot v2: {shards} shards x {count} cubes (budget {budget} conflicts), measure {measure} warm @ {cap_s}s cap")
    res = []
    for r in shard.starmap(jobs):
        res.append(r)
        a, b = r["roundA"], r["roundB"]
        to = sum(1 for x in b["status"] if x == "TIMEOUT")
        sat = sum(1 for x in b["status"] if x == "SAT") + a["sat"]
        print(f"  shard@{r['start']}: A unsat={a['unsat']} undet={a['undet']} mean={a['mean_s']}s | "
              f"B median={b['median']}s p90={b['p90']}s timeouts={to}/{len(b['status'])}"
              + ("  *** SAT ***" if sat else ""), flush=True)
    json.dump(res, open("pilot_v2_results.json", "w"), indent=1)
    print("=== PILOT V2 DONE ===")
