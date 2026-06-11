"""
Experiment 1 driver: incremental cube-and-conquer over canonical diagonal cubes.

Loads the emitted CNF once into CaDiCaL (PySAT Cadical195), then solves each cube as an
ASSUMPTION set (all 108 diagonal-bit literals fixed), retaining the learned-clause database
across cubes (the mechanism under test: cross-row conflicts learned once should prune later
cubes). Records per-cube wall time and status; reports rolling warm median / p90 / trend.

Usage: pilot_drive.py <base> <cubefile> <start> <count> [--budget-conflicts N]
Cube file: lines of 9 coset-bitmask ints (cyclo_dfs --dump-diag order; same R, same pinning).
GO signal (fleet): warm median <= 1s and falling.  KILL: > 10s flat.
"""
import sys, json, time
from pysat.formula import CNF
from pysat.solvers import Cadical195

def main():
    base, cubefile = sys.argv[1], sys.argv[2]
    start, count = int(sys.argv[3]), int(sys.argv[4])
    mp = json.load(open(base + ".map"))
    nc = len(mp["cosets"])
    norb = 9
    diagvars = [[mp["y"][f"{i},{i},{c}"] for c in range(nc)] for i in range(norb)]
    cubes = []
    with open(cubefile) as f:
        for ln, line in enumerate(f):
            if ln < start: continue
            if ln >= start + count: break
            cubes.append([int(t) for t in line.split()])
    print(f"loading {base}.cnf ...", flush=True)
    t0 = time.time()
    cnf = CNF(from_file=base + ".cnf")
    s = Cadical195(bootstrap_with=cnf.clauses)
    print(f"loaded in {time.time()-t0:.1f}s; {len(cubes)} cubes [{start}..{start+count})", flush=True)
    times, statuses = [], []
    for idx, masks in enumerate(cubes):
        assume = []
        for i in range(norb):
            for c in range(nc):
                v = diagvars[i][c]
                assume.append(v if (masks[i] >> c) & 1 else -v)
        t1 = time.time()
        res = s.solve(assumptions=assume)
        dt = time.time() - t1
        times.append(dt); statuses.append(res)
        if res:
            print(f"CUBE {start+idx} SAT !!! dumping model", flush=True)
            model = s.get_model()
            json.dump({"cube": start + idx, "masks": masks, "model": model},
                      open(f"SAT_HIT_{start+idx}.json", "w"))
            break
        if idx < 10:
            print(f"  cube {start+idx}: {'SAT' if res else 'UNSAT'} {dt:.2f}s", flush=True)
        if (idx + 1) % 25 == 0:
            w = sorted(times[max(0, idx - 99):])
            med = w[len(w)//2]; p90 = w[int(len(w)*0.9)]
            print(f"  cube {start+idx+1}: warm median(100)={med:.3f}s p90={p90:.3f}s last={dt:.3f}s "
                  f"total={sum(times):.0f}s", flush=True)
    w = sorted(times)
    out = {"start": start, "n": len(times), "sat": int(any(statuses)),
           "median": w[len(w)//2] if w else None,
           "p90": w[int(len(w)*0.9)] if w else None,
           "mean": sum(times)/len(times) if times else None,
           "first25_mean": sum(times[:25])/min(25, len(times)) if times else None,
           "last25_mean": sum(times[-25:])/min(25, len(times)) if times else None,
           "total_s": sum(times)}
    print("STATS " + json.dumps(out), flush=True)

if __name__ == "__main__":
    main()
