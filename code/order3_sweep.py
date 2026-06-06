"""
COMPLETE order-3 lift sweep. orbit3_enum.cpp enumerates EVERY order-3 coset-compatible 9x9
orbit matrix R (diagonal in {12,18,24}, off-diagonal == 0,1 mod 3) -- the full set of skeletons
that an order-3-invariant (H=<10>={1,10,26}) srg(333,166,82,83) could possibly have.

For each such R we run the exact lift restricted to <10>-invariant connection sets (R fixed,
gauge-fixed). If EVERY skeleton is UNSAT, then no order-3-invariant srg(333,166,82,83) exists --
COMPLETE (no sampling), strengthening Theorem 1 to 7 of 8. Any SAT would be an actual H(668).

Run:  ../.venv/bin/python order3_sweep.py [time_s_per_skel] [n_workers]
"""
import sys, os, time
from multiprocessing import Pool
import exact_lift as el

HERE = os.path.dirname(os.path.abspath(__file__))
P, NORB, K, LAM, MU = 37, 9, 166, 82, 83
MULT = 10  # generator of the order-3 multiplier subgroup H = {1,10,26} of Z_37*


def load(path):
    import itertools
    skels, seen = [], set()
    for line in open(path):
        v = [int(t) for t in line.split()]
        if len(v) != 81:
            continue
        R = [v[9 * i:9 * i + 9] for i in range(9)]
        key = tuple(itertools.chain.from_iterable(R))
        if key not in seen:
            seen.add(key); skels.append(R)
    return skels


def work(args):
    idx, R, time_s = args
    st, D = el.solve_lift(P, NORB, K, LAM, MU, R, time_s=time_s, workers=1, mult=MULT, gauge=True)
    return idx, st, (D is not None)


def main():
    time_s = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    skels = load(os.path.join(HERE, "orbit3_compatible.txt"))
    print(f"order-3 compatible skeletons: {len(skels)};  {time_s}s/skel, {nproc} workers, mult={MULT}")
    t0 = time.time()
    tally = {}; sat = []
    with Pool(nproc) as pool:
        for i, (idx, st, found) in enumerate(pool.imap_unordered(work, [(j, R, time_s) for j, R in enumerate(skels)]), 1):
            tally[st] = tally.get(st, 0) + 1
            if found:
                sat.append(idx)
                print(f"  *** SKELETON {idx}: SAT -> candidate order-3-invariant srg -> H(668)! ***")
            if i % 100 == 0 or i == len(skels):
                print(f"  {i}/{len(skels)}  tally={tally}  {time.time()-t0:.0f}s", flush=True)
    print(f"\n=== DONE === tally={tally}  SAT={sat}")
    unresolved = tally.get("UNKNOWN", 0) + tally.get("MODEL_INVALID", 0)
    if not sat and unresolved == 0:
        print("RESULT: every order-3 compatible skeleton is UNSAT => NO order-3-invariant "
              "srg(333,166,82,83) exists. Order-3 case CLOSED (COMPLETE). Theorem 1 -> 7 of 8.")
    elif sat:
        print("RESULT: SAT found -> verify immediately (possible H(668)).")
    else:
        print(f"RESULT: {unresolved} skeleton(s) timed out (UNKNOWN) -> re-run those with more time.")


if __name__ == "__main__":
    main()
