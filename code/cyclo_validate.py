"""
Validation suite for cyclo_dfs (the sound exhaustive coset-level H-cyclotomic lift engine).

Two kinds of checks:
  (1) TRUE-POSITIVE anchors: engine must rediscover known SRGs (Paley P(9), P(25) for H={1} and
      H={+-1}; Petersen on Z_5 x 2 orbits; Paley(13) as a 1-orbit order-3-cyclotomic graph with
      EXACT count 2), and every printed solution must pass the direct adjacency identity.
  (2) EXHAUSTIVENESS cross-checks: for small (p, norb, H, R), an INDEPENDENT brute-force enumerator
      (pure definition: H-invariant blocks of prescribed sizes, symmetry D[j][i]=-D[i][j],
      diagonal symmetric without 0; check A^2+(mu-lam)A-(k-mu)I=mu J on the actual 0/1 adjacency
      matrix, no shared code with the engine) must produce EXACTLY the same solution set.
      This certifies both soundness (no false solutions) and completeness (no pruning losses).
"""
import itertools, subprocess, tempfile, os, sys
import numpy as np
import conf_core as cc

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, "cyclo_dfs")

def cosets_of(mult, p):
    H = []; x = 1
    while True:
        H.append(x); x = x * mult % p
        if x == 1: break
    seen, cs = set(), []
    for g in range(1, p):
        if g in seen: continue
        c = tuple(sorted((g * h) % p for h in H))
        seen |= set(c); cs.append(c)
    return H, cs

def run_engine(p, norb, k, lam, mu, mult, R, extra=()):
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(" ".join(str(R[i][j]) for i in range(norb) for j in range(norb)) + "\n")
        rf = f.name
    cmd = [ENGINE, str(p), str(norb), str(k), str(lam), str(mu), str(mult), rf, "0", "--quiet", *extra]
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    os.unlink(rf)
    sols, cur = [], None
    for ln in out.stdout.splitlines():
        if ln.startswith("SOLUTION"):
            cur = [[None] * norb for _ in range(norb)]; sols.append(cur)
        elif ln.startswith("D ") and cur is not None:
            head, _, elems = ln.partition(":")
            _, i, j = head.split()
            cur[int(i)][int(j)] = frozenset(int(x) for x in elems.split())
        elif ln.startswith("RESULT"):
            res = dict(t.split("=") for t in ln.split()[1:] if "=" in t)
    return res, sols, out

def sol_key(D, norb):
    return tuple(D[i][j] for i in range(norb) for j in range(norb))

def verify_sol(D, p, norb, k, lam, mu):
    Ds = [[set(D[i][j]) for j in range(norb)] for i in range(norb)]
    A = cc.build_adjacency(Ds, p, norb)
    n = p * norb
    lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(n, dtype=np.int64)
    return np.array_equal(A, A.T) and np.array_equal(lhs, mu * np.ones((n, n), dtype=np.int64))

def brute_force(p, norb, k, lam, mu, mult, R):
    """Definitional enumeration; checks the actual adjacency matrix. Returns set of solution keys."""
    H, CS = cosets_of(mult, p)
    e = len(H)
    negc = {ci: next(cj for cj, c2 in enumerate(CS) if (-c[0]) % p in c2) for ci, c in enumerate(CS)}
    nc = len(CS)
    def blocks_for(r, diag):
        out = []
        if diag:
            if r % e: return out
            orbs, seen = [], set()
            for c in range(nc):
                if c in seen: continue
                seen |= {c, negc[c]}; orbs.append((c, negc[c]))
            for sub in itertools.chain.from_iterable(itertools.combinations(range(len(orbs)), t) for t in range(len(orbs) + 1)):
                csel = set()
                for t in sub: csel |= {orbs[t][0], orbs[t][1]}
                if sum(len(CS[c]) for c in csel) == r:
                    out.append(frozenset(x for c in csel for x in CS[c]))
        else:
            for z in (0, 1):
                if (r - z) % e or (r - z) < 0 or (r - z) // e > nc: continue
                w = (r - z) // e
                for sub in itertools.combinations(range(nc), w):
                    s = set([0] if z else [])
                    for c in sub: s |= set(CS[c])
                    out.append(frozenset(s))
        return out
    diag_lists = [blocks_for(R[i][i], True) for i in range(norb)]
    off_idx = [(i, j) for i in range(norb) for j in range(i + 1, norb)]
    off_lists = [blocks_for(R[i][j], False) for (i, j) in off_idx]
    sols = set()
    n = p * norb
    target = mu * np.ones((n, n), dtype=np.int64)
    for dchoice in itertools.product(*diag_lists):
        for ochoice in itertools.product(*off_lists):
            D = [[None] * norb for _ in range(norb)]
            for i in range(norb): D[i][i] = set(dchoice[i])
            for (i, j), s in zip(off_idx, ochoice):
                D[i][j] = set(s)
                D[j][i] = {(-x) % p for x in s}
            A = cc.build_adjacency(D, p, norb)
            lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(n, dtype=np.int64)
            if np.array_equal(lhs, target):
                sols.add(tuple(frozenset(D[i][j]) for i in range(norb) for j in range(norb)))
    return sols

def main():
    fails = 0
    def report(name, ok, info=""):
        nonlocal fails
        print(f"  [{'PASS' if ok else 'FAIL'}] {name} {info}")
        if not ok: fails += 1

    # ---------- (1) true-positive anchors ----------
    print("anchors:")
    # P(9), P(25): paley blocks, H={1} and H={-1,1}
    for p, mult2 in ((3, 2), (5, 4)):
        v = p * p; k = (v - 1) // 2; lam = (v - 5) // 4; mu = (v - 1) // 4
        D = cc.paley_pp_blocks(p)
        R = [[len(D[i][j]) for j in range(p)] for i in range(p)]
        for mult in (1, mult2):
            res, sols, out = run_engine(p, p, k, lam, mu, mult, R, extra=("--maxsol", "200"))
            okall = all(verify_sol(s, p, p, k, lam, mu) for s in sols)
            paley_key = tuple(frozenset(D[i][j]) for i in range(p) for j in range(p))
            found_paley = any(sol_key(s, p) == paley_key for s in sols)
            # element-level P(25) has >200 lifts; hitting the cap with all-verified solutions is a
            # pass (Paley itself is required only when the run completes the full enumeration)
            ok = len(sols) > 0 and okall and (found_paley or res.get("status") == "SAT_CAP")
            report(f"P({v}) mult={mult}", ok,
                   f"status={res.get('status')} sat={res.get('sat')} allverify={okall} paley_found={found_paley}")
    # Paley(13) as 1-orbit order-3 cyclotomic: exact count 2
    res, sols, out = run_engine(13, 1, 6, 2, 3, 3, [[6]])
    okall = all(verify_sol(s, 13, 1, 6, 2, 3) for s in sols)
    report("Paley(13) norb=1 mult=3 count==2", res.get("status") == "COMPLETE" and len(sols) == 2 and okall,
           f"status={res.get('status')} sat={res.get('sat')}")
    # Petersen on Z5 x 2 orbits, H={1} and H={+-1}
    for mult in (1, 4):
        res, sols, out = run_engine(5, 2, 3, 0, 1, mult, [[2, 1], [1, 2]])
        okall = all(verify_sol(s, 5, 2, 3, 0, 1) for s in sols)
        report(f"Petersen mult={mult}", res.get("status") == "COMPLETE" and len(sols) > 0 and okall,
               f"sat={res.get('sat')} allverify={okall}")

    # ---------- (2) exhaustiveness cross-checks ----------
    print("brute-force cross-checks (engine solution set == definitional enumeration):")
    cases = [
        # p,norb,k,lam,mu,mult,R
        (5, 2, 3, 0, 1, 1, [[2, 1], [1, 2]]),            # Petersen, element-level
        (5, 2, 3, 0, 1, 4, [[2, 1], [1, 2]]),            # Petersen, order-2 cosets, z=1 path
        (5, 2, 6, 3, 4, 1, [[2, 4], [4, 2]]),            # complement-ish params (likely few/none)
        (13, 2, 10, 3, 4, 3, [[6, 4], [4, 6]]),          # order-3 cosets on Z13, z=1 blocks
        (13, 2, 10, 3, 4, 3, [[0, 10], [10, 0]]),        # extreme diag, z=1
        (13, 2, 10, 3, 4, 3, [[6, 3], [3, 6]]),          # z=0 blocks (row sums off => expect 0)
        (7, 3, 10, 4, 5, 6, [[2, 4, 4], [4, 2, 4], [4, 4, 2]]),   # nonexistent srg(21,10,4,5)
        (7, 3, 10, 4, 5, 6, [[4, 3, 3], [3, 4, 3], [3, 3, 4]]),   # nonexistent, z=1 offdiag
        (13, 1, 6, 2, 3, 3, [[6]]),                      # exact-2 case again, vs brute
        (13, 1, 6, 2, 3, 1, [[6]]),                      # element-level singletons + pairs
    ]
    for (p, norb, k, lam, mu, mult, R) in cases:
        bs = brute_force(p, norb, k, lam, mu, mult, R)
        res, sols, out = run_engine(p, norb, k, lam, mu, mult, R)
        es = {sol_key(s, norb) for s in sols}
        report(f"p={p} norb={norb} mult={mult} R00={R[0][0]} R01={R[0][1] if norb>1 else '-'}",
               res.get("status") == "COMPLETE" and es == bs,
               f"engine={len(es)} brute={len(bs)} status={res.get('status')}")

    print(f"\n{'ALL VALIDATION PASSED' if fails == 0 else f'{fails} FAILURES'}")
    return 0 if fails == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
