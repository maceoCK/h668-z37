"""
Modal app: parallel LIFT attack on H(668) via the Z_37 conference route.
Fans out, over hundreds of orbit-matrix skeletons (skeletons.txt from conference_orbit_all),
a per-skeleton search for Z_37 connection sets making the 9x9 block-circulant a genuine
srg(333,166,82,83) == conference matrix C(334) == H(668).

Two tracks (choose with --mode):
  exact : CP-SAT model of the convolution/SRG equations (validated reconstructing Paley P(9)/P(25)).
          Per skeleton returns SAT (-> H(668)), UNSAT (skeleton ruled out), or UNKNOWN (timeout).
  sa    : cheap simulated-annealing probe in the character domain (objective validated =0 on Paley).

Setup (you, once):  modal token new
Run (exact, the real attempt):  modal run modal_lift.py --mode exact --time-s 3600
Run (cheap probe):              modal run modal_lift.py --mode sa --seeds 2000 --iters 4000000

Any SAT solution is verified in-container (A^2+A-83I==83J) before being returned, then saved to
LIFT_SOLUTION.json for local conversion to the 668x668 Hadamard matrix.
"""
import modal

app = modal.App("h668-z37-lift")
image = modal.Image.debian_slim().pip_install("numpy", "ortools")

P, NORB, K, LAM, MU = 37, 9, 166, 82, 83

# ---------- shared verifier (builds the 333x333 adjacency, checks the SRG identity) ----------
def _verify_srg(D, p, norb, k, lam, mu):
    import numpy as np
    V = p * norb
    A = np.zeros((V, V), dtype=np.int64)
    for i in range(norb):
        for j in range(norb):
            for a in range(p):
                base = i * p + a
                for d in D[i][j]:
                    A[base, j * p + (a + d) % p] = 1
    lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(V, dtype=np.int64)
    return bool(np.array_equal(lhs, mu * np.ones((V, V), dtype=np.int64))) and bool(np.array_equal(A, A.T))


@app.function(image=image, cpu=8.0, timeout=86400, retries=0)
def exact_lift(R, time_s=3600, seed=0, mult=None):
    """CP-SAT exact lift for one skeleton R. mult!=None => cyclotomic ansatz (blocks <mult>-invariant,
    fast+definitive). Returns dict(status, solution|None)."""
    from ortools.sat.python import cp_model
    p, norb, k, lam, mu = P, NORB, K, LAM, MU
    m = cp_model.CpModel()
    x = {(i, j, g): m.NewBoolVar(f"x{i}_{j}_{g}")
         for i in range(norb) for j in range(norb) for g in range(p)}
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                m.Add(x[j, i, g] == x[i, j, (-g) % p])
    for i in range(norb):
        m.Add(x[i, i, 0] == 0)
    for i in range(norb):
        for j in range(norb):
            m.Add(sum(x[i, j, g] for g in range(p)) == int(R[i][j]))
    # GAUGE FIX: per-orbit base-point freedom is a p^(norb-1) symmetry; force 0 in D[0][i].
    for i in range(1, norb):
        if R[0][i] > 0:
            m.Add(x[0, i, 0] == 1)
    # CYCLOTOMIC ANSATZ: force each block invariant under <mult> (union of cosets).
    if mult is not None:
        for i in range(norb):
            for j in range(norb):
                for g in range(1, p):
                    m.Add(x[i, j, g] == x[i, j, (mult * g) % p])
    prod = {}
    def AND(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in prod:
            c = m.NewBoolVar(f"p{len(prod)}")
            m.AddMultiplicationEquality(c, [x[a], x[b]])
            prod[key] = c
        return prod[key]
    for i in range(norb):
        for j in range(norb):
            for g in range(p):
                terms = [AND((i, ip, d), (ip, j, (g - d) % p)) for ip in range(norb) for d in range(p)]
                if i == j and g == 0:
                    m.Add(sum(terms) == k)
                else:
                    m.Add(sum(terms) == mu + (lam - mu) * x[i, j, g])
    sv = cp_model.CpSolver()
    sv.parameters.max_time_in_seconds = float(time_s)
    sv.parameters.num_search_workers = 8
    sv.parameters.random_seed = int(seed)
    st = sv.Solve(m)
    name = sv.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        D = [[sorted(g for g in range(p) if sv.Value(x[i, j, g])) for j in range(norb)] for i in range(norb)]
        verified = _verify_srg(D, p, norb, k, lam, mu)
        if verified:
            try:
                import modal as _m
                _m.Dict.from_name("h668-hits", create_if_missing=True)["solution"] = {"R": R, "D": D}
            except Exception:
                pass
        return {"status": "SAT", "verified": verified, "solution": D, "R": R}
    return {"status": name, "verified": False, "solution": None, "R": R}


@app.function(image=image, cpu=2.0, timeout=14400, retries=0)
def sa_lift(R, seed, iters):
    """Cheap SA probe in the character domain. Returns dict(best, solution|None)."""
    import numpy as np
    p, norb = P, NORB
    rng = np.random.default_rng(seed)
    d = np.arange(p); W = np.exp(2j * np.pi / p) ** (np.outer(d, d) % p); I = np.eye(norb)
    def Mt(D):
        M = np.zeros((p, norb, norb), dtype=complex)
        for i in range(norb):
            for j in range(norb):
                if D[i][j]:
                    M[:, i, j] = W[:, np.fromiter(D[i][j], dtype=int)].sum(1)
        return M
    def obj(M):
        return float(np.sum(np.abs((M @ M + (MU - LAM) * M - (K - MU) * I[None])[1:]) ** 2))
    D = [[set() for _ in range(norb)] for _ in range(norb)]; nz = list(range(1, p))
    for i in range(norb):
        pr = list(range(1, p // 2 + 1)); rng.shuffle(pr); S = set()
        for v in pr[:R[i][i] // 2]: S |= {v, (p - v) % p}
        D[i][i] = S
        for j in range(i + 1, norb):
            S = set(int(v) for v in rng.choice(nz, size=R[i][j], replace=False))
            D[i][j] = S; D[j][i] = set((-v) % p for v in S)
    M = Mt(D); F = obj(M); best = F; Thi = max(1.0, F / 3000.0); T = Thi
    blocks = [(i, j) for i in range(norb) for j in range(norb) if i < j] + [(i, i) for i in range(norb)]
    for it in range(iters):
        T = max(0.5, T * 0.9999985)
        i, j = blocks[rng.integers(len(blocks))]; Dij = D[i][j]; oldM = M[:, i, j].copy()
        if i == j:
            ins = [v for v in Dij if v <= p // 2]; out = [v for v in range(1, p // 2 + 1) if v not in Dij]
            if not ins or not out: continue
            dr = ins[rng.integers(len(ins))]; da = out[rng.integers(len(out))]
            newD = (Dij - {dr, (p - dr) % p}) | {da, (p - da) % p}
        else:
            ins = list(Dij); out = [v for v in range(1, p) if v not in Dij]
            if not out: continue
            dr = ins[rng.integers(len(ins))]; da = out[rng.integers(len(out))]
            newD = (Dij - {dr}) | {da}
        M[:, i, j] = W[:, np.fromiter(newD, dtype=int)].sum(1)
        if i != j: M[:, j, i] = np.conj(M[:, i, j])
        nF = obj(M)
        if nF <= F or rng.random() < np.exp((F - nF) / T):
            D[i][j] = newD
            if i != j: D[j][i] = set((-v) % p for v in newD)
            F = nF
            if F < best:
                best = F
                if best < 1e-6:
                    return {"best": 0.0, "solution": [[sorted(D[a][b]) for b in range(norb)] for a in range(norb)]}
        else:
            M[:, i, j] = oldM
            if i != j: M[:, j, i] = np.conj(oldM)
    return {"best": best, "solution": None}


def _load_skeletons(path):
    import itertools
    skels, seen = [], set()
    for line in open(path):
        v = [int(t) for t in line.split()]
        if len(v) != 81: continue
        R = [v[9 * i:9 * i + 9] for i in range(9)]
        key = tuple(itertools.chain.from_iterable(R))
        if key not in seen: seen.add(key); skels.append(R)
    return skels


@app.local_entrypoint()
def main(mode: str = "exact", time_s: int = 3600, seeds: int = 1000,
         iters: int = 4000000, skeleton_file: str = "skeletons.txt", max_skeletons: int = 0):
    import json
    skels = _load_skeletons(skeleton_file)
    if max_skeletons: skels = skels[:max_skeletons]
    print(f"{len(skels)} unique skeletons, mode={mode}")
    found = None
    if mode in ("exact", "cyclo"):
        if mode == "cyclo":
            # one multiplier generator per nontrivial subgroup of Z_37* (cheap, definitive track)
            def primroot(p):
                for g in range(2, p):
                    x, seen = 1, set()
                    for _ in range(p - 1):
                        x = x * g % p; seen.add(x)
                    if len(seen) == p - 1: return g
            g = primroot(37)
            mults = sorted(set(pow(g, 36 // e, 37) for e in (2, 3, 4, 6, 9, 12, 18, 36)))
            jobs = [(R, time_s, 0, mm) for R in skels for mm in mults]
            print(f"cyclo: {len(skels)} skeletons x {len(mults)} multipliers = {len(jobs)} fast jobs")
        else:
            jobs = [(R, time_s, 0, None) for R in skels]
        tally = {}
        for res in exact_lift.starmap(jobs):
            tally[res["status"]] = tally.get(res["status"], 0) + 1
            if res["status"] == "SAT" and res["verified"]:
                found = res; print("*** VERIFIED LIFT -> H(668)! ***"); break
            if sum(tally.values()) % 50 == 0:
                print(f"  progress: {tally}")
        print(f"=== {mode} tally: {tally} ===  (UNSAT/INFEASIBLE=ruled out; UNKNOWN=timeout)")
    else:
        jobs = [(R, s, iters) for R in skels for s in range(seeds)]
        best = 1e18
        for res in sa_lift.starmap(jobs):
            if res["best"] < best: best = res["best"]; print(f"  best objective: {best:.1f}")
            if res["solution"] is not None: found = {"solution": res["solution"], "verified": None}; break
        print(f"=== SA best objective: {best:.2f} (0 => H(668)) ===")
    if found:
        json.dump(found, open("LIFT_SOLUTION.json", "w"))
        print("saved LIFT_SOLUTION.json")
