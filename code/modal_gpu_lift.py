"""
Modal GPU app: batched parallel-tempering LIFT search for H(668) on the Z_37 conference route.
Each GPU container runs B~8192+ chains as tensor ops -> 100 containers ~ 10^6 chains, sidestepping
the 100-container concurrency cap. Validated end-to-end (gpu_lift.py): finds a real lift on P(25).

Run (after `modal token new`):
  modal run modal_gpu_lift.py --b 8192 --steps 20000 --seeds 1
Any chain reaching F=0 -> verified srg(333) -> H(668); written to the h668-hits Modal Dict, then
emit locally with `python emit_h668.py emit`.
"""
import modal

app = modal.App("h668-gpu-lift")
image = (modal.Image.debian_slim(python_version="3.11")
         .pip_install("numpy", "torch"))
P, NORB, K, LAM, MU = 37, 9, 166, 82, 83

def _verify_srg(D, p, norb, k, lam, mu):
    import numpy as np
    V = p * norb
    A = np.zeros((V, V), dtype=np.int64)
    for i in range(norb):
        for j in range(norb):
            for a in range(p):
                for d in D[i][j]:
                    A[i * p + a, j * p + (a + d) % p] = 1
    lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(V, dtype=np.int64)
    return bool(np.array_equal(lhs, mu * np.ones((V, V), dtype=np.int64))) and bool(np.array_equal(A, A.T))


@app.function(image=image, gpu="A10G", timeout=3600, retries=0)
def gpu_attempt(R, seed=0, B=8192, n_levels=16, steps=20000):
    import torch, numpy as np, time
    p, norb, k, lam, mu = P, NORB, K, LAM, MU
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    dt = torch.complex64
    rng = np.random.default_rng(seed)
    t = torch.arange(p, device=dev)
    ang = 2 * np.pi * (torch.outer(t, t) % p).to(torch.float64) / p
    W = torch.complex(torch.cos(ang), torch.sin(ang)).to(dt)
    I = torch.eye(norb, dtype=dt, device=dev)

    def obj(X):
        M = torch.einsum('bijg,tg->btij', X.to(dt), W)
        Rr = (M @ M + (mu - lam) * M - (k - mu) * I)[:, 1:]
        return (Rr.real ** 2 + Rr.imag ** 2).sum(dim=(1, 2, 3))

    X = torch.zeros(B, norb, norb, p, device=dev)
    nz = np.arange(1, p)
    for b in range(B):
        for i in range(norb):
            pr = list(range(1, p // 2 + 1)); rng.shuffle(pr)
            for v in pr[:R[i][i] // 2]:
                X[b, i, i, v] = 1.0; X[b, i, i, (p - v) % p] = 1.0
            for j in range(i + 1, norb):
                for v in rng.choice(nz, size=int(R[i][j]), replace=False):
                    X[b, i, j, int(v)] = 1.0; X[b, j, i, (p - int(v)) % p] = 1.0
    F = obj(X); best = F.min().item()
    Tmax = max(1.0, best / 2000.0)
    levels = torch.logspace(np.log10(0.3), np.log10(Tmax), n_levels, device=dev)
    temp = levels.repeat(B // n_levels + 1)[:B]
    blocks = [(i, j) for i in range(norb) for j in range(norb) if i < j and R[i][j] > 0] + \
             [(i, i) for i in range(norb) if R[i][i] > 0]
    ar = torch.arange(B, device=dev); half = p // 2
    t0 = time.time()
    for s in range(steps):
        i, j = blocks[s % len(blocks)]
        Xp = X.clone()
        if i == j:
            mem = X[:, i, i, 1:half + 1]
            do = torch.multinomial(mem + 1e-9, 1).squeeze(1) + 1
            di = torch.multinomial((1 - mem) + 1e-9, 1).squeeze(1) + 1
            Xp[ar, i, i, do] = 0.0; Xp[ar, i, i, (p - do) % p] = 0.0
            Xp[ar, i, i, di] = 1.0; Xp[ar, i, i, (p - di) % p] = 1.0
        else:
            mem = X[:, i, j, :]
            do = torch.multinomial(mem + 1e-9, 1).squeeze(1)
            di = torch.multinomial((1 - mem) + 1e-9, 1).squeeze(1)
            Xp[ar, i, j, do] = 0.0; Xp[ar, i, j, di] = 1.0
            Xp[ar, j, i, (p - do) % p] = 0.0; Xp[ar, j, i, (p - di) % p] = 1.0
        Fp = obj(Xp)
        acc = (Fp <= F) | (torch.rand(B, device=dev) < torch.exp((F - Fp) / temp))
        X = torch.where(acc.view(-1, 1, 1, 1), Xp, X); F = torch.where(acc, Fp, F)
        if s % 20 == 0 and s:
            Xr = X.view(-1, n_levels, norb, norb, p); Fr = F.view(-1, n_levels); Tr = temp.view(-1, n_levels)
            for l in range(0, n_levels - 1, 2):
                d = (Fr[:, l] - Fr[:, l + 1]) * (1.0 / Tr[:, l] - 1.0 / Tr[:, l + 1])
                sw = torch.rand(Xr.shape[0], device=dev) < torch.exp(torch.clamp(d, max=0))
                a, bb = Xr[:, l].clone(), Xr[:, l + 1].clone()
                Xr[:, l] = torch.where(sw.view(-1, 1, 1, 1), bb, a); Xr[:, l + 1] = torch.where(sw.view(-1, 1, 1, 1), a, bb)
                fa, fb = Fr[:, l].clone(), Fr[:, l + 1].clone()
                Fr[:, l] = torch.where(sw, fb, fa); Fr[:, l + 1] = torch.where(sw, fa, fb)
            X = Xr.view(B, norb, norb, p); F = Fr.view(B)
        c = F.min().item()
        if c < best: best = c
        if best < 1e-6:
            bi = int(torch.argmin(F).item())
            D = [[sorted(g for g in range(p) if X[bi, a, b2, g] > 0.5) for b2 in range(norb)] for a in range(norb)]
            if _verify_srg(D, p, norb, k, lam, mu):
                try:
                    import modal as _m
                    _m.Dict.from_name("h668-hits", create_if_missing=True)["solution"] = {"R": R, "D": D}
                except Exception:
                    pass
                return {"best": 0.0, "verified": True, "solution": D, "step": s}
    # persist per-skeleton best F so a dropped driver never loses coverage data
    try:
        import modal as _m
        _m.Dict.from_name("h668-gpu-results", create_if_missing=True)[str(seed) + "_" + str(hash(tuple(map(tuple, R))) % 10_000_000)] = best
    except Exception:
        pass
    return {"best": best, "verified": False, "solution": None, "secs": time.time() - t0}


def _load(path):
    skels, seen = [], set()
    for line in open(path):
        v = [int(t) for t in line.split()]
        if len(v) != 81: continue
        Rr = [v[9 * i:9 * i + 9] for i in range(9)]
        key = tuple(v)
        if key not in seen: seen.add(key); skels.append(Rr)
    return skels


@app.local_entrypoint()
def main(b: int = 8192, steps: int = 20000, seeds: int = 1, skeleton_file: str = "skeletons.txt",
         max_skeletons: int = 0):
    import json
    skels = _load(skeleton_file)
    if max_skeletons: skels = skels[:max_skeletons]
    jobs = [(R, s, b, 16, steps) for R in skels for s in range(seeds)]
    print(f"GPU sweep: {len(skels)} skeletons x {seeds} seeds = {len(jobs)} jobs, B={b}, steps={steps}")
    best = 1e18; found = None; n = 0
    for res in gpu_attempt.starmap(jobs):
        n += 1
        if res["best"] < best: best = res["best"]; print(f"  [{n}] new best F = {best:.1f}")
        if res.get("verified"):
            found = res; print(f"*** VERIFIED LIFT -> H(668)! step={res['step']} ***"); break
        if n % 50 == 0: print(f"  {n}/{len(jobs)} done, best F={best:.1f}")
    print(f"=== GPU sweep best F = {best:.2f}  (0 => H(668)) ===")
    if found:
        json.dump(found, open("LIFT_SOLUTION.json", "w")); print("saved LIFT_SOLUTION.json")
