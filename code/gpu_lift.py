"""
GPU batched parallel-tempering LIFT search for H(668) (Z_p block-circulant srg).
Runs B chains at once as tensor ops. Objective F = sum_{t=1}^{p-1} ||M_t^2+(mu-lam)M_t-(k-mu)I||_F^2,
M_t[i][j]=sum_{d in D[i][j]} w^{td}.  F=0 <=> srg(v,k,lam,mu).  One GPU runs B~10^4 chains -> sidesteps
the 100-container concurrency cap (each container = thousands of chains).

Validated: (a) F==0 on true Paley P(9)/P(25) blocks; (b) PT search FINDS a lift (F->0) on the P(25)
skeleton (which is constructible) — end-to-end check before any GPU spend on srg(333).
"""
import torch, numpy as np, itertools, time
import conf_core as cc

def make_W(p, device, dtype=torch.complex64):
    t = torch.arange(p, device=device)
    ang = 2 * np.pi * (torch.outer(t, t) % p).to(torch.float64) / p
    return torch.complex(torch.cos(ang), torch.sin(ang)).to(dtype)

def objective(X, W, k, lam, mu):
    # X: [B,norb,norb,p] real(0/1);  returns F: [B]
    M = torch.einsum('bijg,tg->btij', X.to(W.dtype), W)         # [B,p,norb,norb]
    norb = X.shape[1]
    I = torch.eye(norb, dtype=W.dtype, device=X.device)
    Rr = (M @ M + (mu - lam) * M - (k - mu) * I)[:, 1:]          # drop t=0
    return (Rr.real ** 2 + Rr.imag ** 2).sum(dim=(1, 2, 3))

def X_from_D(D, p, norb, device='cpu'):
    X = torch.zeros(1, norb, norb, p, device=device)
    for i in range(norb):
        for j in range(norb):
            for g in D[i][j]:
                X[0, i, j, g] = 1.0
    return X

def random_states(R, p, norb, B, device, rng):
    X = torch.zeros(B, norb, norb, p, device=device)
    nz = np.arange(1, p)
    for b in range(B):
        for i in range(norb):
            pr = list(range(1, p // 2 + 1)); rng.shuffle(pr)
            for v in pr[:R[i][i] // 2]:
                X[b, i, i, v] = 1.0; X[b, i, i, (p - v) % p] = 1.0
            for j in range(i + 1, norb):
                S = rng.choice(nz, size=int(R[i][j]), replace=False)
                for v in S:
                    X[b, i, j, v] = 1.0; X[b, j, i, (p - int(v)) % p] = 1.0
    return X

def pt_lift(R, p, norb, k, lam, mu, B=4096, n_levels=16, steps=4000, exch_every=20,
            device='cpu', seed=0, dtype=torch.complex64):
    rng = np.random.default_rng(seed)
    W = make_W(p, device, dtype)
    X = random_states(R, p, norb, B, device, rng)
    F = objective(X, W, k, lam, mu)
    best = F.min().item()
    # temperature ladder tiled across chains
    Tmax = max(1.0, best / 2000.0)
    levels = torch.logspace(np.log10(0.3), np.log10(Tmax), n_levels, device=device)
    temp = levels.repeat(B // n_levels + 1)[:B]
    blocks = [(i, j) for i in range(norb) for j in range(norb) if i < j and R[i][j] > 0] + \
             [(i, i) for i in range(norb) if R[i][i] > 0]
    ar = torch.arange(B, device=device)
    half = p // 2
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
        Fp = objective(Xp, W, k, lam, mu)
        acc = (Fp <= F) | (torch.rand(B, device=device) < torch.exp((F - Fp) / temp))
        X = torch.where(acc.view(-1, 1, 1, 1), Xp, X)
        F = torch.where(acc, Fp, F)
        # replica exchange between adjacent temperature levels
        if s % exch_every == 0 and s:
            Xr = X.view(-1, n_levels, norb, norb, p); Fr = F.view(-1, n_levels); Tr = temp.view(-1, n_levels)
            for l in range(0, n_levels - 1, 2):
                d = (Fr[:, l] - Fr[:, l + 1]) * (1.0 / Tr[:, l] - 1.0 / Tr[:, l + 1])
                sw = torch.rand(Xr.shape[0], device=device) < torch.exp(torch.clamp(d, max=0))
                a, bb = Xr[:, l].clone(), Xr[:, l + 1].clone()
                Xr[:, l] = torch.where(sw.view(-1, 1, 1, 1), bb, a)
                Xr[:, l + 1] = torch.where(sw.view(-1, 1, 1, 1), a, bb)
                fa, fb = Fr[:, l].clone(), Fr[:, l + 1].clone()
                Fr[:, l] = torch.where(sw, fb, fa); Fr[:, l + 1] = torch.where(sw, fa, fb)
            X = Xr.view(B, norb, norb, p); F = Fr.view(B)
        cur = F.min().item()
        if cur < best: best = cur
        if best < 1e-6:
            bidx = int(torch.argmin(F).item())
            D = [[sorted(g for g in range(p) if X[bidx, a, b, g] > 0.5) for b in range(norb)] for a in range(norb)]
            return 0.0, D, s
    return best, None, steps

def validate():
    print("(a) objective on true Paley blocks (expect ~0):")
    for p in (3, 5, 7):
        v = p * p; k = (v - 1) // 2; lam = (v - 5) // 4; mu = (v - 1) // 4
        D = cc.paley_pp_blocks(p); W = make_W(p, 'cpu', torch.complex128)
        F = objective(X_from_D(D, p, p), W, k, lam, mu)
        print(f"   P({v}): F={F.item():.2e}")
    print("(b) PT search FINDS a lift on the P(25) skeleton (constructible):")
    p = 5; v = 25; k, lam, mu = 12, 5, 6
    R = cc.orbit_matrix(cc.paley_pp_blocks(p), p)
    t0 = time.time()
    best, D, s = pt_lift(R, p, p, k, lam, mu, B=2048, n_levels=16, steps=3000, device='cpu',
                         dtype=torch.complex128, seed=1)
    ok = False
    if D is not None:
        ok = cc.srg_identity_direct(cc.build_adjacency(D, p, p), k, lam, mu)
    print(f"   P(25) PT: best={best:.2e}  found_lift={D is not None}  verified_srg={ok}  step={s}  ({time.time()-t0:.0f}s)")

if __name__ == "__main__":
    validate()
