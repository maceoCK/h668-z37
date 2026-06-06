"""
Lift search: given the 9x9 orbit matrix R (block sizes), find Z_37 connection sets
D[i][j] (|D[i][j]|=R[i][j], D[j][i]=-D[i][j], D[i][i] symmetric) so that the full
block-circulant A is srg(333,166,82,83), i.e. for every nonzero character t,
   M_t^2 + M_t - 83 I = 0,   where M_t[i][j] = sum_{d in D[i][j]} w^{t d}, w=exp(2pi i/37).
(The t=0 condition M_0=R is automatic from the sizes.)
Objective F = sum_{t=1}^{36} || M_t^2 + M_t - 83 I ||_F^2.  F=0  <=>  srg(333) <=> H(668).

Validated first on Paley P(25) (objective must be 0 for its true connection sets).
"""
import numpy as np, sys, time, os
import conf_core as cc
HERE = os.path.dirname(os.path.abspath(__file__))

def make_W(p):
    w = np.exp(2j*np.pi/p)
    d = np.arange(p)
    return w**(np.outer(d, d) % p)   # W[t,d] = w^{t d}

def Mt_from_D(D, p, norb, W):
    M = np.zeros((p, norb, norb), dtype=complex)
    for i in range(norb):
        for j in range(norb):
            if D[i][j]:
                idx = np.fromiter(D[i][j], dtype=int)
                M[:, i, j] = W[:, idx].sum(axis=1)
    return M

def objective(M, k, lam, mu, p):
    I = np.eye(M.shape[1])
    res = M @ M + (mu-lam)*M - (k-mu)*I[None, :, :]   # should be 0 for t!=0
    res = res[1:]                                      # skip t=0
    return float(np.sum(np.abs(res)**2))

def validate():
    for p in (3, 5, 7):
        v=p*p; k=(v-1)//2; lam=(v-5)//4; mu=(v-1)//4
        D = cc.paley_pp_blocks(p)
        W = make_W(p)
        M = Mt_from_D(D, p, p, W)
        f = objective(M, k, lam, mu, p)
        print(f"  validate P({v})=srg({v},{k},{lam},{mu}): objective={f:.3e} (expect ~0)")

# ---- local search on srg(333) ----
def random_init(R, p, norb, rng):
    D = [[set() for _ in range(norb)] for _ in range(norb)]
    nz = list(range(1, p))
    for i in range(norb):
        # diagonal: symmetric subset of {1..p-1}, size R[i][i] (even); pick R[i][i]/2 pairs {d,-d}
        pairs = list(range(1, p//2+1))
        rng.shuffle(pairs)
        need = R[i][i]//2
        S = set()
        for d in pairs[:need]:
            S |= {d, (p-d) % p}
        D[i][i] = S
        for j in range(i+1, norb):
            S = set(int(x) for x in rng.choice(nz, size=R[i][j], replace=False))
            D[i][j] = S
            D[j][i] = set((-d) % p for d in S)
    return D

def lift_search(R, p, norb, k, lam, mu, seed, iters, log_every=20000):
    rng = np.random.default_rng(seed)
    W = make_W(p)
    D = random_init(R, p, norb, rng)
    M = Mt_from_D(D, p, norb, W)
    F = objective(M, k, lam, mu, p)
    best = F
    blocks = [(i, j) for i in range(norb) for j in range(norb) if i < j] + [(i, i) for i in range(norb)]
    T = max(1.0, F/3000)
    Thi = T
    t0 = time.time()
    stall = 0; lastbest = best
    for it in range(iters):
        T *= 0.9999985
        if T < 0.5: T = 0.5
        # reheat if stalled
        if it % 50000 == 0:
            if best >= lastbest - 1e-9:
                stall += 1
                if stall >= 2:
                    T = Thi * 0.6; stall = 0   # reheat
            else:
                stall = 0
            lastbest = best
        i, j = blocks[rng.integers(len(blocks))]
        Dij = D[i][j]
        if i == j:
            # swap a pair {d,-d} for another pair to keep symmetric+size
            inset = sorted(x for x in Dij if x <= p//2)
            outset = [d for d in range(1, p//2+1) if d not in Dij]
            if not inset or not outset: continue
            dr = inset[rng.integers(len(inset))]
            da = outset[rng.integers(len(outset))]
            newD = (Dij - {dr, (p-dr) % p}) | {da, (p-da) % p}
        else:
            inset = list(Dij); outset = [d for d in range(1, p) if d not in Dij]
            if not outset: continue
            dr = inset[rng.integers(len(inset))]
            da = outset[rng.integers(len(outset))]
            newD = (Dij - {dr}) | {da}
        # recompute M for blocks (i,j),(j,i)  (and diagonal just (i,i))
        oldM_ij = M[:, i, j].copy()
        if i == j:
            idx = np.fromiter(newD, dtype=int); M[:, i, i] = W[:, idx].sum(1)
        else:
            idx = np.fromiter(newD, dtype=int); M[:, i, j] = W[:, idx].sum(1)
            M[:, j, i] = np.conj(M[:, i, j])
        newF = objective(M, k, lam, mu, p)
        if newF <= F or rng.random() < np.exp((F-newF)/T):
            D[i][j] = newD
            if i != j:
                D[j][i] = set((-d) % p for d in newD)
            F = newF
            if F < best:
                best = F
                if best < 1e-6:
                    return D, 0.0, it
        else:
            # revert
            M[:, i, j] = oldM_ij
            if i != j: M[:, j, i] = np.conj(oldM_ij)
        if (it+1) % log_every == 0:
            print(f"   seed{seed} it{it+1}: F={F:.1f} best={best:.1f} T={T:.2f} ({time.time()-t0:.0f}s)")
            sys.stdout.flush()
    return D, best, iters

if __name__ == "__main__":
    print("=== validate objective on Paley graphs ===")
    validate()
    print("=== lift search on srg(333,166,82,83) with found orbit matrix ===")
    R = np.load(os.path.join(HERE,"orbitR_found.npy"))
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    best_overall = 1e18
    for s in range(seeds):
        D, best, it = lift_search(R, 37, 9, 166, 82, 83, s, iters)
        print(f"  seed {s}: best objective = {best:.2f}")
        best_overall = min(best_overall, best)
        if best < 1e-6:
            print("  *** LIFT FOUND -> srg(333) -> H(668)! verifying... ***")
            A = cc.build_adjacency(D, 37, 9)
            ok = cc.srg_identity_direct(A, 166, 82, 83)
            print(f"  direct SRG verification: {ok}")
            if ok:
                np.save(os.path.join(HERE,"srg333_adjacency.npy"), A)
                break
    print(f"=== best objective across seeds: {best_overall:.2f} (0 => H(668)) ===")
