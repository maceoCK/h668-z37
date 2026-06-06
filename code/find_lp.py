"""
Fast incremental tabu/SA to FIND Legendre pairs at small lengths.
Maintains PAF_u, PAF_v incrementally; move = single flip (row sum kept near +-1 by penalty).
Objective F = sum_{k=1}^{L-1} (PAF_u(k)+PAF_v(k)+2)^2 + Wr*((sumu)^2-1)^2 + Wr*((sumv)^2-1)^2.
F == 0 with row sums +-1  <=>  Legendre pair.
"""
import numpy as np, sys

def search_lp(L, seed, iters, Wr=8):
    rng = np.random.default_rng(seed)
    u = rng.choice([-1,1], size=L).astype(np.int64)
    v = rng.choice([-1,1], size=L).astype(np.int64)
    def full_paf(x):
        return np.array([int(np.dot(x, np.roll(x,-k))) for k in range(L)], dtype=np.int64)
    pu = full_paf(u); pv = full_paf(v)
    def objective(pu,pv,su,sv):
        d = pu[1:]+pv[1:]+2
        return int(np.dot(d,d)) + Wr*(su*su-1)**2 + Wr*(sv*sv-1)**2
    su = int(u.sum()); sv = int(v.sum())
    F = objective(pu,pv,su,sv); best=F
    # precompute helper for incremental PAF on flip of position p in sequence x with paf px:
    def delta_paf(x, p):
        # returns array d where new_px(k) = px(k) + d(k); d(0)=0
        d = np.zeros(L, dtype=np.int64)
        xp = x[p]
        for k in range(1,L):
            d[k] = -2*xp*(x[(p+k)%L] + x[(p-k)%L])
        return d
    T = max(1.0, F/ (4*L))
    for it in range(iters):
        T *= 0.99998
        if T < 0.02: T = 0.02
        which = rng.integers(2)
        x, px, s = (u,pu,su) if which==0 else (v,pv,sv)
        p = int(rng.integers(L))
        d = delta_paf(x, p)
        npx = px + d
        ns = s - 2*x[p]
        if which==0:
            newF = objective(npx, pv, ns, sv)
        else:
            newF = objective(pu, npx, su, ns)
        if newF <= F or rng.random() < np.exp((F-newF)/T):
            x[p] = -x[p]
            px[:] = npx
            if which==0: su = ns
            else: sv = ns
            F = newF
            if F < best:
                best = F
                if best == 0:
                    return u.copy(), v.copy(), 0, it
    return u, v, best, iters

if __name__ == "__main__":
    lengths = [int(a) for a in sys.argv[1:]] or [45,63]
    import json
    for L in lengths:
        found=False
        for seed in range(40):
            u,v,best,it = search_lp(L, seed, iters=120000)
            if best==0:
                print(f"LP({L}) FOUND  seed={seed} iters={it}")
                np.save(f"deepdive_lp333/lp{L}_u.npy", u)
                np.save(f"deepdive_lp333/lp{L}_v.npy", v)
                # verify
                pu=np.array([int(np.dot(u,np.roll(u,-k))) for k in range(L)])
                pv=np.array([int(np.dot(v,np.roll(v,-k))) for k in range(L)])
                ok = all(pu[k]+pv[k]==-2 for k in range(1,L)) and abs(int(u.sum()))==1 and abs(int(v.sum()))==1
                print(f"   verified LP: {ok}  rowsums=({int(u.sum())},{int(v.sum())})")
                found=True; break
        if not found:
            print(f"LP({L}) not found in budget (best stalled)")
