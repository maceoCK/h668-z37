"""
Probe: F-statistics inside the trace-pinned subspace (U = QR, per-orbit coupling)
vs the unrestricted sigma-equivariant space, class 221_0. Random sampling + greedy descent.
"""
import random, sys
import numpy as np
sys.path.insert(0, "/Users/maceocardinalekwik/git/math/release/h668-z37/code")
from sc_floor_analysis import Tpl, CODE

P = 37
QRs = sorted({(k*k) % P for k in range(1, P)})
def pinned_state(tpl, rng):
    repBit = np.zeros((tpl.NREP, P), dtype=np.int64)
    # orbits with QR-half reps
    seen = [False]*P; orbs = []
    for x in range(1, P):
        if seen[x]: continue
        d = x if x in QRs else (6*x) % P
        ob = (d, P-d, (6*d) % P, P-(6*d) % P)   # (QR half, QNR half)
        for z in ob: seen[z] = True
        orbs.append(ob)
    # rep indices: 0=(0,0) S card12 type1; 8=(2,2) T card18 type1; 12=(7,7) U type2
    # case U=QR, eps=-1: per orbit type A (s=(0,1), t1=t2) or B (t=(0,1), s1=s2)
    nA = rng.choice([0,2,4,6])
    beta = (6-nA)//2; alpha = nA//2
    idx = list(range(9)); rng.shuffle(idx)
    Aorbs = set(idx[:nA]); Borbs = [i for i in idx[nA:]]
    bothS = set(rng.sample(Borbs, beta))
    bothT = set(rng.sample(sorted(Aorbs), alpha))
    for i, (d, nd, d6, nd6) in enumerate(orbs):
        repBit[12, d] = 1; repBit[12, nd] = 1                     # U = QR
        if i in Aorbs:
            repBit[0, d6] = 1; repBit[0, nd6] = 1                 # s=(0,1)
            if i in bothT:
                for z in (d, nd, d6, nd6): repBit[8, z] = 1       # t both
        else:
            repBit[8, d6] = 1; repBit[8, nd6] = 1                 # t=(0,1)
            if i in bothS:
                for z in (d, nd, d6, nd6): repBit[0, z] = 1       # s both
    assert repBit[0].sum() == 12 and repBit[8].sum() == 18 and repBit[12].sum() == 18
    # off-diagonal type-0 reps random
    for r, (i, j, typ, card) in enumerate(tpl.reps):
        if typ == 0:
            el = rng.sample(range(P), card)
            repBit[r, el] = 1
    return repBit

def greedy(tpl, repBit, rng, iters=4000):
    best = tpl.obj_of(tpl.residuals(repBit))
    t0reps = [r for r, (i,j,typ,c) in enumerate(tpl.reps) if typ == 0]
    for _ in range(iters):
        r = rng.choice(t0reps)
        ins = np.nonzero(repBit[r])[0]; outs = np.nonzero(repBit[r]==0)[0]
        a = int(rng.choice(ins)); b = int(rng.choice(outs))
        repBit[r,a] = 0; repBit[r,b] = 1
        o = tpl.obj_of(tpl.residuals(repBit))
        if o <= best: best = o
        else: repBit[r,a] = 1; repBit[r,b] = 0
    return best

def main():
    tpl = Tpl(f"{CODE}/sc_tpl_221_0.txt")
    rng = random.Random(2024)
    Fp = [tpl.obj_of(tpl.residuals(pinned_state(tpl, rng))) for _ in range(150)]
    Fu = [tpl.obj_of(tpl.residuals(tpl.random_state(rng))) for _ in range(150)]
    print(f"random F  pinned subspace: min={min(Fp)} mean={np.mean(Fp):.0f}")
    print(f"random F  full space     : min={min(Fu)} mean={np.mean(Fu):.0f}")
    bp = [greedy(tpl, pinned_state(tpl, rng), rng, 2500) for _ in range(6)]
    bu = [greedy(tpl, tpl.random_state(rng), rng, 2500) for _ in range(6)]
    print(f"greedy-descent floors (6 restarts, off-diag moves only):")
    print(f"  pinned : {sorted(bp)}")
    print(f"  full   : {sorted(bu)}")

if __name__ == "__main__":
    main()
