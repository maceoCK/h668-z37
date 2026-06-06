"""
Validate multiplier_search against ground truth:
 (1) brute-force ALL Legendre pairs at small ell (via PAF bucketing), count them;
 (2) check which are invariant under SOME cyclic multiplier subgroup;
 (3) confirm multiplier_search finds those that are;
 (4) explicit Legendre (quadratic-residue) construction at prime ell -> guaranteed symmetric LP.
"""
import numpy as np
from math import gcd
from itertools import combinations, product
import multiplier_search as ms

def all_lps_bruteforce(ell):
    half = ell//2
    # enumerate u with row sum +-1, bucket by PAF profile (lags 1..half)
    from collections import defaultdict
    buckets = defaultdict(list)
    seqs = []
    idxs = list(range(ell))
    for a in ((ell-1)//2, (ell+1)//2):   # number of -1 positions so sum = +-1
        for neg in combinations(idxs, a):
            u = np.ones(ell, dtype=np.int64);
            for i in neg: u[i] = -1
            pf,_ = ms.paf_fft(u, ell)
            prof = tuple(int(pf[k]) for k in range(1, half+1))
            buckets[prof].append(u.copy())
    lps = []
    for prof, us in buckets.items():
        comp = tuple(-2-x for x in prof)
        if comp in buckets:
            # any u in prof, v in comp is an LP (avoid double count by prof<=comp)
            if prof <= comp:
                for u in us:
                    for v in buckets[comp]:
                        lps.append((u, v))
    return lps

def is_multiplier_invariant(u, ell):
    """return a generator t such that u is constant on <t>-orbits (t in Z_ell*, t!=1), else None."""
    for t in ms.units(ell):
        if t == 1: continue
        ok = all(u[(i*t) % ell] == u[i] for i in range(ell))
        if ok: return t
    return None

def legendre_pair(p):
    qr = set((i*i) % p for i in range(1, p))
    base = np.array([1] + [1 if i in qr else -1 for i in range(1, p)], dtype=np.int64)
    u = base.copy()
    v = base.copy(); v[0] = -1
    return u, v

if __name__ == "__main__":
    print("=== (4) explicit Legendre construction at prime lengths ===")
    for p in (7, 11, 13, 17, 19, 23):
        u, v = legendre_pair(p)
        print(f"  p={p}: legendre_pair is LP? {ms.verify_lp(u,v,p)}  "
              f"u multiplier-invariant under t={is_multiplier_invariant(u,p)}")

    print("\n=== (1)-(3) brute force + symmetry check + does multiplier_search find it? ===")
    for ell in (7, 9, 11, 13, 15):
        lps = all_lps_bruteforce(ell)
        sym = [(u,v) for (u,v) in lps if is_multiplier_invariant(u,ell) and is_multiplier_invariant(v,ell)]
        print(f"  ell={ell}: total LP pairs (rowsum+-1) = {len(lps)};  "
              f"both-multiplier-invariant = {len(sym)}")
        # run multiplier_search
        found = ms.attack(ell, max_log=18, time_budget=120)
        print(f"     multiplier_search found one? {found}  "
              f"(expected: {'YES' if len(sym)>0 else 'maybe-no (none symmetric)'})")
        print()
