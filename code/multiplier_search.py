"""
Multiplier-orbit + PSD search for Legendre pairs (the Kotsireas-Koutschan method).
A sequence invariant under a multiplier subgroup M <= Z_ell* is constant on the orbits
of M acting on Z_ell by multiplication; search 2^(#orbits) sign assignments.

LP condition: PAF_u(k)+PAF_v(k) = -2 for all k!=0  (row sums +-1 are then IMPLIED).
Match u with v by complementary PAF profile.  All arithmetic verified exactly.
"""
import numpy as np
from math import gcd
from itertools import product
import sys, time, os

OUTDIR = os.path.dirname(os.path.abspath(__file__))

def units(ell): return [a for a in range(1, ell) if gcd(a, ell) == 1]

def mult_order(t, ell):
    o, x = 1, t % ell
    while x != 1:
        x = (x*t) % ell; o += 1
    return o

def group_orbits(ell, gens):
    seen = [False]*ell; orbits = []
    for i in range(ell):
        if seen[i]: continue
        orb = [i]; seen[i] = True; fr = [i]
        while fr:
            x = fr.pop()
            for g in gens:
                y = (x*g) % ell
                if not seen[y]: seen[y] = True; orb.append(y); fr.append(y)
        orbits.append(orb)
    return orbits

def paf_fft(u, ell):
    F = np.fft.fft(u)
    psd = (F*np.conj(F)).real
    ac = np.fft.ifft(psd).real
    return np.rint(ac).astype(np.int64), np.rint(psd).astype(np.int64)

def verify_lp(u, v, ell):
    pu = np.array([int(np.dot(u, np.roll(u, -k))) for k in range(ell)])
    pv = np.array([int(np.dot(v, np.roll(v, -k))) for k in range(ell)])
    return all(pu[k]+pv[k] == -2 for k in range(1, ell)) and abs(int(u.sum())) == 1 and abs(int(v.sum())) == 1

def search_family(ell, gens, max_log=22):
    orbits = group_orbits(ell, gens)
    n = len(orbits)
    if n > max_log or n == 0:
        return None, n
    sizes = np.array([len(o) for o in orbits], dtype=np.int64)
    idx_lists = [np.array(o, dtype=np.int64) for o in orbits]
    bound = 2*ell + 2
    half = ell//2
    buckets = {}
    u = np.empty(ell, dtype=np.int64)
    for bits in product((1, -1), repeat=n):
        rs = int(sizes @ np.array(bits))
        if rs*rs != 1:           # row sum must be +-1
            continue
        for val, ol in zip(bits, idx_lists):
            u[ol] = val
        pf, psd = paf_fft(u, ell)
        if psd.max() > bound:    # PSD test: PSD_u(s) <= 2ell+2
            continue
        prof = tuple(int(pf[k]) for k in range(1, half+1))
        comp = tuple(-2 - x for x in prof)
        if prof == comp:                       # self-complementary: u=v is an LP
            return (u.copy(), u.copy()), n
        if comp in buckets:
            return (u.copy(), buckets[comp].copy()), n
        if prof not in buckets:
            buckets[prof] = u.copy()
    return None, n

def candidate_groups(ell, max_log=22, min_log=1):
    """yield distinct multiplier subgroups <t> and <t,-1> with #orbits in [min_log,max_log]."""
    seen_part = set()
    out = []
    us = units(ell)
    for t in us:
        for gens in ([t], [t, ell-1]):
            orbits = group_orbits(ell, gens)
            n = len(orbits)
            if not (min_log <= n <= max_log):
                continue
            key = frozenset(frozenset(o) for o in orbits)
            if key in seen_part: continue
            seen_part.add(key)
            out.append((gens, n))
    out.sort(key=lambda x: x[1])   # fewest orbits first (cheapest)
    return out

def attack(ell, max_log=22, time_budget=600):
    t0 = time.time()
    groups = candidate_groups(ell, max_log=max_log)
    print(f"ell={ell}: {len(groups)} multiplier symmetry classes with <= {max_log} orbits "
          f"(orbit counts: {sorted(set(n for _,n in groups))})")
    tried = 0
    for gens, n in groups:
        if time.time()-t0 > time_budget:
            print(f"  [time budget {time_budget}s hit after {tried} classes]"); break
        res, _ = search_family(ell, gens, max_log=max_log)
        tried += 1
        if res is not None:
            u, v = res
            ok = verify_lp(u, v, ell)
            print(f"  *** LP({ell}) FOUND  gens={gens} (#orbits={n})  verified={ok} ***")
            if ok:
                np.save(os.path.join(OUTDIR, f"LP{ell}_u.npy"), u)
                np.save(os.path.join(OUTDIR, f"LP{ell}_v.npy"), v)
                print(f"     u={list(map(int,u))}")
                print(f"     v={list(map(int,v))}")
                return True
    print(f"  no multiplier-invariant LP({ell}) in {tried} symmetry classes searched "
          f"({time.time()-t0:.1f}s)")
    return False

if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]]
    lens = args[:-1] if len(args) >= 2 else args
    ml = args[-1] if len(args) >= 2 else 20
    if not lens: lens = [45]
    for ell in lens:
        attack(ell, max_log=ml)
        print()
