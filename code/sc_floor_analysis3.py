"""
Stage 2b: extended nonlinear-statistic invariant hunt (cubic, cross-block correlations,
full autocorrelation profile) mod m, plus exact verification of the complement-invariance
identity E(A-bar) = E(A) at the 333-vertex level on a reconstructed dump state.
"""
import json, random, sys
import numpy as np
sys.path.insert(0, "/Users/maceocardinalekwik/git/math/release/h668-z37/code")
from sc_floor_analysis import Tpl, parse_dump, CODE, RESID

def main():
    data = json.load(open(RESID))
    tpl = Tpl(f"{CODE}/sc_tpl_221_0.txt")
    P, N = tpl.P, tpl.NORB
    eqs = tpl.rep_eq_list()
    eqidx = {e:i for i,e in enumerate(eqs)}
    reppairs = sorted(set((a,b) for (a,b,g) in eqs))
    NS = 800
    rng = random.Random(31337)
    R = np.zeros((NS, len(eqs)), dtype=np.int64)
    for s in range(NS):
        res = tpl.residuals(tpl.random_state(rng))
        R[s] = [res[a,b,g] for (a,b,g) in eqs]

    print("== cubic / cross-block / full-AC statistic hunt mod m ==")
    blocks = {}
    for (a,b) in reppairs:
        blocks[(a,b)] = R[:, [eqidx[(a,b,g)] for g in range(P)]]
    stats = {}
    for (a,b) in reppairs:
        Rb = blocks[(a,b)]
        stats[f"S3({a},{b})"] = (Rb**3).sum(axis=1)
        stats[f"S4({a},{b})"] = (Rb**4).sum(axis=1)
        for s in range(1, P//2+1):
            stats[f"AC{s}({a},{b})"] = (Rb*np.roll(Rb,s,axis=1)).sum(axis=1)
    for i,(a,b) in enumerate(reppairs):
        for (a2,b2) in reppairs[i+1:]:
            stats[f"X({a},{b};{a2},{b2})"] = (blocks[(a,b)]*blocks[(a2,b2)]).sum(axis=1)
    hits = []
    for name, vals in stats.items():
        if len(set(vals.tolist()))==1:
            hits.append((name, "Z", int(vals[0])))
            continue
        for m in (2,3,4,5,7,8,9,37):
            if len(set((vals%m).tolist()))==1:
                v0 = int(vals[0]%m)
                hits.append((name, f"mod {m}", v0))
    nonzero = [h for h in hits if h[2]!=0]
    zero = [h for h in hits if h[2]==0]
    print(f"  constant statistics found: {len(hits)} ({len(zero)} with value 0)")
    for h in nonzero:
        print(f"  *** NONZERO CONSTANT: {h} ***")
    if not nonzero: print("  no nonzero-constant statistic — no obstruction from this family")
    # summarize the zero ones briefly
    from collections import Counter
    cnt = Counter(h[1] for h in zero)
    print(f"  zero-constant breakdown: {dict(cnt)}")

    print("\n== full 333-vertex check: E(A) row sums, symmetry, complement invariance ==")
    obj_claim, repBit, _ = parse_dump(data[0]["dump"], tpl)
    mem = tpl.expand(repBit)
    # build 333x333 adjacency: vertex (i, t) -> index 37*i+t ; edge (i,s)~(j,t) iff (t-s) in D_ij
    A = np.zeros((333,333), dtype=np.int64)
    for i in range(N):
        for j in range(N):
            for s in range(P):
                # row (i,s): neighbor (j, s+d) for d in D_ij
                ds = np.nonzero(mem[i,j])[0]
                A[37*i+s, 37*j + ((s+ds)%P)] = 1
    print(f"  A symmetric: {(A==A.T).all()}, diag zero: {(np.diag(A)==0).all()}, "
          f"row sums: {set(A.sum(axis=1).tolist())}")
    J = np.ones((333,333),dtype=np.int64); I = np.eye(333,dtype=np.int64)
    E = A@A + A - 83*I - 83*J
    print(f"  E row sums all zero: {(E.sum(axis=1)==0).all()}")
    Ab = J - I - A
    Eb = Ab@Ab + Ab - 83*I - 83*J
    print(f"  complement invariance E(J-I-A) == E(A): {(Eb==E).all()}")
    print(f"  ||E||_F^2 = {int((E**2).sum())}  (= 37 * F_full_ordered = {37*int((tpl.residuals(repBit)**2).sum())})")
    # objective F is always even: verify on dumps
    objs = []
    for ent in data:
        t2 = Tpl(f"{CODE}/{ent['tpl']}") if ent['tpl']!="sc_tpl_221_0.txt" else tpl
        _, rb, _ = parse_dump(ent["dump"], t2)
        objs.append(t2.obj_of(t2.residuals(rb)))
    print(f"  dump objectives: {objs} ; all even: {all(o%2==0 for o in objs)}")

if __name__=="__main__":
    main()
