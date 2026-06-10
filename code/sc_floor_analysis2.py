"""
Stage 2: identify the extra affine invariants beyond the 70 structural ones;
per-equation residual statistics; quadratic-statistic invariant hunt mod m;
Fourier/defect identity check (F_full vs character matrices).
"""
import json, random, sys
import numpy as np
sys.path.insert(0, "/Users/maceocardinalekwik/git/math/release/h668-z37/code")
from sc_floor_analysis import Tpl, parse_dump, CODE, RESID, gfp_nullspace

PBIG = 2147483647  # prime for exact mod-p linear algebra

def modp_nullspace(D, p=PBIG):
    return gfp_nullspace(D % p, p)

def main():
    data = json.load(open(RESID))
    tpl = Tpl(f"{CODE}/sc_tpl_221_0.txt")
    eqs = tpl.rep_eq_list()
    eqidx = {e:i for i,e in enumerate(eqs)}
    NS = 1200
    rng = random.Random(4242)
    R = np.zeros((NS, len(eqs)), dtype=np.int64)
    for s in range(NS):
        res = tpl.residuals(tpl.random_state(rng))
        R[s] = [res[a,b,g] for (a,b,g) in eqs]
    D = R[1:] - R[0]

    # ---- known structural invariant space (70-dim) ----
    known = []
    P = tpl.P
    reppairs = sorted(set((a,b) for (a,b,g) in eqs))
    for (a,b) in reppairs:
        v = np.zeros(len(eqs));
        for g in range(P): v[eqidx[(a,b,g)]] = 1
        known.append(v)                                  # row sums (13)
        if a==b:
            v = np.zeros(len(eqs)); v[eqidx[(a,a,0)]] = 1
            known.append(v)                              # delta at 0 (3)
            for g in range(1, P//2+1):
                v = np.zeros(len(eqs))
                v[eqidx[(a,a,g)]] = 1; v[eqidx[(a,a,P-g)]] -= 1
                known.append(v)                          # evenness (3*18)
    known = np.array(known)
    print(f"known structural invariants: {known.shape[0]}, "
          f"rank={np.linalg.matrix_rank(known)}")

    # ---- full invariant space via exact mod-p nullspace of D ----
    ns = modp_nullspace(D.astype(np.int64), PBIG)        # basis vectors over F_p
    print(f"nullspace dim over F_p (= #affine invariants): {ns.shape[0]}")
    # verify each on independent fresh samples over the integers? (they're mod-p vectors;
    # lift small-entry representatives: entries are in [0,p); map to signed)
    nss = ns.astype(np.int64).copy()
    nss[nss > PBIG//2] -= PBIG
    # which are independent of 'known'? compute rank growth over Q (floats fine, entries small)
    base = np.linalg.matrix_rank(known)
    extra_vecs = []
    cur = known.astype(float)
    for v in nss:
        test = np.vstack([cur, v.astype(float)])
        if np.linalg.matrix_rank(test) > cur.shape[0] - (cur.shape[0]-np.linalg.matrix_rank(cur)):
            pass
        if np.linalg.matrix_rank(test) > np.linalg.matrix_rank(cur):
            extra_vecs.append(v)
            cur = test
    print(f"extra invariants beyond structural: {len(extra_vecs)}")
    # verify extras exactly on fresh integer samples and print supports
    rng2 = random.Random(777)
    Rf = np.zeros((40, len(eqs)), dtype=np.int64)
    for s in range(40):
        res = tpl.residuals(tpl.random_state(rng2))
        Rf[s] = [res[a,b,g] for (a,b,g) in eqs]
    for t,v in enumerate(extra_vecs):
        vals = Rf @ v
        const = set(vals.tolist())
        sup = [(eqs[i], int(v[i])) for i in np.nonzero(v)[0]]
        pairs_in = sorted(set((e[0][0],e[0][1]) for e in sup))
        print(f"  extra[{t}]: support size={len(sup)}, pairs={pairs_in}, "
              f"verified-constant={'YES val='+str(vals[0]) if len(const)==1 else 'NO (mod-p artifact)'}")
        if len(sup) <= 40:
            print(f"    {sup}")

    # ---- per-equation residual stats over samples + best states ----
    print("\n== per-equation |residual| floor over random samples ==")
    absR = np.abs(R)
    minabs = absR.min(axis=0)
    nevzero = np.nonzero(minabs > 0)[0]
    print(f"equations never hitting 0 in {NS} random samples: {len(nevzero)}")
    if len(nevzero):
        for i in nevzero[:20]:
            print(f"  eq {eqs[i]}: min|r|={minabs[i]}, values seen mod 2: {sorted(set((R[:,i]%2).tolist()))}")

    # ---- quadratic statistic invariant hunt ----
    print("\n== quadratic statistics mod m (constant => candidate invariant) ==")
    stats = {}
    for (a,b) in reppairs:
        idx = [eqidx[(a,b,g)] for g in range(P)]
        Rb = R[:, idx]
        stats[f"S2({a},{b})"] = (Rb**2).sum(axis=1)
        for s in (1,2,3):
            stats[f"AC{s}({a},{b})"] = (Rb * np.roll(Rb, s, axis=1)).sum(axis=1)
    stats["F"] = sum(stats[f"S2({a},{b})"] for (a,b) in reppairs)
    found = 0
    for name, vals in stats.items():
        if len(set(vals.tolist())) == 1:
            print(f"  {name}: CONSTANT = {vals[0]} over {NS} samples")
            found += 1
            continue
        for m in (2,3,4,5,8,9,37):
            vm = set((vals % m).tolist())
            if len(vm) == 1:
                v0 = vals[0] % m
                tag = " *** NONZERO ***" if v0 != 0 else ""
                print(f"  {name} mod {m}: CONSTANT = {v0}{tag}")
                found += 1
    if not found: print("  none found")

    # ---- defect / Fourier identity check on a best-state dump ----
    print("\n== Fourier defect identity: 37*F_full(ordered) = sum_chi!=1 ||M(chi)^2+M(chi)-83I||_F^2 ==")
    ent = data[0]
    obj_claim, repBit, _ = parse_dump(ent["dump"], tpl)
    mem = tpl.expand(repBit)
    res = tpl.residuals(repBit)
    F_ord = float((res**2).sum())            # ordered all-pairs
    w = np.exp(2j*np.pi/37)
    tot = 0.0
    eigsq = []
    for j in range(1,37):
        ch = w**(np.arange(37)*j)
        M = mem @ ch                           # (9,9) complex: sum_d mem[a,b,d] ch[d]
        Eb = M@M + M - 83*np.eye(9)
        tot += float(np.linalg.norm(Eb,'fro')**2)
        th = np.linalg.eigvalsh(M)            # M Hermitian
        eigsq.append(th)
    print(f"  37*F_full_ordered = {37*F_ord:.1f}   sum_chi ||.||^2 = {tot:.1f}   "
          f"match={abs(37*F_ord-tot)<1e-3}")
    th = np.array(eigsq).ravel()
    q = th**2+th-83
    print(f"  eigenvalue displacement at best state: max|q(theta)| = {np.abs(q).max():.3f}, "
          f"mean|q| = {np.abs(q).mean():.3f}; theta range [{th.min():.3f},{th.max():.3f}]; "
          f"targets {(-1+3*np.sqrt(37))/2:.3f}/{(-1-3*np.sqrt(37))/2:.3f}")
    # sanity: F restricted to rep pairs * orbit multiplicities == F_full? check sigma-tie claim:
    # count ordered F vs 4x rep F
    F_rep = tpl.obj_of(res)
    print(f"  F_rep(dump) = {F_rep}, F_full_ordered = {F_ord:.0f}, ratio = {F_ord/F_rep:.3f}")

if __name__=="__main__":
    main()
