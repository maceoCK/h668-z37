"""
Experiment 1 (incremental cube-and-conquer pilot): CNF emitter for the H-cyclotomic coset system.

Emits a DIMACS CNF for the order-3 (or general H) cyclotomic lift of one skeleton:
  vars: y[i][j][c] for i<=j (coset c in D[i][j]); product vars p = y1 AND y2 (3 clauses each);
  totalizer cardinality for block sizes; the 13-class convolution equations as exact PB-equalities
  encoded via pysat ITotalizer-free adder... we use CardEnc.equals on weighted-expanded literals
  (weights <= 3 expanded by literal duplication through auxiliary equality vars).
Pinning a=4 and z-flags forced, mirroring the VALIDATED cyclo_cpsat model (same tables, same
equation set; the emitter is regression-tested by anchors: P(9)/P(25) instances must be SAT and
decode to verified srgs; brute-force-confirmed-empty instances must be UNSAT).

Outputs: <out>.cnf, <out>.map (json: var assignments), and the cube variable list (the 54 diagonal
pair-bits in a canonical order) for the assumption-based driver.
"""
import sys, json
import numpy as np
from pysat.formula import CNF, IDPool
from pysat.card import CardEnc, EncType
from cyclo_cpsat import cosets_of

def emit(p, norb, k, lam, mu, mult, R, a_pin, out):
    H, CS = cosets_of(mult, p)
    nc, e = len(CS), len(H)
    cidx = {g: ci for ci, c in enumerate(CS) for g in c}
    negc = [cidx[(-c[0]) % p] for c in CS]
    QR = {(x * x) % p for x in range(1, p)}
    chi = [1 if c[0] in QR else -1 for c in CS]
    cnt0 = [[e if negc[ci] == cj else 0 for cj in range(nc)] for ci in range(nc)]
    cntC = [[[0] * nc for _ in range(nc)] for _ in range(nc)]
    for ci in range(nc):
        for cj in range(nc):
            t = [0] * p
            for x in CS[ci]:
                for y in CS[cj]:
                    t[(x + y) % p] += 1
            for d in range(nc):
                cntC[ci][cj][d] = t[CS[d][0]]

    pool = IDPool()
    cnf = CNF()
    Y = {}
    Z = {}
    for i in range(norb):
        for j in range(i, norb):
            r = int(R[i][j])
            if i == j:
                need = 2 * e if e % 2 == 1 else e
                assert r % need == 0, "incompatible diagonal"
                Z[i, j] = 0
            else:
                assert r % e in (0, 1), "incompatible block"
                Z[i, j] = r % e
            for c in range(nc):
                Y[i, j, c] = pool.id(f"y{i}_{j}_{c}")
            lits = [Y[i, j, c] for c in range(nc)]
            want = (r - Z[i, j]) // e
            cc = CardEnc.equals(lits=lits, bound=want, vpool=pool, encoding=EncType.totalizer)
            cnf.extend(cc.clauses)
    def y(i, j, c):
        return Y[i, j, c] if i <= j else Y[j, i, negc[c]]
    def z(i, j):
        return Z[(i, j) if i <= j else (j, i)]
    # diagonal symmetry + pinning
    for i in range(norb):
        for c in range(nc):
            a, b = Y[i, i, c], Y[i, i, negc[c]]
            if a != b:
                cnf.append([-a, b]); cnf.append([a, -b])
    if a_pin is not None:
        for c in range(nc):
            t = 3 * a_pin - 9 if chi[c] == 1 else 18 - 3 * a_pin
            lits = [Y[i, i, c] for i in range(norb)]
            cc = CardEnc.equals(lits=lits, bound=t, vpool=pool, encoding=EncType.totalizer)
            cnf.extend(cc.clauses)
    # product vars
    prod = {}
    def AND(va, vb):
        if va > vb: va, vb = vb, va
        if (va, vb) not in prod:
            t = pool.id(f"p{va}_{vb}")
            cnf.append([-t, va]); cnf.append([-t, vb]); cnf.append([t, -va, -vb])
            prod[(va, vb)] = t
        return prod[(va, vb)]
    # equations: weighted sums == rhs; expand weight-w term as w copies (w <= e = 3)
    neq = 0
    for i in range(norb):
        for j in range(i, norb):
            for d in list(range(nc)) + ["zero"]:
                wlits = []                      # (lit, weight)
                const = 0
                for l in range(norb):
                    zil, zlj = z(i, l), z(l, j)
                    for c in range(nc):
                        for cq in range(nc):
                            w = cnt0[c][cq] if d == "zero" else cntC[c][cq][d]
                            if w:
                                wlits.append((AND(y(i, l, c), y(l, j, cq)), w))
                    if d == "zero":
                        const += zil * zlj
                    else:
                        if zil: wlits.append((y(l, j, d), 1))
                        if zlj: wlits.append((y(i, l, d), 1))
                rhs = mu + ((k - mu) if (i == j and d == "zero") else 0) - const
                if d == "zero":
                    rhs += (lam - mu) * z(i, j)   # the indicator [0 in D] is the constant z
                else:
                    wlits.append(((mu - lam) * 0 + y(i, j, d), (mu - lam)))
                total = sum(w for _, w in wlits)
                if rhs < 0 or rhs > total:
                    f = pool.id("false_aux")    # trivially false, loader-safe
                    cnf.append([f]); cnf.append([-f])
                    continue
                from pysat.pb import PBEnc, EncType as PBT
                cc = PBEnc.equals(lits=[l for l, _ in wlits], weights=[w for _, w in wlits],
                                  bound=rhs, vpool=pool, encoding=PBT.adder)
                cnf.extend(cc.clauses)
                neq += 1
    cubevars = [Y[i, i, c] for i in range(norb) for c in range(nc)]
    cnf.to_file(out + ".cnf")
    json.dump({"y": {f"{i},{j},{c}": Y[i, j, c] for (i, j, c) in Y},
               "z": {f"{i},{j}": Z[i, j] for (i, j) in Z},
               "cosets": CS, "cubevars": cubevars, "nvars": pool.top, "neq": neq},
              open(out + ".map", "w"))
    print(f"emitted {out}.cnf: {pool.top} vars, {len(cnf.clauses)} clauses, {neq} equations, {len(prod)} products")

if __name__ == "__main__":
    which = sys.argv[1]
    if which == "anchor9":
        import conf_core as cc
        D = cc.paley_pp_blocks(3)
        R = [[len(D[i][j]) for j in range(3)] for i in range(3)]
        emit(3, 3, 4, 1, 2, 2, R, None, "anchor9")
    elif which == "anchor25":
        import conf_core as cc
        D = cc.paley_pp_blocks(5)
        R = [[len(D[i][j]) for j in range(5)] for i in range(5)]
        emit(5, 5, 12, 5, 6, 4, R, None, "anchor25")
    elif which == "anchor_unsat":
        emit(7, 3, 10, 4, 5, 6, [[2, 4, 4], [4, 2, 4], [4, 4, 2]], None, "anchor_unsat")
    else:
        idx = int(which)
        mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
                for ln in open("order3_skeletons.txt")]
        R = [[int(x) for x in row] for row in mats[idx]]
        emit(37, 9, 166, 82, 83, 10, R, 4, f"order3_sk{idx}")
