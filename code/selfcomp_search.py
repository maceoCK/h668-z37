"""
Self-complementary Z_p-symmetric SRG search (sigma-equivariant CP-SAT model).

Structure: norb orbits x Z_p; adjacency blocks D[i][j] (conf_core conventions); a complementing
isomorphism sigma(i,g) = (pi(i), c*g) with c of multiplicative order 4 mod p (b = 0 WLOG, verified
by selfcomp_twist.py: the coordinate gauge acts transitively on twist solutions).  Block relations:

    D[pi(i)][pi(j)] = co( c * D[i][j] ),       co = complement in Z_p (off-diag) / Z_p^* (diag),

so blocks along a pi-orbit of unordered pairs are alternating scaled complements of the orbit
representative; the fixed-orbit diagonal satisfies the self-referential x(c d) = NOT x(d).

Free variables: one block per pair-orbit (plus symmetry/self-reference constraints), sizes pinned
by the (complement-symmetric) orbit matrix R.  SRG equations at element level, reduced to one
representative per sigma-orbit of pairs (sound: for self-complementary parameters the complement's
equations are literally equivalent).

A solution => srg(v,k,lam,mu) that is self-complementary with a normalized semiregular Z_p.
For production (p=37, norb=9): a hit IS a conference graph for H(668).

Validation (run with --validate):
  V1: p=13, norb=1, c=5: must find Paley(13) (QR is <5>-anticlosed); every solution verifies.
  V2: p=17, norb=1, c=4: Paley(17) likewise.
  V3: brute-force cross-check, p=5, norb=5, pi=(0123)(4), c=2: engine solution set == definitional
      enumeration over ALL sigma-equivariant assignments (18 free bits), for several R choices.
"""
import sys, itertools
import numpy as np
from ortools.sat.python import cp_model

def build_orbits(p, norb, pi, c):
    """pair-orbits of (i,j) i<=j under (i,j)->(pi i, pi j); returns resolution data."""
    cinv = pow(c, -1, p)
    seen = {}
    orbits = []                      # list of reps (i,j) with i<=j
    member = {}                      # ordered-sorted (a,b) -> (rep, k, swapped)
    for i in range(norb):
        for j in range(i, norb):
            if (i, j) in member: continue
            rep = (i, j)
            orbits.append(rep)
            u, v, k = i, j, 0
            while True:
                a, b = (u, v) if u <= v else (v, u)
                if (a, b) not in member:
                    member[(a, b)] = (rep, k, u > v)
                u, v, k = pi[u], pi[v], k + 1
                if (u, v) == (i, j): break
    return orbits, member, cinv

class Model:
    def __init__(self, p, norb, k, lam, mu, R, pi, c):
        self.p, self.norb, self.k, self.lam, self.mu = p, norb, k, lam, mu
        self.R, self.pi, self.c = R, pi, c
        self.orbits, self.member, self.cinv = build_orbits(p, norb, pi, c)
        self.m = cp_model.CpModel()
        self.X = {}                  # (rep, d) -> BoolVar for free rep blocks
        self.prod = {}
        self.ok = self._build()

    def lit(self, a, b, d):
        """membership literal of d in D[a][b]; returns (var|None, neg, const) with const in {0,1,None}."""
        p = self.p
        if a > b:
            return self.lit(b, a, (-d) % p)
        rep, kk, swapped = self.member[(a, b)]
        e = d
        if swapped: e = (-e) % p
        e = (e * pow(self.cinv, kk, p)) % p
        neg = (kk % 2) == 1
        i, j = rep
        if i == j and e == 0:
            # 0 not in any diagonal block: membership False; negated -> True
            return (None, False, 1 if neg else 0)
        v = self.X[(rep, e)]
        return (v, neg, None)

    def product(self, l1, l2):
        """linear expression for AND of two literals (var,neg,const)."""
        v1, n1, c1 = l1
        v2, n2, c2 = l2
        if c1 is not None:
            if c1 == 0: return 0
            return self.expr(l2)
        if c2 is not None:
            if c2 == 0: return 0
            return self.expr(l1)
        # expand (a1 + s1 x)(a2 + s2 y) with literal = x or (1-x)
        a1, s1 = (1, -1) if n1 else (0, 1)
        a2, s2 = (1, -1) if n2 else (0, 1)
        if v1.Index() == v2.Index():
            xy = v1                          # x*x = x
        else:
            key = (min(v1.Index(), v2.Index()), max(v1.Index(), v2.Index()))
            if key not in self.prod:
                t = self.m.NewBoolVar(f"p{len(self.prod)}")
                self.m.Add(t <= v1); self.m.Add(t <= v2); self.m.Add(t >= v1 + v2 - 1)
                self.prod[key] = t
            xy = self.prod[key]
        return a1 * a2 + a1 * s2 * v2 + a2 * s1 * v1 + s1 * s2 * xy

    def expr(self, l):
        v, n, c = l
        if c is not None: return c
        return (1 - v) if n else v

    def _build(self):
        p, norb, R, pi, c = self.p, self.norb, self.R, self.pi, self.c
        m = self.m
        # consistency of sizes with the complement relation
        for i in range(norb):
            for j in range(norb):
                want = (36 if i == j else 37) if p == 37 else ((p - 1) if i == j else p)
                if int(R[pi[i]][pi[j]]) != want - int(R[i][j]):
                    return False
        # free variables per orbit rep
        for (i, j) in self.orbits:
            for d in range(p):
                if i == j and d == 0: continue
                self.X[((i, j), d)] = m.NewBoolVar(f"x{i}_{j}_{d}")
            if i == j:
                # symmetric, no 0
                for d in range(1, p):
                    m.Add(self.X[((i, j), d)] == self.X[((i, j), (p - d) % p)])
                # fixed orbit: self-referential anticlosure x(c d) = NOT x(d) when pi fixes i
                if pi[i] == i:
                    for d in range(1, p):
                        m.Add(self.X[((i, j), (c * d) % p)] + self.X[((i, j), d)] == 1)
                m.Add(sum(self.X[((i, j), d)] for d in range(1, p)) == int(R[i][i]))
            else:
                m.Add(sum(self.X[((i, j), d)] for d in range(p)) == int(R[i][j]))
        # SRG equations, one representative per sigma-orbit of pairs
        # eq (i,j,g): sum_l sum_d [d in D[i][l]][g-d in D[l][j]] = mu + (k-mu)[i=j,g=0] + (lam-mu)[g in D[i][j]]
        for (i, j) in self.orbits:
            for g in range(p):
                terms = []
                for l in range(norb):
                    for d in range(p):
                        t = self.product(self.lit(i, l, d), self.lit(l, j, (g - d) % p))
                        if not (isinstance(t, int) and t == 0):
                            terms.append(t)
                rhs = self.mu + ((self.k - self.mu) if (i == j and g == 0) else 0)
                ind = self.expr(self.lit(i, j, g))
                self.m.Add(sum(terms) == rhs + (self.lam - self.mu) * ind)
        return True

    def solve(self, time_s=600, workers=8, enumerate_all=False):
        sols = []
        sv = cp_model.CpSolver()
        sv.parameters.max_time_in_seconds = float(time_s)
        sv.parameters.num_search_workers = 1 if enumerate_all else workers
        if enumerate_all:
            sv.parameters.enumerate_all_solutions = True
            class Cb(cp_model.CpSolverSolutionCallback):
                def __init__(cb):
                    super().__init__()
                def on_solution_callback(cb):
                    sols.append(self.extract(cb.Value))
            st = sv.Solve(self.m, Cb())
        else:
            st = sv.Solve(self.m)
            if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                sols.append(self.extract(sv.Value))
        return sv.StatusName(st), sols

    def extract(self, val):
        p, norb = self.p, self.norb
        D = [[set() for _ in range(norb)] for _ in range(norb)]
        for a in range(norb):
            for b in range(norb):
                for d in range(p):
                    l = self.lit(a, b, d) if a <= b else self.lit(b, a, (-d) % p)
                    v, n, cst = l
                    mem = cst if cst is not None else (1 - val(v) if n else val(v))
                    if mem: D[a][b].add(d)
        return D

def verify_srg(D, p, norb, k, lam, mu):
    import conf_core as cc
    A = cc.build_adjacency(D, p, norb)
    n = p * norb
    lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(n, dtype=np.int64)
    return bool(np.array_equal(A, A.T) and np.array_equal(lhs, mu * np.ones((n, n), dtype=np.int64)))

def verify_selfcomp(D, p, norb, pi, c):
    for i in range(norb):
        for j in range(norb):
            img = {(c * d) % p for d in D[i][j]}
            uni = set(range(p)) if i != j else set(range(1, p))
            if D[pi[i]][pi[j]] != uni - img: return False
    return True

# ------------------------------------------------------------------ validation
def validate():
    ok_all = True
    def chk(name, cond, info=""):
        nonlocal ok_all
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {info}")
        if not cond: ok_all = False

    # V1/V2: Paley anchors, norb=1 (pi = identity on 1 orbit, c a NONSQUARE of even order:
    # p=13: c=5 (order 4, nonsquare); p=17: c=3 (order 16, nonsquare; the order-4 elements
    # mod 17 are +-4, both squares, so they correctly admit NO solution))
    for p, c, (k, lam, mu) in ((13, 5, (6, 2, 3)), (17, 3, (8, 3, 4))):
        R = [[(p - 1) // 2]]
        M = Model(p, 1, k, lam, mu, R, (0,), c)
        st, sols = M.solve(time_s=60, enumerate_all=True)
        good = all(verify_srg(D, p, 1, k, lam, mu) and verify_selfcomp(D, p, 1, (0,), c) for D in sols)
        qr = {(x * x) % p for x in range(1, p)}
        found_paley = any(D[0][0] == qr for D in sols)
        chk(f"Paley({p}) c={c}", st == "OPTIMAL" and len(sols) >= 2 and good and found_paley,
            f"status={st} count={len(sols)} allverify={good} paley={found_paley}")

    # V3: brute-force cross-check on p=5, norb=5, pi=(0 1 2 3)(4), c=2
    p, norb, c = 5, 2 + 3, 2
    pi = (1, 2, 3, 0, 4)
    co = lambda S, dia: (set(range(1, p)) if dia else set(range(p))) - S
    def brute(R, k, lam, mu):
        import conf_core as cc
        orbits, member, cinv = build_orbits(p, norb, pi, c)
        reps = orbits
        # free choices per rep
        def choices(rep):
            i, j = rep
            r = int(R[i][j])
            out = []
            if i == j:
                univ = [d for d in range(1, p)]
                for sub in itertools.chain.from_iterable(itertools.combinations(univ, t) for t in range(p)):
                    S = set(sub)
                    if len(S) != r: continue
                    if S != {(-d) % p for d in S}: continue
                    if pi[i] == i and any(((c * d) % p in S) == (d in S) for d in range(1, p)):
                        continue                      # fixed orbit: anticlosure x(cd) = NOT x(d)
                    out.append(S)
            else:
                for sub in itertools.chain.from_iterable(itertools.combinations(range(p), t) for t in range(p + 1)):
                    if len(sub) == r: out.append(set(sub))
            return out
        sols = set()
        n = p * norb
        target = mu * np.ones((n, n), dtype=np.int64)
        for combo in itertools.product(*[choices(rep) for rep in reps]):
            assign = dict(zip(reps, combo))
            D = [[None] * norb for _ in range(norb)]
            consistent = True
            for a in range(norb):
                for b in range(norb):
                    aa, bb = (a, b) if a <= b else (b, a)
                    rep, kk, sw = member[(aa, bb)]
                    S = assign[rep]
                    for _ in range(kk):                    # walk forward k steps
                        S = co({(c * d) % p for d in S}, rep[0] == rep[1])
                    if sw: S = {(-d) % p for d in S}
                    if a > b: S = {(-d) % p for d in S}
                    D[a][b] = S
            for i in range(norb):
                if 0 in D[i][i] or D[i][i] != {(-d) % p for d in D[i][i]}: consistent = False
                for j in range(norb):
                    if D[j][i] != {(-d) % p for d in D[i][j]}: consistent = False
                    if len(D[i][j]) != int(R[i][j]): consistent = False
            if not consistent: continue
            A = cc.build_adjacency(D, p, norb)
            lhs = A @ A + (mu - lam) * A - (k - mu) * np.eye(n, dtype=np.int64)
            if np.array_equal(lhs, target):
                sols.add(tuple(frozenset(D[i][j]) for i in range(norb) for j in range(norb)))
        return sols
    # complement-symmetric R choices built from free parameters (d_a0, R01, R02, R40):
    # diag pattern (d, 4-d, d, 4-d, 2); R[12]=5-R[01], R[23]=R[01], R[03]=5-R[01];
    # R[13]=5-R[02]; R[14]=5-R[04], R[24]=R[04], R[34]=5-R[04]; symmetric.
    def mkR(d, r01, r02, r04):
        R = [[0] * 5 for _ in range(5)]
        diag = [d, 4 - d, d, 4 - d, 2]
        for i in range(5): R[i][i] = diag[i]
        def setp(i, j, v): R[i][j] = v; R[j][i] = v
        setp(0, 1, r01); setp(1, 2, 5 - r01); setp(2, 3, r01); setp(0, 3, 5 - r01)
        setp(0, 2, r02); setp(1, 3, 5 - r02)
        setp(0, 4, r04); setp(1, 4, 5 - r04); setp(2, 4, r04); setp(3, 4, 5 - r04)
        return R
    Rs = [mkR(2, 2, 1, 3), mkR(0, 4, 2, 1), mkR(2, 3, 4, 2)]
    for Ridx, R in enumerate(Rs):
        Rok = all(int(R[pi[i]][pi[j]]) == ((p - 1 if i == j else p) - int(R[i][j])) for i in range(norb) for j in range(norb))
        if not Rok:
            chk(f"V3 R{Ridx} size-consistency", False, "bad test matrix"); continue
        bs = brute(R, 12, 5, 6)
        M = Model(p, norb, 12, 5, 6, R, pi, c)
        if not M.ok:
            es = set()
            st = "INCOMPATIBLE"
        else:
            st, sols = M.solve(time_s=120, enumerate_all=True)
            es = {tuple(frozenset(D[i][j]) for i in range(norb) for j in range(norb)) for D in sols}
        chk(f"V3 R{Ridx} engine==brute", es == bs, f"engine={len(es)} brute={len(bs)} status={st}")

    # V4: POSITIVE control for the determination chains: sizes-only counts (no SRG equations)
    # must match the definitional enumeration exactly (nonzero counts).
    import conf_core as cc
    def brute_sizes(R):
        orbits, member, cinv = build_orbits(p, norb, pi, c)
        co2 = lambda S, dia: (set(range(1, p)) if dia else set(range(p))) - S
        def choices(rep):
            i, j = rep
            r = int(R[i][j]); out = []
            if i == j:
                for sub in itertools.chain.from_iterable(itertools.combinations(range(1, p), t) for t in range(p)):
                    S = set(sub)
                    if len(S) != r or S != {(-d) % p for d in S}: continue
                    if pi[i] == i and any(((c * d) % p in S) == (d in S) for d in range(1, p)): continue
                    out.append(S)
            else:
                out = [set(s) for s in itertools.combinations(range(p), r)]
            return out
        cnt = 0
        for combo in itertools.product(*[choices(rep) for rep in orbits]):
            assign = dict(zip(orbits, combo))
            D = [[None] * norb for _ in range(norb)]
            ok = True
            for a in range(norb):
                for b in range(norb):
                    aa, bb = (a, b) if a <= b else (b, a)
                    rep, kk, sw = member[(aa, bb)]
                    S = assign[rep]
                    for _ in range(kk):
                        S = co2({(c * d) % p for d in S}, rep[0] == rep[1])
                    if sw: S = {(-d) % p for d in S}
                    if a > b: S = {(-d) % p for d in S}
                    D[a][b] = S
            for i in range(norb):
                if 0 in D[i][i] or D[i][i] != {(-d) % p for d in D[i][i]}: ok = False
                for j in range(norb):
                    if D[j][i] != {(-d) % p for d in D[i][j]} or len(D[i][j]) != int(R[i][j]): ok = False
            if ok: cnt += 1
        return cnt
    class SizeOnly(Model):
        def _build(self2):
            # same as Model._build but without the SRG equations
            pp, nn, R2, pi2, c2 = self2.p, self2.norb, self2.R, self2.pi, self2.c
            for i in range(nn):
                for j in range(nn):
                    if int(R2[pi2[i]][pi2[j]]) != ((pp - 1 if i == j else pp) - int(R2[i][j])): return False
            for (i, j) in self2.orbits:
                for d in range(pp):
                    if i == j and d == 0: continue
                    self2.X[((i, j), d)] = self2.m.NewBoolVar(f"x{i}_{j}_{d}")
                if i == j:
                    for d in range(1, pp):
                        self2.m.Add(self2.X[((i, j), d)] == self2.X[((i, j), (pp - d) % pp)])
                    if pi2[i] == i:
                        for d in range(1, pp):
                            self2.m.Add(self2.X[((i, j), (c2 * d) % pp)] + self2.X[((i, j), d)] == 1)
                    self2.m.Add(sum(self2.X[((i, j), d)] for d in range(1, pp)) == int(R2[i][i]))
                else:
                    self2.m.Add(sum(self2.X[((i, j), d)] for d in range(pp)) == int(R2[i][j]))
            return True
    for Ridx, R in enumerate(Rs[:2]):
        bc = brute_sizes(R)
        S = SizeOnly(p, norb, 12, 5, 6, R, pi, c)
        st, sols = S.solve(time_s=120, enumerate_all=True)
        chk(f"V4 R{Ridx} sizes-only count engine==brute (positive control)", len(sols) == bc and bc > 0,
            f"engine={len(sols)} brute={bc}")

    # V5: P(25)'s actual orbit matrix, if complement-symmetric under some pi: full-pipeline run
    Dp = cc.paley_pp_blocks(5)
    Rp = [[len(Dp[i][j]) for j in range(5)] for i in range(5)]
    perms = []
    for q in itertools.permutations(range(5)):
        if all(int(Rp[q[i]][q[j]]) == ((4 if i == j else 5) - int(Rp[i][j])) for i in range(5) for j in range(5)):
            cyc_ok = sorted(len(cyc) for cyc in _cycles(q)) == [1, 4]
            if cyc_ok: perms.append(q)
    if perms:
        M = Model(5, 5, 12, 5, 6, Rp, perms[0], 2)
        st, sols = M.solve(time_s=300)
        good = all(verify_srg(D, 5, 5, 12, 5, 6) and verify_selfcomp(D, 5, 5, perms[0], 2) for D in sols)
        chk("V5 srg(25) on P(25)-skeleton", st in ("OPTIMAL", "FEASIBLE", "INFEASIBLE") and good,
            f"status={st} sols={len(sols)} allverify={good}" + (" *** FOUND self-comp Z5 srg(25)!" if sols else ""))
    else:
        print("  [INFO] V5 skipped: P(25) orbit matrix not complement-symmetric under any (4,1) perm")
    print("SELF-COMP MODEL VALIDATION:", "ALL PASS" if ok_all else "FAILURES")
    return ok_all

def _cycles(q):
    seen, out = set(), []
    for s in range(len(q)):
        if s in seen: continue
        cl, x = [], s
        while x not in seen:
            seen.add(x); cl.append(x); x = q[x]
        out.append(cl)
    return out

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--validate":
        sys.exit(0 if validate() else 1)
    print("use --validate, or import and use Model(...) for production runs")
