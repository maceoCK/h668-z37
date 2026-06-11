"""
INDEPENDENT verification of the Paley-pinning / (C1)(C2) necessary condition
for self-complementary Z37-symmetric srg(333,166,82,83) candidates (classes 221/422).

Everything below is re-implemented from scratch from the raw file formats
(template LIT/REP lines, sc_residuals.json dumps). No reuse of the agent's
analysis code paths.
"""
import json, random, itertools
import numpy as np

CODE = "/Users/maceocardinalekwik/git/math/release/h668-z37/code"
RESID = "/Users/maceocardinalekwik/git/math/deepdive_lp333/sc_residuals.json"
P = 37
K, LAM, MU = 166, 82, 83
QR = sorted({pow(k, 2, P) for k in range(1, P)})
QNR = [d for d in range(1, P) if d not in QR]
assert 6 in QNR and pow(6, 2, P) == P - 1, "c=6 must be a QNR with c^2=-1 mod 37"

# ---------- my own template parser ----------
class T:
    def __init__(self, path):
        toks = open(path).read().split()
        it = iter(toks)
        self.P = int(next(it)); self.N = int(next(it)); self.NR = int(next(it))
        self.k = int(next(it)); self.lam = int(next(it)); self.mu = int(next(it)); self.c = int(next(it))
        self.reps = [None] * self.NR
        for _ in range(self.NR):
            assert next(it) == "REP"
            idx = int(next(it))
            self.reps[idx] = tuple(int(next(it)) for _ in range(4))  # i j type card
        self.lr = np.full((self.N, self.N, self.P), -1, np.int64)
        self.le = np.zeros((self.N, self.N, self.P), np.int64)
        self.ln = np.zeros((self.N, self.N, self.P), np.int64)
        for _ in range(self.N * self.N * self.P):
            assert next(it) == "LIT"
            a, b, d, ri, e, ng = (int(next(it)) for _ in range(6))
            self.lr[a, b, d], self.le[a, b, d], self.ln[a, b, d] = ri, e, ng
        self.isrep = np.zeros((self.N, self.N), bool)
        for (i, j, t, c) in self.reps: self.isrep[i, j] = True

    def expand(self, rb):
        mem = np.zeros((self.N, self.N, self.P), np.int64)
        pos = self.lr >= 0
        mem[pos] = rb[self.lr[pos], self.le[pos]] ^ self.ln[pos]
        mem[~pos] = self.ln[~pos]
        return mem

    def resid(self, mem):
        # conv[a,b,g] = sum_l sum_d mem[a,l,d] mem[l,b,g-d]  (direct, no FFT)
        N, p = self.N, self.P
        idx = (np.arange(p)[:, None] - np.arange(p)[None, :]) % p   # idx[g,d] = g-d
        cv = np.zeros((N, N, p), np.int64)
        for l in range(N):
            # mem[l,b,(g-d)%p] -> shape (b,g,d)
            sh = mem[l][:, idx]                     # (N, p, p) over (b, g, d)
            cv += np.einsum('ad,bgd->abg', mem[:, l, :], sh)
        rhs = np.full((N, N, p), self.mu, np.int64)
        for a in range(N): rhs[a, a, 0] += self.k - self.mu
        return cv - (self.lam - self.mu) * mem - rhs

    def obj(self, res):
        s = 0
        for a in range(self.N):
            for b in range(a, self.N):
                if self.isrep[a, b]: s += int((res[a, b] ** 2).sum())
        return s

    def random_state(self, rng):
        rb = np.zeros((self.NR, self.P), np.int64)
        for r, (i, j, typ, card) in enumerate(self.reps):
            if typ == 0:
                rb[r, rng.sample(range(self.P), card)] = 1
            elif typ == 1:
                for d in rng.sample(range(1, self.P // 2 + 1), card // 2):
                    rb[r, d] = rb[r, self.P - d] = 1
            else:  # type 2 anticlosed: orbit {x,6x,-x,-6x}: take {x,-x} or {6x,-6x}
                seen = [False] * self.P
                for x in range(1, self.P):
                    if seen[x]: continue
                    ob = [x, 6 * x % self.P, self.P - x, (self.P - 6 * x % self.P) % self.P]
                    for z in ob: seen[z] = True
                    pick = (ob[0], ob[2]) if rng.random() < .5 else (ob[1], ob[3])
                    for z in pick: rb[r, z] = 1
        return rb

def parse_dump(dump, t):
    rb = np.zeros((t.NR, t.P), np.int64); obj = None
    for ln in dump.splitlines():
        w = ln.split()
        if not w: continue
        if w[0] == "BESTSTATE": obj = int(w[1].split("=")[1])
        elif w[0] == "REPBITS": rb[int(w[1]), [int(x) for x in w[3:]]] = 1
    return obj, rb

# <±6>-orbits of Z37*, oriented (QR half, QNR half)
orbs = []
seen = [False] * P
for x in range(1, P):
    if seen[x]: continue
    ob = [x, 6 * x % P, P - x, (P - 6 * x % P) % P]
    for z in ob: seen[z] = True
    qrh = tuple(sorted(z for z in ob if z in QR))
    qnh = tuple(sorted(z for z in ob if z not in QR))
    assert len(qrh) == 2 and len(qnh) == 2, "each <±6>-orbit must split 2 QR + 2 QNR"
    orbs.append((qrh, qnh))
assert len(orbs) == 9

w37 = np.exp(2j * np.pi / P)
TH, TA = (-1 + 3 * np.sqrt(37)) / 2, (-1 - 3 * np.sqrt(37)) / 2   # roots of x^2+x-83
LATTICE = [m * TH + (9 - m) * TA for m in range(10)]

def f_of(mem):
    f = np.zeros(P, np.int64)
    for a in range(mem.shape[0]): f += mem[a, a]
    return f

def trace_in_lattice(f, tol=1e-7):
    bad = 0.0
    for j in range(1, P):
        tr = complex(f[1:] @ (w37 ** (np.arange(1, P) * j % P)))
        bad = max(bad, min(abs(tr - c) for c in LATTICE))
    return bad < tol, bad

def qr_const(f):
    uq = {int(f[d]) for d in QR}; un = {int(f[d]) for d in QNR}
    return (len(uq) == 1 and len(un) == 1), sorted(uq), sorted(un)

def l2dev(f):
    fr = f[1:].astype(float)
    m = np.array([d in set(QR) for d in range(1, P)])
    dev = fr.copy(); dev[m] -= fr[m].mean(); dev[~m] -= fr[~m].mean()
    return float((dev ** 2).sum())

print("=" * 70)
print("[1] reconstruct all 12 archived dumps; recompute objective; check f")
data = json.load(open(RESID))
tpls = {}
viol = 0
for ent in data:
    nm = ent["tpl"]
    if nm not in tpls: tpls[nm] = T(f"{CODE}/{nm}")
    t = tpls[nm]
    claimed, rb = parse_dump(ent["dump"], t)
    for r, (i, j, typ, card) in enumerate(t.reps):
        assert rb[r].sum() == card, (nm, r)
    mem = t.expand(rb)
    res = t.resid(mem)
    o = t.obj(res)
    f = f_of(mem)
    qc, uq, un = qr_const(f)
    inlat, bad = trace_in_lattice(f)
    viol += (not qc)
    print(f"  {nm} s={ent['seed']}: obj claimed={claimed} mine={o} match={o==claimed}"
          f"  QR-const={qc} L2dev={l2dev(f):.2f} trace-lattice={inlat}")
print(f"  => archived floor states violating the condition: {viol}/12")

print("=" * 70)
print("[2] gauge structure: complement relation + 4-cycle diagonal return")
# infer pi from each template by checking D_{pi a,pi a} = 6*(Z*\D_aa) elementwise on a dump state
for nm, t in tpls.items():
    _, rb = parse_dump([e for e in data if e["tpl"] == nm][0]["dump"], t)
    mem = t.expand(rb)
    # candidate pi: for each a find unique b with D_bb = 6*(complement of D_aa)
    pi = []
    for a in range(9):
        compl = np.zeros(P, np.int64)
        for d in range(1, P): compl[6 * d % P] = 1 - mem[a, a, d]
        cands = [b for b in range(9) if np.array_equal(mem[b, b], compl)]
        pi.append(cands)
    ok_unique = all(len(c) >= 1 for c in pi)
    pi0 = [c[0] for c in pi]
    # 4-cycle return: D_{pi^2 a} = D_aa  (since 36X = -X = X for symmetric X)
    ret = all(np.array_equal(mem[pi0[pi0[a]], pi0[pi0[a]]], mem[a, a]) for a in range(9))
    fix = [a for a in range(9) if pi0[a] == a]
    print(f"  {nm}: complement-relation pi exists={ok_unique} pi={pi0} fixed={fix} pi^2=id on diag={ret}")

print("=" * 70)
print("[3] derivation core, checked exhaustively over the reduced state space")
# f(d) = 4 + 1_U(d) + 2[gS(d)+gT(d)] ; per orbit: f(QRhalf)=4+uO+2e, f(QNRhalf)=5-uO-2e
# QR-constancy <=> uO+2e identical across orbits; with u+v=9 and trace lattice u-v in {±3,±9}
# u-v=2(uO+2e)-1: ±3 -> (uO,e)=(0,+1)/(1,-1); ±9 -> (uO,e)=(1,+2)/(0,-2) needs |S|=|T|=18.
# Exhaustive check of the per-orbit algebra (all (uO,s1,s2,t1,t2) in {0,1}^5):
rows = set()
for uO, s1, s2, t1, t2 in itertools.product((0, 1), repeat=5):
    e = (s1 - s2) + (t1 - t2)
    rows.add((4 + uO + 2 * e, 5 - uO - 2 * e))
# which (u,v) pairs are simultaneously achievable orbit-wise AND lie on the lattice (u-v in ±3,±9, u+v=9)?
lat_uv = [(u, v) for (u, v) in rows if u + v == 9 and abs(u - v) in (3, 9) and 0 <= u <= 9 and 0 <= v <= 9]
print(f"  per-orbit reachable (u,v) on the lattice: {sorted(set(lat_uv))}")
# realizations per (u,v):
for uv in sorted(set(lat_uv)):
    real = [(uO, s1, s2, t1, t2) for uO, s1, s2, t1, t2 in itertools.product((0, 1), repeat=5)
            if (4 + uO + 2 * ((s1 - s2) + (t1 - t2)), 5 - uO - 2 * ((s1 - s2) + (t1 - t2))) == uv]
    uOs = {r[0] for r in real}
    es = {(r[1] - r[2]) + (r[3] - r[4]) for r in real}
    print(f"   (u,v)={uv}: forced uO in {sorted(uOs)}, forced e in {sorted(es)}, #orbit-realizations={len(real)}")
# cardinality kill of |u-v|=9 for the actual classes:
for nm, t in tpls.items():
    cards = {(i, j): card for (i, j, typ, card) in t.reps}
    diag = sorted((i, card) for (i, j, typ, card) in t.reps if i == j)
    print(f"  {nm}: diagonal rep cards {diag}  (|u-v|=9 needs S,T cards == 18,18)")

print("=" * 70)
print("[4] equivalence on random states: trace-lattice <=> QR-constancy of f")
t = tpls["sc_tpl_221_0.txt"]
rng = random.Random(7)
agree = True
nlat = 0
for _ in range(300):
    mem = t.expand(t.random_state(rng))
    f = f_of(mem)
    qc, _, _ = qr_const(f)
    inlat, _ = trace_in_lattice(f)
    agree &= (qc == inlat)
    nlat += inlat
print(f"  300 random states: equivalence holds everywhere: {agree}; states meeting condition: {nlat}")

print("=" * 70)
print("[5] roundtrip: build states satisfying (C1)+(C2), verify f=(3,6)/(6,3) + lattice")
rng = random.Random(11)
ok_all = True
for branch, eps, expected in (("QR", -1, (3, 6)), ("QNR", +1, (6, 3))):
    for trial in range(25):
        S = set(); Tset = set(); U = set()
        for (qrh, qnh) in orbs:
            U.update(qrh if branch == "QR" else qnh)
            # choose (s1,s2,t1,t2) uniformly among realizations with (s1-s2)+(t1-t2)=eps
            cands = [(s1, s2, t1, t2) for s1, s2, t1, t2 in itertools.product((0, 1), repeat=4)
                     if (s1 - s2) + (t1 - t2) == eps]
            s1, s2, t1, t2 = rng.choice(cands)
            if s1: S.update(qrh)
            if s2: S.update(qnh)
            if t1: Tset.update(qrh)
            if t2: Tset.update(qnh)
        f = np.zeros(P, np.int64)
        for d in range(1, P):
            gS = (d in S) - ((-6 * d) % P in S)     # c^{-1} = -6
            gT = (d in Tset) - ((-6 * d) % P in Tset)
            f[d] = 4 + (d in U) + 2 * (gS + gT)
        qc, uq, un = qr_const(f)
        inlat, _ = trace_in_lattice(f)
        ok_all &= qc and inlat and (uq[0], un[0]) == expected
    print(f"  branch {branch}: 25 random (C1)(C2)-states all give f=({expected[0]},{expected[1]}) on lattice: {ok_all}")

print("=" * 70)
print("[6] Fourier bridge on dump 0: 37*F_full == sum_chi!=1 ||M(chi)^2+M(chi)-83I||_F^2")
t = tpls["sc_tpl_221_0.txt"]
_, rb = parse_dump(data[0]["dump"], t)
mem = t.expand(rb)
res = t.resid(mem)
Ffull = sum(int((res[a, b] ** 2).sum()) for a in range(9) for b in range(9))
tot = 0.0
for j in range(1, P):
    M = np.zeros((9, 9), complex)
    for a in range(9):
        for b in range(9):
            M[a, b] = mem[a, b] @ (w37 ** (np.arange(P) * j % P))
    Q = M @ M + M - 83 * np.eye(9)
    tot += float((abs(Q) ** 2).sum())
    herm_ok = np.allclose(M, M.conj().T)
print(f"  37*F_full = {37*Ffull}, sum = {tot:.1f}, match = {abs(37*Ffull-tot)<1e-3}, M(chi) Hermitian = {herm_ok}")

print("=" * 70)
print("[7] Paley anchors: the analog condition HOLDS on the true solutions")
for p, c, k, lam, mu in ((13, 5, 6, 2, 3), (17, 3, 8, 3, 4)):
    QRp = sorted({pow(x, 2, p) for x in range(1, p)})
    QNp = [d for d in range(1, p) if d not in QRp]
    assert c in QNp
    D = QRp
    # anticlosure: D = c*(Z_p^* \ D)
    anti = sorted({c * d % p for d in range(1, p) if d not in D}) == sorted(D)
    # SRG check on the circulant
    A = np.zeros((p, p), np.int64)
    for g in range(p):
        for d in D: A[g, (g + d) % p] = 1
    srg_ok = np.array_equal(A @ A + (mu - lam) * A - (k - mu) * np.eye(p, dtype=np.int64),
                            mu * np.ones((p, p), dtype=np.int64))
    th, ta = (-1 + np.sqrt(p)) / 2, (-1 - np.sqrt(p)) / 2     # roots of x^2+x-(k-mu) when k-mu=(p-1)/4
    wp = np.exp(2j * np.pi / p)
    fD = np.array([1 if d in D else 0 for d in range(p)])
    lat = all(min(abs(complex(fD[1:] @ (wp ** (np.arange(1, p) * j % p))) - z)
                  for z in (th, ta)) < 1e-9 for j in range(1, p))
    print(f"  Paley({p}) c={c}: anticlosed={anti} SRG={srg_ok} trace-condition holds={lat}")
# exhaustive pinning at p=13, c=5 (3 orbits -> 8 anticlosed symmetric candidates)
p, c = 13, 5
QRp = {pow(x, 2, p) for x in range(1, p)}
orb13 = []
seen = [False] * p
for x in range(1, p):
    if seen[x]: continue
    ob = [x, c * x % p, p - x, (p - c * x % p) % p]
    for z in ob: seen[z] = True
    orb13.append(ob)
th, ta = (-1 + np.sqrt(13)) / 2, (-1 - np.sqrt(13)) / 2
wp = np.exp(2j * np.pi / 13)
sat = []
for pick in itertools.product((0, 1), repeat=len(orb13)):
    D = set()
    for ob, b in zip(orb13, pick):
        D.update((ob[0], ob[2]) if b == 0 else (ob[1], ob[3]))
    fD = np.array([1 if d in D else 0 for d in range(13)])
    if all(min(abs(complex(fD[1:] @ (wp ** (np.arange(1, 13) * j % 13))) - z)
               for z in (th, ta)) < 1e-9 for j in range(1, 13)):
        sat.append(sorted(D))
print(f"  Paley(13) exhaustive over 2^3 anticlosed sets: {len(sat)} satisfy; "
      f"= [QR,QNR]: {sat == sorted([sorted(QRp), sorted(set(range(1,13))-QRp)]) or sorted(map(tuple,sat)) == sorted([tuple(sorted(QRp)), tuple(sorted(set(range(1,13))-QRp))])}")
print("  satisfying sets:", sat)

print("=" * 70)
print("[8] C3 parity: per-orbit (s1+s2+t1+t2) odd => |S|/2+|T|/2 odd (9 orbits)")
for nm, t in tpls.items():
    diag = {i: card for (i, j, typ, card) in t.reps if i == j}
    fixedcard = [card for (i, j, typ, card) in t.reps if typ == 2][0]
    cyc = sorted((card for (i, j, typ, card) in t.reps if typ == 1))
    s2, t2 = cyc[0] // 2, cyc[1] // 2
    print(f"  {nm}: |S|/2+|T|/2 = {s2}+{t2} = {s2+t2} ({'ODD ok' if (s2+t2)%2 else 'EVEN -> class CLOSED'})")
print("done")
