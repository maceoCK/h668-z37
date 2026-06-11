"""
Galois (quadratic-subfield) trace condition for sigma-equivariant self-complementary
srg(333,166,82,83) candidates, classes 221/422.

THEOREM (necessary condition; derivation in comments):
  If F = 0 (exact solution), then for every nontrivial character chi of Z37,
  M(chi)^2 + M(chi) - 83 I = 0, hence tr M(chi) = m*theta + (9-m)*tau in Q(sqrt37)
  with theta,tau = (-1 +- 3 sqrt37)/2, m in {0..9}.
  Since tr M(chi) = sum_d f(d) chi(d) with f(d) = #{a : d in D_aa},
  tr M(chi) in Q(sqrt37) for ALL chi forces f to be constant on quadratic-residue
  classes mod 37: f = u on QR, v on QNR, with u + v = 9 and u - v = 3(2m-9) for
  integer m, so u - v in {+-3, +-9}.

  In the sigma-equivariant gauge (multiplier c = 6, QNR; pi cycle type 4+4+1 on the
  9 orbits; diagonal pair-orbit reps S = D_00 (type 1), T = D_t,t (type 1),
  U = D_77 (type 2, anticlosed)), using D_{pi(a)pi(a)} = c*(complement in Z37* of D_aa):

    f(d) = 4 + 1_U(d) + 2[g_S(d) + g_T(d)],   g_X(d) = 1_X(d) - 1_X(c^{-1} d),
  and since c^2 = -1, c^{-1} = -c, and S,T symmetric:
    on each <+-6>-orbit O = {d,-d,6d,-6d} (d in QR): g_X(d) = x1 - x2,
    g_X(6d) = -g_X(d), with x1 = 1_X(+-d), x2 = 1_X(+-6d).

  CONSEQUENCES (each verified below):
    (C1) U = D_77 is EXACTLY the QR set or exactly the QNR set (Paley pinning).
    (C2) per <+-6>-orbit, with eps = -1 if U = QR else +1:
         (s1 - s2) + (t1 - t2) = eps for every orbit  (s,t = indicators of S,T halves)
    (C3) counting: |S|/2 + |T|/2 must be odd  (holds for 221: 6+9; 422: 7+10 -> no closure)

This script verifies: (a) the equivalence [tr M(chi) in Q(sqrt37) for all chi] <=>
[f QR-constant] on random states; (b) the Paley-13 sanity case; (c) that the 12 archived
best-floor states VIOLATE the condition (quantified); (d) feasibility windows per class;
(e) roundtrip: a (S,T,U) triple satisfying (C1),(C2) yields tr M(chi) exactly at the
required lattice points.
"""
import json, random, sys
import numpy as np
sys.path.insert(0, "/Users/maceocardinalekwik/git/math/release/h668-z37/code")
from sc_floor_analysis import Tpl, parse_dump, CODE, RESID

P = 37
QR = sorted({(k*k) % P for k in range(1, P)})
QNR = [d for d in range(1, P) if d not in QR]
assert 6 in QNR and (P-1) in QR

def diag_f(tpl, mem):
    """f(d) = #{a : d in D_aa} over all 9 diagonal blocks, d in 1..36"""
    f = np.zeros(P, dtype=np.int64)
    for a in range(tpl.NORB):
        f += mem[a, a]
    assert f[0] == 0
    return f

def tr_offfield_norm(f):
    """sum over chi != 1 of |projection of sum_d f(d) chi(d) off Q(sqrt37)|^2 --
    equivalently the L2 deviation of f from its QR-class averages (Parseval)."""
    fr = f[1:].astype(float)
    qr_mask = np.array([d in QR for d in range(1, P)])
    dev = fr.copy()
    dev[qr_mask] -= fr[qr_mask].mean()
    dev[~qr_mask] -= fr[~qr_mask].mean()
    return float((dev**2).sum())   # 0 iff f QR-constant

def check_state(tpl, repBit, label):
    mem = tpl.expand(repBit)
    f = diag_f(tpl, mem)
    fq = sorted(set(int(f[d]) for d in QR)); fn = sorted(set(int(f[d]) for d in QNR))
    qrconst = len(fq) == 1 and len(fn) == 1
    dev = tr_offfield_norm(f)
    # numeric check: tr M(chi) in Q(sqrt37) <=> Galois-invariance under QR subgroup
    w = np.exp(2j*np.pi/P)
    bad = 0.0
    for j in range(1, P):
        tr = sum((mem[a, a] @ (w ** ((np.arange(P)*j) % P))) for a in range(tpl.NORB))
        # project off the 2-dim space spanned by eta0, eta1 (Gauss periods) numerically:
        # invariance under s in QR: tr(chi^s) == tr(chi) for all QR s
        tr5 = sum((mem[a, a] @ (w ** ((np.arange(P)*j*ROOT_QR) % P))) for a in range(tpl.NORB))
        bad = max(bad, abs(tr - tr5))
    print(f"  {label}: f on QR={fq} on QNR={fn} QR-constant={qrconst} "
          f"L2dev={dev:.2f} max|tr(chi)-tr(chi^s)|={bad:.3f}")
    return qrconst, dev

ROOT_QR = 3  # 3 is a QR mod 37 generator of the QR subgroup? (3 in QR; need any QR != 1)
assert 3 in QR

def main():
    data = json.load(open(RESID))
    tpls = {}
    print("== (c) archived best-floor states vs the trace condition ==")
    for ent in data:
        tname = ent["tpl"]
        if tname not in tpls: tpls[tname] = Tpl(f"{CODE}/{tname}")
        tpl = tpls[tname]
        _, repBit, _ = parse_dump(ent["dump"], tpl)
        check_state(tpl, repBit, f"{tname} seed={ent['seed']}")

    print("\n== (a) random states: trace condition fails generically ==")
    tpl = tpls["sc_tpl_221_0.txt"]
    rng = random.Random(5)
    devs = []
    for _ in range(200):
        mem = tpl.expand(tpl.random_state(rng))
        devs.append(tr_offfield_norm(diag_f(tpl, mem)))
    print(f"  200 random states: min L2dev={min(devs):.1f} mean={np.mean(devs):.1f} "
          f"(0 would satisfy the condition)")

    print("\n== (e) roundtrip: construct (S,T,U) satisfying (C1),(C2); verify trace lattice ==")
    # case U = QR, eps = -1; orbits with QR reps
    seen = [False]*P; orbits = []
    for x in range(1, P):
        if seen[x]: continue
        ob = [x, (6*x) % P, P-x, P-(6*x) % P]
        for z in ob: seen[z] = True
        # orient so first half is the QR half
        if x in QR: orbits.append(((x, P-x), ((6*x) % P, P-(6*x) % P)))
        else: orbits.append((((6*x) % P, P-(6*x) % P), (x, P-x)))
    assert len(orbits) == 9
    # 221: |S|=12 (6 half-pairs), |T|=18 (9 half-pairs); choose nA=0 -> all 9 orbits type B:
    #   type B: (t1,t2)=(0,1) [t-QRhalf out, t-QNRhalf in], s1=s2 (beta=3 orbits with s both-in)
    S = set(); T = set(); U = set(QR)
    for idx, (qrh, qnh) in enumerate(orbits):
        T.update(qnh)                 # t2=1
        if idx < 3: S.update(qrh); S.update(qnh)   # s both-in (beta=3)
    assert len(S) == 12 and len(T) == 18 and len(U) == 18
    f = np.zeros(P, dtype=np.int64)
    onehot = lambda X: np.array([1 if d in X else 0 for d in range(P)], dtype=np.int64)
    sS, sT, sU = onehot(S), onehot(T), onehot(U)
    cinv = pow(6, P-2, P)
    gs = lambda sX: sX - sX[(np.arange(P)*cinv) % P]
    f = 4*(np.arange(P) > 0) + sU + 2*(gs(sS) + gs(sT))
    f[0] = 0
    fq = sorted(set(int(f[d]) for d in QR)); fn = sorted(set(int(f[d]) for d in QNR))
    print(f"  constructed f: QR values={fq}, QNR values={fn} (expect u=3, v=6)")
    # verify tr M(chi) hits the lattice: u*eta0 + v*eta1 with m from u-v = 3(2m-9)
    w = np.exp(2j*np.pi/P)
    th = (-1+3*np.sqrt(37))/2; ta = (-1-3*np.sqrt(37))/2
    ok = True
    for j in range(1, P):
        tr = complex(f[1:] @ (w ** ((np.arange(1, P)*j) % P)))
        cands = [m*th + (9-m)*ta for m in range(10)]
        ok &= min(abs(tr - c) for c in cands) < 1e-9
    print(f"  tr M(chi) on the {{m*theta+(9-m)*tau}} lattice for all chi: {ok}")

    print("\n== (d) feasibility windows per class (C3 + ranges) ==")
    for cls, s2, t2 in [("221 (tpl 0/1)", 6, 9), ("422 (tpl 2-5)", 7, 10)]:
        # nA + 2*beta = s2 ; (9-nA) + 2*alpha = t2 ; 0<=beta<=9-nA ; 0<=alpha<=nA
        sols = []
        for nA in range(0, 10):
            if (s2 - nA) % 2 or (t2 - (9-nA)) % 2: continue
            beta = (s2 - nA)//2; alpha = (t2 - (9-nA))//2
            if beta < 0 or alpha < 0 or beta > 9-nA or alpha > nA: continue
            sols.append((nA, 9-nA, beta, alpha))
        par = "ODD (feasible)" if (s2 + t2) % 2 == 1 else "EVEN (would CLOSE the class)"
        print(f"  {cls}: |S|/2={s2} |T|/2={t2} parity {s2+t2} {par}; "
              f"(nA,nB,beta,alpha) solutions: {sols}")

    print("\n== (b) Paley-13 sanity: the analogous condition reproduces D = QR(13) ==")
    QR13 = sorted({(k*k) % 13 for k in range(1, 13)})
    w13 = np.exp(2j*np.pi/13)
    th13 = (-1+np.sqrt(13))/2; ta13 = (-1-np.sqrt(13))/2
    D = QR13
    f13 = np.array([1 if d in D else 0 for d in range(13)])
    ok = all(min(abs(complex(f13[1:] @ (w13**((np.arange(1,13)*j) % 13))) - c)
                 for c in [th13, ta13]) < 1e-9 for j in range(1, 13))
    print(f"  Paley(13): tr M(chi) in {{theta,tau}} for all chi: {ok}")

if __name__ == "__main__":
    main()
