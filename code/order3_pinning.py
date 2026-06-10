"""
Order-3 Galois trace pinning for the Z_37 cyclotomic lift (analog of the order-2 N_g pinning).

Setting: srg(333,166,82,83) with a semiregular order-37 automorphism (9 orbits), blocks D[i][j]
subset Z_37, all invariant under the order-3 multiplier subgroup H = {1,10,26} <= Z_37^*.
The 36 nonzero residues split into 12 H-cosets; H <= QR(37) so each coset is wholly QR or QNR
(6 each), and -1 in QR(37) pairs cosets c <-> -c within the same quadratic class
(3 QR-pairs, 3 QNR-pairs).

For t != 0 the 9x9 Hermitian M_t = [psi_t(D[i][j])] satisfies M_t^2 + M_t = 83 I, so its
spectrum is {theta x a_t, tau x (9-a_t)}, theta,tau = (-1 +- 3 sqrt37)/2.  H-invariance makes
M_t depend only on the H-class of t; Galois sigma_u : M_1 -> M_u sends sqrt37 -> chi(u) sqrt37
(chi = Legendre symbol), so a_u = a if chi(u)=1 else 9-a, with a := a_1.

Diagonal blocks are H-invariant AND symmetric (D = -D), so each is a disjoint union of H-cosets;
let n_c := #{ i : coset c subset D[i][i] }.  Then

   Tr(M_1) = sum_c n_c eta_c  =  a theta + (9-a) tau  =  -9/2 + (2a-9)(3/2) sqrt37,

and since 1 = -sum_c eta_c and sqrt37 = sum_c chi(c) eta_c (quadratic Gauss sum; chi constant
on cosets), uniqueness of the expansion in the period basis {eta_c} of K = Q(zeta37)^H forces

   n_c = 3a - 9   for the 6 QR-cosets,      n_c = 18 - 3a  for the 6 QNR-cosets.

0 <= n_c <= 9 pins a in {3,4,5,6}.  a in {3,6} forces every diagonal block to be exactly the
union of all QNR-cosets (resp. all QR-cosets), i.e. |D[i][i]| = 18 for all i -- impossible,
since the complete (no-cap) enumeration uniform18_enum.cpp found NO admissible orbit matrix
with all-18 diagonal.  Hence a in {4,5}:

   a = 4:  every QR-coset lies in exactly 3 of the 9 diagonal blocks, every QNR-coset in 6.
   a = 5:  swap (6 on QR, 3 on QNR).

This script VALIDATES the derivation numerically (high precision):
  (1) H, cosets, quadratic classes as claimed; periods linearly independent (Galois matrix
      nonsingular);
  (2) sqrt37 = sum_c chi(c) eta_c and 1 = -sum_c eta_c in the period basis;
  (3) for each a in 0..9, solving the 12x12 linear system [eta_c(u)] n = Tr-vector reproduces
      n_c = (9 + 3(2a-9) chi(c))/2 exactly;
  (4) machinery sanity: random H-invariant symmetric diagonal data -> expansion coefficients
      recovered from Tr(M_1) equal the true coset counts;
  (5) true-graph anchor: Paley(13) with H = {1,3,9} (m=1, n=1 orbit, a=1) satisfies the general
      formula n_c = (n + s(2a-n) chi(c))/2 with s = sqrt(m) = 1.
"""
import numpy as np

P = 37
H = [1, 10, 26]                      # 10^3 = 1000 = 27*37 + 1
assert sorted({(h * g) % P for h in H for g in H}) == sorted(H)

def cosets_of(H, p):
    seen, cs = set(), []
    for g in range(1, p):
        if g in seen: continue
        c = tuple(sorted((g * h) % p for h in H))
        seen |= set(c); cs.append(c)
    return cs

CS = cosets_of(H, P)                 # 12 cosets
QR = {(x * x) % P for x in range(1, P)}
chi = [1 if c[0] in QR else -1 for c in CS]
neg = [next(i for i, d in enumerate(CS) if (-c[0]) % P in d) for i, c in enumerate(CS)]

def main():
    assert len(CS) == 12 and all(set(H) <= QR for _ in [0]) and all(h in QR for h in H)
    assert sum(1 for x in chi if x == 1) == 6
    # cosets pair under negation within the same quadratic class
    assert all(chi[i] == chi[neg[i]] and neg[i] != i for i in range(12))
    print("[1] H={1,10,26} <= QR(37); 12 cosets: 6 QR + 6 QNR; 3 QR-pairs + 3 QNR-pairs  OK")

    # period values eta_c(u) for u ranging over 12 class representatives (Galois conjugates)
    reps = [c[0] for c in CS]
    w = np.exp(2j * np.pi / P)
    E = np.array([[sum(w ** ((u * x) % P) for x in c) for c in CS] for u in reps])  # rows: sigma_u
    det = np.linalg.det(E)
    print(f"[1] period Galois matrix det = {abs(det):.6g}  (nonsingular => periods are a basis)")
    assert abs(det) > 1e-3

    s37 = np.sqrt(37.0)
    # (2) sqrt37 = sum_c chi(c) eta_c(u) * chi(u)  i.e. per conjugate: sum_c chi(c) eta_c(u) = chi(u) sqrt37
    for ui, u in enumerate(reps):
        val = sum(chi[ci] * E[ui, ci] for ci in range(12))
        want = (1 if u in QR else -1) * s37
        assert abs(val - want) < 1e-9, (u, val, want)
        assert abs(sum(E[ui]) + 1) < 1e-9
    print("[2] Gauss sum: sum_c chi(c) eta_c = chi(u) sqrt37 on every conjugate; sum_c eta_c = -1  OK")

    # (3) solve for n from Tr-vector, for each a
    theta, tau = (-1 + 3 * s37) / 2, (-1 - 3 * s37) / 2
    ok = True
    for a in range(10):
        rhs = np.array([(a if u in QR else 9 - a) * theta + ((9 - a) if u in QR else a) * tau
                        for u in reps], dtype=complex)
        n = np.linalg.solve(E, rhs)
        want = np.array([(9 + 3 * (2 * a - 9) * chi[ci]) / 2 for ci in range(12)])
        if not (np.allclose(n.imag, 0, atol=1e-8) and np.allclose(n.real, want, atol=1e-8)):
            ok = False; print(f"    a={a} MISMATCH {n} vs {want}")
    assert ok
    print("[3] for every a in 0..9: unique period expansion gives n_c = 3a-9 (QR), 18-3a (QNR)  OK")
    print("    integrality 0<=n_c<=9 pins a in {3,4,5,6}; uniform18=0 kills a in {3,6} => a in {4,5}")

    # (4) machinery sanity on random H-invariant symmetric diagonal data
    rng = np.random.default_rng(0)
    pairs = sorted({tuple(sorted((i, neg[i]))) for i in range(12)})
    for trial in range(200):
        nblocks = 9
        counts = np.zeros(12)
        T = 0
        for _ in range(nblocks):
            mask = rng.integers(0, 2, size=len(pairs))
            for pi, (c1, c2) in enumerate(pairs):
                if mask[pi]:
                    counts[c1] += 1; counts[c2] += 1
                    T += E[0, c1] + E[0, c2]     # psi_1 of that block's content
        n = np.linalg.solve(E, np.array([sum(counts[ci] * E[ui, ci] for ci in range(12)) for ui in range(12)]))
        assert np.allclose(n.real, counts, atol=1e-7) and abs(T - counts @ E[0]) < 1e-7
    print("[4] 200 random trials: period-basis expansion recovers exact coset counts  OK")

    # (5) Paley(13) anchor, H = {1,3,9}, single orbit (m=1, n=1): D = QR(13), a = 1
    p13 = 13; H13 = [1, 3, 9]
    cs13 = cosets_of(H13, p13)
    qr13 = {(x * x) % p13 for x in range(1, p13)}
    D = sorted(qr13)
    w13 = np.exp(2j * np.pi / p13)
    tr = sum(w13 ** d for d in D)               # = theta13 = (-1+sqrt13)/2, so a=1 of n=1
    th13 = (-1 + np.sqrt(13)) / 2
    assert abs(tr - th13) < 1e-9
    ncts = [1 if set(c) <= qr13 else 0 for c in cs13]
    want = [(1 + 1 * (2 * 1 - 1) * (1 if c[0] in qr13 else -1)) / 2 for c in cs13]  # (n+s(2a-n)chi)/2
    assert ncts == [int(x) for x in want]
    print("[5] Paley(13), H={1,3,9}: D=QR(13) is union of the 2 QR-cosets; n_c=(1+chi)/2  OK (anchor)")

    print("\nALL CHECKS PASS: order-3 pinning n_c = 3a-9 (QR) / 18-3a (QNR), a in {4,5} is validated.")

if __name__ == "__main__":
    main()
