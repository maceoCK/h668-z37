"""
Self-complementary Z37-symmetric srg(333,166,82,83): twist-consistency analysis.

Ansatz: srg with fixed-point-free rho of order 37 (9 orbits) and a complementing isomorphism
sigma: Gamma -> complement(Gamma) normalizing <rho>:  sigma(i,g) = (pi(i), c*g + b_i),
with pi a complementing permutation of the orbit matrix (R-complement = P R P^T), c^4 = 1.

Derived (and re-verified here):
  * c must have order exactly 4 (c in {6,31}; c=+-1 contradicts the fixed-orbit diagonal
    D[f][f] = Z37^* \\ c D[f][f]);  WLOG c = 6 (sigma <-> sigma^3 swaps (pi,6) <-> (pi^{-1},31)).
  * Block relations: for all i,j:  D[pi(i)][pi(j)] = co( 6*D[i][j] + (b_j - b_i) ),
    where co() is complement in Z37 (off-diagonal) / Z37^* (diagonal).
  * Going around pi-orbits forces 4 linear conditions on b in F37:
      per 4-cycle (a0 a1 a2 a3):   U := c^3 b_{a0} + c^2 b_{a1} + c b_{a2} + b_{a3} = 0
      per 4-cycle (antipodal):     s := c (b_{a2} - b_{a0}) + (b_{a3} - b_{a1}) = 0
    (the antipodal blocks D[a0][a2] are then FULLY FREE; translation-stabilizer argument
    kills any s != 0 since 37 is prime).
  * Gauge: re-coordinatization g -> g - t_i per orbit gives b_i ~ b_i + (c t_i - t_{pi(i)}).

This script: for each of the 3 self-complementary classes x complementing permutation,
build the twist system, compute solution space and gauge image over F37, and decide whether
b = 0 is WLOG (gauge acts transitively on twist solutions).  Also verifies the forced
diagonal pattern (d, 36-d, d, 36-d) along cycles and diag 18 at the fixed orbit.
"""
import numpy as np
import itertools

P = 37
C = 6                                     # 6^2 = 36 = -1 mod 37, order 4

def load_classes(path="orbit_all_classes.txt", idxs=(109, 221, 422)):
    mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
            for ln in open(path)]
    return {k: mats[k] for k in idxs}

def complementing_perms(R):
    Rc = np.where(np.eye(9, dtype=bool), 36 - R, 37 - R)
    out = []
    for p in itertools.permutations(range(9)):
        if np.array_equal(Rc, R[np.ix_(p, p)]):
            out.append(p)
    return out

def cycles_of(p):
    seen, cyc = set(), []
    for s in range(9):
        if s in seen: continue
        c, x = [], s
        while x not in seen:
            seen.add(x); c.append(x); x = p[x]
        cyc.append(c)
    return cyc

def nullspace_F37(A):
    """row-reduce A over F37; return basis of null space (list of vectors)."""
    A = [row[:] for row in A]
    m, n = len(A), len(A[0])
    piv = []
    r = 0
    for col in range(n):
        sel = next((i for i in range(r, m) if A[i][col] % P), None)
        if sel is None: continue
        A[r], A[sel] = A[sel], A[r]
        inv = pow(A[r][col], P - 2, P)
        A[r] = [(v * inv) % P for v in A[r]]
        for i in range(m):
            if i != r and A[i][col] % P:
                f = A[i][col]
                A[i] = [(A[i][k] - f * A[r][k]) % P for k in range(n)]
        piv.append(col); r += 1
        if r == m: break
    free = [c for c in range(n) if c not in piv]
    basis = []
    for fc in free:
        v = [0] * n; v[fc] = 1
        for ri, pc in enumerate(piv):
            v[pc] = (-A[ri][fc]) % P
        basis.append(v)
    return basis, piv

def colspace_rank(M):
    _, piv = nullspace_F37([list(r) for r in zip(*M)])  # rank via transpose nullspace pivots
    return len(piv)

def main():
    classes = load_classes()
    for k, R in classes.items():
        perms = complementing_perms(R)
        print(f"class {k}: diag {[int(R[i,i]) for i in range(9)]}, {len(perms)} complementing perms")
        for p in perms:
            cyc = cycles_of(p)
            four = [c for c in cyc if len(c) == 4]
            fixed = [c[0] for c in cyc if len(c) == 1]
            assert len(four) == 2 and len(fixed) == 1, f"unexpected cycle type {[len(c) for c in cyc]}"
            f = fixed[0]
            # forced diagonal structure checks
            ok_diag = int(R[f, f]) == 18
            for a0, a1, a2, a3 in four:
                ok_diag &= int(R[a1, a1]) == 36 - int(R[a0, a0])
                ok_diag &= int(R[a2, a2]) == int(R[a0, a0])
                ok_diag &= int(R[a3, a3]) == int(R[a1, a1])
            # twist system: 4 conditions on b in F37^9
            A = []
            for a0, a1, a2, a3 in four:
                u = [0] * 9
                u[a0] = pow(C, 3, P); u[a1] = pow(C, 2, P); u[a2] = C; u[a3] = 1
                A.append(u)
                s = [0] * 9
                s[a0] = (-C) % P; s[a2] = C; s[a1] = P - 1; s[a3] = 1
                A.append(s)
            sol_basis, piv = nullspace_F37([row[:] for row in A])
            sol_dim = len(sol_basis)
            # gauge image: b_i ~ b_i + (c t_i - t_{pi(i)}): matrix G (9x9), columns indexed by t
            G = [[0] * 9 for _ in range(9)]
            for i in range(9):
                G[i][i] = (G[i][i] + C) % P
                G[i][p.index(i) if False else 0] = G[i][0]  # placeholder, fixed below
            # rebuild G correctly: delta b_i = c*t_i - t_{pi(i)}
            G = [[0] * 9 for _ in range(9)]
            for i in range(9):
                G[i][i] = (G[i][i] + C) % P
                G[i][p[i]] = (G[i][p[i]] - 1) % P
            # check gauge image lies inside the solution space (A @ G == 0)
            AG = [[sum(A[r][i] * G[i][t] for i in range(9)) % P for t in range(9)] for r in range(len(A))]
            gauge_in_sol = all(v == 0 for row in AG for v in row)
            g_rank = colspace_rank(G)
            # b = 0 WLOG iff gauge image == solution space (then quotient trivial)
            wlog = gauge_in_sol and (g_rank == sol_dim)
            print(f"  perm {p}: cycles {four} fixed {f}  diag-pattern-ok={ok_diag}  "
                  f"twist-solution-dim={sol_dim}  gauge-rank={g_rank}  gauge-in-solutions={gauge_in_sol}  "
                  f"=> b=0 WLOG: {wlog}")

if __name__ == "__main__":
    main()
