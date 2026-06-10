"""
Emit the sigma-equivariant expansion template for the C++ annealer (sc_anneal.cpp), reusing the
VALIDATED literal-resolution logic of selfcomp_search (chains, parities, orientations).

Template format (text):
  header: p norb nreps
  reps: lines "REP idx i j type card"   (type 0=offdiag free, 1=cycle-diag symmetric, 2=fixed-diag
        anticlosed; card = |D| target from R)
  lits: 9*9*37 lines "LIT a b d repIdx e neg"  (repIdx=-1: constant 'neg' is the value 0/1)
  rhs:  for each ordered pair a<=b and g: "RHS a b g value"  (srg right-hand side with the
        membership indicator moved to the LHS is NOT done here; value = mu + (k-mu)[a==b,g==0]
        and the (lam-mu)*[g in D] term is handled by the annealer via the literal of (a,b,g))
Usage: python sc_emit_template.py <class_idx> <perm_csv> <outfile>
"""
import sys, itertools
import numpy as np
from selfcomp_search import Model

def main():
    if sys.argv[1] == "--paley13":
        out = sys.argv[2]
        R = [[6]]
        p, norb, c = 13, 1, 5
        k, lam, mu = 6, 2, 3
        perm = (0,)
    else:
        cls_idx = int(sys.argv[1])
        perm = tuple(int(t) for t in sys.argv[2].split(","))
        out = sys.argv[3]
        mats = [np.array([int(t) for t in ln.split()], dtype=np.int64).reshape(9, 9)
                for ln in open("orbit_all_classes.txt")]
        R = [[int(x) for x in row] for row in mats[cls_idx]]
        p, norb, c = 37, 9, 6
        k, lam, mu = 166, 82, 83

    class Shell(Model):
        def _build(self):                       # structure only; vars are their own keys
            for (i, j) in self.orbits:
                for d in range(self.p):
                    if i == j and d == 0: continue
                    self.X[((i, j), d)] = ((i, j), d)
            return True
    M = Shell(p, norb, k, lam, mu, R, perm, c)
    reps = M.orbits
    ridx = {rep: t for t, rep in enumerate(reps)}
    with open(out, "w") as f:
        f.write(f"{p} {norb} {len(reps)} {k} {lam} {mu} {c}\n")
        for rep in reps:
            i, j = rep
            if i == j:
                typ = 2 if perm[i] == i else 1
            else:
                typ = 0
            f.write(f"REP {ridx[rep]} {i} {j} {typ} {R[i][j]}\n")
        for a in range(norb):
            for b in range(norb):
                for d in range(p):
                    if a <= b:
                        v, neg, cst = M.lit(a, b, d)
                    else:
                        v, neg, cst = M.lit(b, a, (-d) % p)
                    if cst is not None:
                        f.write(f"LIT {a} {b} {d} -1 0 {cst}\n")
                    else:
                        rep, e = v                       # var IS its key ((i,j), e)
                        f.write(f"LIT {a} {b} {d} {ridx[rep]} {e} {1 if neg else 0}\n")
    print(f"template written: {out} ({len(reps)} reps)")

if __name__ == "__main__":
    main()
