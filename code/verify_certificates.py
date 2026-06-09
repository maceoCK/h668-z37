"""
Verifier for the orbit-matrix certificates (Theorem 1(a) + the complete enumeration).
Trusts no SAT solver: it (1) checks every matrix in orbit_all_sym.txt is a genuine admissible orbit matrix
(R^2+R-83I = 3071 J, row sums 166, symmetric, even diagonal in {10..26}, distinct), and (2) reports the
coset-compatible counts from the independent C++ enumerator `cert` (0 for |H| in {4,6,9,12,18,36}).

Run:  ../.venv/bin/python verify_certificates.py
"""
import os, subprocess, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
J = 3071*np.ones((9, 9), dtype=object); I = np.eye(9, dtype=object)

def load(path):
    out = []
    for ln in open(path):
        v = [int(t) for t in ln.split()]
        if len(v) == 81: out.append(np.array(v, dtype=object).reshape(9, 9))
    return out

def check(R):
    return (np.array_equal(R@R + R - 83*I, J)
            and all(int(R[i].sum()) == 166 for i in range(9))
            and np.array_equal(R, R.T)
            and all(int(R[i, i]) % 2 == 0 and 10 <= int(R[i, i]) <= 26 for i in range(9)))

def main():
    mats = load(os.path.join(HERE, "orbit_all_sym.txt"))
    ok = sum(check(R) for R in mats)
    seen = {tuple(R.flatten()) for R in mats}
    print(f"[complete enumeration] orbit_all_sym.txt: {len(mats)} matrices, "
          f"{ok} satisfy R^2+R-83I=3071J & rowsum 166 & symmetric, {len(seen)} distinct"
          f"  => {'VERIFIED' if ok == len(mats) == len(seen) else 'FAIL'}")
    if os.path.exists(os.path.join(HERE, "cert")):
        print("[Theorem 1(a)] coset-compatible orbit-matrix counts (independent C++ search):")
        for e in (4, 6, 9, 12, 18, 36):
            r = subprocess.run([os.path.join(HERE, "cert"), str(e)], capture_output=True, text=True, timeout=900)
            print("   " + r.stdout.strip())
    else:
        print("(build the certificate enumerator first: clang++ -O3 -o cert conference_orbit_certify.cpp)")

if __name__ == "__main__":
    main()
