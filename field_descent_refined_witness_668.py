#!/usr/bin/env python3
"""
PART 7: parity-REFINED relaxation witnesses.

For LP(333), at a character chi of order d | 333, chi(A) = sum_j c_j zeta_d^j
with c_j = (sum of 333/d entries +-1) = odd for every j. Writing the value in
the power basis (zeta^0..zeta^{d-2}) gives ALL EVEN coefficients, i.e.
chi(A) = 2X, chi(B) = 2Y with X, Y in Z[zeta_d].  Hence LP(333) forces

        X Xbar + Y Ybar = 167   in Z[zeta_d]   for every d | 333, d > 1.

(This is the sharpest 'field-descent-visible' consequence: v_167 = 1 is odd and
167 is self-conjugate mod 3, 9, 37 -- a SINGLE norm 167 would be impossible
there. The two-norm version survives only if explicit solutions exist.)

d = 3 (and 9, 111, 333 by embedding):  167 = N(2) + N(14+3w) = 4 + 163. Verify.
d = 37: search the quartic-period subring (the imaginary quartic CM subfield
K4 of Q(zeta_37); sigma_{-1} is not in the index-4 Galois subgroup since
18 != 0 mod 4, so K4 is CM and x xbar is a genuine 2-norm there).
Vectors u in Z^37 constant on the four quartic cosets + free value at 0.
"""
import numpy as np
from itertools import product
from sympy import symbols, Poly, QQ, cyclotomic_poly

# ---------- d = 3 refined witness, exact
z = symbols('z')
phi3 = Poly(z * z + z + 1, z, domain=QQ)
X, Xb = Poly(2, z, domain=QQ), Poly(2, z, domain=QQ)
Y = Poly(14 + 3 * z, z, domain=QQ)
Yb = Poly(14 + 3 * z * z, z, domain=QQ)
v3 = ((X * Xb) + (Y * Yb)) % phi3
print("d=3 refined: 2*2bar + (14+3w)(14+3w^2) mod Phi_3 =", v3.as_expr(),
      "  [4 + 163 = 167]")
print("  => full witness at d=3: (4)(4) + (28+6w)(28+6w^2) = 16 + 652 = 668,")
print("     all coefficients EVEN, congruent to +-1 mod (1-w) as LP requires.")

# ---------- d = 37: quartic cosets mod 37
p = 37
g = 2
# find a primitive root mod 37
def is_proot(a):
    s, x = set(), 1
    for _ in range(36):
        x = x * a % p
        s.add(x)
    return len(s) == 36
while not is_proot(g):
    g += 1
H = sorted(pow(g, 4 * k, p) for k in range(9))            # quartic residues
cosets = [sorted(pow(g, i, p) * h % p for h in H) for i in range(4)]
print("d=37: primitive root", g, "; quartic cosets sizes", [len(c) for c in cosets])

def vec(c0, a):
    u = np.zeros(p, dtype=np.int64)
    u[0] = c0
    for i, C in enumerate(cosets):
        for t in C:
            u[t] = a[i]
    return u

def acorr(u):
    return np.array([int(np.dot(u, np.roll(u, -k))) for k in range(p)],
                    dtype=np.int64)

# coset index of each k != 0
cidx = {}
for i, C in enumerate(cosets):
    for t in C:
        cidx[t] = i

R = range(-3, 4)
entries = {}     # key: (d1,d2,d3) of coset-values minus coset0-value -> list
cands = []
for c0, a0, a1, a2, a3 in product(R, R, R, R, R):
    u = vec(c0, (a0, a1, a2, a3))
    A = acorr(u)
    cv = [A[cosets[i][0]] for i in range(4)]          # 4 coset values
    # sanity: A constant on cosets
    peak_minus = int(A[0] - cv[0])
    key = (cv[1] - cv[0], cv[2] - cv[0], cv[3] - cv[0])
    entries.setdefault(key, []).append((peak_minus, (c0, a0, a1, a2, a3)))
    cands.append((key, peak_minus, (c0, a0, a1, a2, a3)))

found = []
for key, pm, par in cands:
    negkey = (-key[0], -key[1], -key[2])
    for pm2, par2 in entries.get(negkey, []):
        if pm + pm2 == 167:
            found.append((par, par2))
            break
    if found:
        break

if found:
    (pu, pv) = found[0]
    u = vec(pu[0], pu[1:])
    v = vec(pv[0], pv[1:])
    Au, Av = acorr(u), acorr(v)
    S = Au + Av
    print("WITNESS FOUND for 167 = X Xbar + Y Ybar in Z[zeta_37]:")
    print("  X coeff vector (index 0, then coset values):", pu)
    print("  Y coeff vector:", pv)
    print("  joint autocorrelation: peak", int(S[0]), ", off-peak values",
          sorted(set(int(x) for x in S[1:])))
    print("  => value = peak - offpeak =", int(S[0] - S[1]), " [should be 167]")
    # exact symbolic confirmation
    phi37 = Poly(cyclotomic_poly(37, z), z, domain=QQ)
    Xp = Poly(sum(int(u[i]) * z ** i for i in range(p)), z, domain=QQ)
    Xpb = Poly(sum(int(u[i]) * z ** ((p - i) % p) for i in range(p)), z, domain=QQ)
    Yp = Poly(sum(int(v[i]) * z ** i for i in range(p)), z, domain=QQ)
    Ypb = Poly(sum(int(v[i]) * z ** ((p - i) % p) for i in range(p)), z, domain=QQ)
    tot = ((Xp * Xpb) % phi37 + (Yp * Ypb) % phi37) % phi37
    print("  exact check: X Xbar + Y Ybar mod Phi_37 =", tot.as_expr())
else:
    print("no witness in coset-structured family with |coef|<=3 -- widen search")
