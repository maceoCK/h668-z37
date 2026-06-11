#!/usr/bin/env python3
"""Verification scan for the 'Profile of a hypothetical H(668)' section.

Checks, by direct computation at v=333 / n=668:
  (1) arithmetic gates (BRC sum-of-two-squares, field-descent self-conjugacy);
  (2) the master parity lemma N == v (mod 4) applied to prime-order automorphisms:
      admissible fixed-point counts f for every odd prime p <= 333, the residue-class
      confinement for p == 3 (mod 4), and the full exclusion list;
  (3) the order-111 corollary (Fix(sigma^37) forced to exactly 111);
  (4) the Kotsireas-Koutschan mod-3 spectrum at ell=333 (16 PSD values, Loeschian
      pairs summing to 167).
Pure stdlib; runs in < 1 s.
"""

def n_order(a, n):
    k, x = 1, a % n
    while x != 1:
        x = (x * a) % n
        k += 1
    return k

def primes(a, b):
    out = []
    for n in range(max(a, 2), b):
        i = 2
        while i * i <= n:
            if n % i == 0:
                break
            i += 1
        else:
            out.append(n)
    return out

V = 333

# (1) arithmetic gates
assert n_order(2, 167) == 83            # odd => 2 not self-conjugate mod 167 (Turyn-Schmidt silent)
assert n_order(83, 37) == 9             # odd order
assert all(pow(83, k, 37) != 36 for k in range(9))   # -1 not in <83> => 83 not self-conjugate mod 37
assert 18 ** 2 + 3 ** 2 == V            # BRC: 333 sum of two squares; 334 == 2 (mod 4)
print("gates: ord_167(2)=83 (odd), ord_37(83)=9 (odd, -1 not in <83>), 333=18^2+3^2, 334=2 mod 4  OK")

# (2) parity lemma N == v (mod 4): prime-order automorphism, f fixed points,
#     orbits of size 1 or p, N = f + (333-f)/p
def admissible_f(p):
    return [f for f in range(V) if (V - f) % p == 0 and (f + (V - f) // p) % 4 == V % 4]

excluded = [p for p in primes(3, V + 1) if not admissible_f(p)]
assert excluded == [p for p in primes(167, V + 1) if p % 4 == 3]
print("excluded odd prime orders (parity lemma):", excluded)

f3 = admissible_f(3)
assert all(f % 6 == 3 for f in f3) and min(f3) == 3      # order 3: f == 3 (mod 6), >= 3
f7 = admissible_f(7)
assert all(f % 14 == 11 for f in f7) and min(f7) == 11   # order 7: f == 11 (mod 14)
assert admissible_f(83) == [1, 167]                      # order 83: f in {1,167}
assert admissible_f(167) == []                           # order 167: impossible
assert admissible_f(37) == [37 * c for c in range(9)]    # parity silent at 37 (pincer theorem needed)
for p in primes(3, V + 1):
    if p % 4 == 1:
        assert len(admissible_f(p)) == len([f for f in range(V) if (V - f) % p == 0])  # parity silent
    elif p % 4 == 3 and admissible_f(p):
        a = admissible_f(p)
        assert all(f % (2 * p) == a[0] % (2 * p) for f in a)  # single residue class mod 2p
print("p=3: f==3 mod 6 (min 3); p=7: f==11 mod 14; p=83: f in {1,167}; p==1 mod 4: parity silent  OK")

# (3) order-111 corollary: sigma of order 111 => Fix(sigma^37) is <sigma^3>-invariant
#     (37 | f) and f == 3 (mod 6), f < 333  => f = 111 exactly
assert [f for f in f3 if f % 37 == 0] == [111]
print("order-111 corollary: Fix(sigma^37) = 111 exactly  OK")

# (4) K-K mod-3 spectrum at ell=333: PSD_u(111) = 4*L, L Loeschian, L_u + L_v = 167
loesch = sorted({x * x + x * y + y * y for x in range(15) for y in range(-14, 15)
                 if 0 <= x * x + x * y + y * y <= 167})
pairs = sorted({(a, 167 - a) for a in loesch if (167 - a) in loesch})
psd = sorted({4 * a for ab in pairs for a in ab})
assert len(psd) == 16 and all(v % 12 == 4 for v in psd)
assert len([p for p in pairs if p[0] <= p[1]]) == 8
print("LP(333) cube-root spectrum: 16 PSD values, all ==4 mod 12, 8 Loeschian pairs:", psd)

print("ALL CHECKS PASS")
