#!/usr/bin/env python3
"""
Field-descent / self-conjugacy sweep for H(668) = 4*167 across ALL decompositions.

Machinery instantiated (exact statements):
 - Self-conjugacy (Leung-Schmidt DCC 2012, Def 2.4): p self-conjugate mod m iff
   exists j with p^j = -1 mod m', where m = p^a m', gcd(p,m')=1.
 - Turyn (Leung-Schmidt 2012, Result 2.5): X in Z[zeta_m], X Xbar = 0 mod t^{2b},
   t self-conjugate mod m  =>  X = 0 mod t^b.  (SINGLE-norm tool.)
 - Lal-McFarland-Odoni / Schmidt-Tan (DCC 2013, Result 2.5): p = 3 mod 4 prime,
   zeta = zeta_{p^a}; |X|^2+|Y|^2 = 0 mod p^{2b} for X,Y in Z[zeta]
   => X = Y = 0 mod p^b.  (TWO-norm tool; needs p^2 | RHS and p-power modulus.)
 - Schmidt field descent (JAMS 1999; Leung-Schmidt 2012 Def 2.10 + Results 2.11/2.12):
   X Xbar = n in Z[zeta_m] => X zeta^j in Z[zeta_F(m,n)]; if coefficients in [0,C],
   then n <= C^2 F(m,n)^2 / (4 phi(F(m,n))).
 - Turyn 1965: circulant Hadamard order must be 4u^2, u odd non-prime-power.

Formulations swept (decomposition -> equation -> (prime, modulus) pairs):
 (A) LP(333) [project-equivalent to H(668)]: |chi(A)|^2+|chi(B)|^2 = 668 in
     Z[zeta_d], d | 333, d>1.  Two-norm.  Primes 2, 167.
 (B) Williamson/Goethals-Seidel over Z_167: sum of four norms = 668 in Z[zeta_167].
 (C) Ito/dicyclic Q_{8*167} via Golay transversal of Z_668: two-norm = 668 in
     Z[i] and Z[zeta_668].
 (D) Golay transversal of Z_2^2 x Z_167: two-norm = 668 incl. RATIONAL characters.
 (E) abelian (668,2,668,334)-RDS in E of order 1336: SINGLE norm = 668 in
     Z[zeta_m], m in {2,4,8,334,668,1336}.
 (F) two-circulant blocks of order 334: two-norm = 668 incl. rational character.
 (G) C(334)/conference: group-developed weight 333; two-circulant order 167:
     two-norm = 333 in Z[zeta_167]; BRC gate 333 = sum of two squares.
 (H) circulant H(668): Turyn 4u^2 gate.
"""
from sympy import primefactors, multiplicity, gcd, totient, n_order, isprime
from sympy import symbols, Poly, QQ, cyclotomic_poly, jacobi_symbol, expand

# ---------------------------------------------------------------- helpers
def selfconj(p, m):
    """Leung-Schmidt 2012 Def 2.4. Returns (bool, m', ord_{m'}(p) or None)."""
    mp = m
    while mp % p == 0:
        mp //= p
    if mp == 1:
        return True, 1, None          # vacuously self-conjugate (only p-part)
    x, seen = p % mp, set()
    found = False
    while x not in seen:
        if x == mp - 1:
            found = True
        seen.add(x)
        x = (x * p) % mp
    o = n_order(p, mp)
    return found, mp, o

def F(m, n):
    """Leung-Schmidt 2012 Definition 2.10."""
    Dm, Dn = sorted(primefactors(m)), sorted(primefactors(n))
    def m_q(q):
        if m % 2 == 1 or q == 2:
            r = 1
            for p in Dm:
                if p != q:
                    r *= p
        else:
            r = 4
            for p in Dm:
                if p not in (2, q):
                    r *= p
        return r
    def ordmod(q, mq):
        if mq == 1:
            return 1
        return n_order(q % mq, mq)
    def b(r):
        qs = [q for q in Dn if q != r]
        if not qs:
            return 2 if r == 2 else 1          # conventions
        vals = []
        for q in qs:
            mq = m_q(q)
            o = ordmod(q, mq)
            vo = multiplicity(r, o) if o > 1 else 0
            if r == 2:
                vals.append(multiplicity(2, q * q - 1) + vo - 1)
            else:
                vals.append(multiplicity(r, q ** (r - 1) - 1) + vo)
        return max(vals)
    prod = 1
    for p in Dm:
        prod *= p ** b(p)
    return int(gcd(m, prod))

def two_squares(n):
    out = []
    a = 0
    while a * a * 2 <= n:
        b2 = n - a * a
        b = int(b2 ** 0.5)
        for bb in (b - 1, b, b + 1):
            if bb >= a and bb * bb == b2:
                out.append((a, bb))
        a += 1
    return out

print("=" * 78)
print("PART 0: arithmetic gates at 668 = 4*167, 334 = 2*167, 333 = 9*37")
print("=" * 78)
print("167 prime:", isprime(167), " 167 mod 4 =", 167 % 4, " 167 mod 8 =", 167 % 8)
print("668 = sum of two squares?", two_squares(668) or "NO (167 = 3 mod 4 to odd power)")
print("334 = sum of two squares?", two_squares(334) or "NO")
print("333 = sum of two squares?", two_squares(333), " (BRC gate for C(334): PASSES)")
print("668 a perfect square?", int(668 ** .5) ** 2 == 668,
      " -> group-developed H(668) impossible (row-sum s^2 = 668)")
print("333 a perfect square?", int(333 ** .5) ** 2 == 333,
      " -> group-invariant W(334,333)/conference impossible (row-sum s^2 = 333)")
print("Turyn 1965 circulant-Hadamard gate: 668/4 = 167 a square?",
      int(167 ** .5) ** 2 == 167, "-> circulant H(668) impossible classically;")
print("   Leung-Schmidt 2012/2016 CHM machinery has NO TARGET at 668 (not 4u^2).")
# four odd squares for Williamson row sums
sols = [(a, b, c, d) for a in range(1, 26, 2) for b in range(a, 26, 2)
        for c in range(b, 26, 2) for d in range(c, 26, 2)
        if a*a + b*b + c*c + d*d == 668]
print("668 as sum of FOUR ODD squares (Williamson row-sum gate):", sols[:4], "PASSES")

print()
print("=" * 78)
print("PART 1: complete self-conjugacy sweep (Def 2.4) over every (prime, modulus)")
print("        pair any 668-formulation needs")
print("=" * 78)
rows = []
# (prime t of the norm RHS, modulus m of the character) ; n=668 needs t in {2,167};
# n=333 (conference side) needs t in {3,37}.
pairs = []
for m in [3, 9, 37, 111, 333]:                 # LP(333) moduli
    pairs += [(2, m, "LP(333) two-norm 668"), (167, m, "LP(333) two-norm 668")]
for m in [167, 334]:                            # Williamson/GS/conference circulants
    pairs += [(2, m, "Williamson/GS 4-norm 668"), (167, m, "Williamson/GS 4-norm 668"),
              (3, m, "conference 2-norm 333"), (37, m, "conference 2-norm 333")]
for m in [2, 4, 8, 334, 668, 1336]:             # RDS / cocyclic single-norm moduli
    pairs += [(2, m, "RDS(668,2,668,334) single-norm 668"),
              (167, m, "RDS(668,2,668,334) single-norm 668")]
seen = set()
for t, m, ctx in pairs:
    if (t, m) in seen:
        continue
    seen.add((t, m))
    ok, mp, o = selfconj(t, m)
    rows.append((t, m, ok, mp, o, ctx))
    print(f"  t={t:>3}  mod m={m:>4}  m'={mp:>4}  ord_{{m'}}(t)={str(o):>4}  "
          f"self-conjugate: {'YES' if ok else 'NO ':>3}   [{ctx}]")

print()
print("key orders:  ord_167(2) =", n_order(2, 167), " (odd -> 2 NOT s.c. mod 167)")
print("             ord_37(2)  =", n_order(2, 37), "; 2^18 mod 37 =", pow(2, 18, 37))
print("             ord_9(2)   =", n_order(2, 9), ";  2^3  mod 9  =", pow(2, 3, 9))
print("             167=19 mod 37; 19^18 mod 37 =", pow(19, 18, 37))
print("             167=5  mod 9;  5^3   mod 9  =", pow(5, 3, 9))
print("             ord_167(3) =", n_order(3, 167), "  ord_167(37) =", n_order(37, 167))
print("             ord_37(83 mod 37=9) =", n_order(9, 37), " (odd, -1 not in <9>)")
print("CRT clash certificates (why s.c. fails at composite 111, 333):")
print("  2:  -1 mod 9 needs j=3 mod 6 (ODD);  -1 mod 37 needs j=18 mod 36 (EVEN)")
print("  167: -1 mod 9 needs j=3 mod 6 (ODD); -1 mod 37 needs j=18 mod 36 (EVEN)")

print()
print("=" * 78)
print("PART 2: which firings are even POSSIBLE -- valuation parity at 668")
print("=" * 78)
print("v_2(668) = 2 (EVEN: Turyn 2.5 gives only X=0 mod 2, no contradiction);")
print("v_167(668) = 1 (ODD: single-norm X Xbar = 668 dies wherever 167 is s.c.")
print("                 and unramified, i.e. m in {2,4,8}; ramified m (167|m) has")
print("                 v_P(668) = e*1 = 166 EVEN -> no parity fire).")
print("Lal-McFarland-Odoni two-norm tool (Schmidt-Tan 2.5): needs p=3 mod 4,")
print("  p-power modulus, and p^{2b}|668 with b>=1: p=167 -> 167^2 =", 167 * 167,
      "| 668?", 668 % (167 * 167) == 0, " -> VACUOUS at 668.")
print("  (p=2 excluded by p=3 mod 4; no LMO analog at p=2, witness below.)")
print("  LMO p=2 counterexample: X=1, Y=zeta_3-zeta_3^2: |X|^2+|Y|^2 = 1+3 = 4 = 0")
print("  mod 4 yet X != 0 mod 2  -> two-norm descent provably has no p=2 version.")

print()
print("=" * 78)
print("PART 3: single-norm fires against sufficient-only routes (abelian RDS etc.)")
print("=" * 78)
# Z[i]: 668 = a^2+b^2 impossible -- shown in Part 0.
# Z[zeta_8]: X Xbar = 668: brute force the Diophantine system
#   X = a+b z+c z^2+d z^3, X Xbar = (a^2+b^2+c^2+d^2) + (ab+bc+cd-da) sqrt2
hits = []
R = 27
for a in range(-R, R + 1):
    for b in range(-R, R + 1):
        s2 = a * a + b * b
        if s2 > 668:
            continue
        for c in range(-R, R + 1):
            s3 = s2 + c * c
            if s3 > 668:
                continue
            for d in range(-R, R + 1):
                if s3 + d * d == 668 and a * b + b * c + c * d - d * a == 0:
                    hits.append((a, b, c, d))
print("Z[zeta_8] norm-668 brute force (|coef|<=27, exhaustive since a^2<=668):",
      "NO SOLUTION" if not hits else hits[:3])
print("  (theory: 167 = -1 mod 8 self-conjugate, f=2, primes above 167 fixed,")
print("   v_P(668)=1 odd -> impossible. Brute force CONFIRMS the theory.)")
print("=> abelian E of order 1336 (Z_8, Z_4xZ_2, Z_2^3 x Z_167): every choice of")
print("   forbidden N=Z_2 admits a linear character of order 2, 4, or 8 nontrivial")
print("   on N (values in Z, Z[i], Z[zeta_8]) -> chi(D)chi(D)bar = 668 impossible")
print("   -> NO abelian (668,2,668,334)-RDS. FIRES (sufficient-only route).")
print("=> two ±1 circulants of order 334 with AA^T+BB^T = 668 I: order-2 character")
print("   gives X^2+Y^2 = 668 over Z: impossible. FIRES (sufficient-only route).")
print("=> Golay transversal of Z_2^2 x Z_167 (any g): some order-2 chi has")
print("   chi(g)=-1, X^2+Y^2=668 over Z: impossible. FIRES (sufficient-only).")

print()
print("=" * 78)
print("PART 4: explicit WITNESSES -- the norm relaxations of every NECESSARY")
print("        formulation are SOLVABLE, so no ideal-theoretic method can fire")
print("=" * 78)
z = symbols('z')
# (i) Eisenstein witness for every modulus divisible by 3 (d = 3, 9, 111, 333)
phi3 = Poly(z * z + z + 1, z, domain=QQ)
X = Poly(3 + z, z, domain=QQ)            # 3 + zeta_3
Xb = Poly(3 + z * z, z, domain=QQ)       # conjugate: zeta_3 -> zeta_3^2
Y = Poly(20 - 9 * z, z, domain=QQ)
Yb = Poly(20 - 9 * z * z, z, domain=QQ)
val = ((X * Xb) + (Y * Yb)) % phi3
print("(i) Z[zeta_3] (embeds in Z[zeta_d], 3|d, conj-equivariantly):")
print("    (3+w)(3+w^2) + (20-9w)(20-9w^2) mod Phi_3 =", val.as_expr(),
      " [= 7 + 661 = 668]")
# (ii) Gauss-sum witness for every modulus divisible by 37 (d = 37, 111, 333)
phi37 = Poly(cyclotomic_poly(37, z), z, domain=QQ)
g = Poly(sum(jacobi_symbol(t, 37) * z ** t for t in range(1, 37)), z, domain=QQ)
g2 = (g * g) % phi37
gbar = Poly(sum(jacobi_symbol(t, 37) * z ** ((37 - t) % 37) for t in range(1, 37)),
            z, domain=QQ) % phi37
print("(ii) Z[zeta_37]: quadratic Gauss sum g = sum (t|37) z^t;  g^2 mod Phi_37 =",
      g2.as_expr(), " ; gbar - g = 0:", (gbar - g % phi37).is_zero)
Xw = Poly(1, z, domain=QQ) + 3 * g       # 1 + 3 sqrt(37)
Yw = Poly(1, z, domain=QQ) - 3 * g       # 1 - 3 sqrt(37)
tot = ((Xw * Xw) % phi37 + (Yw * Yw) % phi37) % phi37
print("    (1+3g)^2 + (1-3g)^2 mod Phi_37 =", tot.as_expr(),
      " [= 2 + 18*37 = 668; both real, totally positive: 334 +/- 6*sqrt37 > 0]")
print("    NOTE: 1+3sqrt37 = -2s, 1-3sqrt37 = -2r for r,s = (-1+-3sqrt37)/2 the")
print("    srg(333,166,82,83) eigenvalues: 668 = (2r)^2 + (2s)^2. The relaxation")
print("    is solved by the eigenvalue shadow of the conjectured object itself.")
# (iii) Z[i] witness (dicyclic/Ito Golay transversal, order-4 characters)
print("(iii) Z[i]: |18+4i|^2 + |18+2i|^2 =", (18**2 + 4**2), "+", (18**2 + 2**2),
      "=", 18**2 + 4**2 + 18**2 + 2**2, " [= 668]")
# (iv) rational gates
print("(iv) Z: LP trivial char 1^2+1^2 = 2; Williamson rows 1+1+225+441 = 668;")
print("     conference 333 = 18^2+3^2 =", 18**2 + 3**2)

print()
print("=" * 78)
print("PART 5: Schmidt F(m,n) field-descent values + F-bound (Result 2.12, C=2)")
print("=" * 78)
for (m, n) in [(3, 668), (9, 668), (37, 668), (111, 668), (333, 668),
               (167, 668), (334, 668), (668, 668), (1336, 668),
               (167, 333), (334, 333)]:
    f = F(m, n)
    bound = 4 * f * f / (4 * totient(f))
    fires = bound < n
    print(f"  F({m:>4},{n}) = {f:>4} ; F-bound C^2F^2/(4phi(F)) = {float(bound):9.1f} "
          f"{'<' if fires else '>='} {n} -> single-norm +-1 solution "
          f"{'EXCLUDED' if fires else 'allowed'}")
print("  (the F-bound 'EXCLUDED' lines only kill single +-1-sequence formulations,")
print("   each ALREADY impossible by Parseval counting; no necessary route is hit.)")
# Wieferich-type checks used inside F
print("  nu_37(2^36-1) =", multiplicity(37, 2**36 - 1),
      "; nu_167(2^166-1) =", multiplicity(167, 2**166 - 1),
      "; 2^166 mod 167^2 =", pow(2, 166, 167 * 167), "(!=1 -> no Wieferich lift)")

print()
print("=" * 78)
print("PART 6: Chan / Leung-Schmidt 'without self-conjugacy' (2012 Result 2.6)")
print("=" * 78)
print("Needs X Xbar = u^2 a PERFECT SQUARE coprime to p. 668 = 2^2*167 is not a")
print("square (and the only forced equations are multi-norm) -> INAPPLICABLE.")
print("Anti-field-descent (Leung-Schmidt JCTA 2016) likewise targets X Xbar = n")
print("single-norm solutions with root-of-unity coefficient structure -> no")
print("forced 668 instance exists. Both certified mute at 668.")
