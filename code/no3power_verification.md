# Verification report: no srg(333,166,82,83) admits a fixed-point-free automorphism of 3-power order

**Verdict: PROVEN** (complete, gap-free proof below; machine-checked in `no3power_check.py`,
all checks pass, 0.3 s). The corollary — every **vertex-transitive** srg(333,166,82,83) admits a
fixed-point-free automorphism of order exactly 37 — is **PROVEN modulo CFSG** (via
Fein–Kantor–Schacher); its further consequence (coverage by the order-37 orbit-matrix
enumeration) additionally cites the project's completeness result for the 2901 canonical orbit
matrices (FINDINGS §6).

The earlier round's *claimed mechanism* needed repair. Its route (a) (cube-root-of-unity
eigenspace traces) is **vacuous** — the trace identity reduces to 0 = 0 and constrains nothing
(§6.1, machine-checked). Its "Cayley-lemma-style ζ₃ field obstruction" does **not** apply to a
bare semiregular Z₃ action (no group development, no linear character to evaluate; §6.2). The
*intuition* behind the rationality route is nevertheless vindicated: the proof that works is a
**rationality-forced trace parity** argument. Irreducibility of x² + x − 83 over **Q** (because
333 is not a square) forces the θ- and τ-multiplicities in the orbit-quotient characteristic
polynomial to be *equal*, which pins the quotient trace to the exact integer (333 − N)/2; the
handshake lemma forces that trace to be *even*; for fixed-point-free 3-power actions the two
collide. The same collision, run on σ and σ³ simultaneously, kills all orders 9, 27, 81, 243
(where the naive single-σ count does **not** suffice — e.g. orbit type (t, u) = (108, 1) for
order 9 passes the single constraint; the σ³ step is genuinely needed).

---

## Theorem

Let Γ be a strongly regular graph with parameters (v, k, λ, μ) = (333, 166, 82, 83) and let
σ ∈ Aut(Γ) with σ^(3^e) = id, σ ≠ id, e ≥ 1, such that σ has **no fixed vertex**. Then no such σ
exists. Equivalently: every automorphism of Γ of 3-power order fixes at least one vertex.

*(As with the order-37 theorem, the statement is conditional on the existence of Γ, which is
open — equivalent to the open Hadamard order H(668). The proof is valid for any such Γ.)*

In fact the proof gives, with no extra work:

**Theorem (general form).** Let Γ be a conference graph on v vertices with v **not a perfect
square**, and let p be a prime with **p ≡ 3 (mod 4)**. Then Γ admits no fixed-point-free
automorphism of order p^e for any e ≥ 1.

(333 ≡ 1 (mod 4), 333 = 9·37 nonsquare, p = 3 ≡ 3 (mod 4): the title theorem is the special
case. Paley(9) and Paley(49) show the nonsquare hypothesis is necessary: their translations are
fixed-point-free of orders 3 and 7. See §5.)

## 0. Standing facts

A = adjacency matrix of Γ, J all-ones, I identity.

**(F1)** A² + A − 83I = 83J; Γ is connected (μ = 83 > 0), 166-regular, k = 166 simple
(Perron–Frobenius).

**(F2)** char(A) = (x − 166)(x² + x − 83)^166 exactly. Conference condition
2k + (v−1)(λ−μ) = 332 − 332 = 0 gives equal multiplicities (v−1)/2 = 166 for
θ = (−1+3√37)/2 and τ = (−1−3√37)/2. The quadratic q(x) = x² + x − 83 is **irreducible over Q**
(its discriminant is 333 = 9·37, not a perfect square). Also 166² + 166 − 83 = 27639 ≠ 0, so 166
is not a root of q.

**(F3) Equitability of orbit partitions.** For any H ≤ Aut(Γ), the orbit partition of H is
equitable: if u, w ∈ Cᵢ, pick h ∈ H with w = hu; then N(w) ∩ Cⱼ = h(N(u) ∩ Cⱼ), so
|N(w) ∩ Cⱼ| =: Q_{ij} is well defined. With S the v × N indicator matrix of the N orbits
(full column rank): **AS = SQ**.

**(F4) Quotient divisibility** (Godsil–Royle, Thm 9.3.3): char(Q) divides char(A).
(Independent re-derivation used as a machine cross-check: from AS = SQ, p(A)S = S p(Q) for any
polynomial p; taking p = (x−166)q gives p(Q) = 0 since p(A) = 0 and S is injective; this
annihilator is squarefree, so Q is diagonalizable with eigenvalues in {166, θ, τ}, and
rank(Q − 166I) = N − 1 pins the 166-multiplicity to 1.)

**(F5) Handshake lemma.** A finite graph has an even number of odd-degree vertices; in
particular a d-regular graph on an odd number of vertices has d even.

## 1. Master Parity Lemma

**Lemma.** Let H ≤ Aut(Γ) be any subgroup of **odd order**, and let N be the number of H-orbits
on V(Γ). Then **N ≡ 1 (mod 4)**.

*Proof.* Let Q be the N × N quotient matrix of the H-orbit partition (F3).

(i) *Shape of char(Q).* By (F4), char(Q) divides char(A) = (x−166)(x²+x−83)^166. Since
q(x) = x²+x−83 is irreducible over Q (F2) and char(Q) is monic of degree N, the only possible
divisors are char(Q) = (x−166)^a · q(x)^b with a ∈ {0,1}, a + 2b = N. Since Γ is 166-regular,
every row of Q sums to 166, so Q·**1** = 166·**1** and 166 ∈ spec(Q); as 166 is not a root of q
(F2), a = 1. Hence N = 2b + 1 is odd and

  **tr(Q) = 166 + b(θ + τ) = 166 − b = 166 − (N−1)/2 = (333 − N)/2.**  (1.1)

(This is the rationality step: irreducibility of q forces the θ- and τ-multiplicities in
char(Q) to be *equal* (= b), which is what makes tr(Q) a forced function of N alone.)

(ii) *Parity of the diagonal.* Fix an orbit Cᵢ. Its size divides |H| (orbit–stabilizer), hence
is **odd**. H acts on the induced subgraph Γ[Cᵢ] transitively by graph automorphisms, so Γ[Cᵢ]
is vertex-transitive, hence regular — of degree exactly Q_{ii} = |N(u) ∩ Cᵢ|. By the handshake
lemma (F5), a regular graph of odd order has even degree:

  **every Q_{ii} is even, so tr(Q) = Σᵢ Q_{ii} is even.**  (1.2)

(iii) Combining (1.1) and (1.2): (333 − N)/2 ≡ 0 (mod 2), i.e. N ≡ 333 ≡ **1 (mod 4)**. ∎

*Remarks.* (a) Singleton orbits (fixed points of H) contribute Q_{ii} = 0, which is even; the
lemma needs **no** fixed-point-freeness hypothesis. (b) For a general conference graph on
nonsquare v the same proof gives tr(Q) = (v − N)/2 and **N ≡ v (mod 4)**; conference graphs have
v ≡ 1 (mod 4) always. (c) Both hypotheses are sharp: |H| odd is needed for odd orbits (an even
orbit can carry an odd-degree induced circulant), and v nonsquare is needed for step (i) — see
the Paley(9)/Paley(49) controls, §5.

## 2. Proof of the Theorem

Let σ ∈ Aut(Γ) be fixed-point-free with σ^(3^e) = id, σ ≠ id; WLOG ord(σ) = 3^e, e ≥ 1. Every
⟨σ⟩-orbit has size 3^j with **j ≥ 1** (size 1 is excluded by fixed-point-freeness).

**Case e = 1.** All orbits have size 3, so N(σ) = 333/3 = 111 ≡ **3 (mod 4)**, contradicting the
Master Lemma applied to H = ⟨σ⟩ (order 3, odd). ∎

**Case e ≥ 2.** Consider σ³ (order 3^{e−1} ≥ 3, odd). Let O be a ⟨σ⟩-orbit, |O| = 3^j, j ≥ 1.
Since ⟨σ⟩ is abelian and transitive on O, all point stabilizers in O coincide and equal the
unique index-3^j subgroup ⟨σ^{3^j}⟩. For x ∈ O, the ⟨σ³⟩-orbit of x has size
[⟨σ³⟩ : ⟨σ³⟩ ∩ ⟨σ^{3^j}⟩] = 3^{e−1}/3^{e−j} = 3^{j−1} (here ⟨σ^{3^j}⟩ ⊆ ⟨σ³⟩ because j ≥ 1).
Hence **O splits into exactly 3^j / 3^{j−1} = 3 orbits of ⟨σ³⟩**, and so

  **N(σ³) = 3 · N(σ).**  (2.1)

(Note: σ³ may well have *fixed points* — every size-3 σ-orbit contributes three of them — and
that is fine: the Master Lemma does not care.)

Apply the Master Lemma twice: to H = ⟨σ⟩, N(σ) ≡ 1 (mod 4); to H = ⟨σ³⟩,
N(σ³) = 3N(σ) ≡ 3·1 = 3 (mod 4) — but the lemma demands N(σ³) ≡ 1 (mod 4). Contradiction. ∎

This handles **all** mixed orbit-size profiles at once. The σ³ step is *necessary*: e.g. for
order 9 the orbit type (t, u) = (108 three-orbits, 1 nine-orbit) has N = 109 ≡ 1 (mod 4) and
passes the single-σ constraint; only N(σ³) = 327 ≡ 3 (mod 4) kills it. The script checks all
37 + 222 + 260 + 28 orbit-type vectors for orders 9, 27, 81, 243 exhaustively.

For the general form of the theorem (prime p ≡ 3 (mod 4), v ≡ 1 (mod 4) nonsquare): case e = 1
gives N = v/p ≡ p^{−1} ≡ 3 (mod 4) (since 3·3 = 9 ≡ 1), and case e ≥ 2 gives
N(σ^p) = pN(σ) ≡ 3 (mod 4) by the identical splitting argument. Same contradictions. ∎

**Byproduct (order-3 fixed-point structure).** For *any* order-3 automorphism of Γ with f fixed
points and s′ = (333 − f)/3 three-orbits, N = f + s′ = 333 − 2s′ ≡ 1 (mod 4) forces s′ even,
i.e. **f ≡ 3 (mod 6)**; in particular f ≥ 3. (Consistent with the prompt's observation that the
fpf 111-orbit quotient with char (x−166)(x²+x−83)^55 is "dimension-wise consistent": dimensions
are fine — the kill is the parity of its trace, 111, against the even diagonal.)

## 3. Corollary (vertex-transitive case)

**Corollary.** Assume Γ is vertex-transitive. Then Γ admits a fixed-point-free automorphism of
order **exactly 37** (whose 9 orbits of size 37 realize the Z₃₇ orbit-matrix framework). In
particular, every vertex-transitive srg(333,166,82,83) appears among the lifts of the complete
order-37 orbit-matrix enumeration (the 2901 canonical 9×9 matrices, FINDINGS §6).

*Proof.* G = Aut(Γ) acts transitively on 333 > 1 points. By the theorem of
**Fein–Kantor–Schacher** [B. Fein, W. M. Kantor, M. Schacher, *Relative Brauer groups II*,
J. Reine Angew. Math. 328 (1981), 39–57, Thm 1; proof depends on CFSG], every finite transitive
permutation group of degree > 1 contains a fixed-point-free element g of **prime-power** order
q^m. (Prime-power, not prime, is essential in general — elusive groups like M₁₁ on 12 points
have no fpf element of prime order; here this subtlety evaporates, see below.) Since g is
fixed-point-free, every ⟨g⟩-orbit has size q^i with i ≥ 1; the sizes sum to 333, so
q | 333 = 3²·37, i.e. q ∈ {3, 37}.

- q = 3 is impossible: g would be a fixed-point-free automorphism of 3-power order,
  contradicting the Theorem.
- Hence q = 37. Every cycle of g has 37-power length ≤ 333 < 37² = 1369, so every cycle has
  length exactly 37; thus ord(g) = lcm of cycle lengths = 37, with 333/37 = 9 cycles.

So g is a fixed-point-free automorphism of order exactly 37 with 9 orbits of size 37. Choosing a
base point in each orbit, A becomes a 9×9 array of 37×37 circulants, and the orbit quotient
R (R_{ij} = |N(w) ∩ Oⱼ|, w ∈ Oᵢ) is symmetric (37R_{ij} = 37R_{ji} edge count), has row sums
166, even diagonal (Lemma §1(ii)), and satisfies R² + R − 83I = 83·37·J = 3071J (multiply (F1)
by S and use JS = 37·SJ₉) — precisely the orbit-matrix equations whose solution set was
enumerated completely up to equivalence (2901 canonical representatives, FINDINGS §6 /
`central_z111_skeletons.txt` pipeline). ∎

*Dependency ledger for the corollary:* (i) CFSG via FKS; (ii) the Theorem above (unconditional);
(iii) completeness of the 2901-representative enumeration (established project result, cited not
re-proved here). Note the prior theorem "every order-37 automorphism is fixed-point-free"
(`fpf_verification.md`) is **not needed** for this chain — FKS hands us fixed-point-freeness
directly — but it remains what makes the framework exhaustive for *all* (not necessarily
vertex-transitive) Z₃₇-symmetric candidates.

## 4. Why Paley(9) does not contradict the theorem — the exact special feature

Paley(9) = srg(9,4,1,2) is a conference graph with fixed-point-free order-3 automorphisms
(translations), and Paley(49) = srg(49,24,11,12) likewise has fpf order-7 translations
(7 ≡ 3 mod 4). The obstruction *must* fail there, and it fails at exactly one point:

**v is a perfect square ⇒ x² − (λ−μ)x − (k−μ) is reducible ⇒ step (i) of the Master Lemma
collapses.** The restricted eigenvalues are rational, char(Q) need not allot them equal
multiplicities, and tr(Q) is no longer forced to be (v−N)/2. Concretely (machine-checked):

| graph | fpf σ | N | N mod 4 | forced trace (v−N)/2 | actual tr(Q) | char(Q) |
|---|---|---|---|---|---|---|
| Paley(9), x↦x+1 | order 3 | 3 | 3 | 3 (odd — impossible if forced) | 6 (even ✓) | (x−4)(x−1)² — multiplicities (2,0) ≠ |
| Paley(49), x↦x+1 | order 7 | 7 | 3 | 21 (odd — impossible if forced) | 42 (even ✓) | (x−24)(x−3)⁶ — multiplicities (6,0) ≠ |

Every *other* ingredient (equitability, row sums, annihilator, rank, even diagonal/handshake,
even tr(Q)) still holds on these controls — verified — so the irrational eigenvalue field
**Q(√v) ⊄ Q** is precisely the load-bearing hypothesis, exactly as the earlier round guessed
("the obstruction presumably needs irrational sqrt(v)"). For srg(333,·): √333 = 3√37 irrational.

## 5. Sharpness — the lemma does not over-kill

The Master Lemma must *permit* the symmetries we know are consistent, and it does:

- **fpf order 37 on Γ:** N = 9 ≡ 1 (mod 4) ✓ (the Z₃₇ framework survives; the verified orbit
  matrix of `conf_orbit3.out` has even diagonal (10,…,24) summing to tr = 162 = (333−9)/2 ✓).
- **order 3 with f ≡ 3 (mod 6) fixed points on Γ:** N ≡ 1 (mod 4) ✓ — not excluded (and must
  not be: nothing known forbids such automorphisms).
- Real conference graphs, all machine-verified to satisfy N ≡ 1 (mod 4) with the full forced
  char(Q): Paley(13)/x↦3x (order 3, N=5), Paley(37)/x↦10x (order 3, N=13), Paley(29)/x↦16x
  (order 7, N=5), Paley(53)/x↦10x (order 13, N=5), Paley(73)/x↦2x (**order 9**, N=9; and its
  cube x↦8x with N=25 = 3·9 − 2·1, confirming the splitting bookkeeping N(σ³) = 3N(σ) − 2n₁ and
  the exactly-3-suborbits claim), Paley(125)/Frobenius (order 3, f=5, N=45).

## 6. Autopsy of the earlier round's sketch

**6.1 Route (a) — ζ₃ eigenspace traces — is vacuous.** For σ of order 3 with f fixed points:
the permutation matrix P has eigenvalue multiplicities f + (333−f)/3 at 1 and (333−f)/3 at each
of ζ₃, ζ̄₃. V_θ and V_τ are real σ-invariant subspaces, so within each the ζ₃ and ζ̄₃
multiplicities coincide (b_θ, b_τ say): tr(σ|V_θ) = a − b_θ with a + 2b_θ = 166, likewise for τ;
the Perron line is σ-fixed. Then tr(P) = 1 + (a − b_θ) + (a′ − b_τ) = f reduces, after
substituting the multiplicity bookkeeping 1 + a + a′ = f + (333−f)/3, to **0 = 0 identically**
(machine-checked symbolically, Part 2 of the script). The proposed Galois refinement cannot help:
all the eigenspace traces here are already rational integers, and complex conjugation and
√37 ↦ −√37 act exactly as the bookkeeping above already encodes. No constraint on f survives.

**6.2 The "Cayley lemma" route does not transfer.** The project's Cayley/PDS lemma kills
group-developed graphs because a *linear character* of order 3 evaluated on a difference set
would put a root of x²+x−83 into Z[ζ₃] — impossible as Q(√37) ⊄ Q(ζ₃). A bare semiregular Z₃
action provides no group development and no character sum to evaluate; the eigenvalues of A on
each P-isotypic component are still θ, τ (the components are not A-invariant in any useful
rational way beyond §6.1). The *spirit* survives, though: in the actual proof, rationality
enters as "char(Q) ∈ Z[x] + irreducibility ⇒ equal θ/τ multiplicities", which is the legitimate
descendant of the field-obstruction idea for bare actions.

**6.3 The mod-3 free-module route was not needed.** The suggested F₃[Z₃]-module attack
((A−I)² ≡ 2J mod 3, freeness of F₃^333 over the local ring F₃[Z₃]) was not pursued: the parity
pincer closes everything at order 3, and the σ³ trick lifts it to all 3-powers. It remains a
plausible source of *additional* constraints (e.g. on order-3 automorphisms with fixed points)
but nothing here depends on it.

**6.4 What the final proof rests on** (hypothesis checklist):
- equitability of orbit partitions (F3): two lines, no hypotheses ✓
- char(Q) | char(A) (F4): S full column rank ✓ (orbits nonempty, disjoint)
- char(A) = (x−166)(x²+x−83)^166 (F2): connectivity (μ>0), regularity, conference condition ✓
- irreducibility of x²+x−83 over Q: 333 nonsquare ✓ — **the engine**
- 166 not a root of q: 27639 ≠ 0 ✓ (forces a = 1)
- handshake on odd vertex-transitive induced subgraphs: |H| odd ✓
- orbit splitting (2.1): ⟨σ⟩ abelian, j ≥ 1 (fixed-point-freeness) ✓
- corollary only: FKS (CFSG) + enumeration completeness (cited) ✓
- integrality/interlacing/Cauchy–Schwarz/Benson counts: **not used**.

## 7. Machine verification (`no3power_check.py`)

Exact arithmetic (sympy/Fraction) + numpy, 0.3 s, exit 0:

- **Part 0:** conference condition; 333 nonsquare; 166 not a root; tr(Q) = (333−N)/2 symbolic;
  the equivalence {(333−N)/2 even ⟺ N ≡ 1 (mod 4)} for all odd N ≤ 333; the e = 1 kill
  (111 ≡ 3 mod 4); **exhaustive** kills of all orbit-type vectors for fpf orders 9 (37 types),
  27 (222), 81 (260), 243 (28) via the pair {N ≡ 1, 3N ≡ 1}, reporting how many types pass the
  single-σ test (19/108/132/13 — the σ³ step is necessary); the f ≡ 3 (mod 6) byproduct; the
  FKS arithmetic (q ∈ {3,37}; q = 37 forces order exactly 37, 9 cycles).
- **Part 1:** every lemma validated end-to-end on real conference graphs with real odd-order
  automorphisms (§5 list), including the order-9 example Paley(73)/x↦2x with its cube, the
  splitting identities, and direct char(Q) computation up to N = 25; plus the two **controls**
  Paley(9) and Paley(49) with genuine fpf automorphisms of order ≡ 3 (mod 4), where the script
  verifies the master conclusion fails (N ≡ 3 mod 4) and that the *only* broken hypothesis is
  irreducibility (disc = v a perfect square; unequal eigenvalue multiplicities (2,0) and (6,0)
  exhibited), while equitability, annihilator, rank, even diagonal, even trace all still hold.
- **Part 2:** symbolic sympy proof that route (a) reduces to residual 0.

Result: **all checks pass** (full output reproduced at the end of this file).

## 8. Consequences for the project

1. The Z₃-sector of the symmetry landscape is closed: any 3-group of automorphisms of a putative
   srg(333,166,82,83) has a fixed vertex (a fixed-point-free element would have 3-power order).
   Combined with `fpf_verification.md` (order 37 ⇒ fpf), the two prime divisors of 333 now have
   opposite, fully proven behaviors: **37-elements are always fixed-point-free; 3-elements never
   are.**
2. Modulo CFSG, the vertex-transitive case is *completely* reduced to the existing 2901-orbit-
   matrix enumeration: if none of the 2901 skeletons lifts, **no vertex-transitive
   srg(333,166,82,83) exists** (equivalently, no vertex-transitive conference graph of order
   333, hence no such doubly regular tournament-derived H(668) of that symmetry type).
3. The general form (conference graph, v nonsquare, p ≡ 3 (mod 4) ⇒ no fpf p-power automorphism)
   is a clean standalone statement for the paper; the orbit-count parity N ≡ v (mod 4) for
   odd-order groups sits naturally in (and may be folklore within) the Behbahani–Lam orbit-matrix
   tradition, but we have not located this exact statement in the literature — cite as a lemma
   with the self-contained proof above.

---

### Appendix: script output (full, exit code 0, runtime 0.3 s)

```
PART 0: master-lemma arithmetic for srg(333,166,82,83)
  [PASS] conference condition 2k+(v-1)(lam-mu)=0 (multiplicities 166/166)
  [PASS] disc(x^2+x-83) = 333 = 9*37 is NOT a perfect square
  [PASS] 166 is not a root of x^2+x-83 (so a=1 in char(Q))
  [PASS] tr(Q) = 166-(N-1)/2 = (333-N)/2 (symbolic)
  [PASS] for all odd N in [1,333]: (333-N)/2 even <=> N = 1 (mod 4)
  [PASS] fpf order 3: N = 333/3 = 111 = 3 (mod 4)  ==> KILLED by master lemma
  exhaustive orbit-type check, fpf sigma of order 3^e (orbit sizes 3^j, j>=1):
  [PASS] order 3^2=9: all 37 orbit types killed by {N=1 (4)} AND {3N=1 (4)}  (19 types pass the single sigma-constraint (sigma^3 step genuinely needed))
  [PASS] order 3^3=27: all 222 orbit types killed by {N=1 (4)} AND {3N=1 (4)}  (108 types pass the single sigma-constraint (sigma^3 step genuinely needed))
  [PASS] order 3^4=81: all 260 orbit types killed by {N=1 (4)} AND {3N=1 (4)}  (132 types pass the single sigma-constraint (sigma^3 step genuinely needed))
  [PASS] order 3^5=243: all 28 orbit types killed by {N=1 (4)} AND {3N=1 (4)}  (13 types pass the single sigma-constraint (sigma^3 step genuinely needed))
  [PASS] order-3 automorphism: master lemma <=> f = 3 (mod 6); f=0 excluded
  Fein-Kantor-Schacher corollary arithmetic:
  [PASS] fpf element of prime-power q^m order on 333 points forces q | 333, i.e. q in {3,37}
  [PASS] q=37: 37^2 = 1369 > 333, so every cycle has length exactly 37 ==> order exactly 37, 9 cycles
  [PASS] q=3 excluded by the theorem ==> vertex-transitive case has fpf order-37 automorphism  (conditional on FKS (CFSG))

PART 1: framework validation on real conference graphs
  --- Paley(13)/x->3x: srg(13, 6, 2, 3), automorphism of order 3 ---
  [PASS] Paley(13)/x->3x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(13)/x->3x: perm is an automorphism of odd order 3
  [PASS] Paley(13)/x->3x: all 5 orbit sizes odd and divide 3
  [PASS] Paley(13)/x->3x: A S = S Q (orbit partition equitable)
  [PASS] Paley(13)/x->3x: Q row sums = k = 6
  [PASS] Paley(13)/x->3x: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(13)/x->3x: tr(Q) = 4 is even
  [PASS] Paley(13)/x->3x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(13)/x->3x: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(13)/x->3x: disc = 13 = v non-square (quadratic irreducible)
  [PASS] Paley(13)/x->3x: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^2; tr(Q) = (v-N)/2 = 4
  [PASS] Paley(13)/x->3x: b = 2 EVEN, N = 5 = 1 (mod 4)  [MASTER LEMMA]
  [PASS] Paley(13)/x->3x: char(Q) verified directly
  --- Paley(37)/x->10x: srg(37, 18, 8, 9), automorphism of order 3 ---
  [PASS] Paley(37)/x->10x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(37)/x->10x: perm is an automorphism of odd order 3
  [PASS] Paley(37)/x->10x: all 13 orbit sizes odd and divide 3
  [PASS] Paley(37)/x->10x: A S = S Q (orbit partition equitable)
  [PASS] Paley(37)/x->10x: Q row sums = k = 18
  [PASS] Paley(37)/x->10x: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(37)/x->10x: tr(Q) = 12 is even
  [PASS] Paley(37)/x->10x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(37)/x->10x: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(37)/x->10x: disc = 37 = v non-square (quadratic irreducible)
  [PASS] Paley(37)/x->10x: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^6; tr(Q) = (v-N)/2 = 12
  [PASS] Paley(37)/x->10x: b = 6 EVEN, N = 13 = 1 (mod 4)  [MASTER LEMMA]
  [PASS] Paley(37)/x->10x: char(Q) verified directly
  --- Paley(29)/x->16x: srg(29, 14, 6, 7), automorphism of order 7 ---
  [PASS] Paley(29)/x->16x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(29)/x->16x: perm is an automorphism of odd order 7
  [PASS] Paley(29)/x->16x: all 5 orbit sizes odd and divide 7
  [PASS] Paley(29)/x->16x: A S = S Q (orbit partition equitable)
  [PASS] Paley(29)/x->16x: Q row sums = k = 14
  [PASS] Paley(29)/x->16x: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(29)/x->16x: tr(Q) = 12 is even
  [PASS] Paley(29)/x->16x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(29)/x->16x: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(29)/x->16x: disc = 29 = v non-square (quadratic irreducible)
  [PASS] Paley(29)/x->16x: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^2; tr(Q) = (v-N)/2 = 12
  [PASS] Paley(29)/x->16x: b = 2 EVEN, N = 5 = 1 (mod 4)  [MASTER LEMMA]
  [PASS] Paley(29)/x->16x: char(Q) verified directly
  --- Paley(53)/x->10x: srg(53, 26, 12, 13), automorphism of order 13 ---
  [PASS] Paley(53)/x->10x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(53)/x->10x: perm is an automorphism of odd order 13
  [PASS] Paley(53)/x->10x: all 5 orbit sizes odd and divide 13
  [PASS] Paley(53)/x->10x: A S = S Q (orbit partition equitable)
  [PASS] Paley(53)/x->10x: Q row sums = k = 26
  [PASS] Paley(53)/x->10x: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(53)/x->10x: tr(Q) = 24 is even
  [PASS] Paley(53)/x->10x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(53)/x->10x: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(53)/x->10x: disc = 53 = v non-square (quadratic irreducible)
  [PASS] Paley(53)/x->10x: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^2; tr(Q) = (v-N)/2 = 24
  [PASS] Paley(53)/x->10x: b = 2 EVEN, N = 5 = 1 (mod 4)  [MASTER LEMMA]
  [PASS] Paley(53)/x->10x: char(Q) verified directly
  --- Paley(73)/x->2x: srg(73, 36, 17, 18), automorphism of order 9 ---
  [PASS] Paley(73)/x->2x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(73)/x->2x: perm is an automorphism of odd order 9
  [PASS] Paley(73)/x->2x: all 9 orbit sizes odd and divide 9
  [PASS] Paley(73)/x->2x: A S = S Q (orbit partition equitable)
  [PASS] Paley(73)/x->2x: Q row sums = k = 36
  [PASS] Paley(73)/x->2x: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(73)/x->2x: tr(Q) = 32 is even
  [PASS] Paley(73)/x->2x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(73)/x->2x: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(73)/x->2x: disc = 73 = v non-square (quadratic irreducible)
  [PASS] Paley(73)/x->2x: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^4; tr(Q) = (v-N)/2 = 32
  [PASS] Paley(73)/x->2x: b = 4 EVEN, N = 9 = 1 (mod 4)  [MASTER LEMMA]
  [PASS] Paley(73)/x->2x: char(Q) verified directly
  --- Paley(73)/x->8x (=sigma^3): srg(73, 36, 17, 18), automorphism of order 3 ---
  [PASS] Paley(73)/x->8x (=sigma^3): srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(73)/x->8x (=sigma^3): perm is an automorphism of odd order 3
  [PASS] Paley(73)/x->8x (=sigma^3): all 25 orbit sizes odd and divide 3
  [PASS] Paley(73)/x->8x (=sigma^3): A S = S Q (orbit partition equitable)
  [PASS] Paley(73)/x->8x (=sigma^3): Q row sums = k = 36
  [PASS] Paley(73)/x->8x (=sigma^3): induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(73)/x->8x (=sigma^3): tr(Q) = 24 is even
  [PASS] Paley(73)/x->8x (=sigma^3): (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(73)/x->8x (=sigma^3): rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(73)/x->8x (=sigma^3): disc = 73 = v non-square (quadratic irreducible)
  [PASS] Paley(73)/x->8x (=sigma^3): forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^12; tr(Q) = (v-N)/2 = 24
  [PASS] Paley(73)/x->8x (=sigma^3): b = 12 EVEN, N = 25 = 1 (mod 4)  [MASTER LEMMA]
  [PASS] Paley(73)/x->8x (=sigma^3): char(Q) verified directly
  [PASS] Paley(73): orbit-splitting N(sigma^3) = 3 N(sigma) - 2 n_1  (25 = 27-2)
  [PASS] Paley(73): every sigma-orbit of size 9 splits into EXACTLY 3 sigma^3-orbits
  --- Paley(125)/Frobenius: srg(125, 62, 30, 31), automorphism of order 3 ---
  [PASS] Paley(125)/Frobenius: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(125)/Frobenius: perm is an automorphism of odd order 3
  [PASS] Paley(125)/Frobenius: all 45 orbit sizes odd and divide 3
  [PASS] Paley(125)/Frobenius: A S = S Q (orbit partition equitable)
  [PASS] Paley(125)/Frobenius: Q row sums = k = 62
  [PASS] Paley(125)/Frobenius: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(125)/Frobenius: tr(Q) = 40 is even
  [PASS] Paley(125)/Frobenius: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(125)/Frobenius: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(125)/Frobenius: disc = 125 = v non-square (quadratic irreducible)
  [PASS] Paley(125)/Frobenius: forced char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^22; tr(Q) = (v-N)/2 = 40
  [PASS] Paley(125)/Frobenius: b = 22 EVEN, N = 45 = 1 (mod 4)  [MASTER LEMMA]
  --- Paley(9)/x->x+1 [CONTROL]: srg(9, 4, 1, 2), automorphism of order 3 ---
  [PASS] Paley(9)/x->x+1 [CONTROL]: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(9)/x->x+1 [CONTROL]: perm is an automorphism of odd order 3
  [PASS] Paley(9)/x->x+1 [CONTROL]: all 3 orbit sizes odd and divide 3
  [PASS] Paley(9)/x->x+1 [CONTROL]: A S = S Q (orbit partition equitable)
  [PASS] Paley(9)/x->x+1 [CONTROL]: Q row sums = k = 4
  [PASS] Paley(9)/x->x+1 [CONTROL]: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(9)/x->x+1 [CONTROL]: tr(Q) = 6 is even
  [PASS] Paley(9)/x->x+1 [CONTROL]: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(9)/x->x+1 [CONTROL]: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(9)/x->x+1 [CONTROL]: CONTROL -- fpf (no fixed points)
  [PASS] Paley(9)/x->x+1 [CONTROL]: CONTROL -- N = 3 = 3 (mod 4): master CONCLUSION fails
  [PASS] Paley(9)/x->x+1 [CONTROL]: CONTROL -- disc = 9 = v IS a perfect square (irreducibility hypothesis fails; nothing else does)
  [PASS] Paley(9)/x->x+1 [CONTROL]: CONTROL -- forced trace (v-N)/2 = 3 is ODD, actual tr(Q) = 6 is even: equal-multiplicity forcing absent
  [PASS] Paley(9)/x->x+1 [CONTROL]: CONTROL -- char(Q) multiplicities of r=1, s=-2 are (2,0) UNEQUAL (rationality escape)
  --- Paley(49)/x->x+1 [CONTROL]: srg(49, 24, 11, 12), automorphism of order 7 ---
  [PASS] Paley(49)/x->x+1 [CONTROL]: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(49)/x->x+1 [CONTROL]: perm is an automorphism of odd order 7
  [PASS] Paley(49)/x->x+1 [CONTROL]: all 7 orbit sizes odd and divide 7
  [PASS] Paley(49)/x->x+1 [CONTROL]: A S = S Q (orbit partition equitable)
  [PASS] Paley(49)/x->x+1 [CONTROL]: Q row sums = k = 24
  [PASS] Paley(49)/x->x+1 [CONTROL]: induced orbit subgraphs regular of EVEN degree (handshake, odd orbits)
  [PASS] Paley(49)/x->x+1 [CONTROL]: tr(Q) = 42 is even
  [PASS] Paley(49)/x->x+1 [CONTROL]: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0 (annihilator, S injective)
  [PASS] Paley(49)/x->x+1 [CONTROL]: rank(Q-kI) = N-1 (k simple in Q)
  [PASS] Paley(49)/x->x+1 [CONTROL]: CONTROL -- fpf (no fixed points)
  [PASS] Paley(49)/x->x+1 [CONTROL]: CONTROL -- N = 7 = 3 (mod 4): master CONCLUSION fails
  [PASS] Paley(49)/x->x+1 [CONTROL]: CONTROL -- disc = 49 = v IS a perfect square (irreducibility hypothesis fails; nothing else does)
  [PASS] Paley(49)/x->x+1 [CONTROL]: CONTROL -- forced trace (v-N)/2 = 21 is ODD, actual tr(Q) = 42 is even: equal-multiplicity forcing absent
  [PASS] Paley(49)/x->x+1 [CONTROL]: CONTROL -- char(Q) multiplicities of r=3, s=-4 are (6,0) UNEQUAL (rationality escape)

PART 2: route (a) (zeta_3 eigenspace traces) is vacuous -- symbolic
  [PASS] tr(P) = f reduces to 0 = 0 identically (no constraint on f)  (residual = 0)

RESULT: ALL CHECKS PASS -- master parity lemma N=1 (mod 4) verified on all real examples; fpf 3-power orders 3,9,27,81,243 all killed on (333,166,82,83); controls Paley(9)/Paley(49) break exactly at irreducibility (v square). Theorem + corollary arithmetic verified.
```
