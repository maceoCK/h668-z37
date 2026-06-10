# Verification report: every order-37 automorphism of srg(333,166,82,83) is fixed-point-free

**Verdict: PROVEN** (complete, gap-free proof below; machine-checked numerics in `fpf_check.py`).

The earlier agent's sketch needed repair: its part (a) (cyclotomic eigenspace traces) yields
*nothing* beyond the trivial congruence f ≡ 0 (mod 37), and its part (b) as stated (interlacing +
fixed-degree congruence) does **not** close even the smallest case (see §6.2 — for f = 37 the fixed
subgraph is forced to be a putative srg(37,18,8,9), e.g. Paley(37), which passes interlacing
comfortably). The correct closing argument, given in full here, is a **trace pincer on the orbit
quotient matrix**: the two exact trace identities forced by the quotient's characteristic polynomial
collide with Cauchy–Schwarz for every fixed-point count f = 37c, 1 ≤ c ≤ 8. The collision is sharp:
the same inequality is satisfiable precisely in the fixed-point-free case (margin 666), as it must
be, since a valid fixed-point-free orbit matrix for these parameters is already known to exist
(`conf_orbit3.out`, FINDINGS §6).

---

## Theorem

Let Γ be a strongly regular graph with parameters (n, k, λ, μ) = (333, 166, 82, 83) and let
σ ∈ Aut(Γ) with σ^37 = id, σ ≠ id. Then σ has no fixed vertex.

*(The statement is conditional on the existence of Γ, which is open — it is the conference graph
equivalent to the open Hadamard order H(668). The proof below is valid for any such Γ.)*

## 0. Standing facts and the standard theorems invoked

Throughout, A is the adjacency matrix of Γ, J the all-ones matrix, I the identity.

**(F1) SRG identity.** A² = kI + λA + μ(J − I − A) = 166I + 82A + 83(J − I − A), equivalently
A² + A − 83I = 83J.

**(F2) Spectrum.** Eigenvalues of A: k = 166 (Perron) and the roots θ, τ of
x² − (λ−μ)x − (k−μ) = **x² + x − 83**, i.e. θ = (−1+√333)/2 = (−1+3√37)/2 ≈ 8.6241 and
τ = (−1−3√37)/2 ≈ −9.6241. Since 2k + (n−1)(λ−μ) = 332 − 332 = 0 (conference condition), θ and τ
each have multiplicity (n−1)/2 = 166. Useful exact values: θ+τ = −1, θτ = −83, θ²+τ² = 167.
The polynomial x² + x − 83 is irreducible over **Q** (discriminant 333 = 9·37 is not a square).

**(F3) Connectivity / simplicity of k (Perron–Frobenius).** μ = 83 > 0, so Γ is connected of
diameter 2; for a connected k-regular graph the eigenvalue k is simple, with eigenvector **1**.
Hence char(A) = (x − 166)(x² + x − 83)^166 exactly.

**(F4) Orbit partitions are equitable** (Godsil–Royle, *Algebraic Graph Theory*, §9.3): the orbit
partition of any group of automorphisms of a graph is an equitable partition. We re-verify this
from scratch for our situation in Lemma 2/Lemma 3 (it follows in two lines from transitivity of
⟨σ⟩ on each orbit), so the proof is self-contained.

**(F5) Quotient divisibility** (Godsil–Royle, Thm 9.3.3): if π is an equitable partition with
characteristic (indicator) matrix S (full column rank) and quotient matrix Q, then AS = SQ, and
char(Q) divides char(A) (with multiplicity). *Proof recalled:* extend S to an invertible
[S | R]; then A[S | R] = [S | R]·M where M is block upper triangular with top-left block Q
(because AS = SQ lands in the column space of S); so char(A) = char(M) = char(Q)·char(M₂₂). ∎

**(F6) Cauchy–Schwarz.** For reals x₁,…,x_m: Σxᵢ² ≥ (Σxᵢ)²/m.

*Not needed:* Cauchy interlacing, Benson-type counting (tr(AP^i)), Perron bounds on the fixed
subgraph, and integrality of cyclotomic character sums all turn out to be either subsumed or too
weak (§6). The proof below uses only (F1)–(F6) and counting.

## 1. Lemma 1 (fixed-point count)

Let F = Fix(σ), f = |F|. Since σ has prime order 37, every ⟨σ⟩-orbit on V(Γ) has size 1 or 37.
Hence f ≡ n = 333 = 9·37 ≡ 0 (mod 37). Write

- **f = 37c** with c ∈ {0, 1, …, 9}, and
- **s = 9 − c** = number of 37-orbits (so n − f = 37s).

c = 9 (f = 333) means σ = id, excluded. The theorem asserts c = 0. Assume from now on, for
contradiction, **1 ≤ c ≤ 8**, i.e. **1 ≤ s ≤ 8** (so both F ≠ ∅ and at least one 37-orbit exists).

Label the 37-orbits O₁, …, O_s.

## 2. Lemma 2 (all-or-nothing adjacency, degrees, and edge counts)

**(a) All-or-nothing.** Let u ∈ F and let O be a 37-orbit. N(u) is σ-invariant (σu = u), so
N(u) ∩ O is a σ-invariant subset of O. Since ⟨σ⟩ acts transitively on O (orbit of prime size 37),
the only invariant subsets of O are ∅ and O. Hence **u is adjacent to all of O or to none of it**.

**(b) Fixed degrees.** For u ∈ F let aᵤ = #{i : u ~ Oᵢ} and dᵤ = its degree in the induced fixed
subgraph Δ = Γ[F]. By (a), 166 = dᵤ + 37aᵤ, i.e. **dᵤ = 166 − 37aᵤ ≡ 18 (mod 37)**, and
0 ≤ aᵤ ≤ min(s, 4) (since 37·5 > 166).

**(c) Constancy along orbits.**
  - For u ∈ F and a 37-orbit Oᵢ, every vertex of Oᵢ has the same adjacency (0 or 1) to u — this is
    exactly (a). Write N ∈ {0,1}^{F×s} for this incidence (N_{u,i} = 1 iff u ~ Oᵢ), and
    γᵢ = Σ_{u∈F} N_{u,i} (number of fixed vertices adjacent to Oᵢ).
  - For two orbits Oᵢ, Oⱼ (i = j allowed), |N(w) ∩ Oⱼ| is the same for every w ∈ Oᵢ: if w′ = σ^t w
    then N(w′) ∩ Oⱼ = σ^t(N(w) ∩ Oⱼ). Call this constant **B_{ij}** ∈ {0,…,37} (B_{ii} ≤ 36).
    Counting edges between Oᵢ and Oⱼ both ways gives 37B_{ij} = 37B_{ji}, so **B is symmetric**
    (used in (3.4) to identify tr(B²) with Σ B_{ij}²).

**(d) Edge counts.** Let m = |E(Δ)| and **G = Σ_{u∈F} aᵤ = Σᵢ γᵢ** (= total number of incidences
between fixed vertices and orbits). Summing (b) over F:

  **2m = 166f − 37G.**  (2.1)

**(e) Row sums on the orbit side.** A vertex w ∈ Oᵢ has degree 166 = γᵢ + Σⱼ B_{ij}, so

  **Σⱼ B_{ij} = 166 − γᵢ** for each i.  (2.2)

## 3. The orbit quotient and its characteristic polynomial

Order the cells of the orbit partition as: the f singletons {u} (u ∈ F), then O₁,…,O_s. By
Lemma 2(c) this partition is **equitable**, with quotient matrix (size (f+s) × (f+s))

    Q = [ A_Δ      37·N ]
        [ Nᵀ        B   ]

where A_Δ is the adjacency matrix of Δ (entry Q_{u,i} = |N(u) ∩ Oᵢ| = 37N_{u,i}, entry
Q_{i,u} = |N(w) ∩ {u}| = N_{u,i} for w ∈ Oᵢ, entry Q_{i,j} = B_{ij}). Every row of Q sums to 166.

**Lemma 3.** char(Q) = (x − 166)·(x² + x − 83)^b with **b = (f + s − 1)/2 = 166 − 18s**.

*Proof.* By (F5), char(Q) divides char(A) = (x − 166)(x² + x − 83)^166 (by (F3)). Since
x² + x − 83 is irreducible over Q (F2) and char(Q) ∈ Z[x] is monic of degree f + s, we get
char(Q) = (x − 166)^a (x² + x − 83)^b with a ∈ {0, 1} and a + 2b = f + s. Q**1** = 166·**1**
(constant row sums), so 166 is an eigenvalue of Q and a = 1 (166 is not a root of x²+x−83).
Finally f + s = 37c + (9 − c) = 36c + 9 is odd, consistent, and
b = (f + s − 1)/2 = (332 − 36s)/2 = 166 − 18s. ∎

**Corollary (trace identities).** Since the eigenvalues of Q are exactly 166 (once) and θ, τ
(b times each), and θ + τ = −1, θ² + τ² = 167:

  **tr(Q) = 166 − b = 18s,**  (3.1)
  **tr(Q²) = 166² + 167b = 27556 + 167(166 − 18s) = 55278 − 3006s.**  (3.2)

(tr(Q²) = Σ eigenvalues² is legitimate even though Q is not symmetric: all roots of char(Q) are
real and tr(Q²) is the second power-sum of the roots of char(Q).)

Computing the same traces from the block form:

- tr(Q) = tr(A_Δ) + tr(B) = 0 + tr(B), so with (3.1):

  **tr(B) = Σᵢ B_{ii} = 18s =: D.**  (3.3)

- tr(Q²) = tr(A_Δ²) + tr((37N)·Nᵀ) + tr(Nᵀ·(37N)) + tr(B²)
  = 2m + 74·Σ_{u,i} N_{u,i}² + Σ_{i,j} B_{ij}B_{ji} = 2m + 74G + Σ_{i,j} B_{ij}²
  (N is 0/1 so ΣN² = ΣN = G; B is symmetric by 2(c)). With (3.2) and (2.1):

  **Σ_{i,j} B_{ij}² = 55278 − 3006s − (166f − 37G) − 74G = 3136s − 37G.**  (3.4)

  (Using 166f = 166(333 − 37s) = 55278 − 6142s.)

- Total sum of B: by (2.2), Σ_{i,j} B_{ij} = 166s − G. Splitting off the diagonal (3.3), the
  off-diagonal sum is **T := Σ_{i≠j} B_{ij} = 166s − G − 18s = 148s − G**, i.e. G = 148s − T.
  Substituting into (3.4):

  **Σ_{i,j} B_{ij}² = 3136s − 37(148s − T) = 37T − 2340s.**  (3.5)

Every step above is an exact identity in the integers (m, G, T, B_{ij} are integers). Equations
(3.3), (3.5) and Cauchy–Schwarz are all we need.

## 4. The pincer

**Case s = 1 (c = 8, f = 296).** B is the 1×1 matrix (B₁₁); T = 0 (empty off-diagonal sum). Then
(3.5) reads B₁₁² = 37·0 − 2340·1 = **−2340 < 0**, absurd. (Concretely: (3.3) forces B₁₁ = 18,
(2.2) then forces γ₁ = 166 − 18 = 148 = G, and (3.4) demands ΣB² = 3136 − 37·148 = −2340, while in
fact ΣB² = 18² = 324.) Contradiction.

**Case 2 ≤ s ≤ 8 (1 ≤ c ≤ 7).** By Cauchy–Schwarz (F6) applied separately to the s diagonal
entries and the s(s−1) off-diagonal entries of B:

  Σ_{i,j} B_{ij}² = Σᵢ B_{ii}² + Σ_{i≠j} B_{ij}² ≥ D²/s + T²/(s(s−1)) = (18s)²/s + T²/(s(s−1))
  = 324s + T²/(s(s−1)).

Combining with (3.5):

  37T − 2340s ≥ 324s + T²/(s(s−1)),

i.e., multiplying by s(s−1) > 0,

  **q(T) := T² − 37s(s−1)·T + 2664·s²(s−1) ≤ 0.**  (4.1)

(2664 = 324 + 2340.) The discriminant of q is

  disc = 37²s²(s−1)² − 4·2664·s²(s−1) = **s²(s−1)·[1369(s−1) − 10656]**.

For s ≤ 8: 1369(s−1) ≤ 1369·7 = 9583 < 10656, so disc < 0 and q(T) > 0 for *every real* T —
no value of T (a fortiori no integer value, and no choice of G, γᵢ, B_{ij}, regardless of any
further constraints) can satisfy (4.1). Contradiction.

Hence no s ∈ {1,…,8} is possible: a non-identity σ of order 37 has f ∈ {37,…,296} ruled out, and
f = 333 means σ = id. Therefore **f = 0**. ∎

### Numerical table (all verified exactly in `fpf_check.py`)

Write Φ(s) = 1369(s−1) − 10656, so disc q = s²(s−1)·Φ(s). The minimum, over all real T, of
[Cauchy–Schwarz lower bound − required ΣB²] equals −s·Φ(s)/4; infeasibility ⟺ this is > 0 ⟺
Φ(s) < 0. The script's exhaustive integer-G scan reproduces exactly these margins.

| s | c | f   | tr(B)=18s | ΣB² required (3.5) | Φ(s) = 1369(s−1)−10656 | min infeasibility margin −sΦ/4 | status |
|---|---|-----|-----------|--------------------|------------------------|-------------------------------|--------|
| 1 | 8 | 296 | 18        | −2340 (< 0!)       | (separate: T = 0 forced) | 324−(−2340) = 2664          | impossible |
| 2 | 7 | 259 | 36        | 37T − 4680         | −9287                  | 9287/2 (at G = 259)           | impossible |
| 3 | 6 | 222 | 54        | 37T − 7020         | −7918                  | 11877/2 (at G = 333)          | impossible |
| 4 | 5 | 185 | 72        | 37T − 9360         | −6549                  | 6549 (at G = 370)             | impossible |
| 5 | 4 | 148 | 90        | 37T − 11700        | −5180                  | 6475 (at G = 370)             | impossible |
| 6 | 3 | 111 | 108       | 37T − 14040        | −3811                  | 11433/2 (at G = 333)          | impossible |
| 7 | 2 | 74  | 126       | 37T − 16380        | −2442                  | 8547/2 (at G = 259)           | impossible |
| 8 | 1 | 37  | 144       | 37T − 18720        | −1073                  | 2146 (at G = 148)             | impossible |
| 9 | 0 | 0   | 162       | 28224 (forced)     | **+296**               | −666 (satisfied, slack 666)   | **feasible** |

(For s ≤ 8, Φ(s) ≤ 1369·7 − 10656 = −1073 < 0; the exhaustive integer margins differ from −sΦ/4
only by the integrality rounding of G, as visible above.)

## 5. Sharpness — the pincer fails exactly where it must

For s = 9 (f = 0) there is no F-part: G = 0 and T = 148·9 = 1332 are *forced*, and (3.5) gives
ΣB² = 37·1332 − 21060 = 28224 — automatically consistent with tr(Q²). The Cauchy–Schwarz lower
bound is 324·9 + 1332²/72 = 2916 + 24642 = 27558 ≤ 28224, satisfied with **margin 666**. Indeed
1369·8 − 10656 = +296 > 0: the pincer flips sign exactly between s = 8 and s = 9. This is the
correct calibration of the earlier round's "sharp at c = 9" remark: the inequality chain excludes
*precisely* the cases with fixed points and nothing else.

Cross-check against reality: the project's verified fixed-point-free orbit matrix R
(`conf_orbit3.out`, FINDINGS §6 — diagonal (10,14,14,14,18,20,24,24,24), row sums 166,
R² + R − 83I = 3071J) satisfies tr(R) = 162 = 18·9, Σ R_{ij}² = 28224, and
char(R) = (x−166)(x²+x−83)⁴ — exactly the s = 9 instance of Lemma 3 and (3.3)–(3.5). Verified in
`fpf_check.py` Part 2. If our pincer had also excluded s = 9 it would have been provably wrong;
it doesn't, and it reproduces R's invariants on the nose.

## 6. Autopsy of the earlier sketch

**6.1 Part (a) — cyclotomic eigenspace traces give nothing new.** Write the permutation matrix P
of σ; f = tr(P) = 1 + t_θ + t_τ with t_θ = tr(P|V_θ), t_τ = tr(P|V_τ). The eigenvalues of P|V_θ
are 37th roots of unity; for g ∈ Gal(Q(ζ₃₇)/Q) with g(ζ)=ζ^α, g maps the ζ^j-eigenspace of P
inside V_θ to the ζ^{αj}-eigenspace inside V_{g(θ)}, and g fixes θ ∈ Q(√37) iff α is a quadratic
residue mod 37. Hence the multiplicity profile of P on V_θ is constant on the QR and QNR classes:
m₀ + 18(α+β) = 166 and t_θ = m₀ + α(−1+√37)/2 + β(−1−√37)/2 (quadratic Gauss sum,
37 ≡ 1 mod 4). Then f = 1 + t_θ + t_τ = 1 + 2m₀ − (α+β) = 333 − 37(α+β): exactly the congruence
f ≡ 0 (mod 37) of Lemma 1, already known from orbit counting, with no further restriction.
Adding Benson-type counts tr(AP^i) = 166 + θt_θ(P^i) + τt_τ(P^i) only yields the (consistent)
divisibility 2e = (9−c) + 3(α−β) for the integer e = (1/37)·#{v : v ~ σ(v)} — again no exclusion.
Part (a) is therefore dispensable; it is *not* used in the proof above.

**6.2 Part (b) as sketched is insufficient.** For c = 1 (f = 37): Lemma 2(b) forces dᵤ = 18 for
all u ∈ F (55 > 36 = f−1), so Δ is 18-regular on 37 vertices; the σ-invariant common-neighbor
counts force every adjacent pair of fixed vertices to have ≡ 8 (mod 37) and every non-adjacent
pair ≡ 9 (mod 37) common fixed neighbors, and the counting identity
Σ_w d_w(d_w−1) = Σ_{u≠v} |N_Δ(u)∩N_Δ(v)| = 18·17·37 forces *equality*: Δ would be exactly a
(37,18,8,9) strongly regular graph — e.g. Paley(37), which exists. Its spectrum
{18, (−1±√37)/2} sits strictly inside [τ, θ] ∪ {k}, so Cauchy interlacing with Γ raises no
objection whatsoever. The sketched "interlacing + degree congruence" pincer therefore cannot close
even c = 1; the orbit-quotient trace identities of §3–4 are what actually kill it (and all
other c simultaneously). The candidate fixed subgraph Paley(37) fails (3.4)–(3.5), i.e. its
existence as a *fixed* subgraph is incompatible with how the 8 orbits must attach to it — a global
constraint invisible to interlacing.

**6.3 What the final proof rests on** (checklist of invoked theorems and their hypotheses):
- *Prime-order orbit counting* (Lemma 1): needs only |σ| = 37 prime. ✓
- *Equitability of the orbit partition* (F4/Lemma 2c): re-proved inline from transitivity of ⟨σ⟩
  on each orbit; no further hypotheses. ✓
- *Quotient divisibility* (F5): needs S of full column rank — orbit cells are non-empty. ✓
- *char(A) = (x−166)(x²+x−83)^166* (F2, F3): needs Γ connected (μ = 83 > 0) k-regular and the
  conference multiplicity computation 2k + (n−1)(λ−μ) = 0. ✓
- *Irreducibility of x²+x−83 over Q*: disc = 333 not a square. ✓ (This is what forces θ- and
  τ-multiplicities in char(Q) to be equal — the engine of the whole proof.)
- *Cauchy–Schwarz* (F6): unconditional. ✓
- *Integrality*: m, G, T, B_{ij}, γᵢ, aᵤ are integers by definition; in fact the final
  contradiction for 2 ≤ s ≤ 8 holds for all **real** T, so no integrality is even needed there
  (it is needed nowhere: s = 1 uses only T = 0). The proof is robust to any slack in the
  combinatorial side constraints (0 ≤ B_{ij} ≤ 37, B_{ii} even, γᵢ ≤ f, etc. are never used).

## 7. Machine verification (`fpf_check.py`)

Exact-arithmetic (sympy/Fraction) + numpy checks, running in seconds:

- **Part 0**: all spectral facts of (F2) for (333,166,82,83), including irreducibility,
  conference condition, θ+τ = −1, θ²+τ² = 167, and the congruences 166 ≡ 18, 82 ≡ 8, 83 ≡ 9
  (mod 37).
- **Part 1**: every identity of §3 verified *symbolically in s* (b = 166−18s; tr(B) = 18s;
  tr(Q²) = 55278−3006s; ΣB² = 3136s−37G = 37T−2340s; the pincer quadratic and its discriminant
  s²(s−1)(1369(s−1)−10656)), plus the sign table for s = 1…9.
- **Part 2**: for every s = 1…8, *exhaustive* infeasibility over all integers G ∈ [0, 4f]
  (min over G of [CS lower bound − required ΣB²] is strictly positive); for s = 9 feasibility
  with margin 666; and the actual orbit matrix R from `conf_orbit3.out` reproduces every
  identity (row sums, R²+R−83I = 3071J, tr = 162, ΣR² = 28224,
  char(R) = (x−166)(x²+x−83)⁴ via annihilator + rank).
- **Part 3**: the entire framework — all-or-nothing lemma, equitability (AS = SQ), the char(Q)
  factorization (annihilating polynomial + rank(Q−kI) = f+s−1), both trace identities,
  2m = kf − pG, B row sums, and the direction of the Cauchy–Schwarz pincer — validated
  end-to-end on four *real* conference graphs with *real* automorphisms of odd prime order:
  - Paley(13), σ: x↦9x (p = 3, f = 1, s = 4);
  - Paley(29), σ: x↦16x (p = 7, f = 1, s = 4);
  - Paley(37), σ: x↦x+1 (p = 37, f = 0, s = 1 — degenerate FPF case, char(Q) = x−18);
  - Paley(125) = srg(125,62,30,31), σ = Frobenius x↦x⁵ (p = 3, f = 5, s = 40; fixed subgraph
    Δ = Paley(5) = C₅, a genuinely non-trivial fixed subgraph). On each, the identities hold
    exactly and the pincer inequality is *satisfied* (no false kill), confirming all sign
    conventions.

Result: **all checks pass** (see script output reproduced at the end of this file).

## 8. Consequences for the project

Combined with the standing Z₃₇ machinery this upgrades the working hypothesis of FINDINGS §6 to a
theorem: *any* order-37 automorphism of a putative srg(333,166,82,83) is fixed-point-free, so the
9×9-circulant (9 orbits of 37) orbit-matrix/lift framework is not merely one symmetry ansatz among
several — it is the **only** possible shape of a Z₃₇ symmetry on this parameter set. In particular
the 2901-representative skeleton enumeration covers every order-37 symmetric candidate, with no
fixed-point sectors left unexamined.

---

### Appendix: script output (full, exit code 0, runtime 0.4 s)

```
PART 0: spectral facts for srg(333,166,82,83)
  [PASS] 333 = 9*37
  [PASS] quad = x^2+x-83
  [PASS] theta root of quad
  [PASS] tau root of quad
  [PASS] theta+tau = -1, theta*tau = -83
  [PASS] theta^2+tau^2 = 167
  [PASS] disc 333 = 9*37 not a perfect square (quad irreducible /Q)
  [PASS] conference condition 2k+(n-1)(lam-mu)=0
  [PASS] multiplicities (n-1)/2 = 166
  [PASS] 166 mod 37 = 18 (fixed-degree congruence)
  [PASS] 82 mod 37 = 8, 83 mod 37 = 9

PART 1: symbolic identities (in the number of 37-orbits s)
  [PASS] b = (f+s-1)/2 = 166-18s
  [PASS] tr(Q) = 18s  ==> tr(B) = 18s (since tr(A_Delta)=0)
  [PASS] tr(Q^2) = 55278 - 3006 s
  [PASS] Sum B_ij^2 = tr(Q^2) - 2m - 74G = 3136 s - 37 G
  [PASS] Sum B_ij^2 = 37 T - 2340 s
  [PASS] pincer quadratic T^2 - 37s(s-1)T + 2664 s^2(s-1)
  [PASS] discriminant = s^2 (s-1) (1369(s-1) - 10656)
  discriminant sign across the full range:
    s=1 (f=296): disc =          0  <-- s=1 handled separately (T=0 forced)
    s=2 (f=259): disc =     -37148  <-- infeasible (proof case)
  [PASS] disc<0 at s=2
    s=3 (f=222): disc =    -142524  <-- infeasible (proof case)
  [PASS] disc<0 at s=3
    s=4 (f=185): disc =    -314352  <-- infeasible (proof case)
  [PASS] disc<0 at s=4
    s=5 (f=148): disc =    -518000  <-- infeasible (proof case)
  [PASS] disc<0 at s=5
    s=6 (f=111): disc =    -685980  <-- infeasible (proof case)
  [PASS] disc<0 at s=6
    s=7 (f= 74): disc =    -717948  <-- infeasible (proof case)
  [PASS] disc<0 at s=7
    s=8 (f= 37): disc =    -480704  <-- infeasible (proof case)
  [PASS] disc<0 at s=8
    s=9 (f=  0): disc =     191808  <-- feasible (FPF case, as required)
  [PASS] disc>0 at s=9 (pincer correctly does NOT exclude FPF)

PART 2: exhaustive pincer check over all integer G, s=1..9
  [PASS] s=1 (f=296): required Sum B^2 = -2340 < 0 <= 324 = B11^2  (identity forces SumB^2=-2340, but SumB^2=324)
  [PASS] s=2 (f=259): infeasible for every integer G in [0,1036]  (min (CS_lower - required) = 9287/2 at G=259)
  [PASS] s=3 (f=222): infeasible for every integer G in [0,888]  (min (CS_lower - required) = 11877/2 at G=333)
  [PASS] s=4 (f=185): infeasible for every integer G in [0,740]  (min (CS_lower - required) = 6549 at G=370)
  [PASS] s=5 (f=148): infeasible for every integer G in [0,592]  (min (CS_lower - required) = 6475 at G=370)
  [PASS] s=6 (f=111): infeasible for every integer G in [0,444]  (min (CS_lower - required) = 11433/2 at G=333)
  [PASS] s=7 (f=74): infeasible for every integer G in [0,296]  (min (CS_lower - required) = 8547/2 at G=259)
  [PASS] s=8 (f=37): infeasible for every integer G in [0,148]  (min (CS_lower - required) = 2146 at G=148)
  [PASS] s=9 (f=0): required Sum B^2 = 28224 >= CS lower bound 27558 (margin 666)
  [PASS] R symmetric, row sums 166
  [PASS] R^2 + R - 83 I = 3071 J (= 83*37 J)
  [PASS] tr(R) = 162 = 18 s
  [PASS] Sum R_ij^2 = 28224 (matches tr(Q^2) identity at s=9, G=0)
  [PASS] (R-166I)(R^2+R-83I) = 0 and rank(R-166I)=8  ==> char(R)=(x-166)(x^2+x-83)^4

PART 3: framework validation on real conference graphs
  --- Paley(13)/x->9x: srg(13, 6, 2, 3), automorphism of order 3 ---
  [PASS] Paley(13)/x->9x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(13)/x->9x: perm is an automorphism
  [PASS] Paley(13)/x->9x: perm has order 3
  [PASS] Paley(13)/x->9x: all non-fixed orbits have size 3; f=1 = n mod-3-congruent
  [PASS] Paley(13)/x->9x: all-or-nothing (fixed vertex adjacent to all or none of each orbit)
  [PASS] Paley(13)/x->9x: orbit partition is equitable (constant block counts)
  [PASS] Paley(13)/x->9x: A S = S Q (equitability, matrix form)
  [PASS] Paley(13)/x->9x: Q row sums = k
  [PASS] Paley(13)/x->9x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0
  [PASS] Paley(13)/x->9x: rank(Q-kI) = f+s-1 (k simple in Q)
  [PASS] Paley(13)/x->9x: f+s odd
  [PASS] Paley(13)/x->9x: eigenvalue quadratic irreducible (disc 13 non-square)
  [PASS] Paley(13)/x->9x: char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^2 (direct)
  [PASS] Paley(13)/x->9x: tr(B) = k + b(lam-mu) = 4
  [PASS] Paley(13)/x->9x: tr(Q^2) = k^2 + b[(lam-mu)^2+2(k-mu)] = 50
  [PASS] Paley(13)/x->9x: tr(Q^2) = 2m + 2pG + Sum B^2  (50 = 0+12+38)
  [PASS] Paley(13)/x->9x: 2m = kf - pG  (0 = 6 - 6)
  [PASS] Paley(13)/x->9x: B row sums k - gamma_i
  [PASS] Paley(13)/x->9x: pincer C-S inequality satisfied (SumB^2=38 >= 31)
      f=1, s=4, b=2, G=2, 2m=0, trB=4, T=18, SumB^2=38, CS lower=31
  --- Paley(29)/x->16x: srg(29, 14, 6, 7), automorphism of order 7 ---
  [PASS] Paley(29)/x->16x: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(29)/x->16x: perm is an automorphism
  [PASS] Paley(29)/x->16x: perm has order 7
  [PASS] Paley(29)/x->16x: all non-fixed orbits have size 7; f=1 = n mod-7-congruent
  [PASS] Paley(29)/x->16x: all-or-nothing (fixed vertex adjacent to all or none of each orbit)
  [PASS] Paley(29)/x->16x: orbit partition is equitable (constant block counts)
  [PASS] Paley(29)/x->16x: A S = S Q (equitability, matrix form)
  [PASS] Paley(29)/x->16x: Q row sums = k
  [PASS] Paley(29)/x->16x: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0
  [PASS] Paley(29)/x->16x: rank(Q-kI) = f+s-1 (k simple in Q)
  [PASS] Paley(29)/x->16x: f+s odd
  [PASS] Paley(29)/x->16x: eigenvalue quadratic irreducible (disc 29 non-square)
  [PASS] Paley(29)/x->16x: char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^2 (direct)
  [PASS] Paley(29)/x->16x: tr(B) = k + b(lam-mu) = 12
  [PASS] Paley(29)/x->16x: tr(Q^2) = k^2 + b[(lam-mu)^2+2(k-mu)] = 226
  [PASS] Paley(29)/x->16x: tr(Q^2) = 2m + 2pG + Sum B^2  (226 = 0+28+198)
  [PASS] Paley(29)/x->16x: 2m = kf - pG  (0 = 14 - 14)
  [PASS] Paley(29)/x->16x: B row sums k - gamma_i
  [PASS] Paley(29)/x->16x: pincer C-S inequality satisfied (SumB^2=198 >= 183)
      f=1, s=4, b=2, G=2, 2m=0, trB=12, T=42, SumB^2=198, CS lower=183
  --- Paley(37)/x->x+1: srg(37, 18, 8, 9), automorphism of order 37 ---
  [PASS] Paley(37)/x->x+1: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(37)/x->x+1: perm is an automorphism
  [PASS] Paley(37)/x->x+1: perm has order 37
  [PASS] Paley(37)/x->x+1: all non-fixed orbits have size 37; f=0 = n mod-37-congruent
  [PASS] Paley(37)/x->x+1: all-or-nothing (fixed vertex adjacent to all or none of each orbit)
  [PASS] Paley(37)/x->x+1: orbit partition is equitable (constant block counts)
  [PASS] Paley(37)/x->x+1: A S = S Q (equitability, matrix form)
  [PASS] Paley(37)/x->x+1: Q row sums = k
  [PASS] Paley(37)/x->x+1: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0
  [PASS] Paley(37)/x->x+1: rank(Q-kI) = f+s-1 (k simple in Q)
  [PASS] Paley(37)/x->x+1: f+s odd
  [PASS] Paley(37)/x->x+1: eigenvalue quadratic irreducible (disc 37 non-square)
  [PASS] Paley(37)/x->x+1: char(Q) = (x-k)(x^2-(lam-mu)x-(k-mu))^0 (direct)
  [PASS] Paley(37)/x->x+1: tr(B) = k + b(lam-mu) = 18
  [PASS] Paley(37)/x->x+1: tr(Q^2) = k^2 + b[(lam-mu)^2+2(k-mu)] = 324
  [PASS] Paley(37)/x->x+1: tr(Q^2) = 2m + 2pG + Sum B^2  (324 = 0+0+324)
  [PASS] Paley(37)/x->x+1: 2m = kf - pG  (0 = 0 - 0)
  [PASS] Paley(37)/x->x+1: B row sums k - gamma_i
  [PASS] Paley(37)/x->x+1: pincer C-S inequality satisfied (SumB^2=324 >= 324)
      f=0, s=1, b=0, G=0, 2m=0, trB=18, T=0, SumB^2=324, CS lower=324
  --- Paley(125)/Frobenius: srg(125, 62, 30, 31), automorphism of order 3 ---
  [PASS] Paley(125)/Frobenius: srg identity A^2 = kI + lam A + mu(J-I-A)
  [PASS] Paley(125)/Frobenius: perm is an automorphism
  [PASS] Paley(125)/Frobenius: perm has order 3
  [PASS] Paley(125)/Frobenius: all non-fixed orbits have size 3; f=5 = n mod-3-congruent
  [PASS] Paley(125)/Frobenius: all-or-nothing (fixed vertex adjacent to all or none of each orbit)
  [PASS] Paley(125)/Frobenius: orbit partition is equitable (constant block counts)
  [PASS] Paley(125)/Frobenius: A S = S Q (equitability, matrix form)
  [PASS] Paley(125)/Frobenius: Q row sums = k
  [PASS] Paley(125)/Frobenius: (Q-kI)(Q^2-(lam-mu)Q-(k-mu)I) = 0
  [PASS] Paley(125)/Frobenius: rank(Q-kI) = f+s-1 (k simple in Q)
  [PASS] Paley(125)/Frobenius: f+s odd
  [PASS] Paley(125)/Frobenius: eigenvalue quadratic irreducible (disc 125 non-square)
  [PASS] Paley(125)/Frobenius: tr(B) = k + b(lam-mu) = 40
  [PASS] Paley(125)/Frobenius: tr(Q^2) = k^2 + b[(lam-mu)^2+2(k-mu)] = 5230
  [PASS] Paley(125)/Frobenius: tr(Q^2) = 2m + 2pG + Sum B^2  (5230 = 10+600+4620)
  [PASS] Paley(125)/Frobenius: 2m = kf - pG  (10 = 310 - 300)
  [PASS] Paley(125)/Frobenius: B row sums k - gamma_i
  [PASS] Paley(125)/Frobenius: pincer C-S inequality satisfied (SumB^2=4620 >= 3550)
      f=5, s=40, b=22, G=100, 2m=10, trB=40, T=2340, SumB^2=4620, CS lower=3550

RESULT: ALL CHECKS PASS — the pincer kills f = 37c for all 1<=c<=8; only f=0 (and the identity, f=333) survive. Theorem verified.
```
