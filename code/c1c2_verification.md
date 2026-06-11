# Independent verification: Paley-pinning of D[7][7] + per-orbit coupling (C1)/(C2)

**Date:** 2026-06-10
**Subject:** the "bonus theorem-grade necessary condition" produced by the SA-floor-lower-bound
mechanism agent (workflow `wf_30fa9ce0-897`, transcript `agent-a7fd1d7cf53590ebd.jsonl`).
**Verifier:** independent re-derivation + from-scratch numerical re-implementation
(`c1c2_indep_verify.py`, this directory; runs in <1 s).

**VERDICT: SOUND** (theorem-grade necessary condition, modulo the machine-verified gauge inputs
listed in §6; it is necessary-only and does not by itself close classes 221/422).

---

## 1. Provenance

The floor agent's transcript is
`~/.claude/projects/-Users-maceocardinalekwik-git-math/aee5d83d-c167-4240-8abd-3a078b06c680/subagents/workflows/wf_30fa9ce0-897/agent-a7fd1d7cf53590ebd.jsonl`
(260 lines, 173 assistant messages; final `StructuredOutput` at line 257). Its mechanism:
prove the SA objective floor min F > 0 for the σ-equivariant self-complementary Z37
classes 221/422 of srg(333,166,82,83). Main result was **negative** (floor not proven:
the 79-dim affine-invariant space has all constants 0; the McCormick LP relaxation is
exactly 0 with a fractional witness; trace moments collapse). The bonus result is the
necessary condition verified here. Code it left behind (all read for this verification):

- `~/git/math/release/h668-z37/code/sc_trace_condition.py` — the theorem + its validation
- `~/git/math/release/h668-z37/code/sc_floor_analysis{,2,3}.py` — dump certification, invariant hunt
- `~/git/math/release/h668-z37/code/sc_pinned_cpsat.py` — CP-SAT with (C1)/(C2) added
- `~/git/math/release/h668-z37/code/sc_pinned_floor_probe.py` — pinned-landscape probe

---

## 2. Setting

A Z37-symmetric candidate srg(333,166,82,83): vertices (i,g), i ∈ {0,…,8} (9 orbits),
g ∈ Z37; adjacency (i,g) ~ (j,h) ⟺ h−g ∈ D[i][j] ⊆ Z37, with D[j][i] = −D[i][j],
D[i][i] symmetric, 0 ∉ D[i][i]. Orbit matrix R[i][j] = |D[i][j]| (classes 221, 422 from
the 2901-orbit-matrix census). Parameters k=166, λ=82, μ=83, so the SRG equation is
A² + A − 83·I = 83·J.

**Complementing isomorphism** σ(i,g) = (π(i), 6g), where 6² ≡ −1 (mod 37) and 6 is a QNR
(37 ≡ 5 mod 8, so the order-4 units are QNR). σ maps the graph onto its complement, which
forces (verified elementwise on all 6 templates):

- off-diagonal: D[π(i)][π(j)] = 6·(Z37 ∖ D[i][j]);
- diagonal: **D[π(a)][π(a)] = 6·(Z37\* ∖ D[a][a])**.

For classes 221/422, π has cycle type 4+4+1 on the 9 orbits with **fixed orbit f = 7**.
Writing X = D[a][a] for a on a π-4-cycle, the relation iterates to
D[π²(a)][π²(a)] = 36·X = −X = X (X symmetric), so each 4-cycle's diagonal blocks are
{X, 6·Xᶜ, X, 6·Xᶜ} (Xᶜ = Z37\* ∖ X). The fixed orbit satisfies the **anticlosure**
D[7][7] = 6·(Z37\* ∖ D[7][7]), i.e. d ∈ D77 ⟺ 6d ∉ D77; combined with symmetry, on each
⟨±6⟩-orbit O = {d, −d, 6d, −6d} the set D77 contains exactly one half ({±d} or {±6d}),
hence |D77| = 18 automatically and D77 has 2⁹ = 512 states.

Free representative blocks (template REP lines): S = D[0][0] (cycle-diagonal, |S| = 12
for 221, 14 for 422), T = D[t][t] (t = 2 for 221 with |T| = 18; t = 3 for 422 with
|T| = 20), U = D[7][7] (anticlosed, |U| = 18), plus off-diagonal reps. Total 415 free
bits; SA objective F = Σ over the 13 representative pair-blocks (a,b), g ∈ Z37, of the
squared residual of the SRG equation. σ-orbit ties make F = 0 ⟺ all 81·37 SRG
equations hold ⟺ a true srg(333,166,82,83) in the class exists.

Each ⟨±6⟩-orbit splits 2 QR + 2 QNR elements (−1 ∈ QR, 6 ∈ QNR): QR-half {±d} with
d ∈ QR(37), QNR-half {±6d}.

---

## 3. Exact statements

Let QR = quadratic residues mod 37 (18 elements), QNR the non-residues. For each of the
9 ⟨±6⟩-orbits O (QR-rep d) define the half-indicators

- s₁ = [±d ⊆ S], s₂ = [±6d ⊆ S] (well-defined: S symmetric), likewise t₁, t₂ for T,
- u_O = [±d ⊆ U] ∈ {0,1} (anticlosure: U contains exactly one half).

**THEOREM (necessary condition).** If F = 0 (i.e. the class instance carries a true
srg(333,166,82,83) with this σ-equivariant structure), then exactly one of two branches
holds:

> **(C1) Paley pinning of D[7][7]:** U = D[7][7] is **exactly the QR set** or **exactly
> the QNR set** of Z37 — 2 states out of 512.
>
> **(C2) per-orbit coupling of D[0][0] and D[t][t]:** with ε = −1 if U = QR and ε = +1
> if U = QNR, **every** one of the 9 ⟨±6⟩-orbits satisfies
>
> (s₁ − s₂) + (t₁ − t₂) = ε.
>
> Equivalently, the diagonal multiplicity function f(d) = #{a : d ∈ D[a][a]} is constant
> on QR (value u) and on QNR (value v), with (u,v) = (3,6) on the QR branch and (6,3) on
> the QNR branch.
>
> **(C3) parity corollary:** |S|/2 + |T|/2 must be odd. (Holds for both classes —
> 221: 6+9 = 15; 422: 7+10 = 17 — so counting alone does not close them.)

---

## 4. The proof (reconstructed, every step, independently checked)

**Step 1 (Fourier/character reduction).** For each character χ_j (j ≠ 0) of Z37 define the
9×9 matrix M(χ)[a][b] = Σ_{d∈D[a][b]} χ(d). D[b][a] = −D[a][b] makes M(χ) **Hermitian**.
The block-circulant structure diagonalizes A into ⊕_χ M(χ); since the all-ones J has
zero component at χ ≠ 1, the SRG equation A²+A−83I = 83J is equivalent (given the row-sum
component, see Step 1') to

  M(χ)² + M(χ) − 83·I₉ = 0  for all χ ≠ 1.

*Step 1' (χ = 1 component is free).* Summing the SRG equation over g gives
(R² + R − 83I)[a][b] = 37·83 = 3071 for every pair, which holds identically for every
valid orbit matrix R of the census — so the χ=1 residual component is always 0 and F = 0
is exactly the vanishing of all M(χ) quadratics. (Numerically certified:
37·F_full = Σ_{χ≠1} ‖M(χ)²+M(χ)−83I‖²_F, exact equality 467088 = 467088 on the archived
obj=1852 state — both by the agent and by me.)

**Step 2 (trace lattice).** x² + x − 83 has roots θ, τ = (−1 ± √333)/2 = (−1 ± 3√37)/2.
A Hermitian matrix annihilated by this separable quadratic has all eigenvalues in {θ, τ},
so tr M(χ) = m·θ + (9−m)·τ for some integer m ∈ {0,…,9}; in particular
**tr M(χ) ∈ Q(√37)** for every χ ≠ 1.

**Step 3 (Galois ⇒ QR-class constancy).** tr M(χ_j) = f̂(j) where f(d) = #{a : d ∈ D[a][a]}
(only diagonal blocks contribute to the trace), f(0) = 0. The Galois group of Q(ζ37)/Q is
(Z/37)ˣ acting by σ_s(ζ) = ζ^s; Q(√37) is its unique quadratic subfield, fixed exactly by
the index-2 (QR) subgroup. For s ∈ QR: σ_s(f̂(j)) = f̂(js), and σ_s fixes f̂(j) ∈ Q(√37),
so f̂(js) = f̂(j) — f̂ is constant on the two QR-cosets. Fourier inversion (QR-class
functions form an algebra closed under the transform; −1 ∈ QR) gives **f constant on QR
(= u) and on QNR (= v)**.

**Step 4 (lattice arithmetic).** Σ_d f(d) = Σ_a |D[a][a]| = tr R = 162 for both classes
(verified), so 18(u+v) = 162, i.e. **u + v = 9**. With η₀ = Σ_{d∈QR} ζ^d = (−1+√37)/2,
η₁ = (−1−√37)/2 (Gauss periods, 37 ≡ 1 mod 4):
tr M(χ_j) = u·η₀ + v·η₁ (j ∈ QR) or u·η₁ + v·η₀ (j ∈ QNR) = −9/2 ± (u−v)·√37/2.
Matching m·θ + (9−m)·τ = −9/2 + (2m−9)·3√37/2: rational parts agree automatically
(consistency check: tr R = 162 is exactly what makes −(u+v)/2 = −9/2), and the irrational
parts force u − v = ±3(2m−9) ∈ 3·(odd). Since 0 ≤ u, v ≤ 9: **u − v ∈ {±3, ±9}**.

**Step 5 (gauge formula for f).** From §2 the diagonal blocks are
{S, 6Sᶜ, S, 6Sᶜ} ∪ {T, 6Tᶜ, T, 6Tᶜ} ∪ {U}, so for d ≠ 0:

  f(d) = 2[1_S(d) + 1_{6Sᶜ}(d)] + 2[1_T(d) + 1_{6Tᶜ}(d)] + 1_U(d)
       = 4 + 1_U(d) + 2[g_S(d) + g_T(d)],  g_X(d) := 1_X(d) − 1_X(6⁻¹d),

using 1_{6Xᶜ}(d) = 1 − 1_X(6⁻¹d). Since 6⁻¹ = −6 and S, T are symmetric, on each orbit
(QR-rep d): g_X(d) = x₁ − x₂ and g_X(6d) = −g_X(d). Hence with e(O) := (s₁−s₂)+(t₁−t₂):

  f(QR-half of O) = 4 + u_O + 2e(O),  f(QNR-half of O) = 5 − u_O − 2e(O).

(u + v = 9 per orbit automatically — cross-check of Step 4.)

**Step 6 (pinning).** QR-constancy (Step 3) means u_O + 2e(O) = (u−v+1)/2 is the **same**
for all 9 orbits. Exhaustive check of the 32 per-orbit states (u_O,s₁,s₂,t₁,t₂) ∈ {0,1}⁵
(re-done independently, §7 check [3]) gives the unique decompositions:

| u−v | u_O (forced) | e(O) (forced) | meaning |
|---|---|---|---|
| −3 → (u,v)=(3,6) | 1 on all orbits | −1 | U = **QR** exactly; (C2) with ε = −1 |
| +3 → (u,v)=(6,3) | 0 on all orbits | +1 | U = **QNR** exactly; (C2) with ε = +1 |
| +9 → (9,0) | 1 | +2 | needs s₁−s₂ = t₁−t₂ = 1 ∀O ⇒ \|S\| = \|T\| = 18 — **killed** (\|S\| = 12/14) |
| −9 → (0,9) | 0 | −2 | needs s₁−s₂ = t₁−t₂ = −1 ∀O ⇒ \|S\| = \|T\| = 18 — **killed** |

u_O constant = 1 (resp. 0) on all orbits means U is the union of all QR-halves = QR
(resp. all QNR-halves = QNR): **(C1)**. The surviving rows are **(C2)**. ∎

**(C3):** e(O) = ±1 is odd ⇒ s₁+s₂+t₁+t₂ odd per orbit ⇒ summing over 9 orbits,
|S|/2 + |T|/2 ≡ 9 ≡ 1 (mod 2). ∎

The converse direction also holds for these classes: (C1)∧(C2) ⇒ f = (3,6)/(6,3)
QR-constant ⇒ every tr M(χ) on the required lattice — so (C1)∧(C2) is exactly the full
strength of the trace condition (nothing was lost in the reduction).

**Precedent:** for norb = 1 the same argument is Bridges–Mena 1979 (rational-spectrum
circulant SRGs are Paley): tr M(χ) = Σ_{d∈D} χ(d) ∈ {θ_p, τ_p} pins D ∈ {QR, QNR}. The
condition is the 9-orbit-quotient generalization.

---

## 5. How the agent validated it

1. **Data certification first** (`sc_floor_analysis.py`): all 12 BESTSTATE dumps in
   `~/git/math/deepdive_lp333/sc_residuals.json` re-expanded against the 6 archived
   templates; objectives 1852, 1876, 1714, 2528, 1754, 1876, 1796, 1672, 1824, 2146,
   1674, 1748 recomputed exactly, 0 per-equation residual mismatches.
2. **Gauge lemma machine-checked:** D[π(a)][π(a)] = 6·(Z37\* ∖ D[a][a]) verified
   elementwise on all 6 templates (transcript line 131: all True).
3. **Equivalence check (a):** [tr M(χ) ∈ Q(√37) ∀χ] ⟺ [f QR-constant] verified
   numerically on random states (Galois invariance tested as tr(χ) = tr(χ^s), s = 3 ∈ QR).
4. **Anchor (b):** Paley(13) (the analog single-orbit self-complementary case, c = 5,
   srg(13,6,2,3), template emitted through the same pipeline) — the condition **holds**
   on D = QR(13). The agent caught and fixed an operator-precedence bug in its own first
   check (transcript line 98) — initial False became True after the fix.
5. **Violation (c):** all 12 archived floor states violate QR-constancy, L2 deviation
   19.56–171.56; 200/200 random states violate it (min L2dev 16.9).
6. **Roundtrip (e):** 50 states sampled directly from the (C1)∧(C2) subspace all give
   f = (3,6) and traces exactly on the {mθ+(9−m)τ} lattice (transcript line 246).
7. **Feasibility windows (d):** (C3) parity holds for both classes (15 and 17 odd), and
   the per-orbit type counts admit solutions ((nA,nB,β,α) windows printed), so the
   condition reduces but does not close the classes.
8. **Deployment:** `sc_pinned_cpsat.py` adds (C1)/(C2) to the hostile-validated
   `selfcomp_search.Model`; smoke runs (cls 221, perm (1,5,3,4,8,6,0,7,2), 240 s, both
   branches) returned UNKNOWN — consistent with "reduction, not closure".

**File read for the 12-state claim:** `/Users/maceocardinalekwik/git/math/deepdive_lp333/sc_residuals.json` —
12 entries (6 templates `sc_tpl_{221_0,221_1,422_2,422_3,422_4,422_5}.txt` × seeds
101/102), each a BESTSTATE dump (objective + REPBITS for the 13 representative blocks +
nonzero RES lines). What was checked per state: exact reconstruction, then f(d) computed
from the 9 expanded diagonal blocks and tested for QR-class constancy (violated in all
12 cases).

---

## 6. Independent verification (this review)

Re-implemented from scratch — template parser, expansion, direct (non-FFT) convolution,
residuals, f, traces, Gauss periods — in `c1c2_indep_verify.py` (this directory; total
runtime 0.18 s). Results:

1. **All 12 dumps reconstruct exactly** (my objectives = claimed objectives), and all 12
   **violate** QR-constancy and the trace lattice. My L2 deviations match the agent's to
   the digit: 171.56, 46.22, 35.56, 94.22, 72.00, 78.22, 62.22, 48.00, 19.56, 126.22,
   88.00, 43.56 (range 19.56–171.56; agent reported 19.6–171.6). 12/12 violations confirmed.
2. **Gauge structure confirmed on all 6 templates:** the diagonal complement relation has
   a consistent orbit-matching with fixed orbit {7}, and the 4-cycle return
   D[π²(a)][π²(a)] = D[a][a] holds on every diagonal block.
3. **Derivation core exhaustively checked** over all 32 per-orbit states: lattice-reachable
   (u,v) are exactly {(0,9),(3,6),(6,3),(9,0)}; (3,6) forces (u_O,e) = (1,−1), (6,3)
   forces (0,+1); (9,0)/(0,9) force e = ±2, killed by |S| ∈ {12,14} ≠ 18.
4. **Equivalence** [trace lattice ⟺ f QR-constant] holds on 300/300 random states
   (and 0 random states satisfy it).
5. **Roundtrip both branches:** 25+25 random (C1)∧(C2) states all give f = (3,6) (QR
   branch) / (6,3) (QNR branch) with every trace exactly on the lattice.
6. **Fourier bridge exact:** 37·F_full = Σ_{χ≠1}‖M(χ)²+M(χ)−83I‖²_F = 467088 on dump 0;
   all M(χ) Hermitian.
7. **Paley anchors (the analog condition HOLDS where solutions exist):**
   - Paley(13), complementing multiplier c = 5 (QNR, 5² = −1): D = QR(13) is anticlosed,
     satisfies the SRG equation, and satisfies the trace condition. **Exhaustively**, of
     the 2³ = 8 anticlosed symmetric candidates, exactly 2 satisfy it: QR(13) and
     QNR(13) — the pinning is sharp and reproduces Bridges–Mena.
   - Paley(17), complementing multiplier c = 3 (QNR; note 17 ≡ 1 mod 8 so the order-4
     units are QR and the c²=−1 orbit decomposition does not apply — the trace-lattice
     form of the condition is the right analog): D = QR(17) is anticlosed, SRG holds,
     trace condition holds.
8. **(C3) parity** re-derived and consistent: 221: 15 odd, 422: 17 odd (no closure by
   counting).

**Scrutinized for gaps:**
- Hermiticity of M(χ) (needed for the eigenvalue argument): follows from D[b][a] = −D[a][b]; verified numerically. Even without Hermiticity, annihilation by the separable quadratic gives diagonalizability — the step is robust.
- The χ=1 component: identically zero for every valid census orbit matrix (R²+R−83I = 3071·J), so F = 0 genuinely reduces to the χ≠1 quadratics. Verified.
- f receives contributions only from diagonal blocks; tr R = 162 makes the rational part consistent rather than adding a constraint. Verified.
- The objective F covers only the 13 representative pair-blocks, but σ-orbit residual ties (verified by the agent, line 123: all 45 unordered pair-blocks tie to rep blocks) make F = 0 equivalent to all equations. Sound.
- The (9,0)/(0,9) exclusion is class-specific (uses |S| ≠ 18). For a hypothetical class with |S| = |T| = 18 those branches would also need S = T = QR (or QNR) exactly; irrelevant for 221/422.

**VERDICT: SOUND.** The condition is a genuine theorem for classes 221/422:
F = 0 ⇒ (C1) ∧ (C2) (⇒ (C3)). It is **necessary-only**: it reduces the search space
(D77: 512 → 2; per-orbit (s₁,s₂,t₁,t₂): 16 → 4, i.e. roughly 2⁹⁺⁹ ≈ 2¹⁸-fold raw
reduction before cardinality interaction) but does not prove min F > 0 — indeed the
agent's pinned-landscape probe found comparable glassy floors inside the pinned corner
(greedy floors ≈ 2900), so the hardness lives in the off-diagonal blocks.

---

## 7. Precise CP-SAT constraints to add

The self-complementary model (`selfcomp_search.Model`) is at **element level**: bits
x[i][j][d] := [d ∈ D[i][j]] for the free representative blocks (i,j) ∈ orbit reps,
d ∈ Z37, with all other (a,b,d) memberships resolved to literals (possibly negated) of
rep bits via the σ-chains; cosets are trivial (H = {1}), so y[i][j][c] ≡ x[i][j][d].
Already in the model: symmetry x[a][a][d] = x[a][a][37−d], anticlosure
x[7][7][6d mod 37] = 1 − x[7][7][d], and the cardinality sums.

Let f = 7 (fixed orbit), S-rep = (s₀,s₀) with s₀ = 0 (both classes), T-rep = (t₀,t₀)
with t₀ = 2 (class 221) / t₀ = 3 (class 422). QR = {d² mod 37} = {1,3,4,7,9,10,11,12,16,
21,25,26,27,28,30,33,34,36}. The 9 ⟨±6⟩-orbit QR-reps: for each orbit {d,−d,6d,−6d}
pick d ∈ QR (one per orbit; e.g. iterate x = 1..36, skip seen, rep = x if x ∈ QR else 6x mod 37).

**Two-branch form (run as two separate models; INFEASIBLE on both ⇒ instance closed):**

Branch QR (ε = −1):
```
for d in 1..36:  Add( x[7][7][d] == (1 if d in QR else 0) )
for each ⟨±6⟩-orbit with QR-rep d, d6 = 6*d mod 37:
    Add( x[s0][s0][d] - x[s0][s0][d6] + x[t0][t0][d] - x[t0][t0][d6] == -1 )
```

Branch QNR (ε = +1):
```
for d in 1..36:  Add( x[7][7][d] == (0 if d in QR else 1) )
for each ⟨±6⟩-orbit with QR-rep d, d6 = 6*d mod 37:
    Add( x[s0][s0][d] - x[s0][s0][d6] + x[t0][t0][d] - x[t0][t0][d6] == +1 )
```

**Single-model form (one BoolVar b, b = 1 ⟺ QR branch):**
```
for d in QR\{0}:   Add( x[7][7][d] == b )
for d in QNR:      Add( x[7][7][d] == 1 - b )
for each orbit (QR-rep d, d6 = 6d):
    Add( x[s0][s0][d] - x[s0][s0][d6] + x[t0][t0][d] - x[t0][t0][d6] == 1 - 2*b )
```

(Optional redundant-but-propagation-friendly corollaries: the per-orbit constraint plus
cardinalities imply the (nA, nB, β, α) windows — 221: nA ∈ {0,2,4,6}; 422: nA ∈ {1,3,5,7} —
which can be added as a channel on per-orbit type indicators if branching needs help.
A same-derivation strengthening, not yet done: e₃(M(χ)) ∈ Q(√37), cubic constraints.)

This is exactly what `sc_pinned_cpsat.py` implements (verified against the derivation).
Status so far: cls 221 perm (1,5,3,4,8,6,0,7,2), both branches UNKNOWN at the 240 s smoke
budget; the recommended full-budget run is 6 instances × 2 branches × ≥2.5 h × 8 workers.

---

## 8. Summary table of claims vs. verification

| Claim (agent) | Independent check | Status |
|---|---|---|
| 12 dumps in sc_residuals.json genuine (objs 1852…1748) | re-expanded, recomputed objs match 12/12 | CONFIRMED |
| D[π(a)][π(a)] = 6·(Z37\*∖D[a][a]) on all templates | re-checked elementwise, + π²=id on diagonals | CONFIRMED |
| F=0 ⇒ tr M(χ) ∈ {mθ+(9−m)τ} | re-derived (Hermitian + separable quadratic) | SOUND |
| trace lattice ⟺ f QR-constant | re-derived (Galois) + 300/300 random states | CONFIRMED |
| (u,v) ∈ {(3,6),(6,3)} only | exhaustive 32-state per-orbit algebra + cardinality kill of (9,0),(0,9) | CONFIRMED |
| (C1) U ∈ {QR, QNR} exactly; (C2) per-orbit sum = ε | forced uniquely; roundtrip 50/50 states | CONFIRMED |
| (C3) parity odd, satisfied by both classes (no closure) | 15, 17 odd | CONFIRMED |
| All 12 floor states violate (L2dev 19.6–171.6) | recomputed: 19.56–171.56, 12/12 violate | CONFIRMED |
| Paley anchors hold | Paley(13): holds + sharp (exactly {QR,QNR} of 8); Paley(17) c=3: holds | CONFIRMED |
| Generalizes Bridges–Mena 1979 | norb=1 specialization is exactly that statement | CONFIRMED |
