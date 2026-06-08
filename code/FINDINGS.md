# Deep dive: the cyclotomic structure of LP(333) / H(668)

Goal of this dive: find a publishable insight on the Legendre-pair route to H(668)
(LP(333), ℓ=333=9·37, the smallest open Hadamard order). What follows is rigorous,
computationally verified, and **honest about what is and is not new.**

## 0. Framework (verified exactly, `lp_core.py`)

For x ∈ {±1}^ℓ: PAF_x(k)=Σ_i x_i x_{i+k}; PSD_x(s)=|Σ_j x_j ω^{js}|², ω=e^{2πi/ℓ}.
A **Legendre pair** (u,v): PAF_u(k)+PAF_v(k)=−2 for all k≢0. Equivalent flat-spectrum form
(re-derived, was mis-stated as ℓ+2 in the project journal):

> **PSD_u(s) + PSD_v(s) = 2ℓ+2  for all s≢0**,  and = 2 for s=0 (row sums ±1).

For ℓ=333 the constant is **668**, literally the Hadamard order. Verified on the exhaustive
LP(9) set: e.g. PSD_u(3)+PSD_v(3)=4+16=20=2·9+2. ✓

Compression identity (Đoković–Kotsireas, re-derived): m-compression c_r=Σ_{j≡r (m)}x_j gives
PAF_c(k)=Σ_{i≡k (m)}PAF_x(i). So for ℓ=333 the 9-compression (c,c') is a length-9 "compressed
LP" with PAF_c(k)+PAF_c'(k) = −74 (k≠0), 594 (k=0); the 37-compression gives −18, 650.

## 1. Cube-root frequency → Eisenstein integers (reproduces Kotsireas–Koutschan 2021)

At s=ℓ/3, write the 3-compression γ=(γ_0,γ_1,γ_2) (γ_r = sum of ℓ/3 ±1's, hence **odd**).
PSD_u(ℓ/3) = γ_0²+γ_1²+γ_2² − γ_0γ_1−γ_1γ_2−γ_2γ_0 = **the Eisenstein norm** |γ_0+γ_1ζ_3+γ_2ζ_3²|².
Because the γ_r are all odd and 1+ζ_3+ζ_3²=0, this value is **divisible by 4**: writing
x=(γ_0−γ_1)/2, y=(γ_1−γ_2)/2 ∈ Z, PSD_u(ℓ/3) = 4(x²+xy+y²). Hence the mod-3 necessary condition is

> **(ℓ+1)/2 = L_u + L_v**, where L_u,L_v are **Loeschian numbers** (values of x²+xy+y²).

For ℓ=333: **167 = L_u + L_v.** Computation (`explore.py`) gives exactly the
Kotsireas–Koutschan spectrum: **16 admissible PSD values** {16,64,76,112,172,256,268,304,364,
400,412,496,556,592,604,652}, **all ≡ 4 (mod 12)**, forming 8 complementary Loeschian pairs
(4+163, 16+151, 19+148, 28+139, 43+124, 64+103, 67+100, 76+91, after dividing by 4); 2016 ordered
mod-3 states. This **reproduces the published filter**, confirms the framework, not new.

## 2. NEW question: does the SQUARE factor (ninth-root, frequency ℓ/9) give a stronger filter?

The literature (verified by scout) derives spectral filters only at frequencies ℓ/p for a single
**prime** p (mod-3, mod-5, mod-p, quaternary-Galois). **No published filter uses a prime-power /
square-divisor frequency ℓ/p².** Since 333=9·37, the natural new idea is a ninth-root filter.

Structure (derived): at s=ℓ/9 the 9-compression c gives PSD_u(ℓ/9)=|Σ_{r=0}^8 c_r ζ_9^r|², again
divisible by 4 (same parity argument, Σζ_9^r=0). This value is a **totally-non-negative algebraic
integer in the cubic field K=Q(ζ_9)⁺=Q(cos 2π/9)**, and PSD_u(ℓ/9)+PSD_v(ℓ/9)=2ℓ+2 holds at all
three conjugate frequencies ⇒ **(ℓ+1)/2 = N_u + N_v in O_K**, both totally non-negative. This is a
genuine cubic-field refinement, and the formulation appears to be new.

### THE TEST (and the honest result, `ninth_root_multi.py`)
Does this cubic-field condition rule out any mod-3-admissible state? Exhaustively, at every
square-divisor length whose 9-compression alphabet is enumerable:

| ℓ | =9·m | (ℓ+1)/2 | #valid 9-compression pairs | #mod-3 states | mod-3 states NOT lifting | ninth-root stronger? |
|---|---|---|---|---|---|---|
| 27 | 9·3 | 14 | 250,128 | 288 | **0** | **No** |
| 45 | 9·5 | 23 (prime, ≡2 mod 3) | 1,163,808 | 288 | **0** | **No** |
| 63 | 9·7 | 32 | 9,385,956 | 756 | **0** | **No** |

**Every mod-3-admissible (γ,γ′) lifts to a valid 9-compressed LP.** The ninth-root / square-divisor
condition is **implied by** the mod-3 condition, it does **not** strengthen the existence filter.

Crucially, **ℓ=45 has (ℓ+1)/2 = 23, prime and ≡ 2 (mod 3)**, the *same arithmetic type as
167 = (333+1)/2*. So the negative result is tested on a case arithmetically representative of 333,
not merely on small toys. Structural reason: the 9-compression space is enormous (millions of valid
pairs) and is a *weak* filter; the cube-root "norm-down to Q" already captures all the spectral
content that the existence question sees at the 3-part.

> **Finding (new, negative):** The natural square-divisor (ℓ/p²) refinement of the
> Kotsireas–Koutschan cyclotomic filter gives **no additional pruning** of the Legendre-pair
> spectrum. The hoped-for "9|ℓ filter" (flagged as a promising lead in the lit review) is not a
> stronger existence filter. Effort on LP(333) should NOT go into square-divisor spectral filters.

## 3. Where the difficulty actually lives (corroborated)

The bottleneck is the **lift** (compressed → full ±1 sequence), exactly as the literature reports.
Evidence from this dive: (a) the spectral filters (mod-3, mod-9) are weak, millions of admissible
compressed states; (b) a naive incremental SA fails to find even LP(45)/LP(63) (`find_lp.py`),
matching the project's "9·prime is hard" benchmark; (c) the open 9|ℓ lengths nearest 333, **207,
225, 279, 297**, have **no known Legendre pair at all** (largest computational non-series LP in the
literature is only 147). The 3²-structure makes the lift hard (project §16: 18 row-autocorrelations
must sum to −2 ≢ 0 mod 9 ⇒ rows cannot be uniform ⇒ no multiplier-invariant solution), and this dive
confirms the spectral side gives no compensating leverage.

## 4. Honest assessment

No publishable *breakthrough* (no construction of LP(333), no route closed). Genuine outputs:
- A clean cubic-field formulation of the ℓ/9 condition (modest, new).
- A rigorous **negative result** killing the square-divisor-filter direction (verified at 27/45/63,
  including the 333-representative arithmetic at ℓ=45), saves future effort, publishable only as a
  note / part of a larger paper.
- Sharpened localization of the frontier: the lift, on lengths where even 207/225 are unconstructed.

The realistic publishable targets remain compute-bound: (i) **construct a first-ever LP at an open
9|ℓ length (207, 225, 279, 297, or 333)** via coprime 9·37 fixed-marginal compression (Turner) +
multi-level 3→3 uncompression (Lumsden) at cluster scale; (ii) a genuinely new construction
exploiting 9·37. Neither is a laptop-session result.

## 5. Built + validated the multiplier-orbit search; calibration verdict (`multiplier_search.py`)

Built the Kotsireas–Koutschan core method (sequence = union of orbits of a multiplier subgroup
M ≤ Z_ℓ*; PSD filter; complementary-PAF matching). **Validated correct:** it re-finds the
quadratic-residue (Legendre) LP at prime lengths 13, 17, 19 via the 3-orbit QR symmetry, exactly.

Then ran it as a calibration probe (laptop, orbit-count ≤ 20–22):
- **Known composite LPs are NOT low-symmetry.** No cyclic-multiplier-invariant LP exists (≤22
  orbits) at 45, 85, 87, **117, 129, 133, 147**, even though all are known to exist (45 by
  exhaustion <48; 85/87/117/129/133/147 by K–K). So K–K's solutions live in **large orbit-unions
  (2²⁵–2³⁰)** reached only via their spectral enumeration on a cluster, not by laptop brute force.
- **All ten open lengths <200** (115,145,159,161,169,175,177,185,187,195): no low-symmetry LP either.

## 6. Z₃₇ conference attack on srg(333,166,82,83) ≡ C(334) ≡ H(668)  (`conf_core.py`, `orbit_*.py`)

The fresh angle the LP/sequence incumbents have not worked. Assume a **fixed-point-free order-37
automorphism** of the conference graph (333=9·37 ⇒ 9 orbits of 37); the adjacency matrix becomes a
**9×9 array of circulants over Z₃₇**. Why order 37 specifically: the eigenvalue-field theorem (the
project's §14, which kills *group-developed* graphs via an order-3 character) **permits** order 37,
because √37 ∈ Q(ζ₃₇) (quadratic Gauss sum). So this symmetry is not pre-obstructed.

**Built + validated machinery.** SRG identity `A²+A−83I = 83J` ⟺ for every nonzero character
t∈Z₃₇ the 9×9 Hermitian `M_t² + M_t = 83 I`. Verified as a true positive on Paley graphs
**P(9), P(25), P(49)** (both the direct check and the character/Fourier check reconstruct them).

**Orbit-matrix gate (necessary first step).** Any such graph forces a 9×9 integer orbit matrix R
(=|D[i][j]|), symmetric, row sums 166, even diagonal, with eigenvalues {166, r×4, s×4}
(r,s=(−1±√333)/2) ⟺ `R² + R − 83I = 3071 J`. Derived facts:
  - Cauchy–Schwarz pins each diagonal entry to dᵢ∈{10,…,26} even (dᵢ=0 impossible); diagonal sums 162.
  - The sum-of-squares constraint forces each row **near-uniform** (variance ≈0.25; e.g. row-0 must be
    a permutation of four 20s and four 19s).
  - **No spectral obstruction:** an integer matrix with M²=333I exists (333=18²+3² ⇒ [[18,3],[3,−18]]),
    so any nonexistence must be combinatorial, not from the spectrum.

**Status of the gate (honest, UNDECIDED):** CP-SAT did NOT resolve it (UNKNOWN at 240s and at 600s, *not* infeasible). The dedicated DFS (`orbit_dfs.py`, 910 candidate diagonals) did not converge in
900s either, it stalls materializing the (large) completion list for rows whose sum-of-squares
target permits looser variance (an implementation gap: no time-check inside the per-row enumerator).
**Net: orbit-matrix existence is UNRESOLVED. No nonexistence claim is made.** The gate is a clean,
bounded 9×9 integer-feasibility problem that a commercial solver (Gurobi) or a compiled/streaming
backtracker would settle definitively, a *small* compute purchase that returns a definite answer.

**RESOLVED (`conference_orbit.cpp`, compiled DFS):** the orbit-matrix gate is **PASSED**, a valid R
EXISTS (found in 47s; the CP-SAT timeouts were a solver-fit problem, not infeasibility). Verified R:
diagonal (10,14,14,14,18,20,24,24,24), off-diagonals in {15,17,18,19,20,21,22,24}, row/col sums 166,
`R²+R−83I = 3071J` exactly, spectrum {166, r×4, s×4}. So a Z₃₇ skeleton for srg(333) provably exists, the route is NOT obstructed at the gate (unlike the proven-empty group-developed and single-DS routes).
This is a genuinely new result for this parameter set (no prior non-regular Z₃₇ analysis).

**The lift (the open wall), `conference_lift.py`:** realize R as actual Z₃₇ connection sets so every
character condition `M_t²+M_t=83I` holds. Objective validated =0 on Paley P(9)/P(25)/P(49). On the
srg(333) skeleton, local search (SA, slow-cool + reheat, 12M iters, multiple seeds) **plateaus
robustly at objective ≈ 397k–402k (target 0)**, the same route-independent barrier the project
documents everywhere. Also: the found R has entries {17,19,20,22} ∤ 3, so it cannot carry uniform
cyclotomic-coset blocks, the natural "structured" shortcut does not fit this skeleton. So the lift is
genuinely the open problem (≡ constructing srg(333) ≡ H(668)); laptop search does not crack it.

**Net for the $4k:** the gate is passed cheaply (done). The lift is the real sink and is effectively
the open problem, it requires cluster-scale **structured SAT+CAS** over the now-concrete skeleton (a
far better-scoped target than "search 2³³³", but still uncertain, success would BE H(668)). Honest
recommendation: the $4k buys a *real attempt* at the lift with a verified skeleton, but it is a
research-frontier gamble, not a likely win. The LP route (below) is not competitive for $4k at all.

## 7. LP-finding calibration (the multiplier-orbit method)

**Calibration verdict (this is the answer to "should we spend the compute"):** the LP-finding
bottleneck is **competitive compute, not dollars.** The incumbent group (Kotsireas et al.) has the
full spectral+orbit pipeline AND GPU clusters (2026 HPC paper, 6500× speedup), exhaustive search
only reaches ℓ≈55, and they have targeted these exact open lengths without cracking them. A few
thousand dollars running the *same* method cannot surpass that. **Do not spend $4k on the LP brute
search**, it is not a $4k problem; it is a "more compute than the world leaders, with the same
tools" problem. The only LP opening is a *different structure they didn't try* (low-probability
gamble), or a new method/theorem (not purchasable with raw cores).

## 8. Modal cyclo sweep result (2026-06-03), structured Z_37 lift is EMPTY
Ran modal_lift.py --mode cyclo on Modal: 634 enumerated orbit-matrix skeletons x 8 multiplier
subgroups of Z_37* = 5072 CP-SAT jobs (gauge-fixed + cyclotomic-coset-invariant ansatz).
RESULT: **5072/5072 INFEASIBLE, 0 SAT, 0 timeout**, every cyclotomic-structured Z_37 lift
is definitively ruled out (each in seconds). So NO Z_37-cyclotomically-structured srg(333,166,82,83)
exists over the 634 enumerated skeletons. Corroborates the project's pervasive "structured solutions
are empty for 333" pattern. CAVEAT: rules out only the cyclotomic-coset-structured class over the
SAVED skeleton sample (true skeleton count higher); NOT a full nonexistence. The remaining shot is the
exact (unstructured, gauge-fixed, complete) sweep, but the cyclo negative makes it a LONGER shot
(structured solutions, the likeliest to exist, are gone).
