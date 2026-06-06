# H(668) attack pipeline — Z₃₇ conference route (ready to run on Modal)

Goal: construct a Hadamard matrix of order 668 (smallest open order) via a **fixed-point-free
order-37 automorphism** of the conference graph srg(333,166,82,83) ≡ C(334) ≡ H(668).
333 = 9·37 → 9 orbits of 37 → adjacency is a 9×9 array of Z₃₇-circulants.

## The pipeline (all stages validated on small cases)

1. **`conf_core.py`** — block-circulant machinery. Validated: rebuilds Paley P(9)/P(25)/P(49),
   both the direct SRG check and the per-character (`M_t²+M_t=83I`) check agree.
2. **`conference_orbit_all.cpp`** — enumerates the 9×9 orbit-matrix *skeletons*
   (`R²+R−83I=3071J`, spectrum {166, r⁴, s⁴}). Already run → **`skeletons.txt` (634 skeletons)**.
   `clang++ -O2 -o conference_orbit_all conference_orbit_all.cpp && ./conference_orbit_all`
3. **`exact_lift.py`** — CP-SAT *exact* lift for one skeleton: encodes the convolution/SRG
   equations `Σ_{i'}(D[i][i']⊛D[i'][j])(g) = k·[i=j,g=0]+λ·[g∈D[i][j]]+μ·[else]`.
   **Validated** reconstructing Paley P(9)/P(25). **Full pruning:**
   - **Gauge-fixing** (`0∈D[0][i]`): the per-orbit base-point freedom is a **37⁸≈3·10¹²** symmetry;
     fixing it (sound, complete — loses no solutions) is the dominant speedup.
   - **Cyclotomic ansatz** (`mult=m`): force each block `<m>`-invariant (union of cyclotomic cosets).
     Incomplete but **sound + fast + definitive** — proves INFEASIBLE in seconds, or finds a structured H(668).
   Builds the srg(333) model in ~5s; returns **SAT (→H(668))**, **UNSAT/INFEASIBLE (ruled out)**, or **UNKNOWN (timeout)**.
4. **`modal_lift.py`** — Modal app. Fans the lift out over all 634 skeletons across containers.
   - `--mode cyclo` : **run this first** — (skeleton × 8 multipliers) cheap definitive jobs; finds structured H(668) fast or rules each out.
   - `--mode exact` : gauge-fixed CP-SAT per skeleton (complete attempt; in-container SRG verification of any hit).
   - `--mode sa`    : cheap character-domain SA probe (validated objective =0 on Paley).
5. **`emit_h668.py`** — turns a found lift into the matrix: D → A → Seidel S → conference C(334)
   → `H=[[C+I,C−I],[C−I,−C−I]]`, verifies `HHᵀ=668I`. **Validated** P(5)→H(12), P(13)→H(28), P(17)→H(36).

## Run it (after `modal token new`)

```bash
modal token new                                          # one-time auth (your Modal account)
modal run modal_lift.py --mode cyclo --time-s 300        # CHEAPEST+DEFINITIVE: structured search, run first
modal run modal_lift.py --mode exact --time-s 7200       # complete gauge-fixed attack (the long shot)
modal run modal_lift.py --mode sa --seeds 500 --iters 3000000   # optional metaheuristic probe
# if a hit:
python emit_h668.py emit                                 # LIFT_SOLUTION.json -> hadamard_668.csv (verified)
```

**Value-per-dollar order:** `cyclo` first (seconds/job, definitive, finds structured H(668) or rules out
the structured space cheaply) → then `exact` (gauge-fixed, complete, compute-heavy long shot). The gauge-fix
makes every `exact` second productive (no wasted 37⁸ gauge copies).

## Honest status
- The lift **is** the open problem (constructing srg(333) ≡ H(668)). The pipeline is correct
  (validated end-to-end on small cases), so a hit is a genuine H(668). But H(668) has resisted
  experts for 20 years; success is uncertain. This is the best-built, principled, parallel attempt
  available — exact (can find or rule out per skeleton), validated, ready to scale on Modal.
- Skeletons are abundant (634+), so a full "prove no Z₃₇-srg exists" is out of reach; the realistic
  win is a **construction** (one skeleton that lifts), or partial per-skeleton UNSAT certificates.
