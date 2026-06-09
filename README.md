# A Z₃₇ reduction and cyclotomic nonexistence results for the conference graph of H(668)

Validated machinery and the accompanying paper for an attack on the order-668 Hadamard
matrix problem (the smallest order for which existence is currently open) via a
fixed-point-free automorphism of order 37 of the strongly regular (conference) graph
**srg(333, 166, 82, 83)** ⟺ symmetric conference matrix **C(334)** ⟺ **H(668)**.

> **Status.** H(668) remains open. This repository contributes a reduction, a complete enumeration of
> the 2,901 orbit-matrix skeletons, nonexistence of cyclotomically-structured Z₃₇-symmetric solutions
> for **6 of the 8** multiplier subgroups (orders 2 and 3 are open), and the reproducible code +
> certificates behind them. It does **not** claim a construction of H(668).

The paper is in [`paper/h668_z37.pdf`](paper/h668_z37.pdf) (source `paper/h668_z37.tex`).

## Results in one paragraph

1. **Cyclotomic filter (negative).** The natural *square-divisor* refinement of the
   Kotsireas–Koutschan cube-root filter for Legendre pairs of length 333 = 9·37
   (evaluating the power spectrum at a primitive 9th root of unity) adds **no** pruning
   beyond the mod-3 (cube-root) condition.
2. **Z₃₇ reduction.** Assuming a fixed-point-free order-37 automorphism (which the
   eigenvalue-field obstruction *permits*, unlike the group-developed/Cayley case),
   srg(333,166,82,83) reduces to a 9×9 array of Z₃₇-circulants, with strong regularity
   equivalent to a per-character condition `M_t² + M_t − 83I = 0` over Z[ζ₃₇].
3. **Cyclotomic nonexistence (6 of 8, by lift-free arguments).** No `H`-cyclotomic
   (Paley-type) Z₃₇-symmetric srg(333,166,82,83) exists for the **6** proper nontrivial
   multiplier subgroups of order ≥ 4, `|H| ∈ {4,6,9,12,18,36}`: no coset-compatible orbit
   matrix exists at all (SAT-certified UNSAT, independent of any orbit matrix), with an
   independent orbit-matrix-free confirmation for `{12,18,36}`. Both arguments are
   **lift-free** (no base-point gauge), hence unconditional. We also **completely enumerate
   the 2,901 admissible orbit matrices** (the necessary skeletons; sound symmetry reduction,
   every diagonal exhausted).

   > **The two finest multipliers `|H| ∈ {2,3}` are OPEN.** They reduce to a per-skeleton
   > cyclotomic *lift* over the 2,901 orbit matrices, which we do not resolve. A base-point
   > gauge-fixing that accelerates the general lift is **unsound** under a multiplier
   > constraint (it wrongly rejects the genuine negation-invariant Paley graph P(25)), and
   > without it no sound SAT encoding decides a single `|H|∈{2,3}` instance within 30 min.
   > An earlier version of this repo reported these cases as closed; that relied on the
   > unsound gauge-fixing and has been **retracted**. All evidence (below) still points to
   > nonexistence, but it is not proven. Large-scale GPU parallel-tempering finds no
   > unstructured lift either.
4. **The order-2 residue: a sharper ceiling** (paper §4, Prop. order2). The remaining
   order-2 case is *narrowed and explained*:
   - **Narrowing:** a Gauss-sum/Galois argument pins the θ-multiplicity of each
     residue-character to `a ∈ {3,4,5,6}`; `a ∈ {3,6}` force a constant-18 orbit-matrix
     diagonal, which an exhaustive search (sound `S₉` symmetry break; 1.9×10⁸ nodes,
     corroborated to 9×10¹⁰) proves **does not exist**, so `a ∈ {4,5}`.
   - **Why it resists (two meta-theorems):** the group-ring RHS lives entirely on the
     principal character (no character/Gauss-sum/multiplier obstruction can fire; `83`
     isn't self-conjugate mod 37), and the negation involution collapses *every*
     trace-linear functional (spectral, Krein, scheme, counting) to five fixed values.
     Any genuine obstruction must be **non-linear** in the {0,1} indicators, i.e. the
     lift feasibility itself, which is the order-668 open problem. (Order 2 stays open.)

## Repository layout

```
paper/        h668_z37.tex, h668_z37.pdf   the paper
lean/         Lean 4 / Mathlib formalization of the collapse theorem (machine-checked, no sorry)
code/         all validated machinery + data + reference notes
  conf_core.py              block-circulant SRG machinery (validated on Paley graphs)
  conference_z37.py         the Z₃₇ reduction
  conference_orbit.cpp       C++ DFS that found the worked-example orbit matrix
  conference_orbit_all.cpp   full orbit-matrix enumerator (streams skeletons)
  conference_orbit_sym.cpp  ★ COMPLETE orbit-matrix enumerator (sound symmetry break) → 2901 reps
  conference_orbit_certify.cpp ★ certificate: coset-compatible count = 0 for |H|∈{4,6,9,12,18,36} (Thm 1a)
  verify_certificates.py    ★ verifies the 2901 enumeration + reruns the Thm 1(a) certificates
  cyclo_lift_coset.py       ★ SOUND coset-level cyclotomic lift (no gauge; validates on P(9)/P(25))
  orbit3_enum.cpp           order-3-compatible enumerator → 960 skeletons (|H|=3 substrate)
  order3_sweep.py           RETRACTED: order-3 lift sweep via unsound gauge-fix (kept as a pitfall record)
  order2_spectral.py        order-2 (negation-invariant) lift with N_g/a∈{4,5} pruning
  uniform18_enum.cpp        ★ proves no constant-18 orbit matrix exists → a∈{3,6} eliminated
  uniform18_maxbreak.cpp    independent re-check (different sound symmetry break)
  order2_combinatorial.py   the involution-collapse / flat-RHS meta-theorem checks
  modal_order2.py           Modal driver for the a∈{4,5} order-2 frontier search
  orbit_matrix_search.py    CP-SAT model for admissible 9×9 orbit matrices
  orbit_compat.py           ★ COMPLETE cyclotomic closure (Theorem 1(a)), the key result
  cyclo_complete.py         orbit-matrix-free cyclotomic check (Theorem 1(c))
  exact_lift.py             CP-SAT exact lift (convolution encoding + gauge-fixing)
  gpu_lift.py               GPU parallel-tempering lift search (PyTorch)
  emit_h668.py              lift → C(334) → H(668), with end-to-end validation
  multiplier_search.py      Legendre-pair multiplier search
  ninth_root.py, ninth_root_multi.py   the square-divisor cyclotomic filter (Prop 1)
  validate_search.py        ground-truth validation of the multiplier search
  modal_cyclo.py, modal_lift.py, modal_gpu_lift.py   Modal (serverless) drivers for scale
  skeletons.txt (634), skeletons_full.txt (≥1678)    enumerated orbit-matrix skeletons
  orbit_all_sym.txt (2901)      COMPLETE set of admissible orbit matrices (all necessary skeletons)
  orbit3_compatible.txt (960)   order-3-compatible orbit matrices (the |H|=3 substrate)
  orbitR_found.npy          the verified worked-example 9×9 orbit matrix
  LP{13,17,19}_{u,v}.npy    small Legendre pairs used to validate the pipeline
  FINDINGS.md, README_H668_ATTACK.md   detailed research notes
```

## Quick start

```bash
cd code
python -m pip install -r ../requirements.txt    # numpy, ortools
```

**Validate the pipeline (fast, local):**

```bash
python emit_h668.py            # validates lift → C(334) → H(668) on H(12), H(28), H(36)
python exact_lift.py validate  # reconstructs Paley P(9), P(25) via the exact CP-SAT lift
python cyclo_complete.py validate   # reconstructs Paley P(13), P(17) (cyclotomic model)
python validate_search.py      # multiplier search vs brute-force ground truth
python gpu_lift.py             # (needs PyTorch) finds the P(25) lift by parallel tempering
```

**Reproduce the complete (lift-free) nonexistence results — orders ≥ 4 (6 of 8):**

```bash
python orbit_compat.py 300
#   order  9,6,4,12,18,36: INFEASIBLE  -> CLOSED (no such srg, any orbit matrix)
#   order  2,3:            UNKNOWN     -> orbit-matrix argument does not apply; OPEN (see note)
```

**Complete enumeration of the 2,901 admissible orbit matrices (the necessary skeletons):**

```bash
clang++ -O3 -o conference_orbit_sym conference_orbit_sym.cpp && ./conference_orbit_sym
#   -> 2901 representatives, every diagonal exhausted (sound symmetry reduction)
```

**Verifiable certificates (trust no solver):**

```bash
clang++ -O3 -o cert conference_orbit_certify.cpp
python verify_certificates.py
#   [complete enumeration] orbit_all_sym.txt: 2901 matrices ... VERIFIED
#   [Theorem 1(a)] coset-compatible orbit-matrix counts (independent C++ search):
#     e=4,6,9,12,18,36: count 0 -> NO H-cyclotomic srg CERTIFIED
```

**Formal verification (Lean 4).** The paper's central *trace-linearity collapse* theorem, the
moment-collapse it generalizes, and the supporting arithmetic (`a ∈ {3,4,5,6}`; `83` not self-conjugate
mod 37) are machine-checked in [`lean/`](lean/) against Mathlib `v4.30.0`, with **no `sorry`** and depending
only on the standard axioms (`propext`, `Classical.choice`, `Quot.sound`):

```bash
cd lean && lake exe cache get && lake build       # builds H668lean.Collapse
```

**Orders 2 and 3 (OPEN).** These reduce to a per-skeleton cyclotomic *lift*. The sound solver
`cyclo_lift_coset.py` reconstructs Paley P(9)/P(25) (validates) but does **not** decide a single
`|H|∈{2,3}` srg(333) instance within 30 min. The gauge-accelerated `order3_sweep.py` / `order2_*` scripts
are **retracted**: their gauge-fixing is unsound under a multiplier constraint (it rejects the real
negation-invariant P(25)); they are kept only to document the pitfall.

**Orbit-matrix-free confirmation (Theorem 1(b), `{12,18,36}`) and evidence runs on Modal:**

```bash
modal run modal_cyclo.py --time-s 7200     # SOUND orbit-matrix-free check (no gauge) -> {12,18,36} UNSAT
modal run modal_gpu_lift.py                # GPU parallel-tempering, unstructured lift (heuristic evidence)
# NOTE: modal_lift.py --mode cyclo and modal_order2*.py use the unsound multiplier gauge-fix -> RETRACTED.
```

**Re-enumerate orbit-matrix skeletons (C++):**

```bash
g++ -O2 -o conference_orbit conference_orbit.cpp && ./conference_orbit
g++ -O2 -o conference_orbit_all conference_orbit_all.cpp && ./conference_orbit_all
```

## Reproducibility note

All small-case validations above were run locally; the full-scale sweeps were run on
[Modal](https://modal.com). Every encoding is validated against ground truth before being
applied to srg(333); see the `validate` entry points and `FINDINGS.md`. We deliberately
do **not** ship a `hadamard_668.csv`: none was found, and none is claimed.

## Citation

```bibtex
@misc{kwik2026h668z37,
  author = {Kwik, Maceo Cardinale},
  title  = {A {$\mathbb{Z}_{37}$} reduction and cyclotomic nonexistence results
            for the conference graph of {H(668)}},
  year   = {2026},
  note   = {arXiv preprint; code at https://github.com/maceoCK/h668-z37},
}
```

A Zenodo DOI will be minted from the first tagged release.

## License

Code is released under the MIT License (see [`LICENSE`](LICENSE)). The paper text is
© 2026 Maceo Cardinale Kwik.
