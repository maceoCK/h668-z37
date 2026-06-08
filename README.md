# A Z₃₇ reduction and cyclotomic nonexistence results for the conference graph of H(668)

Validated machinery and the accompanying paper for an attack on the order-668 Hadamard
matrix problem — the smallest order for which existence is currently open — via a
fixed-point-free automorphism of order 37 of the strongly regular (conference) graph
**srg(333, 166, 82, 83)** ⟺ symmetric conference matrix **C(334)** ⟺ **H(668)**.

> **Status.** H(668) remains open. This repository contributes a reduction, complete
> nonexistence results for cyclotomically-structured Z₃₇-symmetric solutions, and the
> reproducible code behind them. It does **not** claim a construction of H(668).

The paper is in [`paper/h668_z37.pdf`](paper/h668_z37.pdf) (source `paper/h668_z37.tex`).

## Results in one paragraph

1. **Cyclotomic filter (negative).** The natural *square-divisor* refinement of the
   Kotsireas–Koutschan cube-root filter for Legendre pairs of length 333 = 9·37 —
   evaluating the power spectrum at a primitive 9th root of unity — adds **no** pruning
   beyond the mod-3 (cube-root) condition.
2. **Z₃₇ reduction.** Assuming a fixed-point-free order-37 automorphism (which the
   eigenvalue-field obstruction *permits*, unlike the group-developed/Cayley case),
   srg(333,166,82,83) reduces to a 9×9 array of Z₃₇-circulants, with strong regularity
   equivalent to a per-character condition `M_t² + M_t − 83I = 0` over Z[ζ₃₇].
3. **Cyclotomic nonexistence.** No `H`-cyclotomic (Paley-type) Z₃₇-symmetric
   srg(333,166,82,83) exists for **7 of the 8** proper nontrivial multiplier subgroups —
   *every* subgroup of order ≥ 3 — **completely**:
   - `|H| ∈ {4,6,9,12,18,36}`: no coset-compatible orbit matrix exists at all
     (SAT-certified UNSAT), independent of the orbit matrix;
   - `|H| = 3`: order-3-compatible orbit matrices *do* exist — a complete backtracking
     enumeration (no node cap, 1.99×10¹⁰ nodes) finds **exactly 960** — but *none* of the
     960 admits an order-3-invariant lift (all UNSAT), so no order-3-invariant graph exists.

   Only `|H| = 2` (negation-invariant blocks) remains a *sample* rather than a complete
   proof: order 2 imposes no congruence on the orbit matrix, so it is ruled out over the
   634 enumerated skeletons. Large-scale GPU parallel-tempering finds no unstructured lift
   either.
4. **The order-2 residue — a sharper ceiling** (paper §4, Prop. order2). The remaining
   order-2 case is *narrowed and explained*:
   - **Narrowing:** a Gauss-sum/Galois argument pins the θ-multiplicity of each
     residue-character to `a ∈ {3,4,5,6}`; `a ∈ {3,6}` force a constant-18 orbit-matrix
     diagonal, which an exhaustive search (sound `S₉` symmetry break; 1.9×10⁸ nodes,
     corroborated to 9×10¹⁰) proves **does not exist** — so `a ∈ {4,5}`.
   - **Why it resists (two meta-theorems):** the group-ring RHS lives entirely on the
     principal character (no character/Gauss-sum/multiplier obstruction can fire; `83`
     isn't self-conjugate mod 37), and the negation involution collapses *every*
     trace-linear functional (spectral, Krein, scheme, counting) to five fixed values.
     Any genuine obstruction must be **non-linear** in the {0,1} indicators — i.e. the
     lift feasibility itself, which is the order-668 open problem. (Order 2 stays open.)

## Repository layout

```
paper/        h668_z37.tex, h668_z37.pdf   — the paper
code/         all validated machinery + data + reference notes
  conf_core.py              block-circulant SRG machinery (validated on Paley graphs)
  conference_z37.py         the Z₃₇ reduction
  conference_orbit.cpp       C++ DFS that found the worked-example orbit matrix
  conference_orbit_all.cpp   full orbit-matrix enumerator (streams skeletons)
  orbit3_enum.cpp           ★ order-3-compatible enumerator → exactly 960 skeletons (Thm 1(b))
  order3_sweep.py           ★ exhaustive order-3-invariant lift sweep over the 960 (all UNSAT)
  order2_spectral.py        order-2 (negation-invariant) lift with N_g/a∈{4,5} pruning
  uniform18_enum.cpp        ★ proves no constant-18 orbit matrix exists → a∈{3,6} eliminated
  uniform18_maxbreak.cpp    independent re-check (different sound symmetry break)
  order2_combinatorial.py   the involution-collapse / flat-RHS meta-theorem checks
  modal_order2.py           Modal driver for the a∈{4,5} order-2 frontier search
  orbit_matrix_search.py    CP-SAT model for admissible 9×9 orbit matrices
  orbit_compat.py           ★ COMPLETE cyclotomic closure (Theorem 1(a)) — the key result
  cyclo_complete.py         orbit-matrix-free cyclotomic check (Theorem 1(c))
  exact_lift.py             CP-SAT exact lift (convolution encoding + gauge-fixing)
  gpu_lift.py               GPU parallel-tempering lift search (PyTorch)
  emit_h668.py              lift → C(334) → H(668), with end-to-end validation
  multiplier_search.py      Legendre-pair multiplier search
  ninth_root.py, ninth_root_multi.py   the square-divisor cyclotomic filter (Prop 1)
  validate_search.py        ground-truth validation of the multiplier search
  modal_cyclo.py, modal_lift.py, modal_gpu_lift.py   Modal (serverless) drivers for scale
  skeletons.txt (634), skeletons_full.txt (≥1678)    enumerated orbit-matrix skeletons
  orbit3_compatible.txt (960)   complete set of order-3-compatible orbit matrices (Thm 1(b))
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

**Reproduce the complete nonexistence results (Theorem 1, orders ≥ 3):**

```bash
# (a) orders {4,6,9,12,18,36}: no coset-compatible orbit matrix exists at all (~minutes)
python orbit_compat.py 300
#   order  9,6,4,12,18,36: INFEASIBLE  -> CLOSED (no such srg, any orbit matrix)
#   order  2,3:            UNKNOWN     -> orbit-matrix argument does not apply (see below)

# (b) order 3: enumerate the COMPLETE order-3-compatible set, then sweep lifts (~30 min total)
clang++ -O3 -o orbit3_enum orbit3_enum.cpp && ./orbit3_enum   # -> exactly 960 skeletons
python order3_sweep.py 120 8                                  # -> all 960 UNSAT => order 3 CLOSED
```

**Orbit-matrix-free confirmation (Theorem 1(b)) and the full skeleton sweep run at scale on Modal:**

```bash
modal run modal_cyclo.py --time-s 7200     # orbit-matrix-free per subgroup
modal run modal_lift.py  --mode cyclo      # per-skeleton cyclotomic sweep over 634 skeletons
modal run modal_gpu_lift.py                # GPU parallel-tempering, unstructured lift
```

**Re-enumerate orbit-matrix skeletons (C++):**

```bash
g++ -O2 -o conference_orbit conference_orbit.cpp && ./conference_orbit
g++ -O2 -o conference_orbit_all conference_orbit_all.cpp && ./conference_orbit_all
```

## Reproducibility note

All small-case validations above were run locally; the full-scale sweeps were run on
[Modal](https://modal.com). Every encoding is validated against ground truth before being
applied to srg(333) — see the `validate` entry points and `FINDINGS.md`. We deliberately
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
