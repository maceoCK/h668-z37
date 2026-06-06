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
   srg(333,166,82,83) exists for **6 of the 8** proper nontrivial multiplier subgroups
   `|H| ∈ {4,6,9,12,18,36}` — *completely*, independent of the orbit matrix — via a
   coset-compatibility argument on the 9×9 orbit matrix, SAT-certified UNSAT. The two
   finest subgroups (orders 2, 3) are ruled out over 634 enumerated orbit-matrix
   skeletons, and large-scale GPU parallel-tempering finds no unstructured lift.

## Repository layout

```
paper/        h668_z37.tex, h668_z37.pdf   — the paper
code/         all validated machinery + data + reference notes
  conf_core.py              block-circulant SRG machinery (validated on Paley graphs)
  conference_z37.py         the Z₃₇ reduction
  conference_orbit.cpp       C++ DFS that found the worked-example orbit matrix
  conference_orbit_all.cpp   full orbit-matrix enumerator (streams skeletons)
  orbit_matrix_search.py    CP-SAT model for admissible 9×9 orbit matrices
  orbit_compat.py           ★ COMPLETE cyclotomic closure (Theorem 1(a)) — the key result
  cyclo_complete.py         orbit-matrix-free cyclotomic check (Theorem 1(b))
  exact_lift.py             CP-SAT exact lift (convolution encoding + gauge-fixing)
  gpu_lift.py               GPU parallel-tempering lift search (PyTorch)
  emit_h668.py              lift → C(334) → H(668), with end-to-end validation
  multiplier_search.py      Legendre-pair multiplier search
  ninth_root.py, ninth_root_multi.py   the square-divisor cyclotomic filter (Prop 1)
  validate_search.py        ground-truth validation of the multiplier search
  modal_cyclo.py, modal_lift.py, modal_gpu_lift.py   Modal (serverless) drivers for scale
  skeletons.txt (634), skeletons_full.txt (≥1678)    enumerated orbit-matrix skeletons
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

**Reproduce the main nonexistence result (Theorem 1(a) — complete, ~minutes):**

```bash
python orbit_compat.py 300
#   order  9,6,4,12,18,36: INFEASIBLE  -> CLOSED (no such srg, any orbit matrix)
#   order  2,3:            UNKNOWN     -> handled by the skeleton sweep instead
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
