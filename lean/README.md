# Lean 4 formalization

Machine-checked proofs of the paper's central *rigorous* algebraic results (the large
exhaustive certificates are finite searches verified by the C++/Python code in `../code/`;
what is formalized here is the algebra that Lean can certify in-kernel).

`H668lean/Collapse.lean` contains:

- **`moment_collapse`** — if `u² = c·1 + d·u` in a ℚ-algebra then `uⁿ ∈ span{1, u}` for all `n`
  (the spectral-moment hierarchy collapses).
- **`trace_linearity_collapse`** — the paper's Theorem `thm:collapse`. Elements `A, φ, J`
  satisfying the order-2-cyclotomic relations generate a subalgebra contained in the
  ℚ-span of `{1, A, J, φ, φA}`; hence every trace-linear invariant (spectral, Krein,
  scheme, counting) is pinned to five fixed values, and any genuine obstruction must be
  **non-linear** in the {0,1} connection-set indicators.
- **`a_in_three_to_six`** — the Gauss-sum/Galois pinning `a ∈ {3,4,5,6}`.
- **`mult83_not_selfconjugate`** — `83 ≡ 9 (mod 37)` has odd multiplicative order 9 and
  `-1 ∉ ⟨83⟩`, so `83` is not self-conjugate mod 37 (the multiplier divisibility bounds
  are vacuous).

All four are proved with **no `sorry`** and depend only on the standard Mathlib axioms
`propext`, `Classical.choice`, `Quot.sound` (verifiable with `#print axioms`).

## Build

Requires [`elan`](https://github.com/leanprover/elan); the pinned toolchain
(`leanprover/lean4:v4.30.0`) and Mathlib revision (`v4.30.0`) are declared in
`lean-toolchain` and `lakefile.toml`.

```bash
cd lean
lake exe cache get      # fetch prebuilt Mathlib oleans
lake build              # builds H668lean.Collapse
```

To confirm the proofs are axiom-clean:

```lean
import H668lean.Collapse
#print axioms H668.trace_linearity_collapse
-- 'H668.trace_linearity_collapse' depends on axioms: [propext, Classical.choice, Quot.sound]
```
