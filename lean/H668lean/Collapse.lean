import Mathlib

set_option linter.style.header false
set_option linter.style.longLine false
set_option linter.unnecessarySimpa false
set_option linter.unusedVariables false

open scoped Pointwise

/-!
# Lean verification for "A Z₃₇ orbit-matrix reduction and cyclotomic obstructions for H(668)"

This file machine-checks the paper's core *rigorous* results (the giant exhaustive certificates are
finite searches verified externally; the algebra below is what Lean can certify in-kernel):

* `moment_collapse` — the spectral-moment pinning;
* `trace_linearity_collapse` — the central Theorem `thm:collapse`: in any ℚ-algebra, elements `A, φ, J`
  satisfying the order-2-cyclotomic relations generate a subalgebra contained in the ℚ-span of
  `{1, A, J, φ, φA}`, so every trace-linear invariant is determined by five fixed values;
* small arithmetic facts used in the paper (`a ∈ {3,4,5,6}`; `ord₃₇(83)=9` and `83` not self-conjugate).
-/

namespace H668

variable {𝔸 : Type*} [Ring 𝔸] [Algebra ℚ 𝔸]

/-- **Moment collapse.** If `u² = c•1 + d•u`, then `uⁿ ∈ span{1,u}` for all `n`; hence every trace
`tr(uⁿ)` is determined by `tr 1` and `tr u` (the spectral-moment hierarchy collapses). -/
theorem moment_collapse (u : 𝔸) (c d : ℚ) (hu : u * u = c • (1 : 𝔸) + d • u) (n : ℕ) :
    u ^ n ∈ Submodule.span ℚ ({1, u} : Set 𝔸) := by
  set p := Submodule.span ℚ ({1, u} : Set 𝔸) with hp
  have h1 : (1 : 𝔸) ∈ p := Submodule.subset_span (by simp)
  have hu' : u ∈ p := Submodule.subset_span (by simp)
  have hclosed : ∀ x ∈ p, u * x ∈ p := by
    intro x hx
    induction hx using Submodule.span_induction with
    | mem y hy =>
        rcases hy with rfl | rfl
        · simpa using hu'
        · rw [hu]
          exact Submodule.add_mem _ (Submodule.smul_mem _ _ h1) (Submodule.smul_mem _ _ hu')
    | zero => simpa using Submodule.zero_mem p
    | add x y _ _ hx hy => rw [mul_add]; exact Submodule.add_mem _ hx hy
    | smul a x _ hx => rw [mul_smul_comm]; exact Submodule.smul_mem _ _ hx
  induction n with
  | zero => simpa using h1
  | succ k ih => rw [pow_succ']; exact hclosed _ ih

/-- **Trace-linearity collapse (paper Theorem `thm:collapse`).** Elements `A, φ, J` satisfying the
order-2-cyclotomic relations generate a subalgebra contained in the ℚ-span of `{1, A, J, φ, φA}`. Hence
every trace-linear invariant of the graph is determined by the five fixed traces, and none can obstruct
existence; any obstruction must be non-linear in the `{0,1}` connection-set indicators. -/
theorem trace_linearity_collapse (A φ J : 𝔸)
    (hAA : A * A = (83 : ℚ) • (1 : 𝔸) - A + (83 : ℚ) • J)
    (hφφ : φ * φ = 1) (hφA : φ * A = A * φ)
    (hφJ : φ * J = J) (hJφ : J * φ = J)
    (hJJ : J * J = (333 : ℚ) • J) (hAJ : A * J = (166 : ℚ) • J) (hJA : J * A = (166 : ℚ) • J) :
    (Algebra.adjoin ℚ ({A, φ, J} : Set 𝔸)).toSubmodule ≤
      Submodule.span ℚ ({1, A, J, φ, φ * A} : Set 𝔸) := by
  set p := Submodule.span ℚ ({1, A, J, φ, φ * A} : Set 𝔸) with hp
  have m1 : (1 : 𝔸) ∈ p := Submodule.subset_span (by simp)
  have mA : A ∈ p := Submodule.subset_span (by simp)
  have mJ : J ∈ p := Submodule.subset_span (by simp)
  have mφ : φ ∈ p := Submodule.subset_span (by simp)
  have mφA : φ * A ∈ p := Submodule.subset_span (by simp)
  -- the nonobvious product reduces into the five-dimensional span
  have eAφA : A * (φ * A) = (83 : ℚ) • φ - φ * A + (83 : ℚ) • J := by
    rw [← mul_assoc, ← hφA, mul_assoc, hAA, mul_add, mul_sub, mul_smul_comm, mul_smul_comm, mul_one, hφJ]
  -- left multiplication by each generator preserves `p` (robust `span_induction`)
  have hA : ∀ x ∈ p, A * x ∈ p := by
    intro x hx
    induction hx using Submodule.span_induction with
    | mem g hg =>
        rcases hg with rfl | rfl | rfl | rfl | rfl
        · rw [mul_one]; exact mA
        · rw [hAA]
          exact Submodule.add_mem _ (Submodule.sub_mem _ (Submodule.smul_mem _ _ m1) mA)
            (Submodule.smul_mem _ _ mJ)
        · rw [hAJ]; exact Submodule.smul_mem _ _ mJ
        · rw [← hφA]; exact mφA
        · rw [eAφA]
          exact Submodule.add_mem _ (Submodule.sub_mem _ (Submodule.smul_mem _ _ mφ) mφA)
            (Submodule.smul_mem _ _ mJ)
    | zero => simpa using Submodule.zero_mem p
    | add a b _ _ ha hb => rw [mul_add]; exact Submodule.add_mem _ ha hb
    | smul c a _ ha => rw [mul_smul_comm]; exact Submodule.smul_mem _ _ ha
  have hφc : ∀ x ∈ p, φ * x ∈ p := by
    intro x hx
    induction hx using Submodule.span_induction with
    | mem g hg =>
        rcases hg with rfl | rfl | rfl | rfl | rfl
        · rw [mul_one]; exact mφ
        · exact mφA
        · rw [hφJ]; exact mJ
        · rw [hφφ]; exact m1
        · rw [← mul_assoc, hφφ, one_mul]; exact mA
    | zero => simpa using Submodule.zero_mem p
    | add a b _ _ ha hb => rw [mul_add]; exact Submodule.add_mem _ ha hb
    | smul c a _ ha => rw [mul_smul_comm]; exact Submodule.smul_mem _ _ ha
  have hJc : ∀ x ∈ p, J * x ∈ p := by
    intro x hx
    induction hx using Submodule.span_induction with
    | mem g hg =>
        rcases hg with rfl | rfl | rfl | rfl | rfl
        · rw [mul_one]; exact mJ
        · rw [hJA]; exact Submodule.smul_mem _ _ mJ
        · rw [hJJ]; exact Submodule.smul_mem _ _ mJ
        · rw [hJφ]; exact mJ
        · rw [← mul_assoc, hJφ, hJA]; exact Submodule.smul_mem _ _ mJ
    | zero => simpa using Submodule.zero_mem p
    | add a b _ _ ha hb => rw [mul_add]; exact Submodule.add_mem _ ha hb
    | smul c a _ ha => rw [mul_smul_comm]; exact Submodule.smul_mem _ _ ha
  -- multiplication closure: for fixed `y ∈ p`, `{x | x*y ∈ p}` is a submodule containing the generators
  have hmul : ∀ (x y : 𝔸), x ∈ p → y ∈ p → x * y ∈ p := by
    intro x y hx hy
    induction hx using Submodule.span_induction with
    | mem g hg =>
        rcases hg with rfl | rfl | rfl | rfl | rfl
        · rw [one_mul]; exact hy
        · exact hA y hy
        · exact hJc y hy
        · exact hφc y hy
        · rw [mul_assoc]; exact hφc _ (hA y hy)
    | zero => rw [zero_mul]; exact Submodule.zero_mem p
    | add a b _ _ ha hb => rw [add_mul]; exact Submodule.add_mem _ ha hb
    | smul c a _ ha => rw [smul_mul_assoc]; exact Submodule.smul_mem _ _ ha
  have hsub : Algebra.adjoin ℚ ({A, φ, J} : Set 𝔸) ≤ p.toSubalgebra m1 hmul := by
    apply Algebra.adjoin_le
    intro x hx
    simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
    rcases hx with rfl | rfl | rfl
    exacts [mA, mφ, mJ]
  intro x hx
  exact hsub hx

/-- The order-2 residue-character θ-multiplicity `a` satisfies `a ∈ {3,4,5,6}`: the Gauss-sum/Galois
identity pins the diagonal pair-counts to `3a-9` (residues) and `18-3a` (non-residues), each in `[0,9]`.
(The further elimination of `a ∈ {3,6}` is the external constant-diagonal enumeration certificate.) -/
theorem a_in_three_to_six (a : ℤ)
    (hQR : 0 ≤ 3 * a - 9 ∧ 3 * a - 9 ≤ 9) (hN : 0 ≤ 18 - 3 * a ∧ 18 - 3 * a ≤ 9) :
    3 ≤ a ∧ a ≤ 6 := by omega

/-- `83 ≡ 9 (mod 37)` has odd multiplicative order `9`, and `-1 = 36 ∉ ⟨83⟩`; hence `83` is not
self-conjugate mod `37` and the Turyn–Schmidt/multiplier divisibility bounds are vacuous. -/
theorem mult83_not_selfconjugate :
    (83 : ZMod 37) = 9 ∧ (9 : ZMod 37) ^ 9 = 1 ∧
    (∀ k, 1 ≤ k → k ≤ 8 → (9 : ZMod 37) ^ k ≠ 1) ∧
    (∀ k, (9 : ZMod 37) ^ k ≠ 36) := by
  refine ⟨by decide, by decide, by decide, ?_⟩
  intro k
  have : (9 : ZMod 37) ^ k = (9 : ZMod 37) ^ (k % 9) := by
    conv_lhs => rw [← Nat.div_add_mod k 9, pow_add, pow_mul, (by decide : (9:ZMod 37)^9 = 1),
      one_pow, one_mul]
  rw [this]; have : k % 9 < 9 := Nat.mod_lt _ (by norm_num)
  interval_cases (k % 9) <;> decide

end H668
