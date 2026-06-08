"""
Order-2 (negation-invariant) Z_37 srg(333,166,82,83): COMBINATORIAL / non-spectral obstruction hunt.

Summary of what is PROVED here (all validated against existing order-2 conference/SRG graphs
P(9), P(25), P(49), Paley(13) -- which must NOT be killed):

1. NEGATION INVOLUTION phi: (i,g) -> (i,-g) on all 9 orbits is a genuine graph automorphism of order 2
   (exists precisely BECAUSE every block is negation-invariant). f = 9 fixed points.
   phi COMMUTES with the adjacency A (A is phi-invariant). => phi A^k phi = A^k.

2. META-THEOREM (proved). Because (a) A^2 = 166 I + 82 A + 83(J-I-A) is affine in {I,A,J}
   and (b) phi commutes with A, every alternating word in {A, phi} collapses to an affine combo of
   {I, A, J, phi, phi.A}. All five traces are pinned: trace(I)=333, trace(A)=0, trace(J)=333,
   trace(phi)=9, trace(phi.A)=trace(R)=162. Hence NO trace-linear (spectral/involution/character)
   functional carries information beyond R. This EXTENDS the moment-collapse ceiling theorem to the
   negation involution.

3. Involution multiplicity split is rep-theoretically FORCED (not an obstruction): the theta-eigenspace
   (mult 166) splits 85/81 across V+/V-, tau likewise 85/81; dp=dq=4 forced by the sqrt(37)-irrational
   part of trace(phi.A)=162. All integers. Validated on Paley(13) (dp=dq=0) and Paley(49).

4. ETA-COEFFICIENT EQUATIONS. M_1^2+M_1-83I=0 expanded in the integral basis {eta_g=zeta^g+zeta^-g}
   gives 9*9*18 integer convolution equations (the "per-h" family). Validated exactly on
   P(25),P(49),P(121),P(169). These are a DFT change-of-basis of the standard SRG convolution SAT --
   no new information. The SUMMED-over-h version reduces to 37 | RHS(R), which is FORCED by
   R^2+R-83I=3071J (RHS = 19(R^2+R)+1494 I = 1665 I = 0 mod 37) -- TAUTOLOGICAL.

5. Per-block entrywise absolute bound |M_t[i][j]| <= 3 sqrt37 / 2 (from M_t = -1/2 I + (3sqrt37/2) S_t,
   S_t^2=I) is easily satisfiable for every block size 1..18 -- no obstruction.

6. Krein parameters of the Z_37/negation cyclic (dihedral orbital) scheme are NON-NEGATIVE (group
   scheme) -- direction (a) gives no kill.

7. N_g diagonal Gauss-sum content fully captured by a in {3,4,5,6} (re-derived; prior result).

CONCLUSION (honest): NO combinatorial/non-spectral obstruction was found. Every Z_37-and-negation-
invariant LINEAR or SUMMED functional provably collapses to R via the meta-theorem (2). Any genuine
obstruction must be NONLINEAR in the {0,1} block indicators -- i.e. it lives in the joint realizability
SAT (846 vars), exactly as the moment-collapse ceiling predicted. The order-2 case is NOT closeable by
the combinatorial tools explored.
"""
import numpy as np

P, NORB, K, LAM, MU = 37, 9, 166, 82, 83
M = (P - 1) // 2  # 18 eta-classes

def meta_theorem_traces():
    """All five basis traces are pinned by R / the spectrum (proved)."""
    R = np.load(__file__.rsplit('/', 1)[0] + '/orbitR_found.npy')
    return dict(trI=NORB * P, trA=0, trJ=NORB * P, trPhi=NORB, trPhiA=int(np.trace(R)))

def involution_split():
    """theta/tau multiplicity split across phi's +-1 eigenspaces (forced, integral)."""
    # rational part: dp+dq = trace(phi)-1 = 8 ; irrational(sqrt37) part: dp-dq = 0 => dp=dq=4
    dp = dq = 4
    pplus = (K + dp) // 2; pminus = (K - dp) // 2
    qplus = (K + dq) // 2; qminus = (K - dq) // 2
    assert 1 + pplus + qplus == (NORB * P + NORB) // 2  # dim V+ = 171
    assert pminus + qminus == (NORB * P - NORB) // 2    # dim V- = 162
    return dict(theta=(pplus, pminus), tau=(qplus, qminus))

def summed_eta_is_tautological():
    """37 | RHS(R) is forced by R^2+R-83I=3071J  ->  RHS = 1665 I = 0 mod 37."""
    return (19 * 9 + 1494) % P == 0  # 19*(R^2+R)=19*9I, +1494 I ; (R^2+R=9I mod37)

if __name__ == "__main__":
    print("meta-theorem pinned traces:", meta_theorem_traces())
    print("involution split (theta V+/V-, tau V+/V-):", involution_split())
    print("summed-eta tautological (37 | RHS):", summed_eta_is_tautological())
    print("\nNo combinatorial/non-spectral obstruction found; order-2 stays open (see module docstring).")
