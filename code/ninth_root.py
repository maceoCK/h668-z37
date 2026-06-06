"""
Is the ninth-root (square-divisor) spectral filter STRICTLY stronger than the
mod-3 (cube-root) filter of Kotsireas-Koutschan?

Test on the smallest square-divisor length ell = 45 = 9*5 (m = ell/9 = 5),
where the 9-compression has entries in {-5,-3,-1,1,3,5} -> fully enumerable.

Compression identity: c_r = sum_{j == r (mod 9)} u_j, so for ell=45 each c_r is a
sum of 5 values in {+-1}, hence odd in [-5,5].
Compressed LP condition (length 9): PAF_c(k)+PAF_c'(k) = -2*m = -10 for k!=0,
and = 2*ell - 2*m + 2 = 82 for k=0.

mod-3 admissible states (gamma,gamma'): gamma = 3-compression (length 3),
gamma_j odd in [-15,15], sum = +-1, with Sigma gamma^2 + Sigma gamma'^2 = 4*(ell/3)+2 = 62.
gamma is also the 3-compression of c: gamma_j = c_j + c_{j+3} + c_{j+6}.

QUESTION: does every mod-3 admissible (gamma,gamma') lift to a valid 9-compressed LP
(c,c')? If some do NOT, the ninth-root filter is STRICTLY stronger.
"""
import itertools
from collections import defaultdict

ELL = 45
M = ELL // 9          # 5
VALS = [v for v in range(-M, M+1) if v % 2 != 0]   # {-5,-3,-1,1,3,5}

def paf9(c):
    return tuple(sum(c[r]*c[(r+k) % 9] for r in range(9)) for k in range(5))  # lags 0..4

# target pair-sum profile (lag 0..4)
TARGET = (2*ELL - 2*M + 2,) + (-2*M,)*4          # (82,-10,-10,-10,-10)

def gamma_of(c):
    return (c[0]+c[3]+c[6], c[1]+c[4]+c[7], c[2]+c[5]+c[8])

# ---- enumerate single length-9 compressions c with sum = +-1, bucket by PAF profile ----
# recursive partial-sum pruning to keep it fast
profile_to_c = defaultdict(list)   # paf-profile -> list of c (as tuples)
c = [0]*9
def rec(pos, ssum, ssq):
    if ssq > TARGET[0]:           # sum of squares can't exceed PAF_c(0) bound for a valid pair (<=82-9)
        return
    if pos == 9:
        if ssum in (1, -1):
            profile_to_c[paf9(c)].append(tuple(c))
        return
    # prune: remaining positions contribute at least +-1 each
    for v in VALS:
        c[pos] = v
        rec(pos+1, ssum+v, ssq+v*v)
rec(0, 0, 0)

n_c = sum(len(v) for v in profile_to_c.values())
print(f"ell={ELL}: # length-9 compressions c (sum=+-1): {n_c:,}")

# ---- count valid (c,c') pairs and collect realized mod-3 states ----
realized_mod3 = set()
n_pairs = 0
profiles = list(profile_to_c.keys())
pset = set(profiles)
for P in profiles:
    Pc = (TARGET[0]-P[0],) + tuple(-2*M - P[k] for k in range(1,5))
    if Pc in pset:
        cs = profile_to_c[P]; cps = profile_to_c[Pc]
        n_pairs += len(cs)*len(cps)
        for cc in cs:
            g = gamma_of(cc)
            for ccp in cps:
                realized_mod3.add((g, gamma_of(ccp)))
print(f"ell={ELL}: # valid 9-compressed LP pairs (c,c'): {n_pairs:,}")

# ---- enumerate ALL mod-3 admissible (gamma,gamma') for ell=45 ----
modd = ELL//3   # 15
gvals = [v for v in range(-modd, modd+1) if v % 2 != 0]
def gammas_by_sq():
    by = defaultdict(list)
    for a in gvals:
        for b in gvals:
            s=a+b
            for cc in (1-s, -1-s):
                if cc % 2 != 0 and -modd <= cc <= modd:
                    by[a*a+b*b+cc*cc].append((a,b,cc))
    return by
gby = gammas_by_sq()
TGT3 = 4*modd + 2   # 62
mod3_states = set()
for s, lst in gby.items():
    s2 = TGT3 - s
    if s2 in gby:
        for ga in lst:
            for gb in gby[s2]:
                mod3_states.add((ga, gb))
print(f"ell={ELL}: # mod-3 admissible states (gamma,gamma'): {len(mod3_states):,}")

# ---- THE COMPARISON ----
not_lifting = mod3_states - realized_mod3
print(f"ell={ELL}: # mod-3 states REALIZED by some 9-compression: {len(realized_mod3 & mod3_states):,}")
print(f"ell={ELL}: # mod-3 states that DO NOT lift to any 9-compression: {len(not_lifting):,}")
print(f"=> ninth-root filter strictly stronger than mod-3?  {len(not_lifting) > 0}")
if not_lifting:
    ex = sorted(not_lifting)[:3]
    print("   examples of non-lifting mod-3 states (gamma, gamma'):")
    for e in ex: print("    ", e)
# sanity: realized should be subset of admissible
print(f"   (sanity) realized subset of admissible? {realized_mod3 <= mod3_states}")
