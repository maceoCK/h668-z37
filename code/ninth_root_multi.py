"""
Test at MULTIPLE square-divisor lengths ell = 9*odd whether the ninth-root condition
prunes the mod-3 spectrum: i.e. does every mod-3 admissible (gamma,gamma') lift to a
valid length-9 compressed LP (c,c')?  Generalizes ninth_root.py to ell in {27,45,63}.
"""
import sys
from collections import defaultdict

def analyze(ELL):
    M = ELL // 9
    VALS = [v for v in range(-M, M+1) if v % 2 != 0]
    def paf9(c):
        return tuple(sum(c[r]*c[(r+k) % 9] for r in range(9)) for k in range(5))
    TARGET = (2*ELL - 2*M + 2,) + (-2*M,)*4
    def gamma_of(c):
        return (c[0]+c[3]+c[6], c[1]+c[4]+c[7], c[2]+c[5]+c[8])
    profile_to_c = defaultdict(list)
    c = [0]*9
    cap = TARGET[0]
    def rec(pos, ssum, ssq):
        if ssq > cap: return
        if pos == 9:
            if ssum in (1, -1):
                profile_to_c[paf9(c)].append(tuple(c))
            return
        for v in VALS:
            c[pos]=v; rec(pos+1, ssum+v, ssq+v*v)
    rec(0,0,0)
    n_c = sum(len(v) for v in profile_to_c.values())
    realized = set(); n_pairs=0
    pset=set(profile_to_c)
    for P in profile_to_c:
        Pc=(TARGET[0]-P[0],)+tuple(-2*M-P[k] for k in range(1,5))
        if Pc in pset:
            cs=profile_to_c[P]; cps=profile_to_c[Pc]
            n_pairs+=len(cs)*len(cps)
            gs=[gamma_of(x) for x in cs]; gps=[gamma_of(x) for x in cps]
            for g in gs:
                for gp in gps:
                    realized.add((g,gp))
    # mod-3 admissible
    modd=ELL//3
    gvals=[v for v in range(-modd,modd+1) if v%2!=0]
    gby=defaultdict(list)
    for a in gvals:
        for b in gvals:
            s=a+b
            for cc in (1-s,-1-s):
                if cc%2!=0 and -modd<=cc<=modd:
                    gby[a*a+b*b+cc*cc].append((a,b,cc))
    TGT3=4*modd+2
    mod3=set()
    for s,lst in gby.items():
        s2=TGT3-s
        if s2 in gby:
            for ga in lst:
                for gb in gby[s2]:
                    mod3.add((ga,gb))
    nonlift=mod3-realized
    return dict(ELL=ELL, n_c=n_c, n_pairs=n_pairs, n_mod3=len(mod3),
                realized=len(realized & mod3), nonlift=len(nonlift),
                stronger=len(nonlift)>0)

if __name__=="__main__":
    lens=[int(a) for a in sys.argv[1:]] or [27,45,63]
    for ELL in lens:
        r=analyze(ELL)
        print(f"ell={ELL:3d}=9*{ELL//9}: #9-comp c={r['n_c']:>9,}  #valid pairs={r['n_pairs']:>12,}  "
              f"#mod3 states={r['n_mod3']:>5}  realized={r['realized']:>5}  NOT-lifting={r['nonlift']:>3}  "
              f"ninth-root stronger? {r['stronger']}")
