"""
sc_floor_analysis.py -- "is the SA floor real?" analysis for self-complementary classes 221/422.

Reconstructs BESTSTATE dumps from deepdive_lp333/sc_residuals.json against the archived
templates, verifies objectives/residuals exactly, then hunts for invariant lower bounds:
  1. exact verification of the archived fingerprints (obj + per-equation residuals)
  2. structural identities: per-(a,b) row sums of residuals; (a,a,0) auto-zero; g <-> -g
     evenness on diagonal blocks; F always even
  3. sampling-based affine-invariant hunt over Q, F2, F3 on the residual map
     (any invariant functional with NONZERO constant ==> min F > 0 theorem)
"""
import json, sys, itertools, random
import numpy as np

CODE = "/Users/maceocardinalekwik/git/math/release/h668-z37/code"
RESID = "/Users/maceocardinalekwik/git/math/deepdive_lp333/sc_residuals.json"

# ---------------- template ----------------
class Tpl:
    def __init__(self, path):
        toks = open(path).read().split()
        it = iter(toks)
        self.P = int(next(it)); self.NORB = int(next(it)); self.NREP = int(next(it))
        self.K = int(next(it)); self.LAM = int(next(it)); self.MU = int(next(it)); self.C = int(next(it))
        self.reps = [None]*self.NREP
        for _ in range(self.NREP):
            tag = next(it); idx = int(next(it))
            i,j,typ,card = int(next(it)),int(next(it)),int(next(it)),int(next(it))
            self.reps[idx] = (i,j,typ,card)
        P,N = self.P, self.NORB
        self.litRep = np.full((N,N,P), -2, dtype=np.int64)
        self.litE   = np.zeros((N,N,P), dtype=np.int64)
        self.litNeg = np.zeros((N,N,P), dtype=np.int64)
        for _ in range(N*N*P):
            tag = next(it)
            a,b,d,ri,e,ng = (int(next(it)) for _ in range(6))
            self.litRep[a,b,d]=ri; self.litE[a,b,d]=e; self.litNeg[a,b,d]=ng
        assert not (self.litRep==-2).any()
        self.isRepPair = np.zeros((N,N),dtype=bool)
        for (i,j,typ,card) in self.reps: self.isRepPair[i,j]=True
        # c-orbits {x, cx, -x, -cx} for type-2 sampling
        seen=[False]*P; self.corbs=[]
        for x in range(1,P):
            if seen[x]: continue
            b=(self.C*x)%P; na=P-x; nb=P-b
            for z in (x,b,na,nb): seen[z]=True
            self.corbs.append((x,b,na,nb))

    def expand(self, repBit):
        """repBit: (NREP,P) 0/1 -> mem (N,N,P)"""
        P,N = self.P,self.NORB
        mem = np.zeros((N,N,P), dtype=np.int64)
        r = self.litRep
        pos = r>=0
        mem[pos] = repBit[r[pos], self.litE[pos]] ^ self.litNeg[pos]
        mem[~pos] = self.litNeg[~pos]
        return mem

    def conv(self, mem):
        P,N = self.P,self.NORB
        # conv[a,b,g] = sum_l sum_d mem[a,l,d] mem[l,b,g-d]  (cyclic) -- use FFT, round
        F = np.fft.rfft(mem, axis=2)
        # conv_hat[a,b] = sum_l F[a,l]*F[l,b]
        C = np.einsum('alk,lbk->abk', F, F)
        out = np.fft.irfft(C, n=P, axis=2)
        return np.rint(out).astype(np.int64)

    def residuals(self, repBit):
        """returns res (N,N,P) for ALL pairs; objective restricted to rep pairs a<=b."""
        P,N = self.P,self.NORB
        mem = self.expand(repBit)
        cv = self.conv(mem)
        rhs = np.full((N,N,P), self.MU, dtype=np.int64)
        for a in range(N): rhs[a,a,0] += self.K - self.MU
        res = cv - (self.LAM-self.MU)*mem - rhs
        return res

    def obj_of(self, res, all_pairs=False):
        s = 0
        for a in range(self.NORB):
            for b in range(a, self.NORB):
                if not all_pairs and not self.isRepPair[a,b]: continue
                s += int((res[a,b]**2).sum())
        return s

    def rep_eq_list(self):
        eqs=[]
        for a in range(self.NORB):
            for b in range(a,self.NORB):
                if self.isRepPair[a,b]:
                    for g in range(self.P): eqs.append((a,b,g))
        return eqs

    def random_state(self, rng):
        P = self.P
        repBit = np.zeros((self.NREP,P), dtype=np.int64)
        for r,(i,j,typ,card) in enumerate(self.reps):
            if typ==0:
                el = rng.sample(range(P), card)
                repBit[r,el]=1
            elif typ==1:
                pr = rng.sample(range(1,P//2+1), card//2)
                for d in pr: repBit[r,d]=1; repBit[r,P-d]=1
            else:
                for ob in self.corbs:
                    if rng.random()<0.5: repBit[r,ob[0]]=1; repBit[r,ob[2]]=1
                    else:               repBit[r,ob[1]]=1; repBit[r,ob[3]]=1
        return repBit

def parse_dump(dump, tpl):
    repBit = np.zeros((tpl.NREP, tpl.P), dtype=np.int64)
    res_expected = {}
    obj = None
    for ln in dump.splitlines():
        t = ln.split()
        if not t: continue
        if t[0]=="BESTSTATE": obj = int(t[1].split("=")[1])
        elif t[0]=="REPBITS":
            r = int(t[1]); els=[int(x) for x in t[3:]]
            repBit[r,els]=1
        elif t[0]=="RES":
            a,b,g,v = int(t[1]),int(t[2]),int(t[3]),int(t[4])
            res_expected[(a,b,g)]=v
    return obj, repBit, res_expected

def main():
    data = json.load(open(RESID))
    print("== STEP 1: exact reconstruction of all 12 archived BESTSTATE dumps ==")
    tpls = {}
    all_ok = True
    for ent in data:
        tname = ent["tpl"]
        if tname not in tpls: tpls[tname] = Tpl(f"{CODE}/{tname}")
        tpl = tpls[tname]
        obj_claim, repBit, res_exp = parse_dump(ent["dump"], tpl)
        # cardinality check
        for r,(i,j,typ,card) in enumerate(tpl.reps):
            assert repBit[r].sum()==card, (tname, r, repBit[r].sum(), card)
        res = tpl.residuals(repBit)
        obj = tpl.obj_of(res)
        # compare nonzero residuals on rep pairs
        mism = 0
        for (a,b,g) in tpl.rep_eq_list():
            v = int(res[a,b,g])
            ve = res_exp.get((a,b,g), 0)
            if v != ve: mism += 1
        ok = (obj==obj_claim and mism==0)
        all_ok &= ok
        print(f"  {tname} seed={ent['seed']}: claimed obj={obj_claim} recomputed={obj} "
              f"res-mismatches={mism} -> {'OK' if ok else 'FAIL'}")
    print("ALL DUMPS VERIFIED" if all_ok else "VERIFICATION FAILED", flush=True)

    print("\n== STEP 2: structural identities (verified on dumps + random states) ==")
    tpl = tpls["sc_tpl_221_0.txt"]
    rng = random.Random(12345)
    states = [parse_dump(e["dump"], tpls[e["tpl"]])[1] for e in data if e["tpl"]=="sc_tpl_221_0.txt"]
    states += [tpl.random_state(rng) for _ in range(50)]
    rowsum_max = 0; diag0_max = 0; even_max = 0; objs = []
    for st in states:
        res = tpl.residuals(st)
        for a in range(tpl.NORB):
            for b in range(a,tpl.NORB):
                if not tpl.isRepPair[a,b]: continue
                rowsum_max = max(rowsum_max, abs(int(res[a,b].sum())))
                if a==b:
                    diag0_max = max(diag0_max, abs(int(res[a,a,0])))
                    even_max = max(even_max, int(np.abs(res[a,a,1:]-res[a,a,1:][::-1]).max()))
        objs.append(tpl.obj_of(res))
    print(f"  (I1) per-(a,b) sum_g residual == 0 : max |sum| over {len(states)} states = {rowsum_max}")
    print(f"  (I2) residual(a,a,0) == 0 always   : max |res(a,a,0)| = {diag0_max}")
    print(f"  (I3) diag blocks even in g         : max |res(a,a,g)-res(a,a,-g)| = {even_max}")
    print(f"  (I4) objective always even         : objs mod 2 = {sorted(set(o%2 for o in objs))}")
    print(f"       random-state objective range: [{min(objs)}, {max(objs)}]")

    print("\n== STEP 3: affine-invariant hunt on the residual map (sampling) ==", flush=True)
    # residual vector over rep equations, N samples; find functionals constant across samples
    for tname in ["sc_tpl_221_0.txt", "sc_tpl_422_2.txt"]:
        tpl = tpls[tname]
        eqs = tpl.rep_eq_list()
        NS = 1400
        rng = random.Random(999)
        R = np.zeros((NS, len(eqs)), dtype=np.int64)
        for s in range(NS):
            res = tpl.residuals(tpl.random_state(rng))
            R[s] = [res[a,b,g] for (a,b,g) in eqs]
        D = R[1:] - R[0]                       # differences: invariants = null space (left) of D^T
        # over Q: rank of D
        rank = np.linalg.matrix_rank(D.astype(float))
        ninv = len(eqs) - rank
        print(f"  [{tname}] eqs={len(eqs)} samples={NS} rank(diff)={rank} -> "
              f"#affine invariants over Q = {ninv}")
        # identify the invariant space: null space of D (column null space)
        # via SVD on floats then verify integrally
        u,sv,vt = np.linalg.svd(D.astype(float), full_matrices=True)
        null = vt[rank:]                      # (ninv, neq)
        # check constants: lambda . r for sample 0
        consts = null @ R[0].astype(float)
        nz = np.abs(consts) > 1e-6
        print(f"    invariant functionals with NONZERO constant: {int(nz.sum())} "
              f"(max |const|={np.abs(consts).max():.3e})")
        # mod 2
        R2 = (R % 2).astype(np.uint8)
        D2 = (R2[1:] ^ R2[0])
        r2 = gf2_rank(D2.copy())
        ninv2 = len(eqs) - r2
        # constants of F2-invariants: need null space basis mod 2
        nsb = gf2_nullspace(D2)
        c2 = (nsb @ R2[0]) % 2
        print(f"    mod 2: rank={r2} -> #F2 invariants = {ninv2}; with constant 1: {int(c2.sum())}")
        if int(c2.sum())>0:
            print("    *** F2 INVARIANT WITH CONSTANT 1 FOUND -- POSSIBLE OBSTRUCTION ***")
            np.save(f"/tmp/f2_invariants_{tname}.npy", nsb[c2==1])
        # mod 3
        R3 = (R % 3).astype(np.int64)
        D3 = (R3[1:] - R3[0]) % 3
        r3 = gfp_rank(D3.copy(), 3)
        nsb3 = gfp_nullspace(D3, 3)
        c3 = (nsb3 @ R3[0]) % 3
        print(f"    mod 3: rank={r3} -> #F3 invariants = {len(eqs)-r3}; with constant != 0: {int((c3!=0).sum())}")
        if int((c3!=0).sum())>0:
            print("    *** F3 INVARIANT WITH NONZERO CONSTANT FOUND ***")

def gf2_rank(M):
    M = M.astype(np.uint8) % 2
    rows, cols = M.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i,c]: piv=i; break
        if piv is None: continue
        M[[r,piv]] = M[[piv,r]]
        for i in range(rows):
            if i!=r and M[i,c]: M[i] ^= M[r]
        r += 1
        if r==rows: break
    return r

def gf2_nullspace(M):
    """null space (right) basis of M over F2: vectors v with Mv=0."""
    M = M.astype(np.uint8) % 2
    rows, cols = M.shape
    A = M.copy(); pivcol=[]
    r=0
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i; break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        for i in range(rows):
            if i!=r and A[i,c]: A[i]^=A[r]
        pivcol.append(c); r+=1
        if r==rows: break
    free = [c for c in range(cols) if c not in pivcol]
    basis=[]
    for fc in free:
        v = np.zeros(cols,dtype=np.uint8); v[fc]=1
        for ri,pc in enumerate(pivcol):
            if A[ri,fc]: v[pc]=1
        basis.append(v)
    return np.array(basis, dtype=np.uint8) if basis else np.zeros((0,cols),dtype=np.uint8)

def gfp_rank(M, p):
    M = M % p
    rows, cols = M.shape
    r=0
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if M[i,c]%p: piv=i; break
        if piv is None: continue
        M[[r,piv]]=M[[piv,r]]
        inv = pow(int(M[r,c]), p-2, p)
        M[r] = (M[r]*inv)%p
        for i in range(rows):
            if i!=r and M[i,c]%p: M[i]=(M[i]-M[i,c]*M[r])%p
        r+=1
        if r==rows: break
    return r

def gfp_nullspace(M, p):
    M = M % p
    rows, cols = M.shape
    A = M.copy(); pivcol=[]
    r=0
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]%p: piv=i; break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv = pow(int(A[r,c]), p-2, p)
        A[r]=(A[r]*inv)%p
        for i in range(rows):
            if i!=r and A[i,c]%p: A[i]=(A[i]-A[i,c]*A[r])%p
        pivcol.append(c); r+=1
        if r==rows: break
    free=[c for c in range(cols) if c not in pivcol]
    basis=[]
    for fc in free:
        v=np.zeros(cols,dtype=np.int64); v[fc]=1
        for ri,pc in enumerate(pivcol):
            v[pc]=(-A[ri,fc])%p
        basis.append(v)
    return np.array(basis,dtype=np.int64) if basis else np.zeros((0,cols),dtype=np.int64)

if __name__=="__main__":
    main()
