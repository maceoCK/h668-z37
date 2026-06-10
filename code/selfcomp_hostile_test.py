import numpy as np, itertools, random, sys
sys.path.insert(0,"/Users/maceocardinalekwik/git/math/deepdive_lp333")
from selfcomp_search import Model, build_orbits
mats=[np.array([int(t) for t in ln.split()],dtype=np.int64).reshape(9,9) for ln in open("orbit_all_classes.txt")]
R=[[int(x) for x in row] for row in mats[109]]
Rc=[[(36-R[i][j]) if i==j else (37-R[i][j]) for j in range(9)] for i in range(9)]
perm=next(q for q in itertools.permutations(range(9)) if all(Rc[i][j]==R[q[i]][q[j]] for i in range(9) for j in range(9)))
p, norb, c = 37, 9, 6

class SizeOnly(Model):
    def _build(self):
        pp,nn=self.p,self.norb
        for (i,j) in self.orbits:
            for d in range(pp):
                if i==j and d==0: continue
                self.X[((i,j),d)]=self.m.NewBoolVar(f"x{i}_{j}_{d}")
            if i==j:
                for d in range(1,pp):
                    self.m.Add(self.X[((i,j),d)]==self.X[((i,j),(pp-d)%pp)])
                if self.pi[i]==i:
                    for d in range(1,pp):
                        self.m.Add(self.X[((i,j),(self.c*d)%pp)]+self.X[((i,j),d)]==1)
                self.m.Add(sum(self.X[((i,j),d)] for d in range(1,pp))==int(self.R[i][i]))
            else:
                self.m.Add(sum(self.X[((i,j),d)] for d in range(pp))==int(self.R[i][j]))
        return True
S=SizeOnly(p,norb,166,82,83,R,perm,c)
st,_=S.solve(time_s=60)
print(f"(1) sizes-only class109: {st}")

random.seed(7)
orbits,member,cinv=build_orbits(p,norb,perm,c)
co=lambda S2,dia: (set(range(1,p)) if dia else set(range(p)))-S2
assign={}
for (i,j) in orbits:
    r=R[i][j]
    if i==j:
        if perm[i]==i:
            S2=set(); done=set()
            for d in range(1,p):
                if d in done: continue
                orb={d,(c*d)%p,(-d)%p,(-c*d)%p}; done|=orb
                if random.random()<0.5: S2|={d,(-d)%p}
                else: S2|={(c*d)%p,(-c*d)%p}
            assert len(S2)==18==r
        else:
            pairs=[(d,(p-d)%p) for d in range(1,(p+1)//2)]
            random.shuffle(pairs)
            S2=set()
            for a,b in pairs[:r//2]: S2|={a,b}
        assign[(i,j)]=S2
    else:
        assign[(i,j)]=set(random.sample(range(p),r))
D=[[None]*norb for _ in range(norb)]
for a in range(norb):
    for b in range(norb):
        aa,bb=(a,b) if a<=b else (b,a)
        rep,kk,sw=member[(aa,bb)]
        S2=set(assign[rep])
        for _ in range(kk): S2=co({(c*d)%p for d in S2}, rep[0]==rep[1])
        if sw: S2={(-d)%p for d in S2}
        if a>b: S2={(-d)%p for d in S2}
        D[a][b]=S2
ok=all(D[j][i]=={(-d)%p for d in D[i][j]} for i in range(norb) for j in range(norb))
ok&=all(0 not in D[i][i] and D[i][i]=={(-d)%p for d in D[i][i]} for i in range(norb))
ok&=all(len(D[i][j])==R[i][j] for i in range(norb) for j in range(norb))
ok&=all(D[perm[i]][perm[j]]==co({(c*d)%p for d in D[i][j]}, i==j) for i in range(norb) for j in range(norb))
print(f"(2) plant consistency: {ok}")
def conv_rhs(i,j,g):
    return sum(1 for l in range(norb) for d in D[i][l] if (g-d)%p in D[l][j])
class Planted(Model):
    def _build(self):
        pp,nn=self.p,self.norb
        for (i,j) in self.orbits:
            for d in range(pp):
                if i==j and d==0: continue
                self.X[((i,j),d)]=self.m.NewBoolVar(f"x{i}_{j}_{d}")
            if i==j:
                for d in range(1,pp):
                    self.m.Add(self.X[((i,j),d)]==self.X[((i,j),(pp-d)%pp)])
                if self.pi[i]==i:
                    for d in range(1,pp):
                        self.m.Add(self.X[((i,j),(self.c*d)%pp)]+self.X[((i,j),d)]==1)
                self.m.Add(sum(self.X[((i,j),d)] for d in range(1,pp))==int(self.R[i][i]))
            else:
                self.m.Add(sum(self.X[((i,j),d)] for d in range(pp))==int(self.R[i][j]))
        for (i,j) in self.orbits:
            for g in range(pp):
                terms=[]
                for l in range(nn):
                    for d in range(pp):
                        t=self.product(self.lit(i,l,d), self.lit(l,j,(g-d)%pp))
                        if not (isinstance(t,int) and t==0): terms.append(t)
                self.m.Add(sum(terms)==conv_rhs(i,j,g))
        return True
M=Planted(p,norb,166,82,83,R,perm,c)
st,sols=M.solve(time_s=300)
rec = bool(sols) and all(sols[0][i][j]==D[i][j] for i in range(norb) for j in range(norb))
print(f"(3) planted-RHS engine: {st}  recovered-plant-exactly={rec}")

# (4) DECISIVE: pin all free vars to the plant; model must be OPTIMAL instantly
M2=Planted(p,norb,166,82,83,R,perm,c)
for (rep,d),var in M2.X.items():
    M2.m.Add(var == (1 if d in assign[rep] else 0))
st2,sols2=M2.solve(time_s=120)
rec2 = bool(sols2) and all(sols2[0][i][j]==D[i][j] for i in range(norb) for j in range(norb))
print(f"(4) plant-pinned model: {st2}  extract-matches-plant={rec2}  (OPTIMAL+True required)")
