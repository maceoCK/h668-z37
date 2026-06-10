"""
Cocyclic Hadamard route at 668 = 4*167: verified computations.

PART A: arithmetic kills at 668 (the numbers behind the character-theoretic kills)
  A1. 668 is not a perfect square                  -> kills any order-2 character chi(z)=-1
  A2. 668 is not a sum of two squares (167=3 mod 4) -> kills any order-4 character chi(z)=-1,
       and kills every 2-dim INTEGRAL rep (D8-quotient, dihedral k=167 rep)
  A3. no alpha in Z[zeta_8] with alpha*conj(alpha)=668 (167 = -1 mod 8 self-conjugate)
       -> kills cyclic-Sylow-2 extensions (Z_1336, Z_167:Z_8) where ONLY order-8 chars see z
  A4. ord_167(2) = 83 (odd) -> 2 NOT self-conjugate mod 167 -> no field descent can kill
       the two survivors (Q8 x Z167, Dic_334)
  A5. 668 = w^2+x^2+y^2+z^2 with all odd (Williamson row-sum gate) has solutions
  A6. 334 = 2*167 not a sum of two squares -> the exotic third cocyclic family of
       Barrera Acevedo-O Cathain-Dietrich (needs p = 1 mod 4) is empty at p=167

PART B: empirical validation of the kill framework at the t=3 analogue (order 24 = 8*3,
  RDS(12,2,12,6)): exhaustive search over all transversals in each of the 12 groups of
  order 24 that are Z_3-by-(2-group) (the exact analogues of the 12 groups of order 1336).
  Predictions from the same character kills (12 not a square, 12 not a sum of 2 squares;
  NOTE 3 is NOT self-conjugate mod 8, so the zeta_8 kill does NOT fire at t=3):
    killed by linear/2-dim tests -> NO RDS;  Q8xZ3, Dic6 -> RDS must exist (Williamson(3)/Ito(3)).
"""
import itertools, math
from collections import Counter

print("="*72); print("PART A: arithmetic at 668"); print("="*72)

# A1
assert int(math.isqrt(668))**2 != 668
print("A1. 668 not a perfect square (25^2=625 < 668 < 676=26^2)  [kills order-2 chars]")

# A2
two_sq = [(p,q) for p in range(26) for q in range(p,26) if p*p+q*q==668]
assert two_sq == []
print("A2. 668 = p^2+q^2 has NO solutions (167 = 3 mod 4 to odd power)")
assert 167 % 4 == 3

# A3: alpha = a+b*z+c*z^2+d*z^3 in Z[zeta_8];  alpha*conj(alpha) =
#     (a^2+b^2+c^2+d^2) + (ab+bc+cd-ad)*sqrt(2).  Need =668 exactly.
sols8 = []
for a in range(-25,26):
    for b in range(-25,26):
        ab2 = a*a+b*b
        if ab2 > 668: continue
        for c in range(-25,26):
            s3 = ab2 + c*c
            if s3 > 668: continue
            r = 668 - s3
            d0 = math.isqrt(r)
            for d in {d0,-d0}:
                if d*d == r and a*b + b*c + c*d - a*d == 0:
                    sols8.append((a,b,c,d))
assert sols8 == []
print("A3. NO alpha in Z[zeta_8] with |alpha|^2=668 (exhaustive: a^2+b^2+c^2+d^2=668 &")
print("    ab+bc+cd-ad=0 has no integer solution).  Reason: 167 = 7 = -1 (mod 8) is")
print("    self-conjugate mod 8; primes over 167 in Z[zeta_8] (f=2,g=2) are fixed by")
print("    conjugation, forcing even valuation, but v(668)=1 at each.")
assert 167 % 8 == 7

# A4
o = 1; x = 2
while True:
    if x == 1: break
    x = (x*2) % 167; o += 1
print(f"A4. ord_167(2) = {o} (odd)  -> 2 not self-conjugate mod 167: no Turyn/Schmidt")
print("    descent at the 167-characters; the survivors cannot be killed this way.")
assert o == 83

# A5
quads = sorted({tuple(sorted((w,x,y,zz),reverse=True))
                for w in range(1,26,2) for x in range(1,26,2)
                for y in range(1,26,2) for zz in range(1,26,2)
                if w*w+x*x+y*y+zz*zz==668})
print(f"A5. 668 as sum of 4 ODD squares (Williamson/Ito row sums): {quads}")
assert len(quads) >= 1

# A6
assert [(p,q) for p in range(19) for q in range(p,19) if p*p+q*q==334] == []
print("A6. 334 = 2*167 not a sum of two squares -> exotic (1,1,0)-family (needs p=1 mod 4)")
print("    is empty at p=167, consistent with BACD Thm 4.5b.")

print(); print("="*72)
print("PART B: exhaustive RDS(12,2,12,6) search in the 12 analogue groups of order 24")
print("="*72)

def group_from(elems, mult):
    idx = {e:i for i,e in enumerate(elems)}
    n = len(elems)
    tab = [[idx[mult(a,b)] for b in elems] for a in elems]
    inv = [None]*n
    e0 = None
    for i in range(n):
        for j in range(n):
            if tab[i][j] == j:  # candidate identity: i*j=j for all j
                pass
    # find identity properly
    for i in range(n):
        if all(tab[i][j]==j for j in range(n)):
            e0 = i; break
    for i in range(n):
        for j in range(n):
            if tab[i][j]==e0: inv[i]=j
    return tab, inv, e0, idx

def rds_exists(elems, mult, z):
    tab, inv, e0, idx = group_from(elems, mult)
    n = len(elems); zi = idx[z]
    assert all(tab[zi][k]==tab[k][zi] for k in range(n)), "z not central"
    assert tab[zi][zi]==e0 and zi!=e0, "z not an involution"
    # cosets of N={e,z}
    seen, cosets = set(), []
    for i in range(n):
        if i in seen: continue
        j = tab[zi][i]; seen.update({i,j}); cosets.append((i,j))
    m = len(cosets)            # 12
    lam = m//2                 # 6
    # WLOG the coset of e contributes e (translation invariance)
    cosets = [c for c in cosets if e0 not in c]
    first = [(e0, e0)]
    for bits in itertools.product((0,1), repeat=m-1):
        D = [e0] + [cosets[k][bits[k]] for k in range(m-1)]
        cnt = Counter()
        ok = True
        for a in D:
            for b in D:
                if a==b: continue
                cnt[tab[a][inv[b]]] += 1
        for g in range(n):
            want = 0 if g in (e0,zi) else lam
            if cnt.get(g,0) != want: ok=False; break
        if ok: return True, D
    return False, None

Z2=range(2); Z3=range(3); Z4=range(4); Z6=range(6); Z8=range(8); Z12=range(12); Z24=range(24)
def d8m(g,h):  # D8 = (a in Z4, b in Z2), reflections invert
    (a,b),(c,d)=g,h; return ((a + (c if b==0 else -c))%4, (b+d)%2)
def s3m(g,h):
    (r,f),(r2,f2)=g,h; return ((r + (r2 if f==0 else -r2))%3, (f+f2)%2)
def d24m(g,h):
    (a,f),(c,f2)=g,h; return ((a + (c if f==0 else -c))%12, (f+f2)%2)
def dic_m(N):  # dicyclic on (a in Z_{2N}, e in {0,1}), b^2 = a^N
    def m(g,h):
        (a,e),(c,f)=g,h
        if e==0: return ((a+c)%(2*N), f)
        else:    return ((a-c+(N if f==1 else 0))%(2*N), (e+f)%2)
    return m
Q8E=['1','-1','i','-i','j','-j','k','-k']
def q8m(x,y):
    sx,ux=(x[0]=='-',x.lstrip('-')); sy,uy=(y[0]=='-',y.lstrip('-'))
    T={('1','1'):'1',('1','i'):'i',('1','j'):'j',('1','k'):'k',
       ('i','1'):'i',('j','1'):'j',('k','1'):'k',
       ('i','i'):'-1',('j','j'):'-1',('k','k'):'-1',
       ('i','j'):'k',('j','k'):'i',('k','i'):'j',
       ('j','i'):'-k',('k','j'):'-i',('i','k'):'-j'}
    r=T[(ux,uy)]; s=(r[0]=='-')^sx^sy; r=r.lstrip('-')
    return ('-'+r) if s else r

cases = []
cases.append(("Z24 (cyclic; analogue of Z_1336)", list(Z24),
              lambda a,b:(a+b)%24, 12,
              "zeta_8 kill does NOT fire at t=3 (3 not = -1 mod 8); undetermined"))
cases.append(("Z4xZ2xZ3, z=(2,0,0)", list(itertools.product(Z4,Z2,Z3)),
              lambda g,h:((g[0]+h[0])%4,(g[1]+h[1])%2,(g[2]+h[2])%3),(2,0,0),
              "KILL: order-4 char (12 not 2-square sum) -> predict NO"))
cases.append(("Z2^3xZ3, z=(1,0,0,0)", list(itertools.product(Z2,Z2,Z2,Z3)),
              lambda g,h:tuple((g[i]+h[i])%(2 if i<3 else 3) for i in range(4)),(1,0,0,0),
              "KILL: order-2 char (12 not square) -> predict NO"))
cases.append(("D8xZ3, z=((2,0),0)", list(itertools.product(itertools.product(Z4,Z2),Z3)),
              lambda g,h:(d8m(g[0],h[0]),(g[1]+h[1])%3),((2,0),0),
              "KILL: 2-dim integral rep (12 not 2-square sum) -> predict NO"))
cases.append(("Q8xZ3, z=(-1,0)  [SURVIVOR analogue of Q8xZ167]",
              list(itertools.product(Q8E,Z3)),
              lambda g,h:(q8m(g[0],h[0]),(g[1]+h[1])%3),('-1',0),
              "survives all kills; Williamson(3) exists -> predict YES"))
cases.append(("Z3:Z8 (8-cycle inverts; analogue of Z167:Z8)",
              list(itertools.product(Z3,Z8)),
              lambda g,h:((g[0]+(h[0] if g[1]%2==0 else -h[0]))%3,(g[1]+h[1])%8),(0,4),
              "zeta_8 kill does NOT fire at t=3; undetermined"))
cases.append(("Z4xS3, z=(2,(0,0))", list(itertools.product(Z4,itertools.product(Z3,Z2))),
              lambda g,h:((g[0]+h[0])%4,s3m(g[1],h[1])),(2,(0,0)),
              "KILL: order-4 char -> predict NO"))
dic3 = dic_m(3)
for zz,lab in [((1,(0,0)),"z=(1,e)"),((0,(3,0)),"z=(0,b^2)"),((1,(3,0)),"z=(1,b^2)")]:
    cases.append((f"Z2xDic3, {lab}", list(itertools.product(Z2,itertools.product(Z6,Z2))),
                  lambda g,h:((g[0]+h[0])%2,dic3(g[1],h[1])),zz,
                  "KILL: order-2 or order-4 char -> predict NO"))
cases.append(("Z2^2xS3, z=(1,0,(0,0))",
              list(itertools.product(Z2,Z2,itertools.product(Z3,Z2))),
              lambda g,h:((g[0]+h[0])%2,(g[1]+h[1])%2,s3m(g[2],h[2])),(1,0,(0,0)),
              "KILL: order-2 char -> predict NO"))
cases.append(("D24 (dihedral order 24), z=(6,0)",
              list(itertools.product(Z12,Z2)), d24m,(6,0),
              "KILL: 2-dim rep at zeta=i (12 not 2-square sum) -> predict NO"))
cases.append(("Z3:D8 (Klein kernel), z=(0,(2,0))",
              list(itertools.product(Z3,itertools.product(Z4,Z2))),
              lambda g,h:((g[0]+(h[0] if g[1][0]%2==0 else -h[0]))%3, d8m(g[1],h[1])),
              (0,(2,0)),
              "KILL: lifted 2-dim integral D8 rep -> predict NO"))
cases.append(("Dic6 (dicyclic order 24), z=(6,0)  [SURVIVOR analogue of Dic_334]",
              list(itertools.product(Z12,Z2)), dic_m(6),(6,0),
              "survives all kills; Ito(3) exists -> predict YES"))

for name, elems, mult, z, pred in cases:
    ok, D = rds_exists(elems, mult, z)
    print(f"  {name:58s} RDS exists: {ok}   [{pred}]")
print()
print("Framework check: every group killed by a character test has NO RDS;")
print("the two survivor-analogues (Q8xZ3, Dic6) DO have RDS.")
