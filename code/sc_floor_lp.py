"""
Stage 3 (attack a): certified LP relaxation lower bound on min sum|residual| over the
McCormick (+ optional RLT-1) linearization of the sigma-equivariant space, class 221_0.

If optimum > 0: a genuine continuous lower bound on the defect (THEOREM candidate).
If optimum = 0: the LP relaxation is provably blind (catalog entry), with an explicit
fractional witness.
"""
import sys, numpy as np
from scipy import sparse
from scipy.optimize import linprog
sys.path.insert(0, "/Users/maceocardinalekwik/git/math/release/h668-z37/code")
from sc_floor_analysis import Tpl, CODE

def build_and_solve(tplname="sc_tpl_221_0.txt", rlt=False):
    tpl = Tpl(f"{CODE}/{tplname}")
    P, N = tpl.P, tpl.NORB
    # ---- referenced bits ----
    ref = set()
    for a in range(N):
        for b in range(N):
            for d in range(P):
                ri = tpl.litRep[a,b,d]
                if ri >= 0: ref.add((int(ri), int(tpl.litE[a,b,d])))
    bits = sorted(ref)
    bid = {b:i for i,b in enumerate(bits)}
    nx = len(bits)
    print(f"[{tplname}] bits referenced: {nx}")
    # ---- lazy y vars ----
    ydict = {}
    def yvar(u, v):
        if u > v: u, v = v, u
        if (u,v) not in ydict: ydict[(u,v)] = len(ydict)
        return ydict[(u,v)]

    def lit_of(a,b,d):
        ri = tpl.litRep[a,b,d]
        if ri < 0: return None, int(tpl.litNeg[a,b,d])         # const
        return (bid[(int(ri), int(tpl.litE[a,b,d]))], int(tpl.litNeg[a,b,d])), None

    eqs = tpl.rep_eq_list()
    neq = len(eqs)
    # residual_i = const_i + sum cx*x + sum cy*y ; assemble per-eq dicts
    eq_const = np.zeros(neq)
    eq_x = [dict() for _ in range(neq)]
    eq_y = [dict() for _ in range(neq)]
    def addx(i,u,c):
        if c: eq_x[i][u] = eq_x[i].get(u,0)+c
    def addy(i,p,c):
        if c: eq_y[i][p] = eq_y[i].get(p,0)+c
    for i,(a,b,g) in enumerate(eqs):
        # conv terms
        for l in range(N):
            for d in range(P):
                l1,c1 = lit_of(a,l,d)
                l2,c2 = lit_of(l,b,(g-d)%P)
                if l1 is None and l2 is None:
                    eq_const[i] += c1*c2; continue
                if l1 is None or l2 is None:
                    c, lit = (c1,l2) if l1 is None else (c2,l1)
                    if c == 0: continue
                    u,neg = lit
                    if neg: eq_const[i]+=1; addx(i,u,-1)
                    else: addx(i,u,1)
                    continue
                (u,n1),(v,n2) = l1,l2
                if u==v:
                    if n1==n2:
                        if n1: eq_const[i]+=1; addx(i,u,-1)
                        else: addx(i,u,1)
                    # else x(1-x)=0 exactly
                    continue
                a1,s1 = (1,-1) if n1 else (0,1)
                a2,s2 = (1,-1) if n2 else (0,1)
                eq_const[i] += a1*a2
                addx(i,v,a1*s2); addx(i,u,a2*s1)
                addy(i,(min(u,v),max(u,v)), s1*s2)
        # + mem(a,b,g)  (lam-mu = -1 moves it to LHS with +1)
        lt,c = lit_of(a,b,g)
        if lt is None: eq_const[i] += c
        else:
            u,neg = lt
            if neg: eq_const[i]+=1; addx(i,u,-1)
            else: addx(i,u,1)
        # - rhs
        eq_const[i] -= tpl.MU + ((tpl.K-tpl.MU) if (a==b and g==0) else 0)

    # ---- equality constraints on x ----
    eq_rows = []   # (dict var->coef incl y offset later, rhs)
    for r,(ii,jj,typ,card) in enumerate(tpl.reps):
        row = { bid[(r,e)]: 1 for e in range(P) if (r,e) in bid }
        eq_rows.append((row, card, {}))
        if typ==1:
            for e in range(1,P//2+1):
                if (r,e) in bid and (r,P-e) in bid:
                    eq_rows.append(({bid[(r,e)]:1, bid[(r,P-e)]:-1}, 0, {}))
        if typ==2:
            for e in range(1,P):
                if e < P-e and (r,e) in bid and (r,P-e) in bid:
                    eq_rows.append(({bid[(r,e)]:1, bid[(r,P-e)]:-1}, 0, {}))
            for e in range(1,P):
                e6 = (tpl.C*e)%P
                if e < e6 and (r,e) in bid and (r,e6) in bid:
                    eq_rows.append(({bid[(r,e)]:1, bid[(r,e6)]:1}, 1, {}))
    base_eq_rows = list(eq_rows)

    # ---- RLT-1: multiply every base x-equality (sum c_u x_u = rhs) by every bit w:
    #      sum c_u y_{uw} = rhs * x_w   (y_{ww} := x_w)
    if rlt:
        for (row, rhs, _) in base_eq_rows:
            for w in range(nx):
                xs = {w: -rhs}
                ys = {}
                ok = True
                for u,c in row.items():
                    if u==w: xs[w] = xs.get(w,0)+c
                    else: ys[(min(u,w),max(u,w))] = ys.get((min(u,w),max(u,w)),0)+c
                eq_rows.append((xs, 0, ys))

    ny = None
    # force creation of y vars used anywhere
    for i in range(neq):
        for pkey in eq_y[i]: yvar(*pkey)
    for (_,_,ys) in eq_rows:
        for pkey in ys: yvar(*pkey)
    ny = len(ydict)
    print(f"  eqs={neq}  y-vars={ny}  x-eq-rows={len(eq_rows)}  rlt={rlt}")

    # variable layout: x (nx), y (ny), p (neq), n (neq)
    NV = nx + ny + 2*neq
    def yix(p): return nx + ydict[p]
    rows_eq = []; cols_eq=[]; vals_eq=[]; rhs_eq=[]
    ridx = 0
    # residual defs: const + sum = p - n  ->  sum cx*x + cy*y - p + n = -const
    for i in range(neq):
        for u,c in eq_x[i].items():
            rows_eq.append(ridx); cols_eq.append(u); vals_eq.append(c)
        for pkey,c in eq_y[i].items():
            rows_eq.append(ridx); cols_eq.append(yix(pkey)); vals_eq.append(c)
        rows_eq.append(ridx); cols_eq.append(nx+ny+i); vals_eq.append(-1.0)
        rows_eq.append(ridx); cols_eq.append(nx+ny+neq+i); vals_eq.append(1.0)
        rhs_eq.append(-eq_const[i]); ridx+=1
    for (row, rhs, ys) in eq_rows:
        for u,c in row.items():
            rows_eq.append(ridx); cols_eq.append(u); vals_eq.append(c)
        for pkey,c in ys.items():
            rows_eq.append(ridx); cols_eq.append(yix(pkey)); vals_eq.append(c)
        rhs_eq.append(rhs); ridx+=1
    A_eq = sparse.coo_matrix((vals_eq,(rows_eq,cols_eq)), shape=(ridx,NV)).tocsr()
    b_eq = np.array(rhs_eq, dtype=float)

    # McCormick: y <= xu ; y <= xv ; y >= xu+xv-1  (y>=0 in bounds)
    rows_ub=[]; cols_ub=[]; vals_ub=[]; rhs_ub=[]
    rr=0
    for (u,v),k in ydict.items():
        yc = nx+k
        rows_ub += [rr,rr]; cols_ub += [yc,u]; vals_ub += [1.0,-1.0]; rhs_ub.append(0.0); rr+=1
        rows_ub += [rr,rr]; cols_ub += [yc,v]; vals_ub += [1.0,-1.0]; rhs_ub.append(0.0); rr+=1
        rows_ub += [rr,rr,rr]; cols_ub += [u,v,yc]; vals_ub += [1.0,1.0,-1.0]; rhs_ub.append(1.0); rr+=1
    A_ub = sparse.coo_matrix((vals_ub,(rows_ub,cols_ub)), shape=(rr,NV)).tocsr()
    b_ub = np.array(rhs_ub, dtype=float)

    c = np.zeros(NV); c[nx+ny:] = 1.0
    bounds = [(0,1)]*(nx+ny) + [(0,None)]*(2*neq)
    print(f"  LP: vars={NV} eq-rows={A_eq.shape[0]} ub-rows={A_ub.shape[0]} "
          f"nnz={A_eq.nnz+A_ub.nnz}", flush=True)
    r = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                method="highs")
    print(f"  status={r.status} ({r.message.strip()})  optimum sum|res| = {r.fun:.6e}")
    if r.status==0 and r.fun is not None and r.fun < 1e-7:
        x = r.x[:nx]
        print(f"  -> RELAXATION BLIND: fractional witness exists; x in [{x.min():.3f},{x.max():.3f}], "
              f"#fractional bits = {(np.minimum(x,1-x)>1e-6).sum()}/{nx}")
    elif r.status==0:
        print("  -> *** POSITIVE LP LOWER BOUND -- CHECK CAREFULLY ***")
    return r

if __name__=="__main__":
    rlt = "--rlt" in sys.argv
    tname = sys.argv[1] if len(sys.argv)>1 and not sys.argv[1].startswith("--") else "sc_tpl_221_0.txt"
    build_and_solve(tname, rlt=rlt)
