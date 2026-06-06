"""CP-SAT search for a Legendre pair of length ell (to obtain real examples)."""
import sys
from ortools.sat.python import cp_model
import numpy as np

def find_lp(ELL, time_s=180, seed=0):
    m = cp_model.CpModel()
    b = [m.NewBoolVar(f"b{i}") for i in range(ELL)]   # u_i = 2b_i-1
    d = [m.NewBoolVar(f"d{i}") for i in range(ELL)]    # v_i = 2d_i-1
    half = ELL // 2
    # row sums: sum u = +-1  => sum b in {(ELL-1)/2,(ELL+1)/2}
    sb = m.NewIntVar((ELL-1)//2, (ELL+1)//2, "sb"); m.Add(sum(b)==sb)
    sd = m.NewIntVar((ELL-1)//2, (ELL+1)//2, "sd"); m.Add(sum(d)==sd)
    # agreements: for k=1..half, sum_i [b_i==b_{i+k}] + [d_i==d_{i+k}] = ELL-1
    for k in range(1, half+1):
        eu=[]; ev=[]
        for i in range(ELL):
            j=(i+k)%ELL
            e=m.NewBoolVar(f"eu_{k}_{i}"); m.Add(b[i]==b[j]).OnlyEnforceIf(e); m.Add(b[i]!=b[j]).OnlyEnforceIf(e.Not()); eu.append(e)
            f=m.NewBoolVar(f"ev_{k}_{i}"); m.Add(d[i]==d[j]).OnlyEnforceIf(f); m.Add(d[i]!=d[j]).OnlyEnforceIf(f.Not()); ev.append(f)
        m.Add(sum(eu)+sum(ev)==ELL-1)
    m.Add(b[0]==1)  # negation symmetry break
    sv=cp_model.CpSolver()
    sv.parameters.max_time_in_seconds=time_s
    sv.parameters.num_search_workers=8
    sv.parameters.random_seed=seed
    st=sv.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        u=np.array([1 if sv.Value(b[i]) else -1 for i in range(ELL)])
        v=np.array([1 if sv.Value(d[i]) else -1 for i in range(ELL)])
        return u,v,sv.StatusName(st)
    return None,None,sv.StatusName(st)

if __name__=="__main__":
    ELL=int(sys.argv[1]) if len(sys.argv)>1 else 45
    t=int(sys.argv[2]) if len(sys.argv)>2 else 180
    u,v,status=find_lp(ELL,t)
    print(f"LP({ELL}) status: {status}")
    if u is not None:
        pu=np.array([int(np.dot(u,np.roll(u,-k))) for k in range(ELL)])
        pv=np.array([int(np.dot(v,np.roll(v,-k))) for k in range(ELL)])
        ok=all(pu[k]+pv[k]==-2 for k in range(1,ELL)) and abs(int(u.sum()))==1 and abs(int(v.sum()))==1
        print(f"  verified LP: {ok}  rowsums=({int(u.sum())},{int(v.sum())})")
        np.save(f"deepdive_lp333/lp{ELL}_u.npy",u); np.save(f"deepdive_lp333/lp{ELL}_v.npy",v)
        print("  u=",list(map(int,u))); print("  v=",list(map(int,v)))
