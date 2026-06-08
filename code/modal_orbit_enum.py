"""
Modal: COMPLETE the admissible-orbit-matrix enumeration for srg(333,166,82,83) in parallel.
Single-threaded C++ DFS (sound symmetry break: row-0 sorted within equal-diagonal groups), partitioned
across many containers by the task index of (diagonal, R[0][1], R[0][2]) so even the hardest symmetric
diagonal is split fine-grained. Each container emits its representative orbit matrices; we concatenate to
orbit_all_sym.txt (the COMPLETE representative set), then order2_full_sweep.py closes order 2.

Run:  modal run modal_orbit_enum.py --num-chunks 240
"""
import modal

app = modal.App("h668-orbit-enum")
image = modal.Image.debian_slim().apt_install("g++")

CPP = r'''
#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;
static const int N=9, ROWSUM=166, MAXENT=37;
int R[N][N];
long long nodes=0, nFound=0, taskIdx=0;
int CHUNK_ID, NUM_CHUNKS;
long long diagStartNodes=0; bool abortDiag=false;
const long long CAP=60000000000LL;   // per-(diagonal,this-chunk) node cap
int dchoice[N];
inline int sq(int x){ return x*x; }
bool fillEntry(int i,int j,int rs2,int rq2,int rd[N]);
void record(){ nFound++; for(int a=0;a<N;a++) for(int b=0;b<N;b++) printf("%d ",R[a][b]); printf("\n"); }
bool fillRow(int i){
    if(abortDiag) return false;
    if(i==N){ record(); return false; }
    int bS=R[i][i],bQ=sq(R[i][i]);
    for(int k=0;k<i;k++){ bS+=R[i][k]; bQ+=sq(R[i][k]); }
    int rS=ROWSUM-bS, rQ=(3154-R[i][i])-bQ; if(rS<0||rQ<0) return false;
    int rd[N];
    for(int jp=0;jp<i;jp++){ int b=R[i][i]*R[jp][i]; for(int k=0;k<i;k++) b+=R[i][k]*R[jp][k]; rd[jp]=(3071-R[i][jp])-b; }
    return fillEntry(i,i+1,rS,rQ,rd);
}
inline void knap(const int* w,int c,int bud,int& lo,int& hi){
    int o[N]; for(int t=0;t<c;t++) o[t]=w[t];
    for(int a=1;a<c;a++){int x=o[a],b=a-1;while(b>=0&&o[b]>x){o[b+1]=o[b];b--;}o[b+1]=x;}
    int b=bud,mn=0; for(int t=0;t<c&&b>0;t++){int tk=b<37?b:37;mn+=tk*o[t];b-=tk;}
    b=bud;int mx=0; for(int t=c-1;t>=0&&b>0;t--){int tk=b<37?b:37;mx+=tk*o[t];b-=tk;}
    lo=mn;hi=mx;
}
bool fillEntry(int i,int j,int remSum,int remSumsq,int remDot[N]){
    if(abortDiag) return false;
    if((nodes-diagStartNodes)>CAP){ abortDiag=true; return false; }
    if(j==N){
        if(remSum||remSumsq) return false;
        for(int jp=0;jp<i;jp++) if(remDot[jp]) return false;
        return fillRow(i+1);
    }
    nodes++;
    int slots=N-1-j;
    int vstart=(i==0 && j>=2 && dchoice[j]==dchoice[j-1]) ? R[0][j-1] : 0;
    for(int v=vstart; v<=MAXENT; v++){
        int rs=remSum-v; if(rs<0) break;
        if(rs>MAXENT*slots) continue;
        int rq=remSumsq-sq(v); if(rq<0) continue;
        if(slots>0){
            if((long long)rq*slots<(long long)rs*rs) continue;
            int f=rs/37,rem=rs%37; long long qm=(long long)f*1369+(long long)rem*rem; if((long long)rq>qm) continue;
        } else { if(rs||rq) continue; }
        bool ok=true; int nD[N];
        for(int jp=0;jp<i;jp++){
            int nd=remDot[jp]-v*R[jp][j]; if(nd<0){ok=false;break;}
            if(slots>0){ int w[N]; for(int t=0;t<slots;t++) w[t]=R[jp][j+1+t]; int lo,hi; knap(w,slots,rs,lo,hi); if(nd<lo||nd>hi){ok=false;break;} }
            else if(nd){ok=false;break;}
            nD[jp]=nd;
        }
        if(!ok) continue;
        // PARTITION: each valid (diagonal, R[0][1], R[0][2]) is a task; this container does its residue class.
        if(i==0 && j==2){ if((taskIdx++ % NUM_CHUNKS) != CHUNK_ID) continue; }
        R[i][j]=v; R[j][i]=v;
        fillEntry(i,j+1,rs,rq,nD);
        if(abortDiag){ R[i][j]=0; R[j][i]=0; return false; }
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}
int diagVals[9]={10,12,14,16,18,20,22,24,26};
void runDiag(int pos,int prev,int ssum){
    if(pos==N){
        if(ssum!=162) return;
        for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=dchoice[i]; }
        diagStartNodes=nodes; abortDiag=false;
        fillRow(0);
        if(abortDiag){ fprintf(stderr,"ABORT "); for(int k=0;k<N;k++) fprintf(stderr,"%d ",dchoice[k]); fprintf(stderr,"\n"); }
        return;
    }
    for(int idx=prev;idx<9;idx++){int v=diagVals[idx],ns=ssum+v,rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*26<162) continue; dchoice[pos]=v; runDiag(pos+1,idx,ns);}
}
int main(int argc,char**argv){
    CHUNK_ID=atoi(argv[1]); NUM_CHUNKS=atoi(argv[2]);
    runDiag(0,0,0);
    fprintf(stderr,"SUMMARY chunk=%d found=%lld nodes=%lld\n",CHUNK_ID,nFound,nodes);
    return 0;
}
'''

@app.function(image=image, cpu=1.0, timeout=43200, retries=1)
def enum_chunk(chunk_id, num_chunks):
    import subprocess, tempfile, os
    d = tempfile.mkdtemp()
    src = os.path.join(d, "e.cpp"); open(src, "w").write(CPP)
    subprocess.run(["g++", "-O3", "-o", os.path.join(d, "e"), src], check=True)
    p = subprocess.run([os.path.join(d, "e"), str(chunk_id), str(num_chunks)],
                       capture_output=True, text=True, timeout=43000)
    mats = [ln for ln in p.stdout.splitlines() if ln.strip()]
    aborts = [ln for ln in p.stderr.splitlines() if ln.startswith("ABORT")]
    summ = [ln for ln in p.stderr.splitlines() if ln.startswith("SUMMARY")]
    return {"chunk": chunk_id, "mats": mats, "aborts": aborts, "summary": summ[0] if summ else ""}

@app.local_entrypoint()
def main(num_chunks: int = 240):
    print(f"orbit-matrix enumeration on Modal: {num_chunks} chunks (single-threaded C++ each), parallel")
    all_mats = set(); aborts = []; done = 0
    for r in enum_chunk.starmap([(c, num_chunks) for c in range(num_chunks)]):
        done += 1
        for m in r["mats"]:
            all_mats.add(m.strip())
        aborts += r["aborts"]
        if done % 20 == 0 or r["aborts"]:
            print(f"  {done}/{num_chunks} chunks; unique mats so far={len(all_mats)}; aborts={len(aborts)}", flush=True)
    with open("orbit_all_sym.txt", "w") as f:
        for m in sorted(all_mats):
            f.write(m + "\n")
    print(f"=== ENUMERATION DONE === unique representative orbit matrices = {len(all_mats)}; aborted (diag,chunk) = {len(aborts)}")
    if aborts:
        print("ABORTED diagonals (need finer split / higher cap):")
        for a in aborts[:40]:
            print("  " + a)
    else:
        print("COMPLETE: every diagonal exhausted in every chunk. orbit_all_sym.txt is the full representative set.")
