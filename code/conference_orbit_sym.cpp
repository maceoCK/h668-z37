// COMPLETE enumeration of admissible 9x9 orbit matrices R (up to orbit relabeling) for srg(333,166,82,83),
// with a SOUND symmetry break that finishes the highly-symmetric diagonals the plain enumerator could not.
//
// Diagonal sorted non-decreasing (uses S_9 to order the diagonal multiset). WITHIN each maximal
// equal-diagonal group, orbits are sorted by their connection R[0][.] to orbit 0 (the residual
// orbit-relabeling freedom). Both are sound: every isomorphism class has a representative satisfying them.
// This is the per-group generalization of the all-18 row-0 sort (which finished in 9s / 1.9e8 nodes).
//
// Goal: exhaust ALL diagonals (nUnres==0). The resulting file is the COMPLETE set of representative
// admissible orbit matrices, over which the order-2 lift sweep then runs -> closes order 2 if all UNSAT.
#include <cstdio>
#include <ctime>
using namespace std;
static const int N = 9, ROWSUM = 166, MAXENT = 37;
int R[N][N];
long long nodes = 0, diagsDone = 0, diagsTotal = 0, nFound = 0;
const long long CAP = 200000000000LL;     // per-diagonal node cap (200e9; raise if a diagonal aborts)
long long diagStartNodes = 0;
bool abortDiag = false;
int unresolved[4000][N]; int nUnres = 0;
FILE* fout = nullptr;
time_t t0;
int dchoice[N];
inline int sq(int x){ return x*x; }
bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]);

void record(){
    nFound++;
    if(fout){ for(int a=0;a<N;a++) for(int b=0;b<N;b++) fprintf(fout,"%d ",R[a][b]); fprintf(fout,"\n"); }
}
bool fillRow(int i){
    if(abortDiag) return false;
    if(i==N){ record(); return false; }
    int baseSum=R[i][i], baseSumsq=sq(R[i][i]);
    for(int k=0;k<i;k++){ baseSum+=R[i][k]; baseSumsq+=sq(R[i][k]); }
    int remSum=ROWSUM-baseSum, remSumsq=(3154-R[i][i])-baseSumsq;
    if(remSum<0||remSumsq<0) return false;
    int remDot[N];
    for(int jp=0;jp<i;jp++){ int b=R[i][i]*R[jp][i]; for(int k=0;k<i;k++) b+=R[i][k]*R[jp][k]; remDot[jp]=(3071-R[i][jp])-b; }
    return fillEntry(i,i+1,remSum,remSumsq,remDot);
}
inline void knap_bounds(const int* w,int cnt,int budget,int& lo,int& hi){
    int o[N]; for(int t=0;t<cnt;t++) o[t]=w[t];
    for(int a=1;a<cnt;a++){int x=o[a],b=a-1;while(b>=0&&o[b]>x){o[b+1]=o[b];b--;}o[b+1]=x;}
    int b=budget,mn=0; for(int t=0;t<cnt&&b>0;t++){int tk=b<37?b:37;mn+=tk*o[t];b-=tk;}
    b=budget;int mx=0; for(int t=cnt-1;t>=0&&b>0;t--){int tk=b<37?b:37;mx+=tk*o[t];b-=tk;}
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
    if((nodes & ((1LL<<32)-1))==0){ printf("  ... nodes=%lld diag=%lld/%lld found=%lld %.0fs\n",
        nodes,diagsDone,diagsTotal,nFound,(double)(time(nullptr)-t0)); fflush(stdout); }
    int slotsAfter=N-1-j;
    // SOUND symmetry break: in row 0, within an equal-diagonal group, sort by R[0][.].
    int vstart = (i==0 && j>=2 && dchoice[j]==dchoice[j-1]) ? R[0][j-1] : 0;
    for(int v=vstart; v<=MAXENT; v++){
        int rs=remSum-v; if(rs<0) break;
        if(rs>MAXENT*slotsAfter) continue;
        int rq=remSumsq-sq(v); if(rq<0) continue;
        if(slotsAfter>0){
            if((long long)rq*slotsAfter<(long long)rs*rs) continue;
            int f=rs/37,rem=rs%37; long long qm=(long long)f*1369+(long long)rem*rem; if((long long)rq>qm) continue;
        } else { if(rs||rq) continue; }
        bool ok=true; int nD[N];
        for(int jp=0;jp<i;jp++){
            int nd=remDot[jp]-v*R[jp][j]; if(nd<0){ok=false;break;}
            if(slotsAfter>0){ int w[N]; for(int t=0;t<slotsAfter;t++) w[t]=R[jp][j+1+t]; int lo,hi; knap_bounds(w,slotsAfter,rs,lo,hi); if(nd<lo||nd>hi){ok=false;break;} }
            else if(nd){ok=false;break;}
            nD[jp]=nd;
        }
        if(!ok) continue;
        R[i][j]=v; R[j][i]=v;
        fillEntry(i,j+1,rs,rq,nD);
        if(abortDiag) return false;
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}
int diagVals[9]={10,12,14,16,18,20,22,24,26};
void recDiag(int pos,int prev,int ssum){
    if(pos==N){ if(ssum==162) diagsTotal++; return; }
    for(int idx=prev;idx<9;idx++){int v=diagVals[idx],ns=ssum+v,rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*26<162) continue; recDiag(pos+1,idx,ns);}
}
void runDiag(int pos,int prev,int ssum){
    if(pos==N){
        if(ssum!=162) return;
        for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=dchoice[i]; }
        diagStartNodes=nodes; abortDiag=false;
        long long before=nFound;
        fillRow(0);
        if(abortDiag){ for(int k=0;k<N;k++) unresolved[nUnres][k]=dchoice[k]; nUnres++; }
        diagsDone++;
        if(diagsDone%40==0 || nFound>before || abortDiag){
            printf("  diag %lld/%lld [",diagsDone,diagsTotal);
            for(int k=0;k<N;k++) printf("%d%s",dchoice[k],k<N-1?",":"");
            printf("] found=%lld(+%lld) unres=%d nodes=%lld %.0fs%s\n",
                nFound,nFound-before,nUnres,nodes,(double)(time(nullptr)-t0),abortDiag?"  <ABORT>":"");
            fflush(stdout);
        }
        return;
    }
    for(int idx=prev;idx<9;idx++){int v=diagVals[idx],ns=ssum+v,rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*26<162) continue; dchoice[pos]=v; runDiag(pos+1,idx,ns);}
}
int main(){
    t0=time(nullptr);
    recDiag(0,0,0);
    printf("candidate diagonals: %lld\n",diagsTotal); fflush(stdout);
    fout=fopen("orbit_all_sym.txt","w");
    runDiag(0,0,0);
    if(fout) fclose(fout);
    printf("\n=== DONE === representatives found=%lld  unresolved diagonals=%d/%lld  nodes=%lld  %.0fs\n",
        nFound,nUnres,diagsTotal,nodes,(double)(time(nullptr)-t0));
    if(nUnres==0) printf("ENUMERATION COMPLETE: all %lld diagonals exhausted; %lld representative orbit matrices.\n",diagsTotal,nFound);
    else { printf("UNRESOLVED diagonals (hit %lld-node cap):\n",CAP);
        for(int u=0;u<nUnres&&u<60;u++){ printf("  ["); for(int k=0;k<N;k++) printf("%d%s",unresolved[u][k],k<N-1?",":""); printf("]\n"); } }
    return 0;
}
