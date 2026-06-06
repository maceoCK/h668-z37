// Decide existence of the 9x9 orbit matrix R for a fixed-point-free order-37 automorphism
// of srg(333,166,82,83).  R symmetric, R[i][i]=d_i even in {10..26} (sum 162),
// R[i][j] in [0,37] (i!=j).  Per row: sum=166, sumsq=3154-d_i.  Cross: dot(i,j)=3071-R[i][j].
// Equivalent to R^2 + R - 83 I = 3071 J.  Streaming DFS with incremental pruning.
#include <cstdio>
#include <cstdint>
#include <vector>
#include <array>
#include <ctime>
#include <cstring>
using namespace std;

static const int N = 9;
static const int ROWSUM = 166;
static const int MAXENT = 37;

int R[N][N];
long long nodes = 0;
long long diagsDone = 0, diagsTotal = 0;
bool found = false;
time_t t0;
const long long CAP = 50000000LL;     // per-diagonal node cap
long long diagStartNodes = 0;
bool abortDiag = false;
int unresolved[2000][N]; int nUnres = 0;

// suffix helpers recomputed per row (small)
inline int sq(int x){ return x*x; }

// fill row i, choosing R[i][j] for j=i+1..N-1, given diagonal set and rows<i complete.
// remSum, remSumsq: remaining target for entries j..N-1 in this row.
// remDot[j'] : remaining target dot contribution from entries k=j..N-1 for prior row j'<i.
bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]);

bool fillRow(int i){
    if(found) return true;
    if(abortDiag) return false;
    if(i==N){ found=true; return true; }
    // base contributions from already-known part of row i (k<i and diagonal)
    int baseSum = R[i][i];
    int baseSumsq = sq(R[i][i]);
    for(int k=0;k<i;k++){ baseSum += R[i][k]; baseSumsq += sq(R[i][k]); }
    int remSum = ROWSUM - baseSum;
    int remSumsq = (3154 - R[i][i]) - baseSumsq;
    if(remSum<0 || remSumsq<0) return false;
    int remDot[N];
    for(int jp=0;jp<i;jp++){
        // target dot(i,jp) = 3071 - R[i][jp]; base part over k<=i
        int base = R[i][i]*R[jp][i];
        for(int k=0;k<i;k++) base += R[i][k]*R[jp][k];
        remDot[jp] = (3071 - R[i][jp]) - base;
    }
    return fillEntry(i, i+1, remSum, remSumsq, remDot);
}

// max/min of sum_k x_k * w[k] over integers x_k in [0,37] with sum x_k = budget.
inline void knap_bounds(const int* w, int cnt, int budget, int& lo, int& hi){
    // sort indices by weight (small arrays); compute greedy fills
    int ord[N]; for(int t=0;t<cnt;t++) ord[t]=w[t];
    // simple insertion sort ascending
    for(int a=1;a<cnt;a++){ int x=ord[a],b=a-1; while(b>=0&&ord[b]>x){ord[b+1]=ord[b];b--;} ord[b+1]=x; }
    // min: fill budget into smallest weights
    int b=budget, mn=0;
    for(int t=0;t<cnt && b>0;t++){ int take=b<37?b:37; mn+=take*ord[t]; b-=take; }
    // max: fill into largest weights
    b=budget; int mx=0;
    for(int t=cnt-1;t>=0 && b>0;t--){ int take=b<37?b:37; mx+=take*ord[t]; b-=take; }
    lo=mn; hi=mx;
}

bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]){
    if(found) return true;
    if(abortDiag) return false;
    if((nodes - diagStartNodes) > CAP){ abortDiag = true; return false; }
    if(j==N){
        if(remSum!=0 || remSumsq!=0) return false;
        for(int jp=0;jp<i;jp++) if(remDot[jp]!=0) return false;
        return fillRow(i+1);
    }
    nodes++;
    int slotsAfter = N-1-j;            // entries k=j+1..N-1 after this one
    for(int v=0; v<=MAXENT; v++){
        int rs = remSum - v;
        if(rs<0) break;                       // sum only grows; larger v worse
        if(rs > MAXENT*slotsAfter) continue;  // can't reach budget with remaining
        int rq = remSumsq - sq(v);
        if(rq<0) continue;
        if(slotsAfter>0){
            // Cauchy-Schwarz lower bound on remaining sumsq
            if((long long)rq*slotsAfter < (long long)rs*rs) continue;
            // tight upper bound on remaining sumsq: concentrate mass (f37 entries at 37)
            int f37 = rs/37, rem = rs%37;
            long long qmax = (long long)f37*1369 + (long long)rem*rem;
            if((long long)rq > qmax) continue;
        } else { if(rs!=0 || rq!=0) continue; }
        // tight dot pruning via knapsack on the remaining slots k=j+1..N-1
        bool okdot=true;
        int newDot[N];
        for(int jp=0;jp<i;jp++){
            int nd = remDot[jp] - v*R[jp][j];
            if(nd<0){ okdot=false; break; }
            if(slotsAfter>0){
                int w[N]; for(int t=0;t<slotsAfter;t++) w[t]=R[jp][j+1+t];
                int lo,hi; knap_bounds(w, slotsAfter, rs, lo, hi);
                if(nd<lo || nd>hi){ okdot=false; break; }
            } else if(nd!=0){ okdot=false; break; }
            newDot[jp]=nd;
        }
        if(!okdot) continue;
        R[i][j]=v; R[j][i]=v;
        if(fillEntry(i,j+1,rs,rq,newDot)) return true;
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}

// enumerate non-decreasing diagonals (vals {10..26} even) summing to 162, run DFS for each
int diagVals[9] = {10,12,14,16,18,20,22,24,26};
int dchoice[N];
void recDiag(int pos, int prevIdx, int ssum){
    if(found) return;
    if(pos==N){
        if(ssum==162){
            diagsTotal++;
        }
        return;
    }
    for(int idx=prevIdx; idx<9; idx++){
        int v=diagVals[idx];
        int ns=ssum+v;
        int rem=N-pos-1;
        if(ns + rem*10 > 162) { /* min future too big? actually monotone */ }
        if(ns + rem*v > 162) break;       // even smallest future (>=v) overshoots
        if(ns + rem*26 < 162) continue;   // largest future undershoots
        recDiag(pos+1, idx, ns);
    }
}
// second pass: actually run
void runDiag(int pos, int prevIdx, int ssum){
    if(found) return;
    if(pos==N){
        if(ssum!=162) return;
        for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=dchoice[i]; }
        diagStartNodes = nodes; abortDiag = false;
        if(fillRow(0)){ found=true; }
        if(abortDiag && !found){
            for(int k=0;k<N;k++) unresolved[nUnres][k]=dchoice[k];
            nUnres++;
        }
        diagsDone++;
        if(diagsDone%50==0 || found){
            printf("  diag %lld/%lld done, nodes=%lld, unresolved=%d, %.1fs%s\n",
                   diagsDone, diagsTotal, nodes, nUnres, (double)(time(nullptr)-t0), found?"  <-- FOUND":"");
            fflush(stdout);
        }
        return;
    }
    for(int idx=prevIdx; idx<9; idx++){
        int v=diagVals[idx]; int ns=ssum+v; int rem=N-pos-1;
        if(ns + rem*v > 162) break;
        if(ns + rem*26 < 162) continue;
        dchoice[pos]=v;
        runDiag(pos+1, idx, ns);
        if(found) return;
    }
}

int main(){
    t0=time(nullptr);
    recDiag(0,0,0);
    printf("candidate diagonals: %lld\n", diagsTotal); fflush(stdout);
    runDiag(0,0,0);
    printf("\n=== DONE ===  nodes=%lld  walltime=%.1fs\n", nodes, (double)(time(nullptr)-t0));
    if(found){
        printf("*** ORBIT MATRIX FOUND ***\n");
        for(int i=0;i<N;i++){ for(int j=0;j<N;j++) printf("%3d", R[i][j]); printf("\n"); }
    } else if(nUnres==0){
        printf("*** SEARCH EXHAUSTED: NO orbit matrix exists. ***\n");
        printf("    => srg(333,166,82,83) has NO fixed-point-free order-37 automorphism (rigorous).\n");
    } else {
        printf("*** %d diagonals hit the node cap (UNRESOLVED); all others exhausted with no solution. ***\n", nUnres);
        printf("Unresolved diagonals (resolve separately with a complete solver):\n");
        for(int u=0;u<nUnres;u++){ printf("  ["); for(int k=0;k<N;k++) printf("%d%s",unresolved[u][k],k<N-1?",":""); printf("]\n"); }
    }
    return 0;
}
