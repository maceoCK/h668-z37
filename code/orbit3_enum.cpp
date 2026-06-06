// COMPLETE decision: does an order-3 coset-compatible 9x9 orbit matrix exist for a
// fixed-point-free order-37 automorphism of srg(333,166,82,83)?
//
// An H-cyclotomic graph with |H|=3 (H = {1,10,26} <= Z_37*, the order-3 multiplier
// subgroup, which is NOT closed under negation) forces its orbit matrix R to satisfy:
//   - diagonal R[i][i] ≡ 0 (mod 6)   (diagonal blocks are symmetric => unions of <H,-1>-cosets)
//   - off-diagonal R[i][j] ≡ 0 or 1 (mod 3)   (off-diag blocks are H-invariant, +1 iff 0 in block)
// together with the orbit-matrix identity R^2 + R - 83 I = 3071 J.
//
// This is the SAME verified DFS as conference_orbit_all.cpp, but (a) restricts the diagonal
// alphabet to {12,18,24} and the off-diagonal alphabet to v with v%3 in {0,1}, and (b) removes
// the per-diagonal node cap so the search runs to completion. If NO such R exists, the order-3
// case of Theorem 1 is CLOSED completely (no enumeration sampling), strengthening it to 7 of 8.
#include <cstdio>
#include <cstdint>
#include <ctime>
using namespace std;

static const int N = 9, ROWSUM = 166, MAXENT = 37;
int R[N][N];
long long nodes = 0, diagsDone = 0, diagsTotal = 0, nFound = 0;
bool stopAll = false;
FILE* fout = nullptr;
time_t t0;
inline int sq(int x){ return x*x; }
inline bool offok(int v){ int r=v%3; return r==0 || r==1; }   // order-3 coset-compatible off-diag
bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]);

void record(){
    nFound++;
    if(fout){ for(int a=0;a<N;a++) for(int b=0;b<N;b++) fprintf(fout,"%d ",R[a][b]); fprintf(fout,"\n"); }
}

bool fillRow(int i){
    if(i==N){ record(); return false; }
    int baseSum = R[i][i], baseSumsq = sq(R[i][i]);
    for(int k=0;k<i;k++){ baseSum += R[i][k]; baseSumsq += sq(R[i][k]); }
    int remSum = ROWSUM - baseSum, remSumsq = (3154 - R[i][i]) - baseSumsq;
    if(remSum<0 || remSumsq<0) return false;
    int remDot[N];
    for(int jp=0;jp<i;jp++){
        int base = R[i][i]*R[jp][i];
        for(int k=0;k<i;k++) base += R[i][k]*R[jp][k];
        remDot[jp] = (3071 - R[i][jp]) - base;
    }
    return fillEntry(i, i+1, remSum, remSumsq, remDot);
}

inline void knap_bounds(const int* w, int cnt, int budget, int& lo, int& hi){
    int ord[N]; for(int t=0;t<cnt;t++) ord[t]=w[t];
    for(int a=1;a<cnt;a++){ int x=ord[a],b=a-1; while(b>=0&&ord[b]>x){ord[b+1]=ord[b];b--;} ord[b+1]=x; }
    int b=budget, mn=0; for(int t=0;t<cnt && b>0;t++){ int take=b<37?b:37; mn+=take*ord[t]; b-=take; }
    b=budget; int mx=0; for(int t=cnt-1;t>=0 && b>0;t--){ int take=b<37?b:37; mx+=take*ord[t]; b-=take; }
    lo=mn; hi=mx;
}

bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]){
    if(j==N){
        if(remSum!=0 || remSumsq!=0) return false;
        for(int jp=0;jp<i;jp++) if(remDot[jp]!=0) return false;
        return fillRow(i+1);
    }
    nodes++;
    int slotsAfter = N-1-j;
    for(int v=0; v<=MAXENT; v++){
        if(!offok(v)) continue;                       // <-- order-3 off-diagonal restriction
        int rs = remSum - v;
        if(rs<0) break;
        if(rs > MAXENT*slotsAfter) continue;
        int rq = remSumsq - sq(v);
        if(rq<0) continue;
        if(slotsAfter>0){
            if((long long)rq*slotsAfter < (long long)rs*rs) continue;
            int f37 = rs/37, rem = rs%37;
            long long qmax = (long long)f37*1369 + (long long)rem*rem;
            if((long long)rq > qmax) continue;
        } else { if(rs!=0 || rq!=0) continue; }
        bool okdot=true; int newDot[N];
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
        fillEntry(i,j+1,rs,rq,newDot);
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}

// order-3 diagonal alphabet: even values in {10..26} that are 0 mod 6.
int diagVals[3] = {12,18,24};
int dchoice[N];
void recDiag(int pos, int prev, int ssum){
    if(pos==N){ if(ssum==162) diagsTotal++; return; }
    for(int idx=prev; idx<3; idx++){ int v=diagVals[idx], ns=ssum+v, rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*24<162) continue; recDiag(pos+1,idx,ns); }
}
void runDiag(int pos, int prev, int ssum){
    if(pos==N){
        if(ssum!=162) return;
        for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=dchoice[i]; }
        long long before = nFound;
        fillRow(0);
        diagsDone++;
        printf("  diag %lld/%lld  [", diagsDone, diagsTotal);
        for(int k=0;k<N;k++) printf("%d%s",dchoice[k],k<N-1?",":"");
        printf("]  found=%lld (+%lld)  nodes=%lld  %.1fs\n",
               nFound, nFound-before, nodes,(double)(time(nullptr)-t0));
        fflush(stdout);
        return;
    }
    for(int idx=prev; idx<3; idx++){ int v=diagVals[idx], ns=ssum+v, rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*24<162) continue; dchoice[pos]=v; runDiag(pos+1,idx,ns); }
}

int main(){
    t0=time(nullptr);
    recDiag(0,0,0);
    printf("order-3 coset-compatible diagonals (values in {12,18,24}, sum 162): %lld\n", diagsTotal);
    fflush(stdout);
    fout = fopen("orbit3_compatible.txt","w");
    runDiag(0,0,0);
    if(fout) fclose(fout);
    printf("\n=== DONE (complete; no node cap) ===  order-3-compatible orbit matrices found = %lld  nodes=%lld  %.1fs\n",
           nFound, nodes,(double)(time(nullptr)-t0));
    if(nFound==0)
        printf("RESULT: NO order-3 coset-compatible orbit matrix exists => order-3 cyclotomic case CLOSED (7 of 8).\n");
    else
        printf("RESULT: %lld order-3-compatible skeleton(s) exist (written to orbit3_compatible.txt); lift-sweep needed.\n", nFound);
    return 0;
}
