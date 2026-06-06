// Enumerate ALL 9x9 orbit matrices R for a fixed-point-free order-37 automorphism of
// srg(333,166,82,83) (sorted-diagonal symmetry break). Same verified DFS as conference_orbit.cpp,
// but records every solution instead of stopping. Caps: per-diagonal node cap + total-found cap.
// Output: count + sample + writes skeletons to skeletons.txt for the Modal lift fan-out.
#include <cstdio>
#include <cstdint>
#include <ctime>
using namespace std;

static const int N = 9, ROWSUM = 166, MAXENT = 37;
int R[N][N];
long long nodes = 0, diagsDone = 0, diagsTotal = 0, nFound = 0;
const long long CAP = 4000000000LL;      // per-diagonal node cap (raised to finish hard diagonals)
const long long TOTALCAP = 1000000000LL; // effectively no total cap
long long diagStartNodes = 0;
bool abortDiag = false, stopAll = false;
int unresolved[4000][N]; int nUnres = 0;
FILE* fout = nullptr;                     // stream every skeleton to disk (unbounded count)
time_t t0;
inline int sq(int x){ return x*x; }
bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]);

void record(){
    nFound++;
    if(fout){ for(int a=0;a<N;a++) for(int b=0;b<N;b++) fprintf(fout,"%d ",R[a][b]); fprintf(fout,"\n"); }
    if(nFound >= TOTALCAP) stopAll = true;
}

bool fillRow(int i){
    if(stopAll) return false;
    if(abortDiag) return false;
    if(i==N){ record(); return false; }   // complete orbit matrix: record, keep enumerating
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
    if(stopAll) return false;
    if(abortDiag) return false;
    if((nodes - diagStartNodes) > CAP){ abortDiag = true; return false; }
    if(j==N){
        if(remSum!=0 || remSumsq!=0) return false;
        for(int jp=0;jp<i;jp++) if(remDot[jp]!=0) return false;
        return fillRow(i+1);
    }
    nodes++;
    int slotsAfter = N-1-j;
    for(int v=0; v<=MAXENT; v++){
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
        if(stopAll) return false;
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}

int diagVals[9] = {10,12,14,16,18,20,22,24,26};
int dchoice[N];
void recDiag(int pos, int prev, int ssum){
    if(pos==N){ if(ssum==162) diagsTotal++; return; }
    for(int idx=prev; idx<9; idx++){ int v=diagVals[idx], ns=ssum+v, rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*26<162) continue; recDiag(pos+1,idx,ns); }
}
void runDiag(int pos, int prev, int ssum){
    if(stopAll) return;
    if(pos==N){
        if(ssum!=162) return;
        for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=dchoice[i]; }
        diagStartNodes = nodes; abortDiag = false;
        long long before = nFound;
        fillRow(0);
        if(abortDiag){ for(int k=0;k<N;k++) unresolved[nUnres][k]=dchoice[k]; nUnres++; }
        diagsDone++;
        if(diagsDone%50==0 || nFound>before){
            printf("  diag %lld/%lld  found=%lld (+%lld here)  unresolved=%d  nodes=%lld  %.1fs\n",
                   diagsDone, diagsTotal, nFound, nFound-before, nUnres, nodes,(double)(time(nullptr)-t0));
            fflush(stdout);
        }
        return;
    }
    for(int idx=prev; idx<9; idx++){ int v=diagVals[idx], ns=ssum+v, rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*26<162) continue; dchoice[pos]=v; runDiag(pos+1,idx,ns);
        if(stopAll) return; }
}

int main(){
    t0=time(nullptr);
    recDiag(0,0,0);
    printf("candidate diagonals: %lld\n", diagsTotal); fflush(stdout);
    fout = fopen("skeletons_full.txt","w");
    runDiag(0,0,0);
    if(fout) fclose(fout);
    printf("\n=== DONE ===  total orbit matrices found=%lld  unresolved diagonals=%d (of %lld)  nodes=%lld  %.1fs\n",
           nFound, nUnres, diagsTotal, nodes,(double)(time(nullptr)-t0));
    if(nUnres==0) printf("ENUMERATION COMPLETE: all %lld diagonals resolved; %lld orbit matrices total.\n", diagsTotal, nFound);
    else { printf("UNRESOLVED diagonals (still hit the %lld-node cap):\n", CAP);
        for(int u=0;u<nUnres && u<60;u++){ printf("  ["); for(int k=0;k<N;k++) printf("%d%s",unresolved[u][k],k<N-1?",":""); printf("]\n"); } }
    return 0;
}
