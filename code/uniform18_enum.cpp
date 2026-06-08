// Decide: does an admissible 9x9 orbit matrix R with diagonal ALL 18 exist?
// (a in {3,6} of the order-2-cyclotomic case forces D_ii = Q^-/Q^+ for every orbit i,
//  hence R[i][i]=18 for all i.  If NO such R exists, a in {3,6} is rigorously eliminated
//  and the open order-2 case narrows to a in {4,5}.)
// Same verified DFS as conference_orbit_all.cpp; diagonal fixed at (18,...,18), off-diagonal
// alphabet UNRESTRICTED (0..37), NO node cap -> runs to completion (definitive).
#include <cstdio>
#include <cstdint>
#include <ctime>
using namespace std;
static const int N = 9, ROWSUM = 166, MAXENT = 37;
int R[N][N];
long long nodes = 0, nFound = 0;
FILE* fout = nullptr;
time_t t0;
inline int sq(int x){ return x*x; }
bool fillEntry(int i, int j, int remSum, int remSumsq, int remDot[N]);

void record(){
    nFound++;
    if(fout){ for(int a=0;a<N;a++) for(int b=0;b<N;b++) fprintf(fout,"%d ",R[a][b]); fprintf(fout,"\n"); }
}
bool fillRow(int i){
    if(i==N){ record(); return false; }
    int baseSum=R[i][i], baseSumsq=sq(R[i][i]);
    for(int k=0;k<i;k++){ baseSum+=R[i][k]; baseSumsq+=sq(R[i][k]); }
    int remSum=ROWSUM-baseSum, remSumsq=(3154-R[i][i])-baseSumsq;
    if(remSum<0||remSumsq<0) return false;
    int remDot[N];
    for(int jp=0;jp<i;jp++){
        int base=R[i][i]*R[jp][i];
        for(int k=0;k<i;k++) base+=R[i][k]*R[jp][k];
        remDot[jp]=(3071-R[i][jp])-base;
    }
    return fillEntry(i,i+1,remSum,remSumsq,remDot);
}
inline void knap_bounds(const int* w,int cnt,int budget,int& lo,int& hi){
    int ord[N]; for(int t=0;t<cnt;t++) ord[t]=w[t];
    for(int a=1;a<cnt;a++){int x=ord[a],b=a-1;while(b>=0&&ord[b]>x){ord[b+1]=ord[b];b--;}ord[b+1]=x;}
    int b=budget,mn=0; for(int t=0;t<cnt&&b>0;t++){int tk=b<37?b:37;mn+=tk*ord[t];b-=tk;}
    b=budget;int mx=0; for(int t=cnt-1;t>=0&&b>0;t--){int tk=b<37?b:37;mx+=tk*ord[t];b-=tk;}
    lo=mn;hi=mx;
}
bool fillEntry(int i,int j,int remSum,int remSumsq,int remDot[N]){
    if(j==N){
        if(remSum!=0||remSumsq!=0) return false;
        for(int jp=0;jp<i;jp++) if(remDot[jp]!=0) return false;
        return fillRow(i+1);
    }
    nodes++;
    if((nodes & ((1LL<<31)-1))==0){ printf("  nodes=%lld found=%lld %.0fs\n",nodes,nFound,(double)(time(nullptr)-t0)); fflush(stdout); }
    int slotsAfter=N-1-j;
    // sound symmetry break: orbits 1..8 are interchangeable (uniform diagonal) -> require row 0
    // off-diagonal non-decreasing (relabel orbits by increasing connection to orbit 0).
    int vstart = (i==0 && j>=2) ? R[0][j-1] : 0;
    for(int v=vstart; v<=MAXENT; v++){
        int rs=remSum-v; if(rs<0) break;
        if(rs>MAXENT*slotsAfter) continue;
        int rq=remSumsq-sq(v); if(rq<0) continue;
        if(slotsAfter>0){
            if((long long)rq*slotsAfter<(long long)rs*rs) continue;
            int f37=rs/37,rem=rs%37; long long qmax=(long long)f37*1369+(long long)rem*rem;
            if((long long)rq>qmax) continue;
        } else { if(rs!=0||rq!=0) continue; }
        bool okdot=true; int newDot[N];
        for(int jp=0;jp<i;jp++){
            int nd=remDot[jp]-v*R[jp][j];
            if(nd<0){okdot=false;break;}
            if(slotsAfter>0){
                int w[N]; for(int t=0;t<slotsAfter;t++) w[t]=R[jp][j+1+t];
                int lo,hi; knap_bounds(w,slotsAfter,rs,lo,hi);
                if(nd<lo||nd>hi){okdot=false;break;}
            } else if(nd!=0){okdot=false;break;}
            newDot[jp]=nd;
        }
        if(!okdot) continue;
        R[i][j]=v; R[j][i]=v;
        fillEntry(i,j+1,rs,rq,newDot);
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}
int main(){
    t0=time(nullptr);
    for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=18; }   // diagonal ALL 18
    fout=fopen("uniform18.txt","w");
    fillRow(0);
    if(fout) fclose(fout);
    printf("=== DONE (complete; no node cap) === uniform-18-diagonal orbit matrices found = %lld  nodes=%lld  %.1fs\n",
           nFound, nodes, (double)(time(nullptr)-t0));
    if(nFound==0) printf("RESULT: NO all-18-diagonal orbit matrix exists => a in {3,6} ELIMINATED => open order-2 narrows to a in {4,5}.\n");
    else printf("RESULT: %lld all-18-diagonal orbit matrix/matrices exist => a in {3,6} NOT eliminated; lift-check needed.\n", nFound);
    return 0;
}
