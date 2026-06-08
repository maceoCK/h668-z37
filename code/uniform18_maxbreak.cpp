// INDEPENDENT re-check of "no all-18-diagonal orbit matrix exists", using a DIFFERENT sound
// symmetry break than uniform18_enum.cpp: require R[0][1] = global maximum off-diagonal entry
// (orbits interchangeable under the uniform diagonal => relabel so a max-pair is (0,1), and every
//  other off-diagonal entry <= R[0][1]). If this also yields 0, the a in {3,6} elimination is confirmed
//  by two independent reductions.
#include <cstdio>
#include <ctime>
using namespace std;
static const int N=9, ROWSUM=166, MAXENT=37;
int R[N][N]; long long nodes=0, nFound=0; int GMAX=37; time_t t0;
inline int sq(int x){ return x*x; }
bool fillEntry(int i,int j,int remSum,int remSumsq,int remDot[N]);
bool fillRow(int i){
    if(i==N){ nFound++; return false; }
    int bS=R[i][i], bQ=sq(R[i][i]);
    for(int k=0;k<i;k++){ bS+=R[i][k]; bQ+=sq(R[i][k]); }
    int rS=ROWSUM-bS, rQ=(3154-R[i][i])-bQ; if(rS<0||rQ<0) return false;
    int rD[N];
    for(int jp=0;jp<i;jp++){ int b=R[i][i]*R[jp][i]; for(int k=0;k<i;k++) b+=R[i][k]*R[jp][k]; rD[jp]=(3071-R[i][jp])-b; }
    return fillEntry(i,i+1,rS,rQ,rD);
}
inline void knap(const int* w,int c,int bud,int& lo,int& hi){
    int o[N]; for(int t=0;t<c;t++) o[t]=w[t];
    for(int a=1;a<c;a++){int x=o[a],b=a-1;while(b>=0&&o[b]>x){o[b+1]=o[b];b--;}o[b+1]=x;}
    int b=bud,mn=0; for(int t=0;t<c&&b>0;t++){int tk=b<37?b:37;mn+=tk*o[t];b-=tk;}
    b=bud;int mx=0; for(int t=c-1;t>=0&&b>0;t--){int tk=b<37?b:37;mx+=tk*o[t];b-=tk;}
    lo=mn;hi=mx;
}
bool fillEntry(int i,int j,int remSum,int remSumsq,int remDot[N]){
    if(j==N){ if(remSum||remSumsq) return false; for(int jp=0;jp<i;jp++) if(remDot[jp]) return false; return fillRow(i+1); }
    nodes++;
    if((nodes & ((1LL<<31)-1))==0){ printf("  nodes=%lld found=%lld %.0fs\n",nodes,nFound,(double)(time(nullptr)-t0)); fflush(stdout); }
    int slots=N-1-j;
    int vhi = (i==0&&j==1) ? MAXENT : GMAX;   // every off-diag <= R[0][1] (the global max); (0,1) sets it
    for(int v=0; v<=vhi; v++){
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
        R[i][j]=v; R[j][i]=v;
        int save=GMAX; if(i==0&&j==1) GMAX=v;   // fix the global max
        fillEntry(i,j+1,rs,rq,nD);
        GMAX=save;
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}
int main(){
    t0=time(nullptr);
    for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=18; }
    fillRow(0);
    printf("=== DONE (independent global-max break) === all-18 orbit matrices = %lld  nodes=%lld  %.1fs\n",
           nFound, nodes,(double)(time(nullptr)-t0));
    printf(nFound==0 ? "CONFIRMED: none exist (a in {3,6} eliminated, 2nd independent reduction).\n"
                     : "DISCREPANCY: %lld found -- investigate!\n", nFound);
    return 0;
}
