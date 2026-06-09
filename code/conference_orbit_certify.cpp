// CERTIFICATE for Theorem 1(a): for a multiplier subgroup H of order e, an H-cyclotomic srg(333,166,82,83)
// forces its 9x9 orbit matrix to have H-coset-compatible cardinalities -- off-diagonal R[i][j] = 0,1 (mod e),
// diagonal R[i][i] = 0 (mod e) if -1 in H (e even) else (mod 2e). This program EXHAUSTIVELY enumerates all
// orbit matrices over that restricted alphabet (sound symmetry break, no node cap) and reports the count.
// For e in {4,6,9,12,18,36} the count is 0 -- an auditable complete-search certificate of nonexistence,
// independent of any SAT solver. (For e=3 it returns 960, the substrate of the open |H|=3 case.)
//   build:  clang++ -O3 -o cert conference_orbit_certify.cpp
//   run:    for e in 4 6 9 12 18 36 3; do ./cert $e; done
#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;
static const int N=9, ROWSUM=166, MAXENT=37;
int R[N][N]; long long nodes=0, nFound=0; int dchoice[N];
int E, DIAGMOD, NDIAG, DIAGV[16];
bool abortRun=false;
inline int sq(int x){ return x*x; }
inline bool offok(int v){ int r=v%E; return r==0||r==1; }
bool fillEntry(int i,int j,int rs2,int rq2,int rd[N]);
bool fillRow(int i){
    if(i==N){ nFound++; return false; }
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
    if(j==N){
        if(remSum||remSumsq) return false;
        for(int jp=0;jp<i;jp++) if(remDot[jp]) return false;
        return fillRow(i+1);
    }
    nodes++;
    int slots=N-1-j;
    int vstart=(i==0 && j>=2 && dchoice[j]==dchoice[j-1]) ? R[0][j-1] : 0;
    for(int v=vstart; v<=MAXENT; v++){
        if(!offok(v)) continue;                       // coset-compatible off-diagonal alphabet
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
        fillEntry(i,j+1,rs,rq,nD);
    }
    R[i][j]=0; R[j][i]=0;
    return false;
}
void runDiag(int pos,int prev,int ssum){
    if(pos==N){
        if(ssum!=162) return;
        for(int i=0;i<N;i++){ for(int k=0;k<N;k++) R[i][k]=0; R[i][i]=dchoice[i]; }
        fillRow(0);
        return;
    }
    for(int idx=prev;idx<NDIAG;idx++){int v=DIAGV[idx],ns=ssum+v,rem=N-pos-1;
        if(ns+rem*v>162) break; if(ns+rem*DIAGV[NDIAG-1]<162) continue; dchoice[pos]=v; runDiag(pos+1,idx,ns);}
}
int main(int argc,char**argv){
    E=atoi(argv[1]);
    DIAGMOD = (E%2==0)? E : 2*E;               // -1 in H iff e even
    NDIAG=0; for(int v=10;v<=26;v+=2) if(v%DIAGMOD==0) DIAGV[NDIAG++]=v;
    time_t t0=time(nullptr);
    if(NDIAG==0){ printf("e=%2d: diag alphabet empty (no even diag in {10..26} = 0 mod %d) -> 0 orbit matrices. CERTIFIED.\n",E,DIAGMOD); return 0; }
    runDiag(0,0,0);
    printf("e=%2d: off-diag = 0,1 mod %d; diag in {", E, E);
    for(int k=0;k<NDIAG;k++) printf("%d%s",DIAGV[k],k<NDIAG-1?",":"");
    printf("}; nodes=%lld %.1fs -> coset-compatible orbit matrices = %lld%s\n",
        nodes,(double)(time(nullptr)-t0),nFound,
        nFound==0?"  => NO H-cyclotomic srg (Thm 1(a)) CERTIFIED":"  (exist; |H|=3 substrate / open)");
    return 0;
}
