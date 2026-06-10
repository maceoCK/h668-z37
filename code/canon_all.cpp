#include <cstdio>
#include <algorithm>
#include <vector>
#include <set>
#include <array>
using namespace std;
static array<int,81> canon(const array<int,81>& R){
    int perm[9]={0,1,2,3,4,5,6,7,8};
    array<int,81> best; best.fill(99);
    do{
        array<int,81> cur; bool less=false, greater=false;
        for(int i=0;i<9&&!greater;i++) for(int j=0;j<9;j++){
            int v=R[perm[i]*9+perm[j]]; cur[i*9+j]=v;
            if(!less){ if(v<best[i*9+j]) less=true; else if(v>best[i*9+j]){greater=true;break;} }
        }
        if(!greater) best=cur;
    } while(next_permutation(perm,perm+9));
    return best;
}
int main(){
    FILE* f=fopen("orbit_all_sym.txt","r");
    vector<array<int,81>> a; array<int,81> r; int t=0,v;
    while(fscanf(f,"%d",&v)==1){ r[t++]=v; if(t==81){a.push_back(r); t=0;} }
    fclose(f);
    set<array<int,81>> c;
    int done=0;
    for(auto& m:a){ c.insert(canon(m)); if(++done%500==0){printf("  %d/%zu classes=%zu\n",done,a.size(),c.size());fflush(stdout);} }
    printf("TOTAL: %zu matrices -> %zu isomorphism classes\n", a.size(), c.size());
    FILE* g=fopen("orbit_all_classes.txt","w");
    for(auto& m:c){ for(int x:m) fprintf(g,"%d ",x); fprintf(g,"\n"); }
    fclose(g);
    return 0;
}
