// Verify: which of the 625 canonical orbit-matrix classes are SELF-COMPLEMENTARY
// (complement R' (off-diag 37-r, diag 36-r) is a relabeling of R)?  For each such class,
// report the complementing permutations and their cycle types.
#include <cstdio>
#include <algorithm>
#include <vector>
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
    FILE* f=fopen("orbit_all_classes.txt","r");
    vector<array<int,81>> a; array<int,81> r; int t=0,v;
    while(fscanf(f,"%d",&v)==1){ r[t++]=v; if(t==81){a.push_back(r); t=0;} }
    fclose(f);
    printf("loaded %zu classes\n", a.size());
    int found=0;
    for(size_t k=0;k<a.size();k++){
        array<int,81> C;
        for(int i=0;i<9;i++) for(int j=0;j<9;j++)
            C[i*9+j] = (i==j) ? 36-a[k][i*9+j] : 37-a[k][i*9+j];
        if(canon(C)==canon(a[k])){
            found++;
            printf("class %zu SELF-COMPLEMENTARY: diag", k);
            for(int i=0;i<9;i++) printf(" %d", a[k][i*9+i]);
            // enumerate complementing perms: P with C = P R P^T
            printf("  perms:");
            int perm[9]={0,1,2,3,4,5,6,7,8}; int np=0;
            do{
                bool ok=true;
                for(int i=0;i<9&&ok;i++) for(int j=0;j<9;j++)
                    if(C[i*9+j]!=a[k][perm[i]*9+perm[j]]){ok=false;break;}
                if(ok){
                    np++;
                    // cycle type
                    bool seen[9]={false}; printf(" (");
                    for(int s=0;s<9;s++){ if(seen[s])continue; int len=0,x=s; while(!seen[x]){seen[x]=true;x=perm[x];len++;} printf("%d", len); }
                    printf(")");
                    // fixed points' diagonal entries
                }
            } while(next_permutation(perm,perm+9));
            printf("  nperms=%d\n", np);
        }
    }
    printf("TOTAL self-complementary classes: %d / %zu\n", found, a.size());
    return 0;
}
