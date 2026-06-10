// Consistency check: the 960 order-3 coset-compatible orbit matrices found by the independent
// mod-3-restricted enumeration (orbit3_compatible.txt) and the 28 order-3-compatible representatives
// of the complete enumeration (order3_skeletons.txt) must coincide up to simultaneous row/column
// permutation (orbit relabeling).  Canonical form = lexicographic minimum of the flattened matrix
// over all 9! permutations (exact, no heuristics).
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <set>
#include <array>
using namespace std;

static array<int,81> canon(const array<int,81>& R){
    int perm[9] = {0,1,2,3,4,5,6,7,8};
    array<int,81> best;
    best.fill(99);
    do {
        array<int,81> cur;
        bool less = false, greater = false;
        for (int i = 0; i < 9 && !greater; i++)
            for (int j = 0; j < 9; j++){
                int v = R[perm[i]*9 + perm[j]];
                cur[i*9+j] = v;
                if (!less){
                    if (v < best[i*9+j]) less = true;
                    else if (v > best[i*9+j]){ greater = true; break; }
                }
            }
        if (!greater && less) best = cur;
        else if (!greater && !less) best = cur;     // equal or first
    } while (next_permutation(perm, perm + 9));
    return best;
}

static vector<array<int,81>> load(const char* path){
    vector<array<int,81>> out;
    FILE* f = fopen(path, "r");
    if (!f){ fprintf(stderr, "cannot open %s\n", path); exit(1); }
    array<int,81> r; int t = 0, v;
    while (fscanf(f, "%d", &v) == 1){
        r[t++] = v;
        if (t == 81){ out.push_back(r); t = 0; }
    }
    fclose(f);
    return out;
}

int main(){
    auto a = load("orbit3_compatible.txt");      // 960 from the independent enumeration
    auto b = load("order3_skeletons.txt");       // 28 representatives of the complete enumeration
    printf("loaded %zu and %zu matrices\n", a.size(), b.size());
    set<array<int,81>> ca, cb;
    for (auto& r : b) cb.insert(canon(r));
    printf("28-set canonical classes: %zu\n", cb.size());
    int done = 0;
    for (auto& r : a){
        ca.insert(canon(r));
        if (++done % 100 == 0){ printf("  ... %d/%zu, classes=%zu\n", done, a.size(), ca.size()); fflush(stdout); }
    }
    printf("960-set canonical classes: %zu\n", ca.size());
    bool sub1 = includes(cb.begin(), cb.end(), ca.begin(), ca.end());
    bool sub2 = includes(ca.begin(), ca.end(), cb.begin(), cb.end());
    printf("960-classes subset of 28-classes: %s ; 28 subset of 960: %s ; EQUAL: %s\n",
           sub1?"yes":"no", sub2?"yes":"no", (ca==cb)?"YES":"NO");
    return 0;
}
