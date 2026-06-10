// sc_anneal: simulated annealing over the sigma-equivariant free space of a self-complementary
// Z37-symmetric srg(333,166,82,83) candidate (classes 221/422), SAT-side lane.
//
// Consumes the literal-resolution template emitted by sc_emit_template.py (which reuses the
// hostile-validated chain logic of selfcomp_search.py).  State = free representative blocks:
//   type 0 (off-diag): any subset of Z37 of fixed cardinality;
//   type 1 (cycle diagonal): symmetric subset (union of +- pairs), fixed cardinality, no 0;
//   type 2 (fixed diagonal): anticlosed under c (x(cd) = NOT x(d)): one binary choice per
//          <c>-orbit {x,cx,-x,-cx}: take {x,-x} or {cx,-cx}; cardinality 18 automatic.
// Membership of any (a,b,d) is rep-bit XOR parity (template).  Objective:
//   sum over ordered pairs a<=b that are ORBIT REPS and g in Z37 of
//     ( sum_l #conv(a,l,g) - [mu + (k-mu)[a==b,g==0] + (lam-mu) x_{ab}(g)] )^2
// (equation orbits are sigma-equivalent, so reps suffice; full objective optional via --all-pairs).
// Objective 0 <=> all SRG equations hold <=> solution (verify externally).
//
// Moves: type 0: swap one element in/out; type 1: swap a +- pair in/out; type 2: flip one orbit
// choice.  Incremental conv updates.  Metropolis, geometric cooling with reheats, multi-restart.
//
// Usage: sc_anneal template.txt seed maxSec [--all-pairs] [--check]
// Output: progress lines "BEST obj=..."; on obj==0: "SOLUTION" + 81 block lines, exit 0.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <random>
using namespace std;
typedef long long ll;

static int P, NORB, NREP, KP, LAM, MU, CC;
struct Rep { int i, j, type, card; };
static vector<Rep> reps;
static int litRep[9][9][37], litE[9][9][37], litNeg[9][9][37]; // repIdx (-1 const), element, parity/const
static bool isRepPair[9][9];

static unsigned char repBit[16][37];       // current state of free reps
static unsigned char mem[9][9][37];        // expanded membership
static int conv[9][9][37];                 // conv(a,b,g) = sum_l sum_d mem[a][l][d] mem[l][b][g-d]
static ll obj = 0;
static bool allPairs = false;

static inline int rhsOf(int a, int b, int g){ return MU + ((a == b && g == 0) ? (KP - MU) : 0); }

static inline ll term(int a, int b, int g){
    int lhs = conv[a][b][g] - (LAM - MU) * mem[a][b][g];
    int d = lhs - rhsOf(a, b, g);
    return (ll)d * d;
}

static void expandAll(){
    for (int a = 0; a < NORB; a++)
        for (int b = 0; b < NORB; b++)
            for (int d = 0; d < P; d++){
                int r = litRep[a][b][d];
                if (r < 0) mem[a][b][d] = (unsigned char)litNeg[a][b][d];
                else mem[a][b][d] = (unsigned char)(repBit[r][litE[a][b][d]] ^ litNeg[a][b][d]);
            }
}
static void convAll(){
    for (int a = 0; a < NORB; a++)
        for (int b = 0; b < NORB; b++)
            for (int g = 0; g < P; g++){
                int s = 0;
                for (int l = 0; l < NORB; l++)
                    for (int d = 0; d < P; d++)
                        if (mem[a][l][d] && mem[l][b][(g - d + P) % P]) s++;
                conv[a][b][g] = s;
            }
}
static ll objAll(){
    ll s = 0;
    for (int a = 0; a < NORB; a++)
        for (int b = a; b < NORB; b++){
            if (!allPairs && !isRepPair[a][b]) continue;
            for (int g = 0; g < P; g++) s += term(a, b, g);
        }
    return s;
}

// ---- incremental machinery: flip membership of element d0 in stored direction (a0,b0).
// Affects conv(a0, b, g) for terms l=b0 (and conv(a, b0, .) via l=a0), all g.
// We flip a SET of (a,b,d) memberships atomically (a rep-bit flip touches its whole sigma-orbit).
struct Delta { int a, b, d; };
static vector<Delta> deltas;

static void collectDeltas(int r, int e){
    deltas.clear();
    for (int a = 0; a < NORB; a++)
        for (int b = 0; b < NORB; b++)
            for (int d = 0; d < P; d++)
                if (litRep[a][b][d] == r && litE[a][b][d] == e)
                    deltas.push_back({a, b, d});
}

static ll applyDeltas(bool flip){
    // recompute objective contribution of affected (pair,g) terms; update conv incrementally
    // affected pairs: any (a,b) with a changed mem in row a or column b; simplest: mark pairs
    static bool touched[9][9];
    memset(touched, 0, sizeof touched);
    ll before = 0;
    for (auto& dl : deltas){
        for (int x = 0; x < NORB; x++){
            touched[dl.a][x] = true; touched[x][dl.a] = true;
            touched[dl.b][x] = true; touched[x][dl.b] = true;
        }
    }
    for (int a = 0; a < NORB; a++)
        for (int b = a; b < NORB; b++){
            if (!allPairs && !isRepPair[a][b]) continue;
            if (!touched[a][b]) continue;
            for (int g = 0; g < P; g++) before += term(a, b, g);
        }
    if (flip){
        // apply membership flips with conv updates
        for (auto& dl : deltas){
            int a = dl.a, b = dl.b, d = dl.d;
            int sgn = mem[a][b][d] ? -1 : +1;
            // conv(a, c, g) terms with l=b: mem[a][b][d] * mem[b][c][g-d]
            for (int c2 = 0; c2 < NORB; c2++)
                for (int g2 = 0; g2 < P; g2++)
                    if (mem[b][c2][(g2 - d + P) % P]) conv[a][c2][g2] += sgn;
            // conv(c, b, g) terms with l=a: mem[c][a][g-d'] * mem[a][b][d]... conv(c,b,g): sum_l mem[c][l][x] mem[l][b][g-x]; l=a, g-x=d: x=g-d
            for (int c2 = 0; c2 < NORB; c2++)
                for (int g2 = 0; g2 < P; g2++)
                    if (mem[c2][a][(g2 - d + P) % P]) conv[c2][b][g2] += sgn;
            // diagonal self-square: the product mem[a][a][d]^2 inside conv(a,a,2d) is hit by BOTH
            // passes with the pre-flip value, giving 2*sgn*x; the true bit-square delta is
            // (x^1)^2 - x^2 = sgn*(2x) + 1, so add the missing +1 (both flip directions).
            if (a == b) conv[a][a][(2 * d) % P] += 1;
            mem[a][b][d] ^= 1;
        }
    }
    ll after = 0;
    for (int a = 0; a < NORB; a++)
        for (int b = a; b < NORB; b++){
            if (!allPairs && !isRepPair[a][b]) continue;
            if (!touched[a][b]) continue;
            for (int g = 0; g < P; g++) after += term(a, b, g);
        }
    return after - before;
}

int main(int argc, char** argv){
    if (argc < 4){ fprintf(stderr, "usage: %s template.txt seed maxSec [--all-pairs]\n", argv[0]); return 2; }
    FILE* f = fopen(argv[1], "r");
    if (!f){ fprintf(stderr, "no template\n"); return 2; }
    unsigned seed = (unsigned)atoi(argv[2]);
    double maxSec = atof(argv[3]);
    bool doCheck = false;
    for (int t = 4; t < argc; t++){
        if (!strcmp(argv[t], "--all-pairs")) allPairs = true;
        if (!strcmp(argv[t], "--check")) doCheck = true;
    }
    int nrep;
    fscanf(f, "%d %d %d %d %d %d %d", &P, &NORB, &nrep, &KP, &LAM, &MU, &CC);
    NREP = nrep;
    reps.resize(NREP);
    memset(isRepPair, 0, sizeof isRepPair);
    char tag[8];
    for (int t = 0; t < NREP; t++){
        int idx; Rep r;
        fscanf(f, "%s %d %d %d %d %d", tag, &idx, &r.i, &r.j, &r.type, &r.card);
        reps[idx] = r;
        isRepPair[r.i][r.j] = true;
    }
    for (int t = 0; t < NORB * NORB * P; t++){
        int a, b, d, ri, e, ng;
        fscanf(f, "%s %d %d %d %d %d %d", tag, &a, &b, &d, &ri, &e, &ng);
        litRep[a][b][d] = ri; litE[a][b][d] = e; litNeg[a][b][d] = ng;
    }
    fclose(f);

    mt19937 rng(seed);
    time_t t0 = time(nullptr);
    ll bestEver = -1;

    // <c>-orbit structure for type-2 reps: orbits {x, cx, -x, -cx}
    vector<vector<int>> corbs;
    {
        vector<char> seen(P, 0);
        for (int x = 1; x < P; x++){
            if (seen[x]) continue;
            int a = x, b = (int)((ll)CC * x % P), na = P - x, nb = P - b;
            seen[a] = seen[b] = seen[na] = seen[nb] = 1;
            corbs.push_back({a, b, na, nb});
        }
    }

    while (difftime(time(nullptr), t0) < maxSec){
        // ---- random init respecting structure
        for (int r = 0; r < NREP; r++){
            memset(repBit[r], 0, sizeof repBit[r]);
            Rep& rp = reps[r];
            if (rp.type == 0){
                vector<int> el(P); for (int d = 0; d < P; d++) el[d] = d;
                shuffle(el.begin(), el.end(), rng);
                for (int t = 0; t < rp.card; t++) repBit[r][el[t]] = 1;
            } else if (rp.type == 1){
                vector<int> pr; for (int d = 1; d <= P / 2; d++) pr.push_back(d);
                shuffle(pr.begin(), pr.end(), rng);
                for (int t = 0; t < rp.card / 2; t++){ repBit[r][pr[t]] = 1; repBit[r][P - pr[t]] = 1; }
            } else {
                for (auto& ob : corbs){
                    if (rng() & 1){ repBit[r][ob[0]] = 1; repBit[r][ob[2]] = 1; }
                    else          { repBit[r][ob[1]] = 1; repBit[r][ob[3]] = 1; }
                }
            }
        }
        expandAll(); convAll();
        obj = objAll();
        double T = 60.0;
        ll best = obj;
        ll stall = 0;
        for (ll it = 0; ; it++){
            if ((it & 1023) == 0){
                if (difftime(time(nullptr), t0) > maxSec) break;
                T *= 0.999;
                if (T < 0.4){ T = 25.0; }                 // reheat
            }
            if (doCheck && (it & ((1 << 18) - 1)) == (1 << 18) - 1){
                ll o1 = obj;
                int cv[9][9][37]; memcpy(cv, conv, sizeof conv);
                convAll();
                ll o2 = objAll();
                if (o1 != o2 || memcmp(cv, conv, sizeof conv)){
                    fprintf(stderr, "DRIFT it=%lld obj %lld vs %lld conv %s\n", it, o1, o2,
                            memcmp(cv, conv, sizeof conv) ? "MISMATCH" : "ok");
                    return 3;
                }
            }
            // ---- propose move
            int r = (int)(rng() % NREP);
            Rep& rp = reps[r];
            int e1 = -1, e2 = -1, e3 = -1, e4 = -1;       // bits to flip in rep r
            if (rp.type == 0){
                if (rp.card == 0 || rp.card == P) continue;
                int in, out;
                do { in = (int)(rng() % P); } while (!repBit[r][in]);
                do { out = (int)(rng() % P); } while (repBit[r][out]);
                e1 = in; e2 = out;
            } else if (rp.type == 1){
                if (rp.card == 0) continue;
                int din, dout;
                do { din = 1 + (int)(rng() % (P - 1)); } while (!repBit[r][din]);
                do { dout = 1 + (int)(rng() % (P - 1)); } while (repBit[r][dout] || dout == P - din);
                e1 = din; e2 = P - din; e3 = dout; e4 = P - dout;
                if (e2 == e1) e2 = -1;
            } else {
                auto& ob = corbs[rng() % corbs.size()];
                e1 = ob[0]; e2 = ob[1]; e3 = ob[2]; e4 = ob[3]; // flip all 4 = swap choice
            }
            // ---- evaluate by full-affected recompute (collect deltas for all flipped bits)
            ll dObj = 0;
            int flips[4] = {e1, e2, e3, e4};
            for (int t = 0; t < 4; t++){
                if (flips[t] < 0) continue;
                collectDeltas(r, flips[t]);
                dObj += applyDeltas(true);
                repBit[r][flips[t]] ^= 1;
            }
            bool accept = dObj <= 0 || (double)rng() / rng.max() < exp(-(double)dObj / T);
            if (!accept){
                for (int t = 3; t >= 0; t--){
                    if (flips[t] < 0) continue;
                    collectDeltas(r, flips[t]);
                    applyDeltas(true);
                    repBit[r][flips[t]] ^= 1;
                }
            } else {
                obj += dObj;
                if (obj < best){ best = obj; stall = 0; } else stall++;
                if (obj == 0){
                    printf("SOLUTION seed=%u\n", seed);
                    for (int a = 0; a < NORB; a++)
                        for (int b = 0; b < NORB; b++){
                            printf("D %d %d :", a, b);
                            for (int d = 0; d < P; d++) if (mem[a][b][d]) printf(" %d", d);
                            printf("\n");
                        }
                    fflush(stdout);
                    return 0;
                }
            }
            if (stall > 2000000) break;                    // restart
        }
        if (bestEver < 0 || best < bestEver){
            bestEver = best;
            printf("BEST obj=%lld seed=%u %.0fs\n", bestEver, seed, difftime(time(nullptr), t0));
            fflush(stdout);
        }
    }
    printf("RESULT best=%lld seed=%u\n", bestEver, seed);
    return 1;
}
