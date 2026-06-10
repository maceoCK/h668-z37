// cyclo_dfs: SOUND, COMPLETE coset-level DFS for H-cyclotomic lifts of a fixed orbit matrix.
//
// Conventions (= conf_core.py): vertices = norb orbits x Z_p, A[(i,a),(j,b)] = 1 iff (b-a) in D[i][j];
// D[j][i] = -D[i][j]; D[i][i] = -D[i][i], 0 not in D[i][i].  SRG(v,k,lam,mu):
//   A^2 + (mu-lam) A - (k-mu) I = mu J
// per (i,j) and g in Z_p (group-ring convolution form):
//   sum_l (D[i][l] * D[l][j])(g) = mu + (k-mu)[i==j && g==0] + (lam-mu)[g in D[i][j]].
// Cyclotomic ansatz: every D[i][j] is invariant under H = <mult> <= Z_p^*.  Convolutions of
// H-invariant sets are constant on H-classes of g, so there are (nc+1) equations per pair i<=j
// (nc = #cosets; class nc = the zero element).  Equations for (j,i) follow by g -> -g (A symmetric).
//
// Search order: all diagonal blocks first ((0,0)..(8,8)), then off-diagonal blocks row-major.
// Performance structure:
//   * for every off-diagonal block (i,j) and every pair of diagonal candidates (ci,cj), the
//     spectrally admissible candidate list (two-sided 2x2 interlacing + row budgets) is PRECOMPUTED,
//     sorted by a hash of the candidate's autocorrelation vector conv(D,-D);
//   * at a complete diagonal assignment the per-block lists are pointer lookups; a leaf dies
//     immediately if any list is empty or the adaptive autocorrelation interval check fails;
//   * within each row the LAST block is resolved by hash-join: the diagonal-pair equation (r,r)
//     pins its entire autocorrelation vector exactly, so matching candidates are found by lookup.
// Soundness of every prune:
//   * exact integer equation propagation (static interval bounds while open, equality when closed);
//   * adaptive autocorrelation interval check per diagonal-pair equation over current lists;
//   * spectral: every M_s is Hermitian with spectrum in {theta,tau}; every assigned principal
//     submatrix must have eigenvalues in [tau-eps, theta+eps] (two Cholesky tests; margin => only
//     clearly-violating branches are cut, no true solution lost);
//   * row norms: partial sum_l |M_s[i][l]|^2 <= (k-mu) - (mu-lam) m_ii + eps;
//   * candidate filters: diagonal m_ii(s) in [tau,theta]; off-diagonal |m_ij(s)|^2 <= global cap
//     and the (ci,cj)-conditional caps min((m_ii-tau)(m_jj-tau), (theta-m_ii)(theta-m_jj),
//     row budgets of both rows).
// Optional (gated to srg(333,166,82,83), p=37, norb=9, H <= QR(37)):
//   --pin a       : validated Galois trace pinning (QR-coset diag counts 3a-9, QNR 18-3a);
//                   a=4 is WLOG complete (a QNR decimation maps a=5 solutions to a=4).
//   --break-decim : diagonal-assignment lexicographic minimality under the (pinning-compatible)
//                   decimation action, prefix-pruned during the diagonal stage (ties kept; sound).
// Chunking: --chunk c --nchunks K  partitions surviving diagonal leaves round-robin (disjoint and
// exhaustive across c).  --count-diag counts surviving diagonal leaves only (chunk planning).
//
// Usage:
//   cyclo_dfs p norb k lam mu mult rfile row [--pin a] [--break-decim] [--maxsol N] [--maxsec S]
//             [--maxnodes N] [--chunk c --nchunks K] [--perm 0,1,..] [--count-diag] [--quiet]
//             [--solfile F]
// Output: SOLUTION blocks ("D i j : elements"), then
//   RESULT status=COMPLETE|STOPPED|SAT_CAP sat=<n> nodes=<n> sec=<t> ... dleaves=<n>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <complex>
#include <vector>
#include <algorithm>
using namespace std;
typedef complex<double> cd;
typedef long long ll;

static int P, NORB, KP, LAM, MU, MULT;
static int HE, NC, NCLS;
static vector<int> Hgrp;
static vector<vector<int>> CS;
static vector<int> negc, chiC;
static vector<vector<vector<int8_t>>> TC;
static double THETA, TAU, ENTCAP;
static const double EPS = 1e-6;

static int Rm[12][12];
static int apin = -1;
static bool breakDecim = false, quiet = false, countDiag = false;
static ll maxSol = -1, maxNodes = -1;
static double maxSec = -1;
static int chunkId = 0, nChunks = 1;
static ll diagTaskIdx = 0, diagLeaves = 0;
static bool permArg = false; static int permV[12];

static ll nodeCount = 0, solCount = 0;
static bool stopFlag = false, hitCap = false;
static time_t t0;

struct Cand {
    uint64_t mask, nmask;
    int z;
    vector<int> idx, nidx;
    vector<cd> chv;
    vector<int> aconv;                    // conv(D,-D)(d), d = 0..NC
    uint64_t akey;                        // hash of aconv
};

static int nBlocks;
static vector<int> blkI, blkJ;            // order: 9 diagonals, then off-diag row-major
static int blkOf[12][12];
static vector<vector<Cand>> cands;
static vector<int> assigned;

static int nPairs;
static int pairOf[12][12];
static vector<int> pI, pJ;
static vector<int> eqRhs, eqSum, eqUnk, eqLo, eqHi;
static vector<int> termLo, termHi;
static vector<int8_t> termDone;
struct TermRef { int pair, l; };
static vector<vector<TermRef>> blkTerms;

struct Undo { int eq, dS, dL, dH, dU; };
static vector<Undo> undoLog;
static vector<int> doneStack;

static vector<vector<cd>> Mcur;
static vector<vector<double>> rowNorm;
static int pinTarget[64], pinCount[64];
static vector<vector<int>> decimPerm;

// precomputed per (off-diag block, diag cand of row i, diag cand of row j):
struct Combo {
    vector<uint32_t> idx;                 // candidate indices, sorted by key
    vector<uint64_t> key;                 // parallel sorted keys
    int aLo[40], aHi[40];                 // aconv ranges over this list
};
static vector<vector<Combo>> combos;      // combos[b][ci * ndiag_j + cj]
static const Combo* curCombo[80];         // active combo per off-diag block (set at diag leaf)

static FILE* fsol = nullptr;

static inline int tid(int pair, int l){ return pair * (NORB + 1) + (l + 1); }

static uint64_t hashAconv(const vector<int>& a){
    uint64_t h = 1469598103934665603ull;
    for (int v : a){ h ^= (uint64_t)(uint8_t)v; h *= 1099511628211ull; }
    return h;
}

static void computeAconv(Cand& c){
    c.aconv.assign(NCLS, 0);
    for (int d = 0; d < NC; d++){
        int v = 0;
        for (int a : c.idx) for (int b : c.nidx) v += TC[a][b][d];
        if (c.z){ v += (int)((c.nmask >> d) & 1); v += (int)((c.mask >> d) & 1); }
        c.aconv[d] = v;
    }
    c.aconv[NC] = HE * (int)c.idx.size() + c.z;
    c.akey = hashAconv(c.aconv);
}

// ---------------------------------------------------------------- setup
static void buildScheme(){
    Hgrp.clear(); { int x = 1; do { Hgrp.push_back(x); x = (int)((ll)x * MULT % P); } while (x != 1); }
    HE = (int)Hgrp.size();
    vector<int> cidx(P, -1);
    CS.clear();
    for (int g = 1; g < P; g++){
        if (cidx[g] >= 0) continue;
        vector<int> c;
        for (int h : Hgrp) c.push_back((int)((ll)g * h % P));
        for (size_t a = 1; a < c.size(); a++){ int x = c[a]; int b = (int)a - 1;
            while (b >= 0 && c[b] > x){ c[b+1] = c[b]; b--; } c[b+1] = x; }
        for (int x : c) cidx[x] = (int)CS.size();
        CS.push_back(c);
    }
    NC = (int)CS.size(); NCLS = NC + 1;
    negc.assign(NC, 0);
    for (int c = 0; c < NC; c++) negc[c] = cidx[(P - CS[c][0]) % P];
    vector<char> isQR(P, 0);
    for (int x = 1; x < P; x++) isQR[(int)((ll)x * x % P)] = 1;
    chiC.assign(NC, 0);
    for (int c = 0; c < NC; c++) chiC[c] = isQR[CS[c][0]] ? 1 : -1;
    TC.assign(NC, vector<vector<int8_t>>(NC, vector<int8_t>(NC, 0)));
    for (int a = 0; a < NC; a++)
        for (int b = 0; b < NC; b++){
            vector<int> t(P, 0);
            for (int x : CS[a]) for (int y : CS[b]) t[(x + y) % P]++;
            for (int d = 0; d < NC; d++){
                int v = t[CS[d][0]];
                for (int x : CS[d]) if (t[x] != v){ fprintf(stderr, "class-constancy violated\n"); exit(2); }
                TC[a][b][d] = (int8_t)v;
            }
        }
    double disc = (double)(MU - LAM) * (MU - LAM) + 4.0 * (KP - MU);
    THETA = ((double)(LAM - MU) + sqrt(disc)) / 2.0;
    TAU   = ((double)(LAM - MU) - sqrt(disc)) / 2.0;
    ENTCAP = (double)(KP - MU) + (double)(MU - LAM) * (MU - LAM) / 4.0 + EPS;
}

static vector<vector<cd>> eta;
static void buildChars(){
    eta.assign(NC, vector<cd>(NC));
    for (int s = 0; s < NC; s++){
        int u = CS[s][0];
        for (int c = 0; c < NC; c++){
            cd v(0, 0);
            for (int x : CS[c]){
                double ang = 2.0 * M_PI * ((ll)u * x % P) / P;
                v += cd(cos(ang), sin(ang));
            }
            eta[s][c] = v;
        }
    }
}

static void genDiagCands(int b, int i){
    int target = Rm[i][i];
    if (target % HE != 0) return;
    int want = target / HE;
    vector<pair<int,int>> orbs;
    vector<char> seen(NC, 0);
    for (int c = 0; c < NC; c++){
        if (seen[c]) continue;
        seen[c] = 1; seen[negc[c]] = 1;
        orbs.push_back({c, negc[c]});
    }
    int m = (int)orbs.size();
    for (uint64_t sub = 0; sub < (1ull << m); sub++){
        int cnt = 0; uint64_t mask = 0;
        for (int t = 0; t < m; t++) if ((sub >> t) & 1){
            mask |= 1ull << orbs[t].first;
            mask |= 1ull << orbs[t].second;
            cnt += (orbs[t].first == orbs[t].second) ? 1 : 2;
        }
        if (cnt != want) continue;
        Cand c; c.mask = c.nmask = mask; c.z = 0;
        for (int t = 0; t < NC; t++) if ((mask >> t) & 1) c.idx.push_back(t);
        c.nidx = c.idx;
        bool ok = true;
        c.chv.resize(NC);
        for (int s = 0; s < NC && ok; s++){
            cd v(0, 0);
            for (int t : c.idx) v += eta[s][t];
            if (fabs(v.imag()) > 1e-7) ok = false;
            if (v.real() < TAU - EPS || v.real() > THETA + EPS) ok = false;
            c.chv[s] = v;
        }
        if (ok){ computeAconv(c); cands[b].push_back(c); }
    }
}

static void genOffCands(int b, int i, int j){
    int r = Rm[i][j];
    int zlo = 0, zhi = 0;
    if (HE >= 2){
        int z = r % HE;
        if (z > 1) return;
        zlo = zhi = z;
    } else { zlo = 0; zhi = 1; }
    for (int z = zlo; z <= zhi; z++){
        if ((r - z) % HE || (r - z) < 0) continue;
        int want = (r - z) / HE;
        if (want > NC) continue;
        if (want == 0){
            Cand c; c.mask = c.nmask = 0; c.z = z;
            c.chv.assign(NC, cd(z, 0));
            computeAconv(c);
            cands[b].push_back(c);
            continue;
        }
        vector<int> pick(want);
        for (int t = 0; t < want; t++) pick[t] = t;
        bool more = true;
        while (more){
            uint64_t mask = 0;
            for (int t = 0; t < want; t++) mask |= 1ull << pick[t];
            Cand c; c.mask = mask; c.z = z; c.nmask = 0;
            for (int t = 0; t < NC; t++) if ((mask >> t) & 1){ c.idx.push_back(t); c.nmask |= 1ull << negc[t]; }
            for (int t = 0; t < NC; t++) if ((c.nmask >> t) & 1) c.nidx.push_back(t);
            bool ok = true;
            c.chv.resize(NC);
            for (int s = 0; s < NC && ok; s++){
                cd v(z, 0);
                for (int t : c.idx) v += eta[s][t];
                if (norm(v) > ENTCAP) ok = false;
                c.chv[s] = v;
            }
            if (ok){ computeAconv(c); cands[b].push_back(c); }
            int t = want - 1;
            while (t >= 0 && pick[t] == NC - want + t) t--;
            if (t < 0) more = false;
            else { pick[t]++; for (int u = t + 1; u < want; u++) pick[u] = pick[u-1] + 1; }
        }
    }
}

static inline int convLoOf(int rA, int rB){ int lo = rA + rB - P; return lo > 0 ? lo : 0; }
static inline int convHiOf(int rA, int rB){ return rA < rB ? rA : rB; }

static void buildBlocks(){
    nBlocks = 0;
    blkI.clear(); blkJ.clear();
    for (int i = 0; i < NORB; i++){ blkOf[i][i] = nBlocks; blkI.push_back(i); blkJ.push_back(i); nBlocks++; }
    for (int i = 0; i < NORB; i++)
        for (int j = i + 1; j < NORB; j++){ blkOf[i][j] = nBlocks; blkI.push_back(i); blkJ.push_back(j); nBlocks++; }
    assigned.assign(nBlocks, -1);
    cands.assign(nBlocks, {});
    for (int b = 0; b < nBlocks; b++){
        if (blkI[b] == blkJ[b]) genDiagCands(b, blkI[b]);
        else genOffCands(b, blkI[b], blkJ[b]);
    }
    nPairs = 0;
    pI.clear(); pJ.clear();
    for (int i = 0; i < NORB; i++)
        for (int j = i; j < NORB; j++){ pairOf[i][j] = nPairs; pI.push_back(i); pJ.push_back(j); nPairs++; }
    eqRhs.assign(nPairs * NCLS, 0);
    eqSum.assign(nPairs * NCLS, 0);
    eqUnk.assign(nPairs * NCLS, 0);
    eqLo.assign(nPairs * NCLS, 0);
    eqHi.assign(nPairs * NCLS, 0);
    termLo.assign(nPairs * (NORB + 1), 0);
    termHi.assign(nPairs * (NORB + 1), 0);
    termDone.assign(nPairs * (NORB + 1), 0);
    blkTerms.assign(nBlocks, {});
    int w = MU - LAM;
    if (w < 0){ fprintf(stderr, "need mu>=lam\n"); exit(2); }
    for (int p = 0; p < nPairs; p++){
        int i = pI[p], j = pJ[p];
        for (int d = 0; d < NCLS; d++)
            eqRhs[p * NCLS + d] = MU + ((i == j && d == NC) ? (KP - MU) : 0);
        int sLo = 0, sHi = 0;
        for (int l = 0; l < NORB; l++){
            int rA = Rm[i][l], rB = Rm[l][j];
            termLo[tid(p, l)] = convLoOf(rA, rB);
            termHi[tid(p, l)] = convHiOf(rA, rB);
            sLo += termLo[tid(p, l)]; sHi += termHi[tid(p, l)];
            int bA = (i <= l) ? blkOf[i][l] : blkOf[l][i];
            int bB = (l <= j) ? blkOf[l][j] : blkOf[j][l];
            blkTerms[bA].push_back({p, l});
            if (bB != bA) blkTerms[bB].push_back({p, l});
        }
        termLo[tid(p, -1)] = 0; termHi[tid(p, -1)] = w;
        sHi += w;
        blkTerms[blkOf[i][j]].push_back({p, -1});
        for (int d = 0; d < NCLS; d++){
            eqUnk[p * NCLS + d] = NORB + 1;
            eqLo[p * NCLS + d] = sLo;
            eqHi[p * NCLS + d] = sHi;
        }
    }
    Mcur.assign(NC, vector<cd>(NORB * NORB, cd(0, 0)));
    rowNorm.assign(NC, vector<double>(NORB, 0.0));
}

// precompute spectrally-filtered candidate lists per (block, diag-cand-i, diag-cand-j)
static void buildCombos(){
    combos.assign(nBlocks, {});
    for (int b = NORB; b < nBlocks; b++){
        int i = blkI[b], j = blkJ[b];
        int ni = (int)cands[blkOf[i][i]].size(), nj = (int)cands[blkOf[j][j]].size();
        combos[b].resize((size_t)ni * nj);
        for (int ci = 0; ci < ni; ci++){
            const Cand& di = cands[blkOf[i][i]][ci];
            for (int cj = 0; cj < nj; cj++){
                const Cand& dj = cands[blkOf[j][j]][cj];
                Combo& C = combos[b][(size_t)ci * nj + cj];
                double caps[64];
                bool any = true;
                for (int s = 0; s < NC; s++){
                    double mii = di.chv[s].real(), mjj = dj.chv[s].real();
                    double bi = (KP - MU) - (MU - LAM) * mii - mii * mii + EPS;
                    double bj = (KP - MU) - (MU - LAM) * mjj - mjj * mjj + EPS;
                    double c2a = (mii - TAU) * (mjj - TAU) + EPS;
                    double c2b = (THETA - mii) * (THETA - mjj) + EPS;
                    double cap = bi;
                    if (bj < cap) cap = bj;
                    if (c2a < cap) cap = c2a;
                    if (c2b < cap) cap = c2b;
                    caps[s] = cap;
                    if (cap < -EPS) any = false;
                }
                for (int d = 0; d < NCLS; d++){ C.aLo[d] = 1 << 28; C.aHi[d] = -(1 << 28); }
                if (any){
                    vector<pair<uint64_t, uint32_t>> tmp;
                    for (uint32_t ck = 0; ck < cands[b].size(); ck++){
                        const Cand& c = cands[b][ck];
                        bool ok = true;
                        for (int s = 0; s < NC && ok; s++)
                            if (norm(c.chv[s]) > caps[s]) ok = false;
                        if (!ok) continue;
                        tmp.push_back({c.akey, ck});
                        for (int d = 0; d < NCLS; d++){
                            int v = c.aconv[d];
                            if (v < C.aLo[d]) C.aLo[d] = v;
                            if (v > C.aHi[d]) C.aHi[d] = v;
                        }
                    }
                    // sort by key
                    sort(tmp.begin(), tmp.end());
                    C.idx.reserve(tmp.size()); C.key.reserve(tmp.size());
                    for (auto& pr : tmp){ C.key.push_back(pr.first); C.idx.push_back(pr.second); }
                }
                if (C.idx.empty()) for (int d = 0; d < NCLS; d++){ C.aLo[d] = 0; C.aHi[d] = 0; }
            }
        }
    }
}

static void buildDecim(){
    decimPerm.clear();
    if (!breakDecim) return;
    vector<int> cidx(P, -1);
    for (int c = 0; c < NC; c++) for (int x : CS[c]) cidx[x] = c;
    for (int s = 1; s < NC; s++){
        if (apin >= 0 && chiC[s] != 1) continue;
        int u = CS[s][0];
        vector<int> perm(NC);
        for (int c = 0; c < NC; c++) perm[c] = cidx[(int)((ll)u * CS[c][0] % P)];
        decimPerm.push_back(perm);
    }
}

// ---------------------------------------------------------------- runtime
static inline bool getD(int a, int b, const vector<int>*& idx, int& z, uint64_t& mask){
    int B = (a <= b) ? blkOf[a][b] : blkOf[b][a];
    int ci = assigned[B];
    if (ci < 0) return false;
    const Cand& c = cands[B][ci];
    if (a <= b){ idx = &c.idx;  mask = c.mask;  }
    else       { idx = &c.nidx; mask = c.nmask; }
    z = c.z;
    return true;
}

static inline bool eqOK(int eq){
    int need = eqRhs[eq];
    if (eqUnk[eq] == 0) return eqSum[eq] == need;
    return eqSum[eq] + eqLo[eq] <= need && need <= eqSum[eq] + eqHi[eq];
}

static inline void eqApply(int eq, int dS, int dL, int dH, int dU){
    eqSum[eq] += dS; eqLo[eq] += dL; eqHi[eq] += dH; eqUnk[eq] += dU;
    undoLog.push_back({eq, dS, dL, dH, dU});
}

static bool completeTerm(int p, int l){
    int i = pI[p], j = pJ[p];
    int T = tid(p, l);
    if (termDone[T]) return true;
    const vector<int> *ia, *ib; int za, zb; uint64_t ma, mb;
    if (l == -1){
        if (!getD(i, j, ia, za, ma)) return true;
        termDone[T] = 1; doneStack.push_back(T);
        int w = MU - LAM;
        bool ok = true;
        for (int d = 0; d < NCLS; d++){
            int v = (d == NC) ? za : (int)((ma >> d) & 1);
            int eq = p * NCLS + d;
            eqApply(eq, w * v, -termLo[T], -termHi[T], -1);
            if (!eqOK(eq)) ok = false;
        }
        return ok;
    }
    if (!getD(i, l, ia, za, ma) || !getD(l, j, ib, zb, mb)) return true;
    termDone[T] = 1; doneStack.push_back(T);
    bool ok = true;
    for (int d = 0; d < NCLS; d++){
        int v;
        if (d == NC){
            v = za * zb;
            for (int a : *ia) v += ((mb >> negc[a]) & 1) ? HE : 0;
        } else {
            v = 0;
            for (int a : *ia) for (int bb : *ib) v += TC[a][bb][d];
            if (za) v += (int)((mb >> d) & 1);
            if (zb) v += (int)((ma >> d) & 1);
        }
        int eq = p * NCLS + d;
        eqApply(eq, v, -termLo[T], -termHi[T], -1);
        if (!eqOK(eq)) ok = false;
    }
    return ok;
}

// adaptive autocorrelation interval check (off-diag stage; uses active combo lists)
static bool diagEqAdaptiveOK(){
    for (int i = 0; i < NORB; i++){
        int p = pairOf[i][i];
        for (int d = 0; d < NCLS; d++){
            int eq = p * NCLS + d;
            if (eqUnk[eq] == 0) continue;
            int need = eqRhs[eq] - eqSum[eq];
            int lo = 0, hi = 0;
            if (!termDone[tid(p, -1)]) hi += MU - LAM;
            for (int l = 0; l < NORB; l++){
                if (termDone[tid(p, l)]) continue;
                int b = (i <= l) ? blkOf[i][l] : blkOf[l][i];
                int dd = (l > i) ? d : (d == NC ? NC : negc[d]);
                const Combo* C = curCombo[b];
                lo += C->aLo[dd]; hi += C->aHi[dd];
            }
            if (need < lo || need > hi) return false;
        }
    }
    return true;
}

static bool cholOK(const cd* A, int n, double sgn, double bound){
    cd B[12][12];
    double margin = EPS * (fabs(bound) + 10.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++){
            cd v = A[i * 12 + j];
            B[i][j] = (i == j) ? cd(sgn * (v.real() - bound) + margin, 0) : sgn * v;
        }
    for (int c = 0; c < n; c++){
        double d = B[c][c].real();
        for (int k = 0; k < c; k++) d -= norm(B[c][k]);
        if (d < 0) return false;
        double rt = sqrt(d);
        if (rt < 1e-14){ for (int r = c + 1; r < n; r++) if (norm(B[r][c]) > 1e-20) return false; B[c][c] = 0; continue; }
        B[c][c] = rt;
        for (int r = c + 1; r < n; r++){
            cd v = B[r][c];
            for (int k = 0; k < c; k++) v -= B[r][k] * conj(B[c][k]);
            B[r][c] = v / rt;
        }
    }
    return true;
}

// principal submatrix on S = {0..i} u {j}
static bool spectralOK(int i, int j){
    static cd A[12 * 12];
    int S[12], ns = 0;
    for (int t = 0; t <= i; t++) S[ns++] = t;
    S[ns++] = j;
    for (int s = 0; s < NC; s++){
        const vector<cd>& M = Mcur[s];
        for (int a = 0; a < ns; a++)
            for (int b = 0; b < ns; b++)
                A[a * 12 + b] = M[S[a] * NORB + S[b]];
        if (!cholOK(A, ns, +1.0, TAU))   return false;
        if (!cholOK(A, ns, -1.0, THETA)) return false;
    }
    return true;
}

static bool diagLexMinPrefixOK(int r){
    uint64_t cur[12];
    for (int i = 0; i <= r; i++) cur[i] = cands[blkOf[i][i]][assigned[blkOf[i][i]]].mask;
    for (auto& perm : decimPerm){
        bool decided = false;
        for (int i = 0; i <= r && !decided; i++){
            uint64_t im = 0, m = cur[i];
            for (int c = 0; c < NC; c++) if ((m >> c) & 1) im |= 1ull << perm[c];
            if (im < m) return false;
            if (im > m) decided = true;
        }
    }
    return true;
}

static void printSolution(){
    FILE* f = fsol ? fsol : stdout;
    fprintf(f, "SOLUTION\n");
    for (int i = 0; i < NORB; i++)
        for (int j = 0; j < NORB; j++){
            const vector<int>* idx; int z; uint64_t m;
            getD(i, j, idx, z, m);
            fprintf(f, "D %d %d :", i, j);
            if (z) fprintf(f, " 0");
            for (int c : *idx) for (int x : CS[c]) fprintf(f, " %d", x);
            fprintf(f, "\n");
        }
    fflush(f);
}

static void dfs(int b);

static void tryAssign(int b, int ci){
    if (stopFlag) return;
    nodeCount++;
    if ((nodeCount & ((1 << 22) - 1)) == 0){
        double el = difftime(time(nullptr), t0);
        if (!quiet && (nodeCount & ((1 << 26) - 1)) == 0){
            fprintf(stderr, "  ... nodes=%lld sols=%lld depth=%d dleaves=%lld %.0fs\n",
                    nodeCount, solCount, b, diagLeaves, el);
        }
        if (maxSec > 0 && el > maxSec){ stopFlag = true; hitCap = true; return; }
    }
    if (maxNodes > 0 && nodeCount > maxNodes){ stopFlag = true; hitCap = true; return; }
    int i = blkI[b], j = blkJ[b];
    const Cand& c = cands[b][ci];
    assigned[b] = ci;
    size_t emark = undoLog.size(), dmark = doneStack.size();
    double oldRN[64][2];
    bool ok = true;

    if (i == j && apin >= 0){
        for (int cc : c.idx){ pinCount[cc]++; if (pinCount[cc] > pinTarget[cc]) ok = false; }
        if (ok){
            int remain = NORB - 1 - i;
            for (int cc = 0; cc < NC; cc++) if (pinCount[cc] + remain < pinTarget[cc]){ ok = false; break; }
        }
    }
    for (int s = 0; s < NC; s++){
        cd v = c.chv[s];
        Mcur[s][i * NORB + j] = v;
        Mcur[s][j * NORB + i] = conj(v);
        oldRN[s][0] = rowNorm[s][i];
        rowNorm[s][i] += norm(v);
        if (i != j){ oldRN[s][1] = rowNorm[s][j]; rowNorm[s][j] += norm(v); }
    }
    if (ok && i != j){
        for (int s = 0; s < NC && ok; s++){
            for (int rr : {i, j}){
                double mrr = Mcur[s][rr * NORB + rr].real();
                double cap = (KP - MU) - (MU - LAM) * mrr + 1e-5;
                if (rowNorm[s][rr] > cap){ ok = false; break; }
            }
        }
    }
    if (ok)
        for (const TermRef& t : blkTerms[b])
            if (!completeTerm(t.pair, t.l)){ ok = false; break; }
    if (ok && i == j && breakDecim && !diagLexMinPrefixOK(i)) ok = false;
    if (ok && i != j){
        if (!diagEqAdaptiveOK()) ok = false;
        if (ok && !spectralOK(i, j)) ok = false;
    }
    if (ok) dfs(b + 1);

    while (undoLog.size() > emark){
        Undo& u = undoLog.back();
        eqSum[u.eq] -= u.dS; eqLo[u.eq] -= u.dL; eqHi[u.eq] -= u.dH; eqUnk[u.eq] -= u.dU;
        undoLog.pop_back();
    }
    while (doneStack.size() > dmark){ termDone[doneStack.back()] = 0; doneStack.pop_back(); }
    if (i == j && apin >= 0) for (int cc : c.idx) pinCount[cc]--;
    for (int s = 0; s < NC; s++){
        Mcur[s][i * NORB + j] = cd(0, 0);
        Mcur[s][j * NORB + i] = cd(0, 0);
        rowNorm[s][i] = oldRN[s][0];
        if (i != j) rowNorm[s][j] = oldRN[s][1];
    }
    assigned[b] = -1;
}

static void dfs(int b){
    if (stopFlag) return;
    if (b == nBlocks){
        solCount++;
        printSolution();
        if (maxSol > 0 && solCount >= maxSol) stopFlag = true;
        return;
    }
    if (b == NORB){                               // diagonal complete: activate combos, leaf checks
        for (int bb = NORB; bb < nBlocks; bb++){
            int ii = blkI[bb], jj = blkJ[bb];
            int nj = (int)cands[blkOf[jj][jj]].size();
            const Combo& C = combos[bb][(size_t)assigned[blkOf[ii][ii]] * nj + assigned[blkOf[jj][jj]]];
            if (C.idx.empty()) return;
            curCombo[bb] = &C;
        }
        if (!diagEqAdaptiveOK()) return;
        diagLeaves++;
        if (countDiag) return;
        if (nChunks > 1 && (diagTaskIdx++ % nChunks) != chunkId) return;
        // fall through to first off-diag block
    }
    int i = blkI[b], j = blkJ[b];
    if (i == j){
        for (int ci = 0; ci < (int)cands[b].size(); ci++){ tryAssign(b, ci); if (stopFlag) return; }
        return;
    }
    const Combo& C = *curCombo[b];
    if (j == NORB - 1 && i < NORB - 1){
        // last block of row i: eq (i,i) pins this block's aconv vector exactly IF it is the only
        // unknown conv term left; verify and hash-join, else fall back to full scan
        int p = pairOf[i][i];
        bool only = !termDone[tid(p, j)];
        for (int l = 0; l < NORB && only; l++)
            if (l != j && !termDone[tid(p, l)]) only = false;
        if (only && termDone[tid(p, -1)]){
            vector<int> need(NCLS);
            for (int d = 0; d < NCLS; d++){
                int dd = (j > i) ? d : (d == NC ? NC : negc[d]);     // j > i always here
                (void)dd;
                need[d] = eqRhs[p * NCLS + d] - eqSum[p * NCLS + d];
            }
            // need[] is conv(D[i][j], D[j][i])(d) = aconv of stored candidate (i<j)
            uint64_t k = hashAconv(need);
            size_t lo = 0, hi = C.key.size();
            while (lo < hi){ size_t mid = (lo + hi) / 2; if (C.key[mid] < k) lo = mid + 1; else hi = mid; }
            for (size_t t = lo; t < C.key.size() && C.key[t] == k; t++){
                const Cand& cc = cands[b][C.idx[t]];
                bool exact = true;
                for (int d = 0; d < NCLS && exact; d++) if (cc.aconv[d] != need[d]) exact = false;
                if (exact){ tryAssign(b, C.idx[t]); if (stopFlag) return; }
            }
            return;
        }
    }
    for (uint32_t ci : C.idx){ tryAssign(b, ci); if (stopFlag) return; }
}

int main(int argc, char** argv){
    if (argc < 9){ fprintf(stderr, "usage: %s p norb k lam mu mult rfile row [options]\n", argv[0]); return 2; }
    P = atoi(argv[1]); NORB = atoi(argv[2]); KP = atoi(argv[3]); LAM = atoi(argv[4]); MU = atoi(argv[5]); MULT = atoi(argv[6]);
    const char* rfile = argv[7]; int row = atoi(argv[8]);
    for (int a = 9; a < argc; a++){
        if (!strcmp(argv[a], "--pin")) apin = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--break-decim")) breakDecim = true;
        else if (!strcmp(argv[a], "--maxsol")) maxSol = atoll(argv[++a]);
        else if (!strcmp(argv[a], "--maxsec")) maxSec = atof(argv[++a]);
        else if (!strcmp(argv[a], "--maxnodes")) maxNodes = atoll(argv[++a]);
        else if (!strcmp(argv[a], "--chunk")) chunkId = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--nchunks")) nChunks = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--quiet")) quiet = true;
        else if (!strcmp(argv[a], "--solfile")) fsol = fopen(argv[++a], "w");
        else if (!strcmp(argv[a], "--count-diag")) countDiag = true;
        else if (!strcmp(argv[a], "--perm")){
            char* s = argv[++a];
            for (int t = 0; t < 12; t++) permV[t] = t;
            int t = 0; char* tok = strtok(s, ",");
            while (tok && t < 12){ permV[t++] = atoi(tok); tok = strtok(nullptr, ","); }
            permArg = true;
        }
        else { fprintf(stderr, "bad arg %s\n", argv[a]); return 2; }
    }
    if (NORB > 12 || P > 64){ fprintf(stderr, "limits: norb<=12, p<=64\n"); return 2; }
    {
        FILE* f = fopen(rfile, "r");
        if (!f){ fprintf(stderr, "cannot open %s\n", rfile); return 2; }
        char line[8192]; int ln = 0; bool got = false;
        while (fgets(line, sizeof line, f)){
            if (ln++ == row){
                char* tok = strtok(line, " \t\n");
                for (int t = 0; t < NORB * NORB; t++){
                    if (!tok){ fprintf(stderr, "short R line\n"); return 2; }
                    Rm[t / NORB][t % NORB] = atoi(tok);
                    tok = strtok(nullptr, " \t\n");
                }
                got = true; break;
            }
        }
        fclose(f);
        if (!got){ fprintf(stderr, "row %d not found\n", row); return 2; }
    }
    if (permArg){
        int T[12][12];
        for (int a2 = 0; a2 < NORB; a2++) for (int b2 = 0; b2 < NORB; b2++) T[a2][b2] = Rm[permV[a2]][permV[b2]];
        for (int a2 = 0; a2 < NORB; a2++) for (int b2 = 0; b2 < NORB; b2++) Rm[a2][b2] = T[a2][b2];
    }
    t0 = time(nullptr);
    buildScheme();
    buildChars();
    if (apin >= 0){
        bool gate = (P == 37 && NORB == 9 && KP == 166 && LAM == 82 && MU == 83);
        vector<char> isQR(P, 0);
        for (int x = 1; x < P; x++) isQR[(int)((ll)x * x % P)] = 1;
        for (int h : Hgrp) if (!isQR[h]) gate = false;
        if (!gate){ fprintf(stderr, "--pin only valid for srg(333) with H<=QR(37)\n"); return 2; }
        for (int c = 0; c < NC; c++) pinTarget[c] = (chiC[c] == 1) ? (3 * apin - 9) : (18 - 3 * apin);
        memset(pinCount, 0, sizeof pinCount);
    }
    buildBlocks();
    buildCombos();
    buildDecim();
    bool feasible = true;
    for (int b = 0; b < nBlocks && feasible; b++) if (cands[b].empty()) feasible = false;
    for (int e = 0; e < nPairs * NCLS && feasible; e++) if (!eqOK(e)) feasible = false;
    if (feasible) dfs(0);
    double el = difftime(time(nullptr), t0);
    const char* status = hitCap ? "STOPPED" : "COMPLETE";
    if (maxSol > 0 && solCount >= maxSol) status = "SAT_CAP";
    printf("RESULT status=%s sat=%lld nodes=%lld sec=%.1f p=%d norb=%d mult=%d row=%d pin=%d chunks=%d/%d dleaves=%lld\n",
           status, solCount, nodeCount, el, P, NORB, MULT, row, apin, chunkId, nChunks, diagLeaves);
    if (fsol) fclose(fsol);
    return 0;
}
