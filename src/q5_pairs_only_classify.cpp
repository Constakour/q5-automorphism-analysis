// q5_pairs_only_classify.cpp
// Enumerate ALL pairs-only partitions {S, S^c} in Q5 (6435 of them),
// and classify whether S can be mapped to S^c by an automorphism of Q5.
//
// Output JSON contains three lists:
//   - symmetric_translation:  swap via XOR only (perm = identity)
//   - symmetric_requires_permutation: swap via XOR + coordinate permutation, but NOT via XOR alone
//   - non_symmetric: no cube automorphism swaps S and S^c
//
// Build: g++ -O3 -std=c++20 -fopenmp q5_pairs_only_classify.cpp -o q5_pairs_only_classify
// Run:   ./q5_pairs_only_classify --threads 12

#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

static constexpr int M = 5;
static constexpr int N = 32;                 // 2^5
static constexpr uint32_t FULL = 0xFFFFFFFFu;

static inline int popcount_u32(uint32_t x) { return __builtin_popcount(x); }
static inline int ctz_u32(uint32_t x)      { return __builtin_ctz(x); }

// Gosper: next integer with same popcount (for <= 64-bit)
static inline uint64_t nextCombPopcount(uint64_t x){
    uint64_t u = x & -x;
    uint64_t v = u + x;
    if (v == 0) return 0;
    return v + (((v ^ x) / u) >> 2);
}

// Exact binom for small n (here 15 choose 8)
static uint64_t binom_u64(int n, int k){
    if (k < 0 || k > n) return 0;
    k = min(k, n-k);
    uint64_t num = 1, den = 1;
    for (int i = 1; i <= k; i++){
        num *= (uint64_t)(n - k + i);
        den *= (uint64_t)i;
        uint64_t g = std::gcd(num, den);
        num /= g; den /= g;
    }
    return num / den;
}

// Build all 120 coordinate permutations of [0,1,2,3,4].
// Convention: p[b] is the old bit-position that moves into new bit b.
// So new_bit[b] = old_bit[p[b]].
static vector< array<int,5> > allPerms() {
    vector< array<int,5> > perms;
    vector<int> base = {0,1,2,3,4};
    do {
        array<int,5> p{};
        for (int i=0;i<5;i++) p[i]=base[i];
        perms.push_back(p);
    } while (next_permutation(base.begin(), base.end()));
    return perms;
}

static inline bool is_identity_perm(const array<int,5>& p){
    for(int i=0;i<5;i++) if(p[i]!=i) return false;
    return true;
}

// Apply affine automorphism to a vertex label v in {0..31}:
//   1) XOR by flip
//   2) permute bits by p: new_bit[b] = old_bit[p[b]]
static inline uint32_t transform_vertex(uint32_t v, uint32_t flip, const array<int,5>& p){
    uint32_t x = (v ^ flip) & 31u;
    uint32_t y = 0;
    for(int b=0;b<5;b++){
        uint32_t oldbit = (x >> p[b]) & 1u;
        y |= (oldbit << b);
    }
    return y & 31u;
}

// Transform an entire set mask S (32-bit) under (flip, p)
static inline uint32_t transform_mask(uint32_t S, uint32_t flip, const array<int,5>& p){
    uint32_t T = 0;
    uint32_t mm = S;
    while(mm){
        uint32_t lsb = mm & -mm;
        int v = ctz_u32(lsb);
        uint32_t u = transform_vertex((uint32_t)v, flip, p);
        T |= (1u << u);
        mm ^= lsb;
    }
    return T;
}

// Check if there exists a pure translation (XOR) sending S -> comp(S).
// If so, return true and output the flip.
static bool find_translation_swap(uint32_t S, uint32_t &out_flip){
    uint32_t comp = FULL ^ S;
    array<int,5> id = {0,1,2,3,4};
    for(uint32_t flip=0; flip<32; flip++){
        if(transform_mask(S, flip, id) == comp){
            out_flip = flip;
            return true;
        }
    }
    return false;
}

// Check if there exists ANY cube automorphism (flip + perm) sending S -> comp(S).
// Return witness (flip, perm) if exists.
static bool find_any_swap(uint32_t S, uint32_t &out_flip, array<int,5> &out_perm,
                          const vector< array<int,5> > &perms){
    uint32_t comp = FULL ^ S;
    for(const auto& p : perms){
        for(uint32_t flip=0; flip<32; flip++){
            if(transform_mask(S, flip, p) == comp){
                out_flip = flip;
                out_perm = p;
                return true;
            }
        }
    }
    return false;
}

struct Pair { int a,b; };

enum class Kind : uint8_t { TRANS=0, REQUIRES_PERM=1, NONE=2 };

struct Result {
    uint32_t S = 0;
    Kind kind = Kind::NONE;
    uint32_t flip = 0;
    array<int,5> perm = {0,1,2,3,4};
};

// Turn mask into sorted vertex list
static vector<int> mask_to_vertices(uint32_t S){
    vector<int> v;
    v.reserve(16);
    for(int i=0;i<N;i++){
        if(S & (1u<<i)) v.push_back(i);
    }
    return v;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int threads = 12; // default to 12 as requested
    string outname = "q5_pairs_only_classification.json";

    for(int i=1;i<argc;i++){
        string s = argv[i];
        if(s=="--threads" && i+1<argc) threads = stoi(argv[++i]);
        else if(s=="--out" && i+1<argc) outname = argv[++i];
        else {
            cerr << "Usage: ./q5_pairs_only_classify [--threads T] [--out file.json]\n";
            return 1;
        }
    }

#ifdef _OPENMP
    if(threads > 0) omp_set_num_threads(threads);
    cerr << "OpenMP threads: " << omp_get_max_threads() << "\n";
#else
    cerr << "Warning: compiled without OpenMP; running single-threaded.\n";
#endif

    // Build antipodal pairs (v, v^31).
    vector<bool> used(N,false);
    vector<Pair> pairsList;
    pairsList.reserve(16);
    for(int i=0;i<N;i++){
        if(!used[i]){
            int j = i ^ 31;
            used[i]=true; used[j]=true;
            pairsList.push_back({min(i,j), max(i,j)});
        }
    }
    if((int)pairsList.size()!=16){
        cerr << "Error: expected 16 antipodal pairs.\n";
        return 1;
    }

    // Find the pair containing vertex 31; fix it to avoid double counting modulo complement.
    int fixPairIdx = -1;
    for(int i=0;i<16;i++){
        if(pairsList[i].b == 31){
            fixPairIdx = i;
            break;
        }
    }
    if(fixPairIdx < 0){
        cerr << "Error: no pair containing 31 found.\n";
        return 1;
    }

    // Mask per antipodal pair
    array<uint32_t,16> pairMask{};
    for(int i=0;i<16;i++){
        pairMask[i] = (1u << pairsList[i].a) | (1u << pairsList[i].b);
    }

    // Indices of the remaining 15 pairs (excluding fixed pair)
    array<int,15> idxList{};
    {
        int t=0;
        for(int i=0;i<16;i++){
            if(i==fixPairIdx) continue;
            idxList[t++] = i;
        }
    }

    // Enumerate all 15-bit combos of popcount 8: total = C(15,8)=6435
    uint64_t total = binom_u64(15,8);
    vector<uint16_t> combs;
    combs.reserve((size_t)total);

    uint64_t comb = (1ULL<<8) - 1ULL;        // lowest popcount-8
    uint64_t endC = comb << (15-8);          // highest popcount-8

    while(true){
        combs.push_back((uint16_t)comb);
        if(comb == endC) break;
        uint64_t nxt = nextCombPopcount(comb);
        if(nxt==0) break;
        comb = nxt;
    }
    if(combs.size() != total){
        cerr << "Warning: expected " << total << " combos, got " << combs.size() << "\n";
    }

    // Precompute all perms
    auto perms = allPerms();
    array<int,5> id = {0,1,2,3,4};

    // Results array (thread-safe: each index written once)
    vector<Result> res(combs.size());

    // Helper: build S from a 15-bit combination mask
    auto build_S = [&](uint16_t c)->uint32_t{
        uint32_t S = 0;
        uint16_t x = c;
        while(x){
            uint16_t lsb = (uint16_t)(x & (uint16_t)(- (int16_t)x));
            int b = __builtin_ctz((unsigned)lsb); // 0..14
            x ^= lsb;
            int pairIdx = idxList[b];
            S |= pairMask[pairIdx];
        }
        // sanity: should be 16 vertices
        // (each chosen pair contributes 2 vertices)
        return S;
    };

    // Parallel classification
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(long long i=0;i<(long long)combs.size(); i++){
        uint32_t S = build_S(combs[(size_t)i]);

        Result r;
        r.S = S;

        // 1) translation only?
        uint32_t flipT = 0;
        if(find_translation_swap(S, flipT)){
            r.kind = Kind::TRANS;
            r.flip = flipT;
            r.perm = id;
            res[(size_t)i] = r;
            continue;
        }

        // 2) any automorphism?
        uint32_t flipA = 0;
        array<int,5> permA = id;
        if(find_any_swap(S, flipA, permA, perms)){
            r.kind = Kind::REQUIRES_PERM;
            r.flip = flipA;
            r.perm = permA;
            res[(size_t)i] = r;
            continue;
        }

        // 3) none
        r.kind = Kind::NONE;
        res[(size_t)i] = r;
    }

    // Count classes
    uint64_t ct_trans=0, ct_req=0, ct_none=0;
    for(const auto& r: res){
        if(r.kind==Kind::TRANS) ct_trans++;
        else if(r.kind==Kind::REQUIRES_PERM) ct_req++;
        else ct_none++;
    }

    // Write JSON
    ofstream fout(outname);
    if(!fout){
        cerr << "Cannot open output file: " << outname << "\n";
        return 1;
    }

    auto dump_entry = [&](const Result& r){
        char buf[16];
        sprintf(buf, "0x%08X", r.S);
        auto verts = mask_to_vertices(r.S);

        fout << "    {\n";
        fout << "      \"mask_hex\": \"" << buf << "\",\n";
        fout << "      \"vertices\": [";
        for(size_t i=0;i<verts.size();i++){
            if(i) fout << ", ";
            fout << verts[i];
        }
        fout << "]";
        if(r.kind != Kind::NONE){
            fout << ",\n";
            fout << "      \"flip_mask\": " << r.flip << ",\n";
            fout << "      \"perm\": [";
            for(int b=0;b<5;b++){
                if(b) fout << ", ";
                fout << r.perm[b];
            }
            fout << "]\n";
        } else {
            fout << "\n";
        }
        fout << "    }";
    };

    fout << "{\n";
    fout << "  \"meta\": {\n";
    fout << "    \"m\": 5,\n";
    fout << "    \"family\": \"pairs-only halves (8 antipodal pairs)\",\n";
    fout << "    \"count_partitions_mod_complement\": " << res.size() << ",\n";
    fout << "    \"notes\": \"Each entry represents a partition {S,S^c} (mod complement) by fixing the antipodal pair containing vertex 31 to lie in S^c.\"\n";
    fout << "  },\n\n";

    // symmetric_translation
    fout << "  \"symmetric_translation\": [\n";
    bool first=true;
    for(const auto& r: res){
        if(r.kind!=Kind::TRANS) continue;
        if(!first) fout << ",\n";
        first=false;
        dump_entry(r);
    }
    fout << "\n  ],\n\n";

    // symmetric_requires_permutation
    fout << "  \"symmetric_requires_permutation\": [\n";
    first=true;
    for(const auto& r: res){
        if(r.kind!=Kind::REQUIRES_PERM) continue;
        if(!first) fout << ",\n";
        first=false;
        dump_entry(r);
    }
    fout << "\n  ],\n\n";

    // non_symmetric
    fout << "  \"non_symmetric\": [\n";
    first=true;
    for(const auto& r: res){
        if(r.kind!=Kind::NONE) continue;
        if(!first) fout << ",\n";
        first=false;
        dump_entry(r);
    }
    fout << "\n  ]\n";
    fout << "}\n";

    fout.close();

    cerr << "DONE\n";
    cerr << "partitions_total=" << res.size() << "\n";
    cerr << "symmetric_translation=" << ct_trans << "\n";
    cerr << "symmetric_requires_permutation=" << ct_req << "\n";
    cerr << "non_symmetric=" << ct_none << "\n";
    cerr << "wrote: " << outname << "\n";

    return 0;
}
