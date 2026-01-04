#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

static constexpr int m = 5;
static constexpr int N = 32;                 // vertices 0..31
static constexpr uint32_t FULL = 0xFFFFFFFF; // all 32 bits

// -------------------- helpers --------------------

static inline int popcount_u32(uint32_t x) { return __builtin_popcount(x); }

static string hexmask(uint32_t x){
    char buf[16];
    sprintf(buf, "0x%08X", x);
    return string(buf);
}

static vector<int> mask_to_vertices(uint32_t S){
    vector<int> verts;
    verts.reserve(16);
    for(int v=0; v<N; v++) if (S & (1u<<v)) verts.push_back(v);
    return verts;
}

// -------------------- induced subgraph on a 16-vertex half (k<=16) --------------------

struct SmallGraph {
    int k = 0;                // <= 16
    vector<uint16_t> adj;     // adjacency bitmask rows (size k)

    int edge_count() const {
        int s = 0;
        for (int i = 0; i < k; i++) s += __builtin_popcount((unsigned)adj[i]);
        return s / 2;
    }
    vector<int> degrees_sorted() const {
        vector<int> d(k);
        for (int i = 0; i < k; i++) d[i] = __builtin_popcount((unsigned)adj[i]);
        sort(d.begin(), d.end());
        return d;
    }
};

static vector<int> wl_colors(const SmallGraph& G) {
    int k = G.k;
    vector<int> col(k, 0);

    vector<int> deg(k);
    for (int i = 0; i < k; i++) deg[i] = __builtin_popcount((unsigned)G.adj[i]);

    vector<int> uniq = deg;
    sort(uniq.begin(), uniq.end());
    uniq.erase(unique(uniq.begin(), uniq.end()), uniq.end());
    for (int i = 0; i < k; i++)
        col[i] = (int)(lower_bound(uniq.begin(), uniq.end(), deg[i]) - uniq.begin());

    for (int it = 0; it < 20; it++) {
        vector<vector<int>> sig(k);
        for (int v = 0; v < k; v++) {
            vector<int> neigh_cols;
            uint16_t nb = G.adj[v];
            while (nb) {
                uint16_t lsb = nb & (uint16_t)(- (int16_t)nb);
                int u = __builtin_ctz((unsigned)lsb);
                neigh_cols.push_back(col[u]);
                nb ^= lsb;
            }
            sort(neigh_cols.begin(), neigh_cols.end());
            sig[v].push_back(col[v]);
            sig[v].insert(sig[v].end(), neigh_cols.begin(), neigh_cols.end());
        }

        vector<int> newcol(k, 0), order(k);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b){ return sig[a] < sig[b]; });

        int c = 0;
        for (int idx = 0; idx < k; idx++) {
            if (idx > 0 && sig[order[idx]] != sig[order[idx-1]]) c++;
            newcol[order[idx]] = c;
        }
        if (newcol == col) break;
        col.swap(newcol);
    }
    return col;
}

static bool is_isomorphic_small(const SmallGraph& A, const SmallGraph& B) {
    if (A.k != B.k) return false;
    int k = A.k;

    if (A.edge_count() != B.edge_count()) return false;
    if (A.degrees_sorted() != B.degrees_sorted()) return false;

    vector<int> ca = wl_colors(A), cb = wl_colors(B);
    {
        auto sa = ca, sb = cb;
        sort(sa.begin(), sa.end());
        sort(sb.begin(), sb.end());
        if (sa != sb) return false;
    }

    int maxc = 0;
    for (int x : ca) maxc = max(maxc, x);
    for (int x : cb) maxc = max(maxc, x);
    int C = maxc + 1;

    vector<vector<int>> classB(C);
    for (int v = 0; v < k; v++) classB[cb[v]].push_back(v);

    vector<int> degA(k), degB(k);
    for (int i = 0; i < k; i++) {
        degA[i] = __builtin_popcount((unsigned)A.adj[i]);
        degB[i] = __builtin_popcount((unsigned)B.adj[i]);
    }

    vector<int> order(k);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int u, int v){
        int cu = ca[u], cv = ca[v];
        if ((int)classB[cu].size() != (int)classB[cv].size())
            return classB[cu].size() < classB[cv].size();
        if (degA[u] != degA[v]) return degA[u] > degA[v];
        return u < v;
    });

    vector<int> mapA(k, -1);
    vector<char> usedB(k, 0);

    function<bool(int)> dfs = [&](int idx) -> bool {
        if (idx == k) return true;
        int v = order[idx];
        int colv = ca[v];

        for (int w : classB[colv]) {
            if (usedB[w]) continue;
            if (degA[v] != degB[w]) continue;

            bool ok = true;
            for (int j = 0; j < idx; j++) {
                int u = order[j];
                int mu = mapA[u];
                if (mu < 0) continue;
                bool a_edge = (A.adj[v] >> u) & 1;
                bool b_edge = (B.adj[w] >> mu) & 1;
                if (a_edge != b_edge) { ok = false; break; }
            }
            if (!ok) continue;

            mapA[v] = w;
            usedB[w] = 1;
            if (dfs(idx + 1)) return true;
            usedB[w] = 0;
            mapA[v] = -1;
        }
        return false;
    };

    return dfs(0);
}

static SmallGraph induced_graph_cube(uint32_t mask) {
    SmallGraph G;
    G.k = 16;
    G.adj.assign(16, 0);

    vector<int> verts;
    verts.reserve(16);
    for (int v = 0; v < N; v++) if ((mask >> v) & 1u) verts.push_back(v);

    int idx[N];
    for (int i = 0; i < N; i++) idx[i] = -1;
    for (int i = 0; i < 16; i++) idx[verts[i]] = i;

    for (int i = 0; i < 16; i++) {
        int v = verts[i];
        uint16_t row = 0;
        for (int bit = 0; bit < m; bit++) {
            int u = v ^ (1 << bit);
            int j = idx[u];
            if (j >= 0) row |= (uint16_t)(1u << j);
        }
        G.adj[i] = row;
    }
    return G;
}

// -------------------- Aut(Q5) as (flip XOR) + (coordinate permutation) --------------------
// perm[b] = oldBit index that becomes new bit b

static vector< array<int,5> > ALL_PERMS;

static void build_all_perms_once(){
    vector<int> base{0,1,2,3,4};
    do {
        array<int,5> p;
        for(int i = 0; i < 5; i++) p[i] = base[i];
        ALL_PERMS.push_back(p);
    } while(next_permutation(base.begin(), base.end()));
}

static inline uint32_t applyXor(uint32_t S, uint32_t flipMask){
    uint32_t T = 0, mm = S;
    while(mm){
        uint32_t lsb = mm & -mm;
        int v = __builtin_ctz(lsb);
        int x = v ^ (int)flipMask;
        T |= (1u << x);
        mm ^= lsb;
    }
    return T;
}

static inline int permute_bits(int v, const array<int,5> &p){
    int newV = 0;
    for (int b = 0; b < 5; b++){
        int oldBit = (v >> p[b]) & 1;
        newV |= (oldBit << b);
    }
    return newV; // in 0..31
}

static inline uint32_t applyPermutation(uint32_t S, const array<int,5> &p){
    uint32_t P = 0, mm = S;
    while(mm){
        uint32_t lsb = mm & -mm;
        int v = __builtin_ctz(lsb);
        int newV = permute_bits(v, p);
        P |= (1u << newV);
        mm ^= lsb;
    }
    return P;
}

// The actual automorphism on vertices is:
// g(v) = permute_bits( v XOR flip, perm )
static inline int apply_automorphism_to_vertex(int v, uint32_t flip, const array<int,5>& perm){
    return permute_bits(v ^ (int)flip, perm);
}

static inline array<int,32> build_vertex_map(uint32_t flip, const array<int,5>& perm){
    array<int,32> mp{};
    for(int v=0; v<32; v++) mp[v] = apply_automorphism_to_vertex(v, flip, perm);
    return mp;
}

// translation-only: exists flip with flip(A)=~A
static bool checkSingleTranslation(uint32_t A, uint32_t &out_flip){
    uint32_t B = (~A) & FULL;
    array<int,5> id{0,1,2,3,4};
    for(uint32_t flip=0; flip<32; flip++){
        // translation-only means perm is identity, so mapping is A^flip
        if (applyXor(A, flip) == B) {
            out_flip = flip;
            return true;
        }
    }
    return false;
}

// full: exists perm,flip with perm(flip(A))=~A
static bool checkFlipPerm(uint32_t A, uint32_t &out_flip, array<int,5> &out_perm){
    uint32_t B = (~A) & FULL;
    for(const auto &p : ALL_PERMS){
        for(uint32_t flip = 0; flip < 32; flip++){
            uint32_t X = applyXor(A, flip);
            uint32_t T = applyPermutation(X, p);
            if(T == B){
                out_flip = flip;
                out_perm = p;
                return true;
            }
        }
    }
    return false;
}

// -------------------- enumerate pairs-only partitions (6435) --------------------

struct Pair { int a,b; };

// Gosper next comb with same popcount
static inline unsigned long long nextCombPopcount(unsigned long long x){
    if(x == 0ULL) return 0ULL;
    unsigned long long u = x & -x;
    unsigned long long v = u + x;
    if(v == 0ULL) return 0ULL;
    return v + (((v ^ x)/u)>>2);
}

static unsigned long long binom(int n,int k){
    long double r = 1.0;
    for(int i=1; i<=k; i++) r *= (n-k+i)/(long double)i;
    return (unsigned long long) llroundl(r);
}

struct Entry {
    // canonical unordered partition representation:
    uint32_t half_lo = 0;
    uint32_t half_hi = 0;

    bool graphic = false;

    // symmetry classification
    bool symmetric = false;
    bool sym_translation = false;

    // witness parameters (ALWAYS maps half_lo -> half_hi)
    uint32_t flip = 0;
    array<int,5> perm{0,1,2,3,4};
};

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int threads = 0;
    string outname = "q5_pairs_only_partitions_split.json";

    for (int i=1;i<argc;i++){
        string s = argv[i];
        if (s=="--threads" && i+1<argc) threads = stoi(argv[++i]);
        else if (s=="--out" && i+1<argc) outname = argv[++i];
        else {
            cerr << "Usage: ./q5_pairs_only_partitions_split [--threads T] [--out FILE]\n";
            return 1;
        }
    }

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    cerr << "OpenMP threads: " << omp_get_max_threads() << "\n";
#else
    cerr << "Warning: compiled without OpenMP.\n";
#endif

    build_all_perms_once();

    // build antipodal pairs in Q5: (v, v^31)
    vector<bool> used(N,false);
    vector<Pair> pairsList;
    for(int i=0;i<N;i++){
        if(!used[i]){
            int j = i ^ 31;
            used[i]=true; used[j]=true;
            pairsList.push_back({ min(i,j), max(i,j) });
        }
    }
    if(pairsList.size()!=16){
        cerr<<"Error building antipodal pairs.\n";
        return 1;
    }

    // Fix one pair (the one containing 31) so each partition is counted once
    int fixPairIdx=-1;
    for(int i=0;i<16;i++){
        if(pairsList[i].b==31){ fixPairIdx=i; break; }
    }
    if(fixPairIdx<0){
        cerr<<"No pair containing vertex31??\n";
        return 1;
    }

    vector<uint32_t> pairMask(16,0U);
    for(int i=0;i<16;i++){
        pairMask[i] = (1u<<pairsList[i].a) | (1u<<pairsList[i].b);
    }

    vector<int> idxList; idxList.reserve(15);
    for(int i=0;i<16;i++) if(i!=fixPairIdx) idxList.push_back(i);

    auto buildSetMask = [&](unsigned long long combBits){
        uint32_t S=0U;
        while(combBits){
            unsigned long long lsb = combBits & -combBits;
            int b = __builtin_ctzll(combBits);
            combBits ^= lsb;
            int pairIdx = idxList[b];
            S |= pairMask[pairIdx];
        }
        return S;
    };

    const unsigned long long total = binom(15,8); // 6435
    unsigned long long startC = (1ULL<<8)-1ULL;
    unsigned long long endC   = startC<<(15-8);

    vector<unsigned long long> combs;
    combs.reserve(total);
    {
        unsigned long long comb = startC;
        while(true){
            combs.push_back(comb);
            if(comb==endC) break;
            unsigned long long nxt = nextCombPopcount(comb);
            if(nxt==0ULL) break;
            comb = nxt;
        }
        if(combs.size()!=total){
            cerr << "Internal error: expected " << total << " combs, got " << combs.size() << "\n";
            return 1;
        }
    }

    vector<Entry> entries(total);

    atomic<long long> cnt_trans(0), cnt_flipperm(0), cnt_nonsym(0);
    atomic<long long> cnt_graphic(0), cnt_nonsym_graphic(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long i=0;i<(long long)total;i++){
        uint32_t S = buildSetMask(combs[(size_t)i]);
        uint32_t comp = (~S) & FULL;

        // canonical unordered pair
        uint32_t A = min(S, comp);
        uint32_t B = max(S, comp);

        Entry e;
        e.half_lo = A;
        e.half_hi = B;

        // symmetry classification + witness ALWAYS from A -> B
        uint32_t flipT=0;
        if (checkSingleTranslation(A, flipT)){
            e.symmetric = true;
            e.sym_translation = true;
            e.flip = flipT;
            e.perm = {0,1,2,3,4};
            cnt_trans.fetch_add(1, memory_order_relaxed);
        } else {
            uint32_t flip=0; array<int,5> perm;
            if (checkFlipPerm(A, flip, perm)){
                e.symmetric = true;
                e.sym_translation = false;
                e.flip = flip;
                e.perm = perm;
                cnt_flipperm.fetch_add(1, memory_order_relaxed);
            } else {
                e.symmetric = false;
                cnt_nonsym.fetch_add(1, memory_order_relaxed);
            }
        }

        // graphic test on A vs B
        SmallGraph GA = induced_graph_cube(A);
        SmallGraph GB = induced_graph_cube(B);

        bool isGraphic = false;
        if (GA.edge_count() == GB.edge_count() && GA.degrees_sorted() == GB.degrees_sorted()) {
            auto ca = wl_colors(GA);
            auto cb = wl_colors(GB);
            sort(ca.begin(), ca.end());
            sort(cb.begin(), cb.end());
            if (ca == cb) isGraphic = is_isomorphic_small(GA,GB);
        }
        e.graphic = isGraphic;

        if (isGraphic) cnt_graphic.fetch_add(1, memory_order_relaxed);
        if (!e.symmetric && isGraphic) cnt_nonsym_graphic.fetch_add(1, memory_order_relaxed);

        entries[(size_t)i] = e;
    }

    // bucket by symmetry type
    vector<Entry> bucket_trans, bucket_flipperm, bucket_nonsym;
    bucket_trans.reserve((size_t)cnt_trans.load());
    bucket_flipperm.reserve((size_t)cnt_flipperm.load());
    bucket_nonsym.reserve((size_t)cnt_nonsym.load());

    for (const auto &e : entries){
        if (e.symmetric && e.sym_translation) bucket_trans.push_back(e);
        else if (e.symmetric) bucket_flipperm.push_back(e);
        else bucket_nonsym.push_back(e);
    }

    ofstream fout(outname);
    if(!fout){
        cerr << "Cannot open output file: " << outname << "\n";
        return 1;
    }

    auto write_entry = [&](const Entry& e, int indent){
        string sp(indent, ' ');
        auto loV = mask_to_vertices(e.half_lo);
        auto hiV = mask_to_vertices(e.half_hi);

        fout << sp << "{\n";
        fout << sp << "  \"partition\": {\n";
        fout << sp << "    \"half_lo_mask_hex\": \"" << hexmask(e.half_lo) << "\",\n";
        fout << sp << "    \"half_hi_mask_hex\": \"" << hexmask(e.half_hi) << "\",\n";
        fout << sp << "    \"half_lo_vertices\": [";
        for (size_t i=0;i<loV.size();i++){ if(i) fout << ", "; fout << loV[i]; }
        fout << "],\n";
        fout << sp << "    \"half_hi_vertices\": [";
        for (size_t i=0;i<hiV.size();i++){ if(i) fout << ", "; fout << hiV[i]; }
        fout << "]\n";
        fout << sp << "  },\n";

        fout << sp << "  \"graphic\": " << (e.graphic ? "true" : "false") << ",\n";

        if (e.symmetric){
            // Explicit vertex map for unambiguous witness: maps half_lo -> half_hi
            auto vmap = build_vertex_map(e.flip, e.perm);

            fout << sp << "  \"witness\": {\n";
            fout << sp << "    \"maps\": \"half_lo -> half_hi\",\n";
            fout << sp << "    \"flip_mask\": " << e.flip << ",\n";
            fout << sp << "    \"perm\": [" << e.perm[0] << ", " << e.perm[1] << ", " << e.perm[2] << ", " << e.perm[3] << ", " << e.perm[4] << "],\n";

            // Full map on all 32 vertices: vertex_map[v] = g(v)
            fout << sp << "    \"vertex_map\": [";
            for(int v=0; v<32; v++){
                if (v) fout << ", ";
                fout << vmap[v];
            }
            fout << "],\n";

            // Restriction to half_lo (pairs list)
            fout << sp << "    \"half_lo_map_pairs\": [";
            bool first = true;
            for(int v=0; v<32; v++){
                if (e.half_lo & (1u<<v)){
                    if (!first) fout << ", ";
                    first = false;
                    fout << "[" << v << ", " << vmap[v] << "]";
                }
            }
            fout << "]\n";

            fout << sp << "  }\n";
        } else {
            fout << sp << "  \"note\": \"no cube automorphism swaps the two halves\"\n";
        }

        fout << sp << "}";
    };

    auto write_array = [&](const string& name, const vector<Entry>& bucket){
        fout << "  \"" << name << "\": [\n";
        for (size_t i=0;i<bucket.size();i++){
            write_entry(bucket[i], 4);
            if (i+1<bucket.size()) fout << ",";
            fout << "\n";
        }
        fout << "  ]";
    };

    long long trans   = cnt_trans.load();
    long long flipperm= cnt_flipperm.load();
    long long nonsym  = cnt_nonsym.load();
    long long graphic = cnt_graphic.load();
    long long nonsym_graphic = cnt_nonsym_graphic.load();

    fout << "{\n";
    fout << "  \"meta\": {\n";
    fout << "    \"graph\": \"Q5\",\n";
    fout << "    \"class\": \"pairs-only partitions (8 antipodal pairs)\",\n";
    fout << "    \"count_partitions\": " << total << "\n";
    fout << "  },\n";
    fout << "  \"summary\": {\n";
    fout << "    \"translation\": " << trans << ",\n";
    fout << "    \"flip_perm\": " << flipperm << ",\n";
    fout << "    \"non_symmetric\": " << nonsym << ",\n";
    fout << "    \"graphic_total\": " << graphic << ",\n";
    fout << "    \"non_symmetric_graphic\": " << nonsym_graphic << "\n";
    fout << "  },\n";

    write_array("translation", bucket_trans); fout << ",\n";
    write_array("flip_perm", bucket_flipperm); fout << ",\n";
    write_array("non_symmetric", bucket_nonsym); fout << "\n";
    fout << "}\n";
    fout.close();

    cout << "DONE\n";
    cout << "partitions_total=" << total << "\n";
    cout << "translation=" << trans << "\n";
    cout << "flip_perm=" << flipperm << "\n";
    cout << "non_symmetric=" << nonsym << "\n";
    cout << "graphic_total=" << graphic << "\n";
    cout << "non_symmetric_graphic=" << nonsym_graphic << "\n";
    cout << "wrote: " << outname << "\n";

    return 0;
}
