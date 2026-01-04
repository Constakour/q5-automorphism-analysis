// cube_halves_exhaustive_orbits.cpp
// Exhaustive orbit-by-orbit enumeration in Q5 (no Monte-Carlo):
//   - Enumerates 16-vertex halves H up to Aut(Q5) using isomorph-free generation
//   - Deduplicates at the unordered partition level {H,H^c}
//   - Counts partition-orbits that are (graphic AND non-symmetric)
//   - Uses the same vertex labeling + bit convention as your code (0..31, LSB=bit0)
//
// Compile:
//   g++ -O3 -std=c++20 -fopenmp cube_halves_exhaustive_orbits.cpp -o cube_halves_exhaustive
//
// Run:
//   ./cube_halves_exhaustive --exhaustive-q5 --threads 12

#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

static inline int popcount_u32(uint32_t x) { return __builtin_popcount(x); }
static inline int ctz_u32(uint32_t x)      { return __builtin_ctz(x); }

static constexpr int MAX_M = 16;

// -------------------- Small induced graph iso (k<=16) --------------------

struct SmallGraph {
    int k = 0;
    vector<uint16_t> adj; // adjacency bitmask rows

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
    for (int i = 0; i < k; i++) {
        col[i] = (int)(lower_bound(uniq.begin(), uniq.end(), deg[i]) - uniq.begin());
    }

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

        vector<int> newcol(k, 0);
        vector<int> order(k);
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

// returns only bool (fast). (You can extend to output mapping if you want certificates.)
static bool is_isomorphic_small(const SmallGraph& A, const SmallGraph& B) {
    if (A.k != B.k) return false;
    int k = A.k;

    if (A.edge_count() != B.edge_count()) return false;
    if (A.degrees_sorted() != B.degrees_sorted()) return false;

    vector<int> ca = wl_colors(A);
    vector<int> cb = wl_colors(B);
    {
        auto sa = ca; auto sb = cb;
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

static SmallGraph induced_graph_cube(int m, uint32_t mask) {
    int n = 1 << m;      // m=5 => 32
    int k = n / 2;       // => 16
    SmallGraph G;
    G.k = k;
    G.adj.assign(k, 0);

    vector<int> verts;
    verts.reserve(k);
    for (int v = 0; v < n; v++) if ((mask >> v) & 1u) verts.push_back(v);

    vector<int> idx(n, -1);
    for (int i = 0; i < k; i++) idx[verts[i]] = i;

    for (int i = 0; i < k; i++) {
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

// -------------------- Aut(Q5) enumeration, with precomputed vertex maps --------------------

static inline array<uint8_t, MAX_M> perm_identity(int m) {
    array<uint8_t, MAX_M> id{};
    for (int i = 0; i < m; i++) id[i] = (uint8_t)i;
    return id;
}

static inline uint32_t apply_P_to_vec(uint32_t v, const array<uint8_t, MAX_M>& p, int m) {
    // (P v)_i = v_{p[i]}    (LSB convention: bit i is coordinate i)
    uint32_t y = 0;
    for (int i = 0; i < m; i++) {
        if ((v >> p[i]) & 1u) y |= (1u << i);
    }
    return y;
}

struct AutoGroup {
    int m = 0;
    int n = 0;
    uint32_t full_mask = 0;

    struct Elem {
        array<uint8_t, MAX_M> p{};
        uint32_t c = 0;

        // precomputed vertex image: vmap[v] = g(v)
        array<uint8_t, 32> vmap{};
    };

    vector<Elem> elems;

    explicit AutoGroup(int m_) : m(m_), n(1<<m_) {
        full_mask = (n == 32) ? 0xFFFFFFFFu : ((1u << n) - 1u);

        vector<int> pvec(m);
        iota(pvec.begin(), pvec.end(), 0);

        do {
            array<uint8_t, MAX_M> p{};
            for (int i = 0; i < m; i++) p[i] = (uint8_t)pvec[i];

            for (uint32_t c = 0; c < (uint32_t)n; c++) {
                Elem e;
                e.p = p;
                e.c = c;

                for (int v = 0; v < n; v++) {
                    uint32_t pv = apply_P_to_vec((uint32_t)v, e.p, m);
                    uint32_t img = pv ^ e.c;
                    e.vmap[v] = (uint8_t)img;
                }

                elems.push_back(e);
            }
        } while (next_permutation(pvec.begin(), pvec.end()));
    }

    uint32_t image_mask(uint32_t mask, const Elem& e) const {
        uint32_t res = 0;
        uint32_t mm = mask;
        while (mm) {
            uint32_t lsb = mm & -mm;
            int v = ctz_u32(lsb);
            int img = e.vmap[v];
            res |= (1u << img);
            mm ^= lsb;
        }
        return res;
    }

    uint32_t orbit_rep(uint32_t mask) const {
        uint32_t best = UINT32_MAX;
        for (const auto& e : elems) best = min(best, image_mask(mask, e));
        return best;
    }

    // exact symmetry test: exists g with g(H) = H^c (no filters here; can add filters later)
    bool has_swapper(uint32_t mask) const {
        uint32_t comp = full_mask ^ mask;
        for (const auto& e : elems) {
            if (image_mask(mask, e) == comp) return true;
        }
        return false;
    }

    // setwise stabilizer indices for current mask: { g : g(mask)=mask }
    vector<int> stabilizer_indices(uint32_t mask) const {
        vector<int> stab;
        stab.reserve(elems.size());
        for (int i = 0; i < (int)elems.size(); i++) {
            if (image_mask(mask, elems[i]) == mask) stab.push_back(i);
        }
        return stab;
    }
};

// -------------------- orbit-by-orbit generator (isomorph-free) --------------------

struct WorkNode {
    uint32_t mask;
    int last;
    int depth;
};

static void build_frontier(const AutoGroup& G,
                           int target_depth,
                           vector<WorkNode>& out)
{
    // We enumerate in increasing vertex order to avoid duplicates from permutations of addition order.
    function<void(uint32_t,int,int)> rec = [&](uint32_t mask, int last, int depth) {
        if (depth == target_depth) {
            out.push_back({mask, last, depth});
            return;
        }
        auto stab = G.stabilizer_indices(mask);

        // Allowed candidates: vertices > last not in mask
        // We'll compute orbits of candidates under stab, and take 1 rep per orbit.
        array<char, 32> used{};
        used.fill(0);

        for (int v = last + 1; v < 32; v++) {
            if ((mask >> v) & 1u) continue;
            if (used[v]) continue;

            int minv = v;
            // mark whole orbit in "used"
            for (int idx : stab) {
                int y = G.elems[idx].vmap[v];
                if (y <= last) continue;
                if ((mask >> y) & 1u) continue;
                used[y] = 1;
                minv = min(minv, y);
            }
            // If v is not the minimum of its orbit among allowed candidates, skip.
            if (minv != v) continue;

            rec(mask | (1u << v), v, depth + 1);
        }
    };

    rec(0u, -1, 0);
}

// DFS from a node to depth 16; collect local counterexample orbit reps.
static void dfs_from_node(const AutoGroup& G,
                          uint32_t mask,
                          int last,
                          int depth,
                          unordered_set<uint32_t>& local_seen_partrep,
                          vector<pair<uint32_t,uint32_t>>& local_counter_orbits,
                          long long& local_leaves)
{
    if (depth == 16) {
        local_leaves++;

        uint32_t comp = 0xFFFFFFFFu ^ mask;

        // Fast invariants + iso test first
        SmallGraph A = induced_graph_cube(5, mask);
        SmallGraph B = induced_graph_cube(5, comp);
        if (A.edge_count() != B.edge_count()) return;
        if (A.degrees_sorted() != B.degrees_sorted()) return;
        {
            auto ca = wl_colors(A);
            auto cb = wl_colors(B);
            sort(ca.begin(), ca.end());
            sort(cb.begin(), cb.end());
            if (ca != cb) return;
        }
        if (!is_isomorphic_small(A, B)) return;

        // Now canonicalize to partition-rep and de-dup at partition level.
        uint32_t repH = G.orbit_rep(mask);
        uint32_t repC = G.orbit_rep(comp);
        uint32_t partrep = min(repH, repC);

        if (!local_seen_partrep.insert(partrep).second) return;

        // Symmetry test only after we know it is graphic
        bool sym = G.has_swapper(mask);
        if (!sym) {
            // record (partrep, witness_mask)
            local_counter_orbits.push_back({partrep, mask});
        }
        return;
    }

    auto stab = G.stabilizer_indices(mask);

    array<char, 32> used{};
    used.fill(0);

    for (int v = last + 1; v < 32; v++) {
        if ((mask >> v) & 1u) continue;
        if (used[v]) continue;

        int minv = v;
        for (int idx : stab) {
            int y = G.elems[idx].vmap[v];
            if (y <= last) continue;
            if ((mask >> y) & 1u) continue;
            used[y] = 1;
            minv = min(minv, y);
        }
        if (minv != v) continue;

        dfs_from_node(G, mask | (1u << v), v, depth + 1,
                      local_seen_partrep, local_counter_orbits, local_leaves);
    }
}

// -------------------- CLI --------------------

static void usage() {
    cerr << "Usage:\n"
         << "  ./cube_halves_exhaustive --exhaustive-q5 [--threads T]\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    bool exhaustive = false;
    int threads = 0;

    for (int i = 1; i < argc; i++) {
        string s = argv[i];
        if (s == "--exhaustive-q5") exhaustive = true;
        else if (s == "--threads" && i+1 < argc) threads = stoi(argv[++i]);
        else {
            cerr << "Unknown arg: " << s << "\n";
            usage();
            return 1;
        }
    }

    if (!exhaustive) {
        usage();
        return 1;
    }

#ifndef _OPENMP
    cerr << "This binary was compiled without OpenMP. Recompile with -fopenmp.\n";
    return 1;
#else
    if (threads > 0) omp_set_num_threads(threads);
    cerr << "OpenMP threads: " << omp_get_max_threads() << "\n";
#endif

    const int m = 5;
    AutoGroup G(m);

    // Split recursion into a frontier for parallel processing.
    // Depth 6 is a decent default (enough tasks, not too many).
    const int SPLIT_DEPTH = 6;

    vector<WorkNode> frontier;
    frontier.reserve(50000);
    build_frontier(G, SPLIT_DEPTH, frontier);
    cerr << "Frontier nodes at depth " << SPLIT_DEPTH << ": " << frontier.size() << "\n";

    // Global aggregation (merge after parallel phase)
    vector< unordered_set<uint32_t> > per_thread_seen;
    vector< vector<pair<uint32_t,uint32_t>> > per_thread_ce;
    vector<long long> per_thread_leaves;

    int T = 1;
#ifdef _OPENMP
    T = omp_get_max_threads();
#endif
    per_thread_seen.resize(T);
    per_thread_ce.resize(T);
    per_thread_leaves.assign(T, 0);

    using clock = chrono::steady_clock;
    auto t0 = clock::now();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < (int)frontier.size(); i++) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        auto& local_seen = per_thread_seen[tid];
        auto& local_ce   = per_thread_ce[tid];
        auto& local_lv   = per_thread_leaves[tid];

        const auto& node = frontier[i];
        dfs_from_node(G, node.mask, node.last, node.depth,
                      local_seen, local_ce, local_lv);
    }

    // Merge results
    unordered_set<uint32_t> global_seen_partrep;
    global_seen_partrep.reserve(200000);

    unordered_map<uint32_t,uint32_t> global_ce; // partrep -> witness
    global_ce.reserve(64);

    long long total_leaves = 0;
    for (int tid = 0; tid < T; tid++) {
        total_leaves += per_thread_leaves[tid];

        // merge seen_partrep
        for (auto x : per_thread_seen[tid]) global_seen_partrep.insert(x);

        // merge counterexamples (dedup by partrep)
        for (auto &kv : per_thread_ce[tid]) {
            uint32_t partrep = kv.first;
            uint32_t wit = kv.second;
            auto it = global_ce.find(partrep);
            if (it == global_ce.end()) global_ce.emplace(partrep, wit);
            else it->second = min(it->second, wit); // keep a deterministic witness
        }
    }

    auto t1 = clock::now();
    double elapsed = chrono::duration<double>(t1 - t0).count();

    // Print catalog (sorted)
    vector<pair<uint32_t,uint32_t>> ce_list;
    ce_list.reserve(global_ce.size());
    for (auto &kv : global_ce) ce_list.push_back(kv);
    sort(ce_list.begin(), ce_list.end(),
         [](auto& a, auto& b){ return a.first < b.first; });

    cout << "\nDONE (exhaustive orbit enumeration)\n";
    cout << "elapsed_seconds=" << fixed << setprecision(3) << elapsed << "\n";
    cout << "frontier_nodes=" << frontier.size() << " (split_depth=" << SPLIT_DEPTH << ")\n";
    cout << "dfs_leaves_reached=" << total_leaves << " (includes duplicates across partition level)\n";
    cout << "unique_partition_orbits_processed=" << global_seen_partrep.size() << "\n";

    cout << "\nCOUNTEREXAMPLE PARTITION-ORBIT REPS (graphic & non-symmetric)\n";
    cout << "count=" << ce_list.size() << "\n";
    for (auto &kv : ce_list) {
        uint32_t partrep = kv.first;
        uint32_t mask = kv.second;
        uint32_t comp = 0xFFFFFFFFu ^ mask;
        cout << "  partrep=0x" << hex << partrep
             << "  witness_mask=0x" << mask
             << "  witness_comp=0x" << comp << dec << "\n";
    }

    return 0;
}
