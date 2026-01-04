import itertools

m = 5
Hmask = 0x7bb49071
full  = (1<<32) - 1
Cmask = full ^ Hmask

def ham(a,b): 
    return (a^b).bit_count()

def induced_graph(mask):
    V = [v for v in range(32) if (mask >> v) & 1]
    idx = {v:i for i,v in enumerate(V)}
    adj = [0]*len(V)
    for v in V:
        i = idx[v]
        for b in range(m):
            u = v ^ (1<<b)
            j = idx.get(u, -1)
            if j >= 0:
                adj[i] |= (1 << j)
    return V, adj

def edge_count(adj):
    return sum(a.bit_count() for a in adj)//2

def degree_seq(adj):
    return sorted(a.bit_count() for a in adj)

def find_iso(adjA, adjB):
    k = len(adjA)
    if edge_count(adjA) != edge_count(adjB): 
        return None
    if degree_seq(adjA) != degree_seq(adjB): 
        return None

    degA = [a.bit_count() for a in adjA]
    degB = [b.bit_count() for b in adjB]

    order = list(range(k))
    order.sort(key=lambda v: (-degA[v], v))

    usedB = [False]*k
    mapA  = [-1]*k

    cand_by_deg = {}
    for w in range(k):
        cand_by_deg.setdefault(degB[w], []).append(w)

    def dfs(pos):
        if pos == k:
            return True
        v = order[pos]
        for w in cand_by_deg.get(degA[v], []):
            if usedB[w]:
                continue
            ok = True
            for j in range(pos):
                u = order[j]
                mu = mapA[u]
                if mu < 0:
                    continue
                a_edge = (adjA[v] >> u) & 1
                b_edge = (adjB[w] >> mu) & 1
                if a_edge != b_edge:
                    ok = False
                    break
            if not ok:
                continue
            mapA[v] = w
            usedB[w] = True
            if dfs(pos+1):
                return True
            usedB[w] = False
            mapA[v] = -1
        return False

    if dfs(0):
        return mapA
    return None

def permute_bits(x, p):
    y = 0
    for i, pi in enumerate(p):
        if (x >> pi) & 1:
            y |= (1 << i)
    return y

def image_mask(mask, p, t):
    res = 0
    mm = mask
    while mm:
        lsb = mm & -mm
        v = (lsb.bit_length() - 1)
        img = permute_bits(v, p) ^ t
        res |= (1 << img)
        mm ^= lsb
    return res

HV, Hadj = induced_graph(Hmask)
CV, Cadj = induced_graph(Cmask)

print("H vertices:", HV)
print("C vertices:", CV)
print("Edges(H) =", edge_count(Hadj))
print("Edges(C) =", edge_count(Cadj))

iso_map = find_iso(Hadj, Cadj)
print("\nInduced-graph isomorphic?", iso_map is not None)

if iso_map is not None:
    phi = {HV[i]: CV[iso_map[i]] for i in range(16)}
    print("\nOne induced-graph isomorphism phi: H -> H^c")
    for a in sorted(phi):
        print(f"  {a:2d} -> {phi[a]:2d}")

    keys = sorted(phi)
    bad = None
    for i in range(len(keys)):
        for j in range(i+1, len(keys)):
            a,b = keys[i], keys[j]
            da = ham(a,b)
            db = ham(phi[a], phi[b])
            if da != db:
                bad = (a,b,da,phi[a],phi[b],db)
                break
        if bad:
            break

    if bad:
        a,b,da,fa,fb,db = bad
        print("\nWitness phi is NOT a cube automorphism (doesn't preserve Hamming distance):")
        print(f"  d({a},{b}) = {da}, but d(phi({a}),phi({b})) = d({fa},{fb}) = {db}")

found = None
for p in itertools.permutations(range(m)):
    for t in range(1<<m):
        if image_mask(Hmask, p, t) == Cmask:
            found = (p,t)
            break
    if found:
        break

print("\nExists cube automorphism swapping H and H^c?", found is not None)
