#!/usr/bin/env python3
"""
verify_witness.py

Given a witness_mask (for a half H in Q5), verify:
  1) H has size 16 (half)
  2) Q5[H] is isomorphic to Q5[H^c]  (graphic)
  3) whether there exists a cube automorphism g(x)=P x XOR c with g(H)=H^c (swappable)
Optionally:
  4) print an explicit isomorphism phi: H -> H^c if graphic.
  5) when --print-phi is used, ALSO print one explicit Hamming-distance mismatch
     (if found) certifying that phi is not a cube automorphism.

This is intended as a "checkable certificate" companion to your C++ output.

Conventions (matches your C++):
- Vertices are integers 0..31, LSB is coordinate 0.
- Neighbor edges in Q5 are v <-> v XOR (1<<bit), bit in {0..4}.
- An automorphism is g(v)=P(v) XOR c where:
    (P(v))_i = v_{perm[i]}   using LSB coordinate indexing.

Usage examples:
  python3 verify_witness.py 0x2695b9f
  python3 verify_witness.py 0x2695b9f --print-phi
  python3 verify_witness.py 0x2695b9f --json-line 'partrep=... witness_mask=0x...'
  python3 verify_witness.py 0x2695b9f --find-swapper   (prints a swapper if one exists)
"""

import argparse
import itertools
from collections import defaultdict

N = 5
NUM_VERTICES = 1 << N
FULL_MASK_Q5 = 0xFFFFFFFF


# ------------------------
# Basic mask / set helpers
# ------------------------

def parse_mask(s: str) -> int:
    s = s.strip().lower()
    if s.startswith("0x"):
        return int(s, 16) & 0xFFFFFFFF
    return int(s, 10) & 0xFFFFFFFF


def extract_hex_from_text(text: str):
    """
    Convenience: pull the first 0x... substring from a line if user pastes
    'partrep=... witness_mask=0x....'
    """
    import re
    m = re.search(r"0x[0-9a-fA-F]+", text)
    if not m:
        raise ValueError("No hex literal like 0x... found in the provided text.")
    return m.group(0)


def mask_to_vertices(mask: int):
    return [v for v in range(NUM_VERTICES) if (mask >> v) & 1]


def popcount32(x: int) -> int:
    return int(x & 0xFFFFFFFF).bit_count()


# ------------------------
# Hamming distance sanity check helpers
# ------------------------

def wt(x: int) -> int:
    """Hamming weight (popcount)."""
    return x.bit_count()


def hamming_dist(u: int, v: int) -> int:
    """Hamming distance in Q5: d(u,v)=wt(u XOR v)."""
    return wt(u ^ v)


def find_distance_mismatch_in_H(H_list, phi_map):
    """
    Look for u,v in H_list such that:
        d(u,v) != d(phi(u), phi(v)).
    If found, return (u, v, d_uv, phi_u, phi_v, d_phi).
    Otherwise return None.
    """
    n = len(H_list)
    for i in range(n):
        u = H_list[i]
        if u not in phi_map:
            continue
        pu = phi_map[u]
        for j in range(i + 1, n):
            v = H_list[j]
            if v not in phi_map:
                continue
            pv = phi_map[v]
            duv = hamming_dist(u, v)
            dp = hamming_dist(pu, pv)
            if duv != dp:
                return (u, v, duv, pu, pv, dp)
    return None


# ------------------------
# Q5 induced graph builder
# ------------------------

def induced_adj_rows(mask: int):
    """
    Return:
      verts: list of the 16 vertices in increasing order
      idx:   dict vertex->0..15
      rows:  list of 16 integers, where bit j of rows[i] indicates adjacency
             between verts[i] and verts[j] inside the induced subgraph.
    """
    verts = mask_to_vertices(mask)
    idx = {v: i for i, v in enumerate(verts)}
    rows = [0] * len(verts)

    for i, v in enumerate(verts):
        row = 0
        for bit in range(N):
            u = v ^ (1 << bit)
            j = idx.get(u, -1)
            if j >= 0:
                row |= (1 << j)
        rows[i] = row
    return verts, idx, rows


def edge_count(rows):
    # sum degrees / 2
    return sum(r.bit_count() for r in rows) // 2


def degree_multiset(rows):
    return sorted(r.bit_count() for r in rows)


# ------------------------
# Explicit iso (backtracking)
# ------------------------

def find_isomorphism(rowsA, rowsB):
    """
    Try to find an explicit isomorphism between two graphs on k vertices.
    Inputs are adjacency rows as bitmasks (length k).

    Returns:
      mapping list mapA_to_B of length k, or None if not isomorphic.
    """
    k = len(rowsA)
    if k != len(rowsB):
        return None
    if edge_count(rowsA) != edge_count(rowsB):
        return None
    if degree_multiset(rowsA) != degree_multiset(rowsB):
        return None

    degA = [r.bit_count() for r in rowsA]
    degB = [r.bit_count() for r in rowsB]

    # Partition vertices by degree for pruning
    bucketB = defaultdict(list)
    for j, d in enumerate(degB):
        bucketB[d].append(j)

    # order A vertices: higher degree first, then smaller index
    order = list(range(k))
    order.sort(key=lambda i: (-degA[i], i))

    usedB = [False] * k
    mapA = [-1] * k

    def compatible(a, b, upto_idx):
        # Check adjacency consistency between a and previously assigned vertices in order[:upto_idx]
        for t in range(upto_idx):
            u = order[t]
            mu = mapA[u]
            if mu < 0:
                continue
            a_edge = (rowsA[a] >> u) & 1
            b_edge = (rowsB[b] >> mu) & 1
            if a_edge != b_edge:
                return False
        return True

    def dfs(pos):
        if pos == k:
            return True
        a = order[pos]
        d = degA[a]
        for b in bucketB[d]:
            if usedB[b]:
                continue
            if not compatible(a, b, pos):
                continue
            mapA[a] = b
            usedB[b] = True
            if dfs(pos + 1):
                return True
            usedB[b] = False
            mapA[a] = -1
        return False

    ok = dfs(0)
    return mapA if ok else None


# ------------------------
# Aut(Q5) and swappers
# ------------------------

def apply_perm_to_vertex(v: int, perm):
    """
    perm[i] means: new bit i gets old bit perm[i]
    i.e. (P v)_i = v_{perm[i]} with LSB indexing.
    """
    out = 0
    for i in range(N):
        bit = (v >> perm[i]) & 1
        out |= (bit << i)
    return out


def apply_affine_to_mask(mask: int, perm, c: int) -> int:
    res = 0
    mm = mask
    while mm:
        lsb = mm & -mm
        v = (lsb.bit_length() - 1)
        img = apply_perm_to_vertex(v, perm) ^ c
        res |= (1 << img)
        mm ^= lsb
    return res & 0xFFFFFFFF


def find_swapper(mask: int):
    """
    Search all 3840 automorphisms g(v)=P(v) XOR c to see if g(H)=H^c.
    Returns (perm, c) if found, else None.
    """
    comp = FULL_MASK_Q5 ^ mask
    for perm in itertools.permutations(range(N)):
        for c in range(NUM_VERTICES):
            img = apply_affine_to_mask(mask, perm, c)
            if img == comp:
                return list(perm), c
    return None


def is_translation_only_swap(mask: int):
    """
    Check if there exists c such that (v -> v XOR c) maps H to H^c.
    Equivalent: translated mask equals complement.
    """
    comp = FULL_MASK_Q5 ^ mask
    ident = list(range(N))
    for c in range(NUM_VERTICES):
        img = apply_affine_to_mask(mask, ident, c)
        if img == comp:
            return c
    return None


# ------------------------
# Main
# ------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mask", nargs="?", help="witness_mask as hex (0x...) or decimal.")
    ap.add_argument("--json-line", help="Alternatively paste a whole output line; script will extract first 0x... as mask.")
    ap.add_argument("--print-phi", action="store_true", help="If graphic, print an explicit isomorphism phi: H -> H^c.")
    ap.add_argument("--find-swapper", action="store_true", help="Search and print a swapper (perm,c) if one exists.")
    ap.add_argument("--quick", action="store_true", help="Skip printing vertex lists (just summary).")
    args = ap.parse_args()

    if args.json_line:
        mask_str = extract_hex_from_text(args.json_line)
        mask = parse_mask(mask_str)
    elif args.mask:
        mask = parse_mask(args.mask)
    else:
        ap.error("Provide a mask (e.g. 0x2695b9f) or use --json-line ...")

    comp = FULL_MASK_Q5 ^ mask

    # decode
    H = mask_to_vertices(mask)
    Hc = mask_to_vertices(comp)

    print("=== verify_witness (Q5) ===")
    print(f"witness_mask = 0x{mask:08x}  popcount={len(H)}")
    print(f"witness_comp = 0x{comp:08x}  popcount={len(Hc)}")

    if not args.quick:
        print("\nH  =", H)
        print("H^c=", Hc)

    # build induced graphs
    vertsA, _, rowsA = induced_adj_rows(mask)
    vertsB, _, rowsB = induced_adj_rows(comp)

    if len(vertsA) != 16 or len(vertsB) != 16:
        print("\nWARNING: mask is not a half (size 16). Induced-graph tests may be meaningless.")
        return

    # graphic test (and optional explicit isomorphism)
    phi_idx = find_isomorphism(rowsA, rowsB)
    is_graphic = (phi_idx is not None)
    print("\nGraphic test: Q5[H] ≅ Q5[H^c] ?  ", "YES" if is_graphic else "NO")

    if is_graphic and args.print_phi:
        print("\nIsomorphism certificate φ : H → H^c (vertex labels):")

        # Build vertex-level map phi_map: vertex in H -> vertex in H^c
        phi_map = {}
        mapping_pairs = []
        for i, j in enumerate(phi_idx):
            a = vertsA[i]
            b = vertsB[j]
            phi_map[a] = b
            mapping_pairs.append((a, b))
        mapping_pairs.sort(key=lambda x: x[0])

        for a, b in mapping_pairs:
            print(f"  {a:2d} -> {b:2d}")

        # NEW: print one human-checkable Hamming-distance mismatch (if any)
        mismatch = find_distance_mismatch_in_H(H, phi_map)
        if mismatch is None:
            print("\nDistance sanity check on pairs in H: no mismatch found.")
            print("Note: this does NOT prove φ extends to a cube automorphism;")
            print("it only means φ happened to preserve distances on H×H.")
        else:
            (u, v, duv, pu, pv, dp) = mismatch
            print("\nDistance sanity check (human-checkable obstruction):")
            print("FOUND distance mismatch, so φ cannot be a cube automorphism:")
            print(f"  u={u}, v={v}")
            print(f"  φ(u)={pu}, φ(v)={pv}")
            print(f"  d(u,v) = {duv}")
            print(f"  d(φ(u),φ(v)) = {dp}")

    # symmetry (swapper) test
    c_trans = is_translation_only_swap(mask)
    if c_trans is not None:
        print("\nSwappability: YES (translation-only)")
        print(f"  swapper: g(v)=v XOR c with c={c_trans} (0x{c_trans:02x})")
    else:
        print("\nSwappability: translation-only?  NO")

    if args.find_swapper:
        res = find_swapper(mask)
        if res is None:
            print("Swappability: YES/NO under full Aut(Q5)?  NO (no (perm,c) exists)")
        else:
            perm, c = res
            print("Swappability: YES under full Aut(Q5)")
            print(f"  c (flipmask) = {c} (0x{c:02x})")
            print(f"  perm         = {perm}")
            print("  note: (P v)_i = v_{perm[i]} in LSB=coord0 convention")

    print("\nDone.")


if __name__ == "__main__":
    main()
