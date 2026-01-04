#!/usr/bin/env python3
"""
q5_pairs_only_partitions_split_verify.py

Verifies partitions from q5_pairs_only_partitions_split.json (NEW schema)
produced by q5_pairs_only_full_classification_split.cpp.

Each entry (in translation[] or flip_perm[]) has the form:
  {
    "partition": {
      "half_lo_mask_hex": "...",
      "half_hi_mask_hex": "...",
      "half_lo_vertices": [...16...],
      "half_hi_vertices": [...16...]
    },
    "graphic": true/false,
    "witness": {
      "maps": "half_lo -> half_hi",
      "flip_mask": 0..31,
      "perm": [p0,p1,p2,p3,p4],
      "vertex_map": [...32...],              # optional
      "half_lo_map_pairs": [[u,g(u)], ...]   # optional
    }
  }

We check:
  (1) half_hi is exactly the complement of half_lo (on {0..31})
  (2) witness really maps half_lo -> half_hi
  (3) if vertex_map exists, it matches flip+perm
  (4) if half_lo_map_pairs exists, it's consistent
  (5) optional: summary counts match array sizes

Run:
  python q5_pairs_only_partitions_split_verify.py q5_pairs_only_partitions_split.json

Optional detailed inspection:
  inspect_entry("flip_perm", 0)  # category in {"translation","flip_perm","non_symmetric"}
"""

import json
import sys

N = 5
NUM_VERTICES = 1 << N  # 32

##############################
# Helper Functions
##############################

def mask_to_set(mask_int):
    """Convert a 32-bit integer mask to a set of vertices {0..31}."""
    S = set()
    for v in range(NUM_VERTICES):
        if (mask_int & (1 << v)) != 0:
            S.add(v)
    return S

def to_binary5(x):
    """Return x as a 5-bit binary string, e.g. x=13 => '01101'."""
    return f"{x:05b}"

def permute_bits(val, perm):
    """
    perm[i] means: oldBit = (val >> perm[i]) & 1 becomes newBit i.
    This matches the C++ applyPermutation.
    """
    newv = 0
    for i in range(N):
        old_bit = (val >> perm[i]) & 1
        newv |= (old_bit << i)
    return newv & 0x1F

def apply_witness_to_vertex(v, flip_mask, perm):
    """g(v) = permute_bits(v XOR flip_mask, perm)."""
    return permute_bits((v ^ flip_mask) & 0x1F, perm)

def apply_witness_to_set(S, flip_mask, perm):
    return {apply_witness_to_vertex(v, flip_mask, perm) for v in S}

def compute_vertex_map(flip_mask, perm):
    return [apply_witness_to_vertex(v, flip_mask, perm) for v in range(NUM_VERTICES)]

def format_set_bitpatterns(num_set):
    arr = sorted(num_set)
    return ", ".join(f"{v}({to_binary5(v)})" for v in arr)

##############################
# Core checks
##############################

def check_partition(part):
    lo_hex = part["half_lo_mask_hex"]
    hi_hex = part["half_hi_mask_hex"]
    lo = int(lo_hex, 16) & 0xFFFFFFFF
    hi = int(hi_hex, 16) & 0xFFFFFFFF

    Slo = mask_to_set(lo)
    Shi = mask_to_set(hi)
    U = set(range(NUM_VERTICES))

    # basic set sanity
    if len(Slo) != 16 or len(Shi) != 16:
        raise AssertionError(f"half sizes not 16/16: {len(Slo)}/{len(Shi)}")
    if (Slo & Shi):
        raise AssertionError("halves are not disjoint")
    if (Slo | Shi) != U:
        raise AssertionError("halves do not cover all vertices")

    # complement sanity (bitwise, on 32 vertices)
    if (lo ^ hi) != 0xFFFFFFFF:
        raise AssertionError(f"half_hi is not ~half_lo: lo={lo_hex}, hi={hi_hex}")

    # optional: cross-check vertex lists
    if "half_lo_vertices" in part:
        if sorted(part["half_lo_vertices"]) != sorted(Slo):
            raise AssertionError("half_lo_vertices disagrees with half_lo_mask_hex")
    if "half_hi_vertices" in part:
        if sorted(part["half_hi_vertices"]) != sorted(Shi):
            raise AssertionError("half_hi_vertices disagrees with half_hi_mask_hex")

    return lo, hi, Slo, Shi

def check_witness(entry):
    part = entry["partition"]
    lo, hi, Slo, Shi = check_partition(part)

    wit = entry.get("witness", None)
    if wit is None:
        raise AssertionError("expected witness, but none found")

    if wit.get("maps", "").strip() != "half_lo -> half_hi":
        raise AssertionError(f'expected witness["maps"]="half_lo -> half_hi", got {wit.get("maps")}')

    flip_mask = int(wit["flip_mask"])
    perm = list(wit["perm"])

    if sorted(perm) != [0,1,2,3,4]:
        raise AssertionError(f"perm not a permutation of 0..4: {perm}")

    image = apply_witness_to_set(Slo, flip_mask, perm)
    if image != Shi:
        missing = Shi - image
        extra = image - Shi
        raise AssertionError(
            "witness does not map half_lo to half_hi\n"
            f"  missing (in half_hi, not in image): {sorted(missing)}\n"
            f"  extra   (in image, not in half_hi): {sorted(extra)}"
        )

    # optional: verify vertex_map
    if "vertex_map" in wit:
        stored = list(wit["vertex_map"])
        calc = compute_vertex_map(flip_mask, perm)
        if stored != calc:
            for v in range(NUM_VERTICES):
                if stored[v] != calc[v]:
                    raise AssertionError(f"vertex_map mismatch at v={v}: stored={stored[v]}, calc={calc[v]}")

    # optional: verify half_lo_map_pairs
    if "half_lo_map_pairs" in wit:
        for (u, fu) in wit["half_lo_map_pairs"]:
            if u not in Slo:
                raise AssertionError(f"half_lo_map_pairs includes u={u} not in half_lo")
            if fu not in Shi:
                raise AssertionError(f"half_lo_map_pairs maps to fu={fu} not in half_hi")
            calc_fu = apply_witness_to_vertex(u, flip_mask, perm)
            if calc_fu != fu:
                raise AssertionError(f"half_lo_map_pairs wrong for u={u}: stored={fu}, calc={calc_fu}")

##############################
# Main
##############################

def main():
    if len(sys.argv) < 2:
        print("Usage: python q5_pairs_only_partitions_split_verify.py q5_pairs_only_partitions_split.json")
        sys.exit(2)

    filename = sys.argv[1]
    try:
        with open(filename, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading JSON from {filename}: {e}")
        sys.exit(1)

    translation = data.get("translation", [])
    flip_perm   = data.get("flip_perm", [])
    non_sym     = data.get("non_symmetric", [])
    summary     = data.get("summary", {})

    # summary sanity (if present)
    def chk(name, arr):
        if name in summary and int(summary[name]) != len(arr):
            raise AssertionError(f"summary[{name}]={summary[name]} but len({name})={len(arr)}")

    chk("translation", translation)
    chk("flip_perm", flip_perm)
    chk("non_symmetric", non_sym)

    # check witnesses in translation and flip_perm
    for i, e in enumerate(translation):
        try:
            check_witness(e)
        except AssertionError as ex:
            raise AssertionError(f"[translation entry {i}] {ex}")

    for i, e in enumerate(flip_perm):
        try:
            check_witness(e)
        except AssertionError as ex:
            raise AssertionError(f"[flip_perm entry {i}] {ex}")

    # non_symmetric: must still be valid partitions; should NOT have witness
    for i, e in enumerate(non_sym):
        check_partition(e["partition"])
        if "witness" in e and e["witness"] is not None:
            raise AssertionError(f"[non_symmetric entry {i}] unexpectedly has a witness")

    # verify the reported non_symmetric_graphic if present
    if "non_symmetric_graphic" in summary:
        claimed = int(summary["non_symmetric_graphic"])
        actual = sum(1 for e in non_sym if bool(e.get("graphic", False)))
        if claimed != actual:
            raise AssertionError(f"non_symmetric_graphic mismatch: claimed={claimed}, actual={actual}")

    print("OK: all partitions and witnesses verified successfully.")
    print(f"Counts: translation={len(translation)}, flip_perm={len(flip_perm)}, non_symmetric={len(non_sym)}")

##############################
# Inspector
##############################

def inspect_entry(category, index, filename="q5_pairs_only_partitions_split.json"):
    """
    category: "translation" | "flip_perm" | "non_symmetric"
    index: integer index into that array
    """
    with open(filename, "r") as f:
        data = json.load(f)

    arr = data.get(category, None)
    if arr is None:
        print(f"No category {category!r} in JSON.")
        return
    if index < 0 or index >= len(arr):
        print(f"Index out of range: 0..{len(arr)-1}")
        return

    e = arr[index]
    part = e["partition"]
    lo, hi, Slo, Shi = check_partition(part)

    print(f"\n=== {category}[{index}] ===")
    print("half_lo_mask_hex:", part["half_lo_mask_hex"])
    print("half_hi_mask_hex:", part["half_hi_mask_hex"])
    print("\nhalf_lo vertices:")
    print(" ", format_set_bitpatterns(Slo))
    print("\nhalf_hi vertices:")
    print(" ", format_set_bitpatterns(Shi))

    wit = e.get("witness", None)
    if wit is None:
        print("\n(no witness in this entry)")
        return

    flip_mask = int(wit["flip_mask"])
    perm = list(wit["perm"])
    print("\nwitness:")
    print("  maps:", wit.get("maps"))
    print("  flip_mask:", flip_mask, f"({flip_mask:05b})")
    print("  perm:", perm)

    X = {((v ^ flip_mask) & 0x1F) for v in Slo}
    T = apply_witness_to_set(Slo, flip_mask, perm)

    print("\nStep 1: flip(half_lo):")
    print(" ", format_set_bitpatterns(X))
    print("\nStep 2: perm(flip(half_lo)):")
    print(" ", format_set_bitpatterns(T))

    print("\nResult:")
    print("  perm(flip(half_lo)) == half_hi ?", T == Shi)

if __name__ == "__main__":
    main()
