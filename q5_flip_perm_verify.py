#!/usr/bin/env python
"""
q5_flip_perm_verify.py

This script verifies that half-sets from the file q5_required_permutation.json
conform to the transformation S -> flip_mask -> perm => complement(S)
in a 5-dimensional hypercube Q_5 (32 vertices labeled 0..31).

Each JSON entry is expected to have:
  {
    "mask_hex": "0x....",     # a 32-bit mask for 16 vertices
    "vertices": [...16 ints...],
    "flip_mask": someInt(0..31),
    "perm": [p0,p1,p2,p3,p4]  # a permutation of [0,1,2,3,4]
  }

The code checks:
  T = perm( flip(S) ) == complement(S) ?

You can run this script via:
  python q5_flip_perm_verify.py

And optionally use the inspect_entry(index) function at the end for a detailed
step-by-step on any single entry in the JSON.
"""

import json
import sys

##############################
# Configuration
##############################
FILENAME = "q5_required_permutation.json"  # JSON file to load
N = 5
NUM_VERTICES = 1 << N  # 32 for Q_5

##############################
# Helper Functions
##############################

def mask_to_set(mask_int):
    """
    Convert an integer mask to a set of vertices 0..31
    by checking bits [0..31].
    """
    S = set()
    for v in range(NUM_VERTICES):  # 0..31
        if (mask_int & (1 << v)) != 0:
            S.add(v)
    return S

def to_binary5(x):
    """Return x as a 5-bit binary string, e.g. x=13 => '01101'."""
    return f"{x:05b}"

def apply_flip(S, flipmask):
    """
    Flip each vertex v in S => (v XOR flipmask).
    Force it to 5 bits by & 0x1F as a safety net.
    """
    return { ((v ^ flipmask) & 0x1F) for v in S }

def permute_bits(val, perm):
    """
    perm[i] means oldBit( val >> perm[i] ) => newBit i.
    We force result to <32 by &0x1F.
    """
    newv = 0
    for i in range(N):
        old_bit = (val >> perm[i]) & 1
        newv |= (old_bit << i)
    newv &= 0x1F
    return newv

def apply_perm(S, perm):
    """
    Apply permute_bits(...) to each vertex in S.
    """
    out = set()
    for v in S:
        newv = permute_bits(v, perm)
        out.add(newv)
    return out

##############################
# Main Verification
##############################

def main():
    # Load JSON
    try:
        with open(FILENAME, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading JSON from {FILENAME}: {e}")
        sys.exit(1)

    if isinstance(data, dict) and "sets" in data:
        # If the JSON has a top-level dict with a 'sets' key
        sets_data = data["sets"]
    elif isinstance(data, list):
        # If the JSON is directly a list of entries
        sets_data = data
    else:
        print("Unexpected JSON structure. Not a list or doesn't contain 'sets' key.")
        sys.exit(1)

    total = len(sets_data)
    if total == 0:
        print("No sets found in JSON.")
        sys.exit(0)

    print(f"Loaded {total} entries from {FILENAME}.")

    mismatches = []
    for idx, entry in enumerate(sets_data):
        # parse mask
        mask_hex_str = entry["mask_hex"]
        mask_int = int(mask_hex_str, 16) & 0xFFFFFFFF  # parse as hex => integer
        S = mask_to_set(mask_int)
        comp = set(range(NUM_VERTICES)) - S

        flip_mask = entry["flip_mask"]
        perm_list = entry["perm"]

        # flip->perm
        X = apply_flip(S, flip_mask)
        T = apply_perm(X, perm_list)

        if T != comp:
            mismatches.append(idx)

    passed = total - len(mismatches)
    print(f"Verification complete: {passed} out of {total} entries matched the complement transformation.")
    if mismatches:
        print(f"Mismatch found in these entries: {mismatches}")
    else:
        print("All entries successfully map S to its complement.")

    # optional: interactively call inspect_entry(...) here if desired.

##############################
# Detailed Inspector
##############################
def inspect_entry(index):
    """
    Inspect a single entry in detail. 
    Print step-by-step transformations:
      - Original set S
      - Flip => X
      - Perm => T
      - Compare T with complement(S)
    """
    try:
        with open(FILENAME, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading JSON from {FILENAME}: {e}")
        return

    # same logic as main for the data structure
    if isinstance(data, dict) and "sets" in data:
        sets_data = data["sets"]
    elif isinstance(data, list):
        sets_data = data
    else:
        print("Unexpected JSON structure. Not a list or missing 'sets' key.")
        return

    if index<0 or index>= len(sets_data):
        print(f"Index {index} out of range [0..{len(sets_data)-1}].")
        return

    entry = sets_data[index]
    mask_hex_str = entry["mask_hex"]
    mask_int = int(mask_hex_str, 16) & 0xFFFFFFFF
    S = mask_to_set(mask_int)
    comp = set(range(NUM_VERTICES)) - S

    flip_mask = entry["flip_mask"]
    perm_list = entry["perm"]

    # flip->perm
    X = apply_flip(S, flip_mask)
    T = apply_perm(X, perm_list)

    def format_set_bitpatterns(num_set):
        arr = sorted(num_set)
        return ", ".join(f"{v}({to_binary5(v)})" for v in arr)

    print(f"\n=== Detailed Entry {index} ===")
    print(f"mask_hex: {mask_hex_str}")
    print(f"flip_mask: {flip_mask} ({flip_mask:05b} in binary)")
    print(f"perm: {perm_list}\n")

    print(f"Original set S (size={len(S)}):")
    print(f"  {format_set_bitpatterns(S)}")

    print(f"\nComplement of S (size={len(comp)}):")
    print(f"  {format_set_bitpatterns(comp)}")

    print(f"\nStep1: Flip => X:")
    print(f"  {format_set_bitpatterns(X)}")

    print("\nStep2: Perm => T:")
    print(f"  {format_set_bitpatterns(T)}")

    print()
    if T == comp:
        print("SUCCESS: T matches the complement of S.")
    else:
        print("MISMATCH: T != complement(S).")
        missing = comp - T
        extra   = T - comp
        if missing:
            print(" In complement but not in T:", format_set_bitpatterns(missing))
        if extra:
            print(" In T but not in complement:", format_set_bitpatterns(extra))

##############################
# If run as main, do main()
##############################
if __name__ == "__main__":
    main()
