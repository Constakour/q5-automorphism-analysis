#!/usr/bin/env python3
"""
q5_is_equitable_check.py

Check equitable 2-partitions in Q5 for a JSON produced by the C++ splitter, e.g.
  q5_pairs_only_partitions_split.json

Expected top-level structure:
{
  "meta": {...},
  "summary": {...},
  "translation": [ { "partition": {...}, "graphic": bool, "witness": {...} }, ... ],
  "flip_perm":   [ ... ],
  "non_symmetric":[ ... ]
}

Each entry has:
  entry["partition"]["half_lo_mask_hex"]
  entry["partition"]["half_hi_mask_hex"]
(or vertices arrays, but we recompute from masks).

We test equitability of the partition {H, H^c} by checking that:
- all v in H have the same #neighbors in H
- all v in H^c have the same #neighbors in H^c

In Q5 (5-regular) that is equivalent to also having constant cross-neighbor counts.

Outputs:
- total partitions tested (across all categories found)
- equitable partitions total
- breakdown by category
- quotient matrix frequencies
- optional: print a few equitable examples
"""

import argparse
import json
from collections import Counter, defaultdict

N = 5
NUM_VERTICES = 1 << N  # 32
FULL_MASK = (1 << NUM_VERTICES) - 1  # 0xFFFFFFFF


def mask_from_hex(s: str) -> int:
    return int(s, 16) & FULL_MASK


def neighbors(v: int):
    # Q5 neighbors: flip one bit
    for b in range(N):
        yield v ^ (1 << b)


def equitable_quotient(half_mask: int):
    """
    Return quotient matrix (a,b,c,d) where:
      a = neighbors in H for v in H
      b = neighbors in H^c for v in H
      c = neighbors in H for v in H^c
      d = neighbors in H^c for v in H^c
    or None if not equitable.
    """
    H = half_mask
    C = FULL_MASK ^ H

    # collect counts for vertices in H and in C
    a_vals = []
    d_vals = []

    for v in range(NUM_VERTICES):
        if (H >> v) & 1:
            cnt_in_H = 0
            for u in neighbors(v):
                if (H >> u) & 1:
                    cnt_in_H += 1
            a_vals.append(cnt_in_H)
        else:
            cnt_in_C = 0
            for u in neighbors(v):
                if (C >> u) & 1:
                    cnt_in_C += 1
            d_vals.append(cnt_in_C)

    # equitability requires all equal within each side
    if len(set(a_vals)) != 1:
        return None
    if len(set(d_vals)) != 1:
        return None

    a = a_vals[0]
    d = d_vals[0]
    b = N - a
    c = N - d
    return (a, b, c, d)


def iter_entries(data: dict):
    """
    Yield (category_name, entry_dict) for each category array present.
    """
    for cat in ("translation", "flip_perm", "non_symmetric"):
        if cat in data and isinstance(data[cat], list):
            for entry in data[cat]:
                yield cat, entry


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("json_file", help="e.g. q5_pairs_only_partitions_split.json")
    ap.add_argument("--examples", type=int, default=0,
                    help="print up to this many equitable examples")
    args = ap.parse_args()

    with open(args.json_file, "r") as f:
        data = json.load(f)

    tested = 0
    equitable_total = 0

    equitable_by_cat = Counter()
    cat_sizes = Counter()
    quotient_freq = Counter()

    examples = []  # (cat, (a,b,c,d), lo_hex, hi_hex)

    for cat, entry in iter_entries(data):
        cat_sizes[cat] += 1
        part = entry.get("partition", {})
        lo_hex = part.get("half_lo_mask_hex")
        hi_hex = part.get("half_hi_mask_hex")
        if lo_hex is None or hi_hex is None:
            continue

        H = mask_from_hex(lo_hex)
        # sanity: hi should be complement
        if mask_from_hex(hi_hex) != (FULL_MASK ^ H):
            # if your file ever changes conventions, this catches it
            raise ValueError(f"Complement mismatch in category={cat}, lo={lo_hex}, hi={hi_hex}")

        tested += 1
        q = equitable_quotient(H)
        if q is not None:
            equitable_total += 1
            equitable_by_cat[cat] += 1
            quotient_freq[q] += 1
            if args.examples and len(examples) < args.examples:
                examples.append((cat, q, lo_hex, hi_hex))

    if args.examples:
        for cat, q, lo_hex, hi_hex in examples:
            a, b, c, d = q
            print("\nEQUITABLE example:")
            print(f"  category: {cat}")
            print(f"  quotient matrix: [[{a},{b}],[{c},{d}]]")
            print(f"  half_lo_mask_hex: {lo_hex}")
            print(f"  half_hi_mask_hex: {hi_hex}")

    print("\n=== Equitable 2-partition check in Q5 ===")
    print(f"file: {args.json_file}")
    print(f"partitions_tested: {tested}")
    print(f"equitable_partitions: {equitable_total}\n")

    print("By category:")
    for cat in ("translation", "flip_perm", "non_symmetric"):
        if cat_sizes[cat]:
            e = equitable_by_cat[cat]
            t = cat_sizes[cat]
            pct = 100.0 * e / t
            print(f"  {cat:12s} equitable={e:5d} / {t:5d}  ({pct:6.3f}%)")
    print()

    if quotient_freq:
        print("Quotient matrices (a b; c d) frequency:")
        for (a, b, c, d), cnt in sorted(quotient_freq.items(), key=lambda kv: (-kv[1], kv[0])):
            print(f"  [[{a},{b}],[{c},{d}]] : {cnt}")
    else:
        print("No equitable partitions found.")

    print("\nNote:")
    print("  This counts *partitions* (unordered {H,H^c}) exactly once, as they appear in your JSON.")


if __name__ == "__main__":
    main()
