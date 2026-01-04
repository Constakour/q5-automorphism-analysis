#!/usr/bin/env python3
"""
decode_mask.py

Decode a 32-bit (Q5) vertex-mask into the corresponding vertex set.

Convention (matches your C++):
- Vertices are integers v in {0,...,31}
- Bit v of the mask indicates membership of vertex v in H
- Complement uses full_mask ^ mask, with full_mask = 0xFFFFFFFF for Q5

Usage:
  python3 decode_mask.py 0x2695b9f
  python3 decode_mask.py 0x2695b9f --binary
  python3 decode_mask.py 0x2695b9f --msb-vectors
  python3 decode_mask.py 0x2695b9f --complement
"""

import argparse

N = 5
NUM_VERTICES = 1 << N
FULL_MASK_Q5 = 0xFFFFFFFF


def parse_mask(s: str) -> int:
    s = s.strip().lower()
    if s.startswith("0x"):
        return int(s, 16) & 0xFFFFFFFF
    # allow plain decimal too
    return int(s, 10) & 0xFFFFFFFF


def mask_to_vertices(mask: int):
    return [v for v in range(NUM_VERTICES) if (mask >> v) & 1]


def to_binary5(v: int) -> str:
    return f"{v:05b}"


def to_msb_vector(v: int):
    # (b4,b3,b2,b1,b0) from integer v
    b = to_binary5(v)
    return tuple(int(ch) for ch in b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mask", help="Mask as hex (0x...) or decimal.")
    ap.add_argument("--binary", action="store_true", help="Also print each vertex in 5-bit binary.")
    ap.add_argument("--msb-vectors", action="store_true", help="Also print each vertex as (b4,b3,b2,b1,b0).")
    ap.add_argument("--complement", action="store_true", help="Decode complement mask instead of the mask.")
    args = ap.parse_args()

    mask = parse_mask(args.mask)
    if args.complement:
        mask = FULL_MASK_Q5 ^ mask

    verts = mask_to_vertices(mask)
    print(f"mask = 0x{mask:08x}  (popcount={len(verts)})")
    print("H = { " + ", ".join(map(str, verts)) + " }")

    if args.binary or args.msb_vectors:
        print("\nDetailed:")
        for v in verts:
            line = f"  {v:2d}"
            if args.binary:
                line += f"  bin={to_binary5(v)}"
            if args.msb_vectors:
                line += f"  (b4..b0)={to_msb_vector(v)}"
            print(line)


if __name__ == "__main__":
    main()
