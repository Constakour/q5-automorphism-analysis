#!/usr/bin/env python3
"""
check_hamming.py

Utilities for Hamming distance in Q5 (vertices 0..31).
Also includes an optional "phi sanity check":
Given H (as a 32-bit mask) and a mapping phi: H -> Hc,
find a pair u,v in H where distances disagree:
    d(u,v) != d(phi(u), phi(v))
This certifies phi is not a cube automorphism.

Examples:
  # basic distance
  python3 check_hamming.py --dist 0 14

  # weight / popcount
  python3 check_hamming.py --wt 0x1e

  # print all distances from a vertex
  python3 check_hamming.py --from 7

  # print full 32x32 distance matrix
  python3 check_hamming.py --matrix

  # check phi-distance preservation on your counterexample witness
  python3 check_hamming.py --mask 0x2695b9f --phi "0:31,1:30,2:27,3:26,4:15,7:10,8:29,9:28,11:24,12:13,14:5,16:23,19:18,21:6,22:17,25:20"
"""

import argparse
from typing import Dict, List, Tuple


def popcount(x: int) -> int:
    return int(x).bit_count()


def hamming_dist(u: int, v: int) -> int:
    return popcount(u ^ v)


def decode_mask(mask: int, n: int = 32) -> List[int]:
    return [v for v in range(n) if (mask >> v) & 1]


def parse_phi(phi_str: str) -> Dict[int, int]:
    """
    Parse phi mapping from a string like:
      "0:31,1:30,2:27,..."
    returns dict {0:31, 1:30, ...}
    """
    phi: Dict[int, int] = {}
    if not phi_str.strip():
        return phi
    parts = [p.strip() for p in phi_str.split(",") if p.strip()]
    for p in parts:
        if ":" not in p:
            raise ValueError(f"Bad phi entry '{p}'. Expected 'u:v'.")
        a, b = p.split(":")
        u = int(a.strip(), 0)
        v = int(b.strip(), 0)
        phi[u] = v
    return phi


def phi_distance_counterexample(H: List[int], phi: Dict[int, int]) -> Tuple[int, int, int, int]:
    """
    Return one witness (u,v, d(u,v), d(phi(u),phi(v))) where distances differ.
    Raises ValueError if phi is incomplete on H.
    Returns (-1,-1,0,0) if no counterexample is found (i.e., phi preserves all distances on H).
    """
    for u in H:
        if u not in phi:
            raise ValueError(f"phi missing value for u={u}")
    for i in range(len(H)):
        for j in range(i + 1, len(H)):
            u, v = H[i], H[j]
            duv = hamming_dist(u, v)
            dimg = hamming_dist(phi[u], phi[v])
            if duv != dimg:
                return (u, v, duv, dimg)
    return (-1, -1, 0, 0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wt", nargs=1, help="popcount / Hamming weight of x (decimal or 0x..)")
    ap.add_argument("--dist", nargs=2, help="Hamming distance d(u,v) for u v (0..31)")
    ap.add_argument("--from", dest="from_v", nargs=1, help="Print distances from a fixed vertex to all 0..31")
    ap.add_argument("--matrix", action="store_true", help="Print full 32x32 distance matrix")

    # Optional phi-check mode:
    ap.add_argument("--mask", help="Half mask (32-bit) defining H (popcount should be 16)")
    ap.add_argument("--phi", help="Mapping string 'u:v,...' for u in H (typically phi: H -> Hc)")
    args = ap.parse_args()

    if args.wt:
        x = int(args.wt[0], 0)
        print(f"wt({args.wt[0]}) = {popcount(x)}")
        return

    if args.dist:
        u = int(args.dist[0], 0)
        v = int(args.dist[1], 0)
        print(f"d({u},{v}) = wt({u} XOR {v}) = wt({u ^ v}) = {hamming_dist(u,v)}")
        return

    if args.from_v:
        u = int(args.from_v[0], 0)
        print(f"Distances from {u}:")
        for v in range(32):
            print(f"  d({u},{v}) = {hamming_dist(u,v)}")
        return

    if args.matrix:
        print("32x32 Hamming distance matrix for Q5 vertices 0..31:")
        header = "    " + " ".join(f"{v:2d}" for v in range(32))
        print(header)
        for u in range(32):
            row = " ".join(f"{hamming_dist(u,v):2d}" for v in range(32))
            print(f"{u:2d}: {row}")
        return

    # phi-check mode
    if args.mask and args.phi:
        mask = int(args.mask, 0)
        H = decode_mask(mask)
        phi = parse_phi(args.phi)

        print("=== phi distance-preservation check (Q5) ===")
        print(f"mask = 0x{mask:08x}  popcount={len(H)}")
        print(f"H   = {H}")

        # sanity: phi(H) should be inside complement (not required for distance test, but helpful)
        comp = ((1 << 32) - 1) ^ mask
        Hc = decode_mask(comp)
        images = [phi[u] for u in H if u in phi]
        bad = [x for x in images if x not in set(Hc)]
        if bad:
            print("Warning: some phi(u) are not in H^c:", bad)

        u, v, d1, d2 = phi_distance_counterexample(H, phi)
        if u == -1:
            print("No pair (u,v) in H found with d(u,v) != d(phi(u),phi(v)).")
            print("So phi preserves all Hamming distances on H (this does NOT prove phi is a cube automorphism,")
            print("but it means this particular distance-based disproof did not fire).")
        else:
            print("FOUND distance mismatch, so phi cannot be a cube automorphism:")
            print(f"  u={u}, v={v}")
            print(f"  phi(u)={phi[u]}, phi(v)={phi[v]}")
            print(f"  d(u,v) = {d1}")
            print(f"  d(phi(u),phi(v)) = {d2}")
        return

    ap.print_help()


if __name__ == "__main__":
    main()
