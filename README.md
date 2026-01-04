# q5-automorphism-analysis

C++ and Python tools for computations on half-sets (size 16) in the 5-dimensional hypercube \(Q_5\), including:
- exhaustive orbit-level enumeration of halves/partitions,
- detection/classification of complement symmetries under cube automorphisms,
- verification scripts and sanity checks,
- (optional) JSON datasets produced by the programs.

Repository structure:

- `src/` — C++ sources (main programs)
- `scripts/` — Python helpers / verifiers
- `data/` — optional JSON outputs (datasets)

## Contents

### C++ programs (`src/`)

- `cube_halves_exhaustive_orbits.cpp` **(paper program)**  
  Isomorph-free exhaustive orbit generation over halves of \(V(Q_5)\), working at the unordered partition level \(\{H,H^c\}\), and classifying orbits (e.g., graphic vs swappable).

- `cube_halves_exhaustive_orbits_classified.cpp`  
  Variant that additionally records swappability-class information (translation / mixed, etc.), as used in the computational census.

- `q5_pairs_only_classify.cpp`  
  Exhaustive classification within the pairs-only family (halves that are unions of antipodal pairs).

- `q5_pairs_only_full_classification_split.cpp`  
  Produces the full split classification dataset for pairs-only partitions.

### Python scripts (`scripts/`)

- `verify_witness.py`  
  Verifies a claimed witness half-set / partition and (optionally) prints an induced isomorphism map and/or distance obstructions.

- `q5_pairs_only_partitions_split_verify.py`  
  Verifies the consistency of the pairs-only split dataset.

- `q5_is_equitable_check.py`  
  Tests equitability of a 2-partition \(\{H,H^c\}\) and reports quotient matrices.

- `decode_mask.py`  
  Decodes hexadecimal bitmasks into vertex sets.

- `check_hamming.py`  
  Small helper for Hamming-distance sanity checks.

- `witness.py`  
  Helper for working with / printing witness instances.

### Data (`data/`)

This folder may contain JSON outputs produced by the programs, e.g.
- `q5_pairs_only_classification.json`
- `q5_pairs_only_partitions_split.json`

If you prefer a lightweight repo, you can delete large JSON files and keep only a small example output.

## Build and run (Ubuntu / WSL)

### Requirements
- `g++` with C++17 support (e.g. `g++ >= 9`)
- Python `>= 3.6`

Install compiler (Ubuntu):
```bash
sudo apt update
sudo apt install -y build-essential

