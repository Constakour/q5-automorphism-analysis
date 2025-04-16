# q5-automorphism-analysis

C++ and Python tools for detecting **pairs-only half-set complement symmetries** in the 5-dimensional hypercube \( Q_5 \) via automorphisms.

This project contains two key programs:

- **`Q5_new_auto.cpp`**: A C++ program that enumerates half-sets and detects those requiring full automorphisms (flip + coordinate permutation) to map to their complements.
- **`q5_flip_perm_verify.py`**: A Python script that verifies each transformation produced by the C++ program.

---

## Installation

### Requirements

- **C++17 or later**
- **Python 3.6+**

For Python:
```bash
pip install --upgrade pip
# No external dependencies required
