# Makefile for q5-automorphism-analysis
# Builds C++ programs into ./build/

CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra
OMPFLAGS := -fopenmp

BUILD := build
SRC   := src

EXE_ORBITS     := $(BUILD)/cube_halves_exhaustive_orbits
EXE_CLASSIFIED := $(BUILD)/cube_halves_exhaustive_orbits_classified
EXE_PAIRS_ONLY := $(BUILD)/q5_pairs_only_classify
EXE_SPLIT      := $(BUILD)/q5_pairs_only_full_classification_split

all: $(EXE_ORBITS) $(EXE_CLASSIFIED) $(EXE_PAIRS_ONLY) $(EXE_SPLIT)

$(BUILD):
	mkdir -p $(BUILD)

# Orbit program (no OpenMP required)
$(EXE_ORBITS): $(SRC)/cube_halves_exhaustive_orbits.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Classified orbit program (uses OpenMP)
$(EXE_CLASSIFIED): $(SRC)/cube_halves_exhaustive_orbits_classified.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -o $@ $<

$(EXE_PAIRS_ONLY): $(SRC)/q5_pairs_only_classify.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(EXE_SPLIT): $(SRC)/q5_pairs_only_full_classification_split.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Convenience targets
run_orbits: $(EXE_ORBITS)
	./$(EXE_ORBITS) --exhaustive-q5

# Example: make run_classified T=14
run_classified: $(EXE_CLASSIFIED)
	./$(EXE_CLASSIFIED) --exhaustive-q5 --threads $(T)

clean:
	rm -rf $(BUILD)

.PHONY: all clean run_orbits run_classified
