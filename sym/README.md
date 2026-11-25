# sym: Kakeya sets over F_3 with symmetry from the field F_{3^n}

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name.
**No human mathematician has checked this.** Code Apache-2.0, text CC BY 4.0.

## Question

The smallest Kakeya sets over F_3 are known only in dimensions up to 4 (3, 7, 13, 27); in
dimension 5 the answer lies between 45 and 53, in dimension 6 between 78 and 105, in
dimension 7 between 137 and 213 (lower bounds: our unchecked moment method, `moments/SD.md`;
upper bounds: Chung, Du and Wesley for 53, our searches for 105 and 213, `k6/RECORDS.md`).
All of the small sets found so far came from unstructured local search. Do good Kakeya sets
have symmetry?

Identify F_3^n with the field F_{3^n}. Multiplication by a nonzero element and the Frobenius
map x -> x^3 are F_3-linear, so they map lines to lines. For each subgroup H_r of order r of
the multiplicative group (r divides 3^n - 1) and each subgroup F_f of order f of the Galois
group (f divides n), `ksym.py` finds the **least size of a Kakeya set invariant under
G(r, f) = H_r F_f**. A G-invariant set is a union of G-orbits, and it is Kakeya as soon as it
covers one direction from each G-orbit of directions, so the problem shrinks to a small 0/1
program over orbits. Example: in F_3^6, H_13 has orbits of size 13 and only 28 orbits of
directions to cover.

## Method

- `gf.py`: the field (first lexicographic primitive polynomial), group generators as point
  permutations, point orbits, direction orbits, and for each representative direction the
  minimal sets of orbits that contain one of its lines.
- `ksym.py`: CP-SAT (OR-Tools) minimises the size over orbit unions; every set found is
  re-checked by `check.py`, a stand-alone checker that knows only the definition. For orbit
  problems with at most 70 orbits, `bb.py` (an exact branch and bound that uses no solver)
  re-proves each optimum. Phase 2 starts the WalkSAT search of `k6/fastwalk.py` from the
  best symmetric sets.
- `oracle.py`: a standard-library pass (branch and bound only) over the groups with at most
  60 orbits.
- `test_ksym.py`: field and group checks, and the orbit model against brute force for n = 2
  and 3 (k_2 = 7, k_3 = 13).

## What would count

A symmetric set below a record (52 or fewer in F_3^5, 104 or fewer in F_3^6, 212 or fewer in
F_3^7) would be a new upper bound with an explicit structure that may generalise. Exact
minima for each symmetry class are also new facts: they say which symmetries good Kakeya sets
can and cannot have. Results will be added to this folder when the runs report.
