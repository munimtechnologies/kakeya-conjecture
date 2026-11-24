# kcut: a solver-free exhaustive search for 25- and 26-point Kakeya sets in F_3^4

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Code Apache-2.0, text CC BY 4.0.
Python standard library only (`kcut.py`, tests in `test_kcut.py`).

## Why

k_4, the smallest size of a Kakeya set in F_3^4, is 25, 26 or 27 (`moments/SD.md` proves 25; a
27-point set has been known since 2009, and Lev conjectured 27 is optimal). The exact searches
`kprof25`, `kprof26` and `kdrat26` use CP-SAT and SAT solvers. This folder is a third,
independent route with **no solver at all**: plain enumeration plus an exact set-cover search,
small enough to read and re-implement.

## The reduction

Cut F_3^4 = F_3^3 x F_3 by the class {s = 0, 1, 2}, with sections A, B, C of F_3^3. K is Kakeya iff

1. for every d in F_3^3 some a in A has a + d in B and a + 2d in C (the 27 crossing directions), and
2. each of the 13 directions of F_3^3 has a line inside A, B or C.

Affine maps of F_3^4 let us put the smallest section at s = 0 as a fixed AGL(3,3)-orbit
representative, the middle one at s = 1 containing the origin (a shear (x, s) -> (x - s v, s)
translates B and fixes A), and the largest at s = 2. For each such A and every B, the program
searches for C exactly. C must meet each set W[d] = {-(a + b) : b - a = d}; the search branches on
the smallest unmet W[d], excludes earlier choices in later branches, prunes with a disjoint-packing
bound, and then completes the parallel directions line by line.

Which cuts must be searched comes from the moment identities (`moments/BAD25.md`,
`moments/min_section.py`):

- **N = 25:** some class has cut (6,6,13), (5,7,13) or (6,7,12).
- **N = 26:** some class has cut (6,6,14), (5,7,14), (6,7,13) or (7,7,12). Otherwise every class
  has an 8-point section and some hyperplane has at most 7 points, which leaves a class with cut
  (5,8,13), (6,8,12) or (7,8,11).

If no cut in a list admits a set, that size is impossible: k_4 >= 26 from the first list, and
k_4 = 27 from the second. Both conclusions rest on the moment lemmas, which have short proofs and
committed tests.

## Controls (in the job and the tests)

- In F_3^3, the same engine searches every cut of 12 and of 13 points. It must find nothing at 12
  and a set at 13, because k_3 = 13.
- In F_3^2 sections, `solve` and `run` are compared with brute force over all C (and all B) on
  random instances. The B count must match exactly.
- The known 27-point set: its cut (6,8,13) must be found by the full enumeration, and no 26-point
  subset of it is Kakeya.
- The cut (6,6,12) must give nothing, since k_4 >= 25 is already proved.
- Every reported set is re-checked against the plain definition of a Kakeya set.

## Status

Job `kcut-v1` is queued. Its results will be added here; no claim is made until then.
