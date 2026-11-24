# A plane-class relaxation for Kakeya sets in F_3^4

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

Code: `plane_lp.py` (standard library only). Tests: `test_plane_lp.py` (sanity) and the readout
tests listed at the end. Builds on `SD.md` (k_4 >= 25 by hyperplane moments and a sum-difference
lemma).

## Why

`SD.md` uses only the 40 hyperplane classes of F_3^4. At N = 25 that relaxation is feasible: an
averaged mix of section profiles such as (8,8,9) and (6,6,13) survives every condition, with between
40 and 47 lines. This note adds the 130 classes of parallel planes, which see much finer structure.

## The arrays

Let K be a Kakeya set in F_3^4 with N points and L >= 40 lines. For a 2-dimensional direction space
W, the 9 planes parallel to W are the cosets of W; indexing them by (u.x, v.x) for two normals u, v
spanning the annihilator of W turns the section sizes into a 3x3 array M with entries 0..9 and total
N. The 12 lines of this AG(2,3) are exactly the 12 hyperplanes of the 4 hyperplane classes that
contain W, so each of the 4 line directions of M gives the profile of one hyperplane class.

Write c_k(M) for the sum over the 9 cells of C(m, k), and h_k(M) for the sum over the 12 lines of
C(line sum, k).

## Exact identities (sums over the 130 classes; C_k = C(N, k))

* (I2) sum c_2 = 13 C_2: a pair of points lies in one cell exactly when its direction is in W, and a
  direction lies in 13 of the spaces W.
* (I3) sum c_3 = C_3 + 12 L: three non-collinear points span one plane direction W; a line of K lies
  in a cell for 13 spaces W.
* (H3) sum h_3 = 13 (4 C_3 + 9 L): this is the hyperplane triple identity of `SD.md`, counted once
  for each of the 13 spaces W inside each hyperplane class.
* (I4) sum (h_4 - 39 c_4) = 13 C_4. Four points of K are never collinear. If they span a 3-flat they
  lie in exactly one hyperplane; if they are coplanar they lie in 4 hyperplanes and in exactly one
  cell of one class W. So the hyperplane 4-moment is C_4 + 3 Q_2 with Q_2 = sum c_4, and each
  hyperplane class is counted 13 times. This identity does not involve L.
* (Ik) sum (h_k - 39 c_k) <= 13 C_k for k = 5..13 (k points may also span all of F_3^4).
* (LW) sum linecap >= 13 L: a line of K whose direction lies in W lies inside one cell, and a cell of
  m points contains at most 0,0,0,1,1,2,3,5,8,12 lines (m = 0..9).
* (LX) sum crosscap >= 117 L: a line whose direction is not in W meets three collinear cells and is
  fixed by its points in any two of them, so a quotient line with cells a, b, c carries at most
  min(ab, ac, bc) such lines; each line of K is counted for 117 spaces W.

## Conditions on a single array

1. Every line direction of M gives an allowed hyperplane profile under (T), (SD) and (HZ) of
   `SD.md`.
2. Each of the 4 line directions has a line whose 3 cells are all non-empty: the image of K in
   F_3^4 / W is a Kakeya set of AG(2,3).
3. The 4 directions of W are carried by lines inside cells, and a cell of m points carries at most
   0,0,0,1,1,2,3,4,4,4 directions (m = 0..9; 7 = k_2), so these numbers sum to at least 4.
4. For each line direction (a hyperplane class U containing W), the 9 directions of U outside W are
   carried by lines that cross three collinear cells inside one hyperplane. Three parallel plane
   sections of sizes a, b, c inside AG(3,3) carry at most min(9, ab, ac, bc) of them, and all 9 only
   if the sum-difference condition holds with T = 9. Summed over the 3 parallel lines this must
   reach 9.

## The relaxation

A Kakeya set of size N gives 130 arrays that satisfy 1-4, together with its L, satisfying all the
identities. `plane_feasible(N, Lmin)` enumerates every valid array up to row permutation, merges
arrays with equal moment vectors, and decides with an exact rational simplex (Bland's rule) whether
nonnegative weights with total 130 and some L >= Lmin satisfy the system. If there are none, there
is no Kakeya set with N points and at least Lmin lines.

## Checks built into the tests

* On the 27-point Kakeya set found by the kakeya-min-v1 integer program (in `plane_lp.py`): the set
  is checked to be Kakeya, all 130 of its arrays satisfy conditions 1-4, every identity holds
  exactly with its true L, and its arrays (rows sorted) are found by the enumeration.
* The relaxation is infeasible at N = 24 (implied by `SD.md`) and feasible at N = 27.
* The exact simplex is checked on small systems.

## Results: read them off the repository

Each readout test is committed only if it passes, so the files present in this folder are the
results.

| file present | meaning |
|---|---|
| `test_plane_k4_at_least_26.py` | no Kakeya set of 25 points: k_4 >= 26 |
| `test_plane_k4_at_least_27.py` | none of 25 or 26 points either: k_4 = 27 |
| `test_plane_n25_lines_below_X.py` | a 25-point Kakeya set would have fewer than X lines |
| `test_k7_at_least_X.py` | `SD.md`'s hyperplane relaxation gives k_7 >= X (Bukh-Chao: 103) |

If none of the plane readouts is present, the plane relaxation does not improve on `SD.md`, and we
will say so.

## Limits

These are relaxation bounds. The conditions are necessary, not sufficient. The code was written and
checked by AI agents and has not been reviewed by a human mathematician.
