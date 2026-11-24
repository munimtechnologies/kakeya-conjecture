# Which cuts a 25-point Kakeya set in F_3^4 must have

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

Read `SD.md` first. There, conditions (T) and (SD) and the two moment identities prove k_4 >= 25,
and the averaged (linear) relaxation is feasible at N = 25, so moments alone stop there. This note
turns that relaxation into a structural statement about any 25-point Kakeya set, which then
narrows the exact search for k_4 >= 26 to three kinds of cut. Code: `bad25.py`. Every claim
below is checked by `test_bad25.py` in exact integer arithmetic and by brute force on random sets.

## The identity

Let K be any set of 25 points of F_3^4 with L lines. For each of the 40 parallel classes, write its
cut (n1, n2, n3) and s2 = sum C(n_i, 2), s3 = sum C(n_i, 3). The identities of `SD.md` (a1 = 13,
a2 = 4) say that the 40 values of s2 add up to 13 C(25, 2) = 3900 and the 40 values of s3 to
4 C(25, 3) + 9 L = 9200 + 9 L. Put

    delta(cut) = 7 s3 - 51 s2 + 3302.

Then, for every 25-point set,

    sum over the 40 classes of delta = 7 (9200 + 9 L) - 51 * 3900 + 40 * 3302 = 63 L - 2420.   (*)

A Kakeya set contains a line in each of the 40 directions, so L >= 40 and the sum is **at least 100**.

## The 14 cuts and their delta

(T) and (SD) leave exactly 14 cuts for a 25-point Kakeya set (`allowed_profiles(25)`):

| cut | s2 | s3 | delta |
|---|---|---|---|
| (6,6,13) | 108 | 326 | **+76** |
| (5,7,13) | 109 | 331 | **+60** |
| (6,7,12) | 102 | 275 | **+25** |
| (7,7,11) | 97 | 235 | 0 |
| (5,8,12) | 104 | 286 | 0 |
| (6,8,11) | 98 | 241 | -9 |
| (7,8,10) | 94 | 211 | -15 |
| (8,8,9) | 92 | 196 | -18 |
| (7,9,9) | 93 | 203 | -20 |
| (6,9,10) | 96 | 224 | -26 |
| (5,9,11) | 101 | 259 | -36 |
| (5,10,10) | 100 | 250 | -48 |
| (4,9,12) | 108 | 308 | -50 |
| (4,10,11) | 106 | 289 | -81 |

The weights 7 and 51 are the slope of the upper hull of the points (s2, s3) through (7,7,11) and
(5,8,12) once the three top cuts are removed.

## Theorem

Let K be a Kakeya set in F_3^4 with 25 points. Then:

1. **At least two** of its 40 classes have cut (6,6,13), (5,7,13) or (6,7,12): the positive deltas
   must add up to at least 100, and one class gives at most 76. Two classes suffice only for the
   pairs {(6,6,13), (6,6,13)}, {(6,6,13), (5,7,13)}, {(5,7,13), (5,7,13)} and {(6,6,13), (6,7,12)}.
2. If K has no hyperplane section of 13 points, **at least four** classes have cut (6,7,12)
   (3 x 25 < 100 <= 4 x 25). Sections never exceed 13 points, and 13-point sections occur only in
   the cuts (6,6,13) and (5,7,13).
3. K has L = (2420 + sum of delta)/63 lines, so L is at most (2420 + the sum of the positive deltas)/63.

*Proof.* (*) with L >= 40, and the table. QED.

Without the three top cuts the relaxation of `SD.md` is infeasible at N = 25; with them it is
feasible. So the theorem is exactly what the moment method still says about 25 points.

## Use: a narrower exact search

An affine map of F_3^4 can move any class to the "last coordinate" class, with its smallest
section at s = 0 and its largest at s = 2 (translate s, then possibly s -> -s). An affine map of
the first three coordinates then moves the s = 0 section to a canonical representative of its
orbit under AGL(3,3). So if a 25-point Kakeya set exists, one exists in one of these positions:

- **A.** the last class has cut (6,6,13), slice 0 a canonical 6-point set;
- **B.** no class has cut (6,6,13), the last class has cut (5,7,13), slice 0 a canonical 5-point set;
- **C.** no section has 13 points, the last class has cut (6,7,12), slice 0 a canonical 6-point set
  (and then at least four classes have cut (6,7,12)).

In each position, all 40 classes must have one of the allowed cuts (minus the ones excluded in B
and C), and the identities hold with L counted exactly. The job `kprof25` (folder `kprof25/`)
decides these cases with CP-SAT and, as an independent check, the HiGHS MILP solver. If every case
is infeasible, k_4 >= 26. It replaces the case split of `kmin43` (smallest section n0 = 1..6, all
other cuts free) by fixed cuts and a table of allowed cuts in every class.

## The same identity at N = 26

For 26 points the identities give sum s2 = 4225 and sum s3 = 10400 + 9 L, so with
delta26(cut) = s3 - 8 s2 + 576,

    sum over the 40 classes of delta26 = 9 (L - 40)   for every 26-point set.

(T) and (SD) leave 17 cuts (`HAND26`). delta26 is positive on four of them: (6,6,14) +12,
(5,7,14) +9, (6,7,13) +5, (7,7,12) +2. It is zero on exactly the five allowed cuts that contain an
8-point section: (5,8,13), (6,8,12), (7,8,11), (8,8,10), (8,9,9). It is negative on the other
eight. Hence **every 26-point Kakeya set in F_3^4 has a class with cut (6,6,14), (5,7,14),
(6,7,13) or (7,7,12), or else every one of its 40 classes has a section of exactly 8 points and it
contains exactly 40 lines (one per direction).** The second alternative is very rigid. It is the
natural first target of an exact search for k_4 = 27 (Lev's conjecture). `test_bad25.py` checks
this part too.

## Limits

- This is a structural lemma, not a new bound: k_4 >= 25 stands, and 25 is not excluded until the
  exact search finishes and its result is checked.
- It relies on (T) and (SD) from `SD.md`. Both are short proofs, but no human mathematician has
  checked them or this note yet.
- Sanity checks in `test_bad25.py`: the 14-cut list is re-derived independently of `sd_bound.py`
  and matches `SD.md` at N = 24 (12 cuts before (HZ)); the identity (*) holds on random 25-point
  sets with brute-force line counts; every cut of the known 27-point Kakeya set passes (T) and (SD).
