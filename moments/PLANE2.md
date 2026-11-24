# A sharper plane-class relaxation for Kakeya sets in F_3^4

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

Code: `plane_lp2.py` (standard library only, builds on `plane_lp.py`). Sanity tests:
`test_plane_lp2.py`. Readout tests: listed at the end. Read `PLANE.md` first.

## Two exact tables for three parallel planes

Inside one hyperplane (a copy of AG(3,3)) take three parallel planes, identify each with F_3^2 and put
them at heights 0, 1, 2. Points x, y, z at heights 0, 1, 2 are collinear exactly when
x + y + z = 0 in F_3^2 (the third point of the line through (x, 0) and (y, 1) is (2y - x, 2), and
2y - x = -x - y). That line has direction (y - x, 1), so distinct lines of this kind carry distinct
directions exactly when their differences y - x are distinct.

For sections A, B, C of sizes a, b, c define

* T3(a, b, c) = the largest number of triples (x, y, z) in A x B x C with x + y + z = 0;
* G(a, b, c) = the largest number of distinct differences y - x over those triples,

maximised over all subsets of F_3^2 of those sizes. `triple_tables()` computes both exactly: A and B
may be translated to contain 0 (this moves C by a translation and shifts all differences by one
vector), and for fixed A, B the best C of size k is found among all 2^9 sets of k sums. Both tables
are symmetric in a, b, c (the direction of a line can be read from any two of its points; checked
in the tests), never exceed the old capacities min(ab, ac, bc) and `cross_cap`, and match hand
values: T3(2,2,2) = 3 (the 4 sums 0, a, b, a + b of A = {0, a}, B = {0, b} take at least 3
values, because a = b and a + b = 0 cannot both hold), G(2,2,2) = 3 (A = B = {0, a}, C = {0, -a}),
T3(3,3,3) = 9 (three copies of one line through 0), G(1,9,9) = 9.

## What changes in the relaxation

1. Condition 4 of `PLANE.md` (the 9 crossing directions of each hyperplane class U containing W)
   now uses G instead of `cross_cap`: summed over the 3 parallel lines of the quotient, G must
   reach 9.
2. The crossing-line capacity in (LX) uses T3 instead of min(ab, ac, bc).
3. **New per-class condition.** For a plane class W with array M, every line of K either has its
   direction in W, and then lies in a cell (at most 0,0,0,1,1,2,3,5,8,12 lines for a cell of
   0..9 points), or crosses three collinear cells (at most T3 of them per quotient line). So every
   one of the 130 arrays satisfies linecap(M) + crosscap3(M) >= L >= Lmin. Arrays with too little
   room for lines are removed before the linear program is solved.

Each change only uses facts that hold for every Kakeya set, so `plane2_feasible(N, Lmin) = False`
still proves that no Kakeya set with N points and at least Lmin lines exists.

## Guards against a wrong bound

* `test_plane_lp2.py` checks the tables (symmetry under all 6 permutations, monotonicity, the hand
  values, never above the old capacities), that all 130 arrays of the real 27-point Kakeya set K27
  pass every new condition with its true number of lines, that K27 solves the sharper linear
  system exactly, that its arrays are found by the enumeration, and that N = 27 stays feasible.
* Every readout test that claims a bound first runs `real_set_guard()` (K27 passes everything).

## Results: read them off the repository

Each readout test is committed only if it passes (and within the 600 s test limit), so the files
present in this folder are the results.

| file present | meaning |
|---|---|
| `test_plane2_k4_at_least_26.py` | no Kakeya set of 25 points: k_4 >= 26 |
| `test_plane2_k4_at_least_27.py` | none of 25 or 26 points either (needs an independent check) |
| `test_plane2_n25_not_excluded.py` | the sharper relaxation still admits 25 points |
| `test_plane2_n25_lines_below_X.py` | a 25-point Kakeya set would have fewer than X lines |
| `test_plane_n25_not_excluded.py` | the first plane relaxation admits 25 points |

If neither `test_plane2_k4_at_least_26.py` nor `test_plane2_n25_not_excluded.py` is present, the
decision did not finish inside the time limit, and the answer will come from a separate computation.

## Limits

These are relaxation bounds: the conditions are necessary, not sufficient. A bound of 27 would settle
Lev's conjecture and would be published only after an independent check. The code was written and
checked by AI agents and has not been reviewed by a human mathematician.
