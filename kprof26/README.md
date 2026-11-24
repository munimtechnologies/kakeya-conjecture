# kprof26: is there a Kakeya set of 26 points in F_3^4?

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. No human
mathematician has checked this. Code licence: Apache-2.0; text CC BY 4.0.

The smallest Kakeya set known in F_3^4 has 27 points (Lev, 2009), and Lev conjectured that 27
is optimal. A Kakeya set with at most 26 points can be enlarged to exactly 26 points, so **if no
26-point Kakeya set exists, then k_4 = 27**. This folder holds an exact search for one. It
imports `../kprof25/kprof25.py` (the 25-point search) for the shared model, group and checks.

## The reduction (proof sketch; the code checks every number)

Let K be any 26-point set in F_3^4. For each of the 40 classes of parallel hyperplanes, K is cut
into three sections of sizes (n1, n2, n3). With s2 = sum C(n_i, 2), s3 = sum C(n_i, 3) and L the
number of lines inside K, double counting gives, summed over the 40 classes,

    sum s2 = 13 C(26,2) = 4225,        sum s3 = 4 C(26,3) + 9 L = 10400 + 9 L.

Put delta26 = s3 - 8 s2 + 576 = g(n1) + g(n2) + g(n3) + 576 with g(n) = n(n-1)(n-26)/6. Then

    sum over the 40 classes of delta26 = 9 (L - 40).

Because g(n) + g(18 - n) = -408 for every n, **every cut with an 8-point section has delta26 = 0**.
Among the cuts allowed for a Kakeya set by the counting rule (T) and the sum-difference lemma (SD)
of `moments/SD.md`:

| delta26 | cuts |
|---|---|
| +12, +9, +5, +2 | (6,6,14), (5,7,14), (6,7,13), (7,7,12) |
| 0 | (5,8,13), (6,8,12), (7,8,11), (8,8,10), (8,9,9) |
| < 0 | every other allowed cut |

A Kakeya set has L >= 40. So a 26-point Kakeya set has a class with one of the four positive
cuts (cases P1-P4, taking the first one in the list that occurs anywhere and forbidding the
earlier ones in every class), or else every class has a zero cut and L = 40 exactly, i.e. exactly
one line of K in each of the 40 directions (cases Z1-Z5, again taking the first zero cut that
occurs).

**Hand remark (cases Z3-Z5 need no search).** With only zero cuts, the s2 identity reads
sum over classes of (s2 - 100) = 225, and s2 - 100 equals 16, 9, 4, 1, 0 for the five zero cuts in
the order above. Without (5,8,13) and (6,8,12) the left side is at most 4 x 40 = 160 < 225, so
Z3, Z4 and Z5 are impossible by arithmetic. Only Z1 and Z2 need the computer. (The search still
runs Z3-Z5 as a check; the solver should refute them at once.)

## The search

For each case, an affine map sends the chosen class to the last coordinate with sections sorted
by size, and the smallest section is replaced by a representative of its orbit under AGL(3,3)
(orbit sizes are checked to add up to C(27, k)). The model has 54 Boolean point variables for
the other two slices, a Boolean for each of the 1080 lines (AND of its points), at least one line
per direction, the allowed-cut table (with s2 and s3) in all 40 classes, and both identities with
exact L. In Z cases it also fixes L = 40 and at most one line per direction (both implied by the
table and identities, added for propagation). CP-SAT decides each representative; HiGHS
re-solves every infeasible one with kprof25's 0/1 profile-choice encoding as an independent
second solver. Every solution found is re-checked by brute force.

Sanity checks in every run: the known 27-point set is accepted with fixed slices, rejected when
one of its own cuts is forbidden, and re-found by a free search (kprof25); the extra L and
lines-per-direction constraints accept the 27-point set's true values and reject values one off;
and for each of the 27 sets "known set minus one point" the fixed-slice model agrees with a
brute-force check.

## Jobs

- `kprof26-p-v1`: cases P1-P4.  `kprof26-z-v1`: cases Z1-Z5.
- Verdicts: `PART-CLOSED` (both solvers, all representatives infeasible), `PART-CLOSED-CPSAT-ONLY`,
  `UNDECIDED` (some representatives timed out), `FOUND-26` (a 26-point Kakeya set: it would
  refute Lev's conjecture), `BUG`.
- Both jobs `PART-CLOSED` would mean k_4 = 27, as a computer-assisted result that still needs an
  independent certificate (for example DRAT proofs from a SAT solver) before it is claimed.
