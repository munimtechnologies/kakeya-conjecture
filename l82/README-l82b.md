# l82b: the open cases of l82sat-v1, with symmetry breaking

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this note.** Text CC BY 4.0; code Apache-2.0.

## Where things stand

The question (see `prereg-l82sat.md`): is the smallest 1-saturating set of PG(7,2) of size 25
(Struik's 1994 lower bound) or 26 (best known construction)? Equivalently: is there a binary linear
[25, 17] code with covering radius 2? As far as we found, this has been open since 1994; the
2026 tables of Davydov, Marcugini, Pambianco and Wu (arXiv:2609.16078v2) still list 26 as an
upper bound only.

`l82sat-v1` (finished 2026-11-24 06:04 UTC, 5.2 hours, 4 cores) passed every self-check: with the
same pipeline it reproduced l(5,2) = 9, l(6,2) = 13 and Struik's l(7,2) = 19. Then for n = 8,
N = 25 it found no set, and CP-SAT proved infeasible:

| min weight a | cases | infeasible | open (time limit) | not run |
|---|---|---|---|---|
| 7  | 14  | 13 | 1  | 0   |
| 8  | 38  | 3  | 35 | 0   |
| 9  | 123 | 89 | 34 | 0   |
| 10 | 499 | 0  | 0  | 499 |

So the verdict is UNDECIDED. Every case at a = 8 with affine rank 3 to 7 except three, and every
a = 9 case of affine rank 6 or 7, hit the 300 s limit. Since Ytrehus and Helleseth (1990) showed
that no [25, 8, 10] code exists, closing the 70 open cases with a = 7, 8, 9 would settle
l(8,2) = 26. The a = 10 cases would only be needed for a proof that does not rely on that result.

## What v2 adds

After the minimum-weight functional f and the set U (the last-coordinate-1 part, |U| = a) are
fixed, the search is for T (|T| = 25 - a) in F_2^7 \ {0}. The maps

    phi_(A,c) : (y, s) -> (A y + s c, s),   A in GL(7,2), c in F_2^7, with A U + c = U,

form a group G_U. Each one sends a solution (T, U) to the solution (A T, U):

- phi is an invertible linear map of F_2^8 that fixes the functional f = last bit, so it maps
  the set S = T u {(u,1)} to S' = A T u {(u,1) : u in A U + c} = A T u {(u,1) : u in U}.
  A linear bijection preserves 1-saturation and the multiset of weights of the other
  functionals, so (C1), (C2) and (C3) hold for (A T, U), and |A T| = |T|.

So the CP-SAT model for fixed U is invariant under T -> A T. For any subset H of G_U, the
lexicographically largest member of each orbit satisfies t >=_lex t o A for every A in H.
Adding these constraints therefore removes only symmetric copies and never turns a feasible case
into an infeasible one. l82b uses two kinds of elements of G_U:

1. **Transvections** e_j -> e_j + e_i with e_j outside the span of U (c = 0). They fix U
   pointwise. When U spans only an r-dimensional subspace, the group they generate has order
   2^(r(7-r)) |GL(7-r,2)|, which is large for the low-rank cases.
2. **Affine maps that permute U**, found by depth-first search over the images of the affine
   basis {0, e_1, ..., e_r} of U (up to 40, those that move the fewest points). For the
   affine-basis case of a = 8 this group is the symmetric group S_8.

The tests (`test_l82b.py`) check that every map found really fixes U, that the solver's
solutions are mapped to solutions, that the lex constraint means exactly t >=_lex t o P (brute
force on 4 variables), and, most importantly, that with and without the new constraints the
model gives the same feasibility answer on every case of l(5,2) = 9 and l(6,2) = 13, at N = l - 1
and N = l. Each job also re-runs the l(5,2) and l(6,2) reproductions (one job also l(7,2) = 19)
with the constraints on.

## Jobs (2026-11-24)

`l82b-a8p0-v1` and `l82b-a8p1-v1` split the open a = 7 case and the 35 open a = 8 cases;
`l82b-a9p0-v1` and `l82b-a9p1-v1` split the 34 open a = 9 cases. Each case gets 120 s, then 900 s,
then an equal share of what is left of about 5 hours. Verdicts per job: FOUND-25 (a verified
set: l(8,2) = 25), PART-CLOSED (all its cases infeasible), PART-OPEN, SANITY-INCOMPLETE or BUG.

**Only if all four jobs end PART-CLOSED** do we have l(8,2) = 26, computer-assisted with CP-SAT
and resting on Ytrehus-Helleseth for a >= 10. CP-SAT infeasibility is not an independently
checkable certificate, so that claim would then be re-proved with a DRAT-checked SAT encoding
before any article calls it a theorem.
