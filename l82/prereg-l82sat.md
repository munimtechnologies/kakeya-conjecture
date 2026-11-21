# Pre-registration: l82sat-v1. Is l(8,2) = 25 or 26? (equivalently, K(8,2) = 26 or 27)

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-21, written before the job runs.
Made with AI agents working in MV2's name. **No human mathematician has checked this note.**
Text CC BY 4.0. Code `l82.py` and tests `test_l82.py` are in this folder, Apache-2.0.

## Question

l(r,2) is the shortest length of a binary linear code with codimension r and covering radius 2.
Equivalently, it is the smallest size of a set S of nonzero vectors of F_2^r such that every
nonzero vector is in S or is a sum of two distinct elements of S, i.e. the smallest 1-saturating
set of PG(r-1,2). Over F_2 a line is a pair {x, x + d}, so a Kakeya set B in F_2^r is exactly a
set with B + B = F_2^r. Translating so that 0 is in B and putting S = B \ {0} gives
B + B = {0} ∪ S ∪ (S ⊕ S), which means B is Kakeya iff S is 1-saturating. Hence
**K(r,2) = l(r,2) + 1**.

Exact values are known for r <= 7: 2, 4, 5, 9, 13, 19. The values 13 and 19 are due to Struik
(IEEE Trans. Inf. Theory 40, 1994), who also proved l(8,2) >= 25. The best upper bound is 26 (the
family phi(r) = 27 * 2^(r/2-4) - 1). The table of Davydov, Marcugini, Pambianco and Wu (arXiv:2609.16078v2,
11 Oct 2026, Table 5.1) lists r = 8 as "26", without the asterisk it uses for exact values. Kéri's
covering-code tables give K(25,2) <= 2^17. That bound is only the trivial doubling of
Etzion-Greenberg's nonlinear K(23,2) <= 2^15, not a linear [25,17] code (checked against the
per-cell sources in github.com/chrispbsmolen/covering-codes). **As far as we found, whether
l(8,2) is 25 or 26 has been open since 1994.**

## Reduction (proved by hand; fixed in advance)

Let S be 1-saturating in F_2^8 with |S| = 25.

- **L1.** S spans F_2^8. Otherwise S and all its sums lie in a proper subspace.
- **L2.** For each nonzero functional g, let a_g = |{s in S : g(s) = 1}|. The 128 vectors with g = 1
  are in S_1 or of the form s_0 + s_1 with g(s_0) = 0 and g(s_1) = 1. So a_g (26 - a_g) >= 128,
  i.e. 7 <= a_g <= 19.
- **L3.** The a_g are the nonzero weights of the binary [25, 8] code {(g(s))_s} (dimension 8 by
  L1). Ytrehus and Helleseth (IEEE Trans. Inf. Theory 36, 1990) proved that no [25,8,10] code
  exists, so min_g a_g <= 9. The Griesmer bound alone already gives min_g a_g <= 10, because
  d = 11 would need length 26.
- **Coordinates.** Choose f of minimum weight a in {7, 8, 9} (or 10) and coordinates in which f is
  the last bit. Then S = T ∪ {(u, 1) : u in U}, with T ⊂ F_2^7 \ {0}, |T| = 25 - a, U ⊂ F_2^7 and
  |U| = a, and:
  - (C1) U + (T ∪ {0}) = F_2^7;
  - (C2) T ∪ (T ⊕ T) ∪ (U ⊕ U) ⊇ F_2^7 \ {0}, where ⊕ means sums of two distinct elements;
  - (C3) every other functional has weight in [a, 19].
- **Symmetry.** The maps (y, s) -> (Ay + sc, s), with A in GL(7,2), fix f and send (T, U) to
  (AT, AU + c). So U may be replaced by any affine image. Choose an affine basis of U inside U:
  then U = {0, e_1, ..., e_r} ∪ Y, where r is the affine rank and Y is a set of vectors of weight
  >= 2 in <e_1..e_r>. Y is taken up to permutations of e_1..e_r.
- For each such U, the existence of T is one CP-SAT model with 127 point variables. **If every
  model is infeasible for a = 7, 8, 9, then l(8,2) = 26 (given Ytrehus-Helleseth); if also for a = 10,
  then l(8,2) = 26 using only the Griesmer bound.** A feasible model gives a 25-point set, which is
  re-checked by brute force on all 255 points.

## Self-checks (a failure = BUG, or SANITY-INCOMPLETE if only time ran out)

- Unit tests:
  - brute-force minima 2, 4, 5 for r = 2, 3, 4;
  - (C1)+(C2) agree with the direct test on 300 random sets;
  - the split weights agree with the dual weights;
  - the normal forms meet every AGL(3,2)-orbit of subsets of F_2^3;
  - the CP-SAT model with T fixed agrees with the direct test on 40 random sets;
  - the pipeline reproduces l(4,2) = 5, l(5,2) = 9 and l(6,2) = 13.
- In the job, the same pipeline (with the averaging bound in place of L3) must reproduce
  l(5,2) = 9, l(6,2) = 13 and **l(7,2) = 19 (Struik)**: every case at N = l - 1 infeasible, and a
  verified set found at N = l. The time limit for this part is 50 minutes.
- Every set a solver returns passes the brute-force check.

## Verdicts

- **L82-EQ-25**: a verified 25-point 1-saturating set in PG(7,2). Then l(8,2) = 25 and K(8,2) = 26.
  This would be a new record: it improves the table entry 26 and meets Struik's 1994 lower bound.
  The set itself is the certificate.
- **L82-EQ-26**: every case for a = 7..10 infeasible. Then l(8,2) = 26 and K(8,2) = 27
  (computer-assisted, CP-SAT, Griesmer bound only).
- **L82-EQ-26-YH**: a = 7..9 all infeasible, a = 10 not finished. Same conclusion, also relying on
  Ytrehus-Helleseth.
- **UNDECIDED**: some cases are open. We report the open (a, U) and file -v2 for those only.
- **SANITY-INCOMPLETE**: no bug, but the reproduction of l(7,2) = 19 did not finish. The main
  counts are reported but not claimed.
- **BUG**: a self-check fails. Nothing is claimed; we fix the code and file -v2.

CP-SAT infeasibility is not an independently checkable certificate. An EQ-26 verdict is therefore
reported as "computer-assisted (CP-SAT)", and a -v2 will re-prove the cases with a second method
(a DRAT-checked SAT encoding).

## Predictions (written before the run)

- Unit tests pass: 80%.
- L82-EQ-25: 15%. L82-EQ-26: 25%. L82-EQ-26-YH: 10%. UNDECIDED: 25%. SANITY-INCOMPLETE: 10%.
  BUG: 15%.
- The a = 7 cases (excess 5 in (C1)) are all infeasible in under 30 s each: 70%.
