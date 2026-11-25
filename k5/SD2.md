# SD2: a sharper sum-difference bound, and a second proof of k_5 >= 45

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24, result added 2026-11-25. Made with AI
agents working in MV2's name. **No human mathematician has checked this yet.** Text CC BY 4.0,
code Apache-2.0.

## Status in one paragraph

k_5 is the least size of a Kakeya set in F_3^5 (a set containing a line in each of the 121
directions). Known: k_5 <= 53 (Chung, Du, Wesley, arXiv:2608.23691) and, from our own
moment + sum-difference method, k_5 >= 44 (`moments/SD.md`). This note proves a sharper form
of the sum-difference bound (SD2). With SD2 the moment linear program at N = 44 survives only
through hyperplane cuts of profile (10, 11, 23), and only if at least 20 of the 121 cuts have
that profile. The exhaustive search `k5cut.py` (compete job k5cut-v1, run 2026-11-25
19:42-19:46 UTC, 108 seconds on 4 cores) found **no** (10, 11, 23) cut: verdict CLOSED
(`results/k5cut-v1.json`). So **k_5 >= 45** follows by this route as well. It is a second
proof of the bound proved by hand in `QUAD45.md` (quadric constraint on section pair counts);
the two routes share only the section-counting identities.

## Setting

A Kakeya set K in F_3^n with |K| = N splits along each of the H = (3^n-1)/2 parallel classes
of hyperplanes into sections of sizes (n1, n2, n3). Put T = 3^(n-1). In one class, identify the
three hyperplanes with G = F_3^(n-1) and call the sections A, B, C. The line with direction (v, 1)
meets them in a, a+v, a+2v. Since K is Kakeya, (*) holds: for every v in G there is a in A
with a+v in B and a+2v in C. Write

- r(v) = |A ∩ (B - v)| (columns: pairs (a, b) with b - a = v), and
- r+(s) = |A ∩ (s - B)| (rows: pairs with a + b = s).

Here sum_v r(v) = sum_s r+(s) = xy (x = |A|, y = |B|), and the energy is
E = sum_v r(v)^2 = sum_s r+(s)^2 (both count solutions of a + b' = a' + b). The third point of
the line through a and b is -(a+b). So (*) says: B - A = G (every column is non-empty), and the
rows s with -s in C cover every column.

## Lemma SD2

If A, B, C satisfy (*) with |A| = x, |B| = y, |C| = z, then

    E >= T^2/z + (xy - T),   so   z >= T^2 / (Emax(x, y) - (xy - T)),

where Emax(x, y) is the largest value of sum of r_v^2 over T integers 1 <= r_v <= min(x, y)
with sum xy (fill columns to min(x, y) greedily; `sd2.emax`).

*Proof.* Let S = -C (|S| = z) and W = sum over s in S of r+(s). The T columns are covered by
rows in S, and different columns use different pairs, so W >= T. By Cauchy-Schwarz the rows in
S contribute at least W^2/z to E. Every other row contributes r+(s)^2 >= r+(s), together at
least xy - W. So E >= W^2/z + xy - W. This is increasing in W for W >= z/2, and T >= z/2,
so E >= T^2/z + xy - T. Finally E <= Emax(x, y), because r(v) >= 1 for all T columns,
r(v) <= min(x, y) and sum r(v) = xy. QED.

The earlier SD lemma used only E >= T^2/z. SD2 adds the xy - T pairs that are not needed.

**Example (the case that matters).** For (x, y) = (10, 10), Emax = 2*100 + 4 + 78 = 282. SD gives
z >= 6561/282 = 23.3, so 24. SD2 gives z >= 6561/(282 - 19) = 24.9, so **25**. So the cut
(10, 10, 24) is impossible.

## The moment LP at N = 44

The identities for n = 5 (T = 81, H = 121) are

    sum over classes of s2 = 40 C(44,2) = 37840,  s2 = C(n1,2)+C(n2,2)+C(n3,2),
    sum over classes of s3 = 13 C(44,3) + 27 L >= 172172 + 27*121 = 175439.

The most balanced profile (14, 15, 15) has s2 = 301 and s3 = 1274. Every other allowed profile
has larger s2. Relative to it the totals need an excess of 37840 - 121*301 = 1419 in s2 and at
least 175439 - 121*1274 = 21285 = 15 * 1419 in s3. So some allowed profile must have
(s3 - 1274) >= 15 (s2 - 301). Among all profiles allowed by (T) and SD2, exactly one does:

| profile | s2 - 301 | s3 - 1274 | ratio | status |
|---|---|---|---|---|
| (10, 11, 23) | 52 | 782 | 15.04 | allowed: SD2 needs only z >= 19 |
| (10, 10, 24) | 65 | 990 | 15.23 | excluded by SD2 (needs 25), allowed by SD |
| (9, 11, 24) | 66 | 999 | 15.14 | excluded (needs 29) |
| (9, 10, 25), (8, 11, 25), (9, 9, 26) | | | > 15 | excluded |
| (11, 11, 22) | 40 | 596 | 14.90 | allowed |
| (9, 12, 23) | 54 | 801 | 14.83 | allowed |

`test_sd2.py` checks all of this by exact enumeration: the profile lists, that (10, 11, 23) is
the only allowed profile with ratio >= 15, and that the LP is feasible with it and infeasible
without it. Since every other allowed profile has ratio <= 14.9, k cuts of profile (10, 11, 23)
give at most 782k + 14.9(1419 - 52k) in s3. This reaches 21285 only for **k >= 20**.

**Theorem.** No three sections of sizes (10, 11, 23) in F_3^4 satisfy (*) (search below), hence
k_5 >= 45. (Adding points keeps a set Kakeya, so excluding N = 44 excludes all N <= 44.)

## Facts about a (10, 11, 23) cut, proved by hand

Let |A| = 10, |B| = 11, |C| = 23 satisfy (*). Then:

1. **E >= 315** (SD2: 6561/23 + 29 = 314.3).
2. **Some column has r(v) >= 8.** With every r(v) <= 7, E <= 4*49 + 36 + 76 = 308.
3. **No column has r(v) = 10.** Suppose r(0) = 10 (translate B), so B = A + {y}. Let rho(d)
   be the number of pairs of A in direction d, D the number of directions with rho >= 1, and
   X = 45 - D = sum (rho - 1)+ the "excess". Coverage G = (A - A) ∪ (y - A) and
   |A - A| = 1 + 2D force the overlap to be exactly 20 - 2X, so X <= 10. Expanding
   r(v) = r_AA(v) + [v in y - A] gives E = 110 + 2T', where
   T' = sum rho^2 + sum over a of r_AA(y - a), and so T' >= 102.5. Look at the lines of
   direction d that meet A in 3 or 2 points. Cross pairs between two such lines spread over
   the three directions w + t d. That costs excess 6, 3 or at least 1 for line types
   (3,3), (3,2) and (2,2). It follows that rho >= 5 forces X >= 11 and rho = 4 forces X >= 6.
   With e = rho - 1 <= 3 this gives T' <= 65 - X + sum e^2 + min(2X, (20 - 2X) e_max) <= 95
   in every case (X = 6..10 for e_max = 3, and 65 + 3X otherwise). That contradicts
   T' >= 102.5.
4. So the largest column is k1 = 8 or 9. Translate B so that it sits at v = 0, and let
   I = A ∩ B. Since B - A lies in (I-I) ∪ (I-X_A) ∪ (Y-I) ∪ (Y-X_A), I spans at least 26
   directions if k1 = 9 (excess <= 10), and at least 17 if k1 = 8 (excess <= 11).

## The search (k5cut.py)

The search takes, for k1 in {8, 9, 10} (10 as an independent check of fact 3), every
AGL(4,3)-orbit representative I of k1-subsets of F_3^4 with the allowed excess. Excess never
increases on subsets, so these are generated level by level with canonical forms. It then
takes every completion A = I + X_A and B = I + Y with B - A = F_3^4 and largest column k1.
Every pair (A, B) with E >= 315 is decided exactly by a set-cover search over rows. Any
solution found is re-checked against (*) directly.

The controls are:

- Orbit-stabiliser sums equal C(81, k) on the complete levels.
- In F_3^2 the whole pipeline agrees with brute force for every (x, y).
- In F_3^3 the filtered orbit counts agree with brute force.
- Lev's 27-point set gives a (6, 8, 13) cut that satisfies (*) and every pruning condition.
- In F_3^3, 25 random instances satisfy every pruning condition at their exact minimum z.

## Result (k5cut-v1, 2026-11-25)

Unit tests (`test_sd2`, `test_k5cut_small`) passed inside the job before the search ran. The
orbit-stabiliser checks passed on levels 1 to 5, and Lev's (6, 8, 13) cut passed every pruning
condition (E = 110 >= 78, largest column 5 allowed, excess of I 2 allowed, minimum cover 13).

| largest column k1 | allowed excess | orbit representatives I | completions X_A | kept (missing differences at most 10 times the size of Y) | pairs with B - A = F_3^4 | solutions |
|---|---|---|---|---|---|---|
| 8 | 11 | 75 | 197,100 | 129,367 | 0 | 0 |
| 9 | 10 | 162 | 11,664 | 3,159 | 0 | 0 |
| 10 | 10 | 20 | 20 | 20 | 0 | 0 |

No pair (A, B) of the required kind even satisfies B - A = F_3^4, so the energy and set-cover
stages were never reached. The k1 = 10 row re-derives fact 3 independently of its hand proof.

## Limits

The search relies on facts 1, 2 and 4 to restrict k1 to 8 and 9 (fact 3 is re-checked by the
search), and on the SD2 lemma and the moment LP to reduce N = 44 to the (10, 11, 23) cut. All
of these were proved and checked by AI agents working in MV2's name; none has been checked by
a human mathematician. The hand proof in `QUAD45.md` reaches the same bound by a different
route.

## Credits

The moment method over hyperplane sections follows the Polymath project's "slices" idea for
finite-field Kakeya sets over F_3 (Polymath wiki, 2009). Lev (2009) found the 27-point set in F_3^4.
Chung, Du and Wesley (arXiv:2608.23691) built the 53-point set in F_3^5. The first SD lemma and
k_5 >= 44 are ours (`moments/SD.md`).
