# A sum-difference lemma for hyperplane sections: k_4 >= 25 and k_5 >= 44

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-23. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

k_n is the least size of a Kakeya set in F_3^n (a set containing a full line in every direction).
Known: k_2 = 7, k_3 = 13, k_4 <= 27 (a 2009 computer search; the Polymath wiki says "it seems likely
that k_4 = 27"), k_5 <= 53 (Chung-Du-Wesley, arXiv:2608.23691). The general lower bound of Bukh and
Chao gives k_4 >= 18 and k_5 >= 32. This folder earlier proved k_4 >= 24 (README.md) and k_5 >= 41
(F35.md) from hyperplane moments.

**New here:** k_4 >= 25 and k_5 >= 44. The code is `sd_bound.py`; `test_sd_bound.py` checks every
claim below in exact integer arithmetic. As a sanity check, the same code gives back the known
values k_2 = 7 and k_3 = 13 exactly, and it does not rule out the known sets of 27 points in F_3^4
or 53 points in F_3^5.

## Setting

K is a Kakeya set in F_3^n, |K| = N, with L >= H = (3^n - 1)/2 lines inside it. The hyperplanes
fall into H parallel classes of three. For each class, the profile (n1, n2, n3) gives the sizes of
K's three sections. With a1 = (3^(n-1) - 1)/2 and a2 = (3^(n-2) - 1)/2:

- sum over classes of C(n1,2) + C(n2,2) + C(n3,2) = a1 C(N,2)
- sum over classes of C(n1,3) + C(n2,3) + C(n3,3) = a2 C(N,3) + (a1 - a2) L

Let T = 3^(n-1), the number of directions that are not parallel to a given class.

**(T)** n_i n_j >= T. Each non-parallel direction has a line of K. That line meets H_i and H_j in
one point each, and different directions give different pairs.

## The lemma (SD)

Fix a class and two of its sections, S_i in H_i and S_j in H_j, with third section S_k. Then

    n_k >= T^2 / Emax(n_i, n_j),

where Emax(x, y) is the largest value of sum_d r(d)^2 over T integers with 1 <= r(d) <= min(x, y)
and sum_d r(d) = x y. It is attained by putting as many r(d) as possible equal to min(x, y), one at
the remainder and the rest at 1.

*Proof.* Choose one line of K in each of the T non-parallel directions, and let R be the set of the
T pairs (p, q) in S_i x S_j where these lines meet H_i and H_j. Two lines through p and q in
different directions have different difference vectors q - p: all such vectors have the same
component across the class, so the direction fixes the vector. So the T differences in R are
distinct, and they are all T possible values. Over F_3 the third point of the line through p and q
is -(p + q), and it lies in S_k. So the sums p + q over R take at most n_k values. Let sigma(x) be
the number of pairs in R with sum x. By Cauchy-Schwarz,

    T^2 = (sum sigma)^2 <= n_k * sum sigma(x)^2.

Also sum sigma^2 <= #{(p,q,p',q') in (S_i x S_j)^2 : p + q = p' + q'}, the additive energy. Since
p + q = p' + q' holds exactly when q - p' = q' - p, the map (p,q,p',q') -> (p',q,p,q') shows that
this equals the difference energy #{q - p = q' - p'} = sum_d r(d)^2. Here r(d) counts the pairs in
S_i x S_j with difference d. Each r(d) is at least 1, because R uses every difference, and at most
min(n_i, n_j); the r(d) add up to n_i n_j. Hence sum sigma^2 <= Emax(n_i, n_j). QED.

Example in F_3^4 (T = 27): sections of 5 and 6 points give Emax = 42 (one difference value
repeated 4 times, the other 26 once), so the third section has at least 729/42 > 17 points. The
profile (5,6,13) is therefore impossible.

**(HZ), n = 4 only** (from N24.md, re-checked here by exhaustive search). Each of the 13 directions
parallel to a class has its line inside one section. A section of m <= 8 points carries at most
D3(m) directions, with D3(5..8) = 2, 3, 4, 4, because the fewest points on lines in k distinct
directions of AG(3,3) are f(1..5) = 3, 5, 6, 7, 9. So the profile (8,8,8) is impossible.

## k_4 >= 25

For N = 24, 17 profiles satisfy (T). (SD) removes (3,9,12), (3,10,11), (4,7,13), (4,8,12) and
(5,6,13), and (HZ) removes (8,8,8). Every one of the 11 profiles that are left satisfies

    11 s3 - 85 s2 <= -5300     (with equality at (7,8,9) and (6,6,12)),

where s2 = sum C(n_i,2) and s3 = sum C(n_i,3). Summing over the 40 classes and using the identities
with a1 = 13 and a2 = 4:

    11 (4 C(24,3) + 9 L) - 85 * 13 C(24,2) <= 40 * (-5300)   =>   99 L <= 3924   =>   L <= 39,

which contradicts L >= 40. Without (SD), the profile (5,6,13) keeps N = 24 possible; without (HZ),
(8,8,8) does. N <= 23 was already excluded (README.md). **So k_4 >= 25.** The same relaxation is
feasible at N = 25, so this method stops at 25. The candidates left for k_4 are 25, 26 and 27.

## k_5 >= 44

Here a1 = 40, a2 = 13, T = 81 and H = 121.

- N = 41: moments and (T) alone already exclude it, with 13 s3 - 185 s2 <= -34918 on every allowed
  profile (equality at (13,14,14) and (9,9,23)), so 351 L <= 41382, i.e. L <= 117 < 121. (F35.md
  recorded only N <= 40.)
- N = 42: (SD) removes (9,9,24) (equal products 81 force a third section of 81 points) and
  (9,10,23) (it needs at least 41). Then 3 s3 - 43 s2 <= -8463 on every allowed profile, with
  equality at (14,14,14), so 81 L <= 9177, i.e. L <= 113.
- N = 43: (SD) removes (9,11,23) (it needs 27), (10,10,23) (it needs 24) and (9,10,24),
  (9,9,25) and (8,11,24). Then 22 s3 - 321 s2 <= -66101 with equality at (14,14,15) and
  (10,11,22), so 594 L <= 66773, i.e. L <= 112.

**So k_5 >= 44**, against the published 32 (Bukh-Chao) and our earlier 41. At N = 44 the relaxation
is feasible, for example with (14,15,15) and (10,10,24).

## Higher dimensions (added 2026-11-24)

The same relaxation (`sd_bound.feasible(n, N, sd=True)`), decided in exact arithmetic, gives the
table below. Each lower bound is a test file that was committed only after it passed; the next
size up was tried the same way and its test failed (the relaxation is feasible there), so each
value is the most this method gives.

| n | Bukh-Chao lower bound | this method | test file | first size the method allows |
|---|---|---|---|---|
| 4 | 18 | 25 | `test_sd_bound.py` | 25 |
| 5 | 32 | 44 | `test_sd_bound.py` | 44 |
| 6 | 57 | 78 | `test_k6_at_least_78.py` | 78 |
| 7 | 103 | 137 | `test_k7_at_least_137.py` | 137 |

For n = 6 and 7 we have not surveyed the best known constructions; the classical one has
2^(n+1) - 1 points (127 and 255). The method's advantage over Bukh-Chao shrinks in relative terms
as n grows, so it is a small-dimension tool, not a new asymptotic bound.

A second fact used by the exact search (`moments/min_section.py`, tested): summing |K n H|^2 over
the 120 hyperplanes of F_3^4 gives 40 N + 13 N (N - 1) for every set, so every 25-point set has a
hyperplane section of at most 6 points and every 26-point set one of at most 7.

## Limits

- These are lower bounds from a linear relaxation; no Kakeya set of 25 points in F_3^4 is claimed,
  and Lev's conjectured value k_4 = 27 is not proved.
- Two finer relaxations over the 130 plane classes (`PLANE.md`, `PLANE2.md`) still admit 25 points
  (readouts `test_plane_n25_not_excluded.py`, `test_plane2_n25_not_excluded.py`), even with 46 or
  more lines.
- The (SD) lemma is a standard additive-energy argument (Cauchy-Schwarz plus E(A,B) = E(A,-B));
  what is new is using it on the sections of a Kakeya set together with the moment identities.
- Checked by the code in this folder and by hand. No human mathematician has checked it yet.
