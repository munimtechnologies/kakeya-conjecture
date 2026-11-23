# Every Kakeya set in F_3^4 has at least 24 points

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-22. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

Let k_4 be the smallest size of a Kakeya set in F_3^4 (a set containing a full line in every
direction). The smallest set known has 27 points (Lev, 2009; the Polymath wiki adds "it seems
likely that k_4 = 27 holds"). The best lower bound we could find in print is 18, from the general
bound of Bukh and Chao (K(n,q) >= q^n / (2 - 1/q)^(n-1), which gives 17.5 for q = 3, n = 4).

**Theorem.** k_4 >= 24. So k_4 is one of 24, 25, 26, 27.

The proof below is short and uses only counting. The same method gives exactly the known
values K(2,3) = 7 and K(3,3) = 13, which is a check that it does not overreach.
`kakeya_moments.py` (standard library, exact integer arithmetic) re-checks every number, and
`test_kakeya_moments.py` asserts them. The same method in F_3^5 gives k_5 >= 41 (`F35.md`), and
`N24.md` narrows the remaining case N = 24.

## Set-up

Let K be a Kakeya set in F_3^4 with N points. An affine hyperplane is a 3-flat; there are 120 of
them, in 40 parallel classes of 3. For a hyperplane H write n_H = |K n H|. For one class, the
sorted triple of its three section sizes is the class's *profile*.

Standard flat counts in AG(4,3): a point lies in 40 hyperplanes (one per class), a line in 13,
an affine plane in 4. (Hyperplanes through a line with direction d correspond to the 2-dimensional
subspaces of F_3^4 / <d>, which number 13; through a plane, to the 1-dimensional subspaces of a
2-dimensional quotient, 4.)

**(T) Transversal condition.** Fix a class. Of the 40 directions, 13 are parallel to it and 27
are not. A line in a non-parallel direction meets each of the three hyperplanes of the class in
exactly one point, and any two of those points determine the line. K contains a line in each of
the 27 directions, so for any two hyperplanes H, H' of one class, n_H * n_H' >= 27. In
particular n_H + n_H' >= 11.

**(M) Moment identities.** Count pairs and triples of points of K by the hyperplanes containing
them. Two points span a line (13 hyperplanes). Three points span a line if collinear (13
hyperplanes) and otherwise an affine plane (4 hyperplanes). Over F_3 a collinear triple is a whole
line. Let L be the number of lines contained in K; L >= 40 because K has a line in each of the 40
directions. Then

- sum over all 120 hyperplanes of C(n_H, 2) = 13 C(N, 2),
- sum over all 120 hyperplanes of C(n_H, 3) = 4 C(N, 3) + 9 L >= 4 C(N, 3) + 360.

For a profile P = (a, b, c) write s2(P) = C(a,2) + C(b,2) + C(c,2) and s3(P) likewise with C(.,3).
The sums above are the sums of s2 and s3 over the 40 classes.

## Proof that N <= 23 is impossible

Each N is excluded separately (a Kakeya set of N points plus any extra point is a Kakeya set of
N + 1 points, so excluding every N <= 23 proves k_4 >= 24).

**N <= 16.** By (T) the three section sizes pairwise sum to at least 11, so N >= 17.

**N = 17, 18, 19.** The only profiles allowed by (T) are (5,6,6) for N = 17;
(4,7,7), (5,6,7), (6,6,6) for N = 18; and (4,7,8), (5,6,8), (5,7,7), (6,6,7) for N = 19.
Their largest s2 values are 40, 48 and 55, so the sum of C(n_H, 2) over the 40 classes is at most
1600, 1920 and 2200. But it must equal 13 C(N,2) = 1768, 1989 and 2223. Contradiction.

**N = 20, 21, 22, 23.** Pick a number c. Summing s3 - c s2 over the 40 classes gives
4 C(N,3) + 9L - 13c C(N,2), and this is at most 40 times the largest value of s3 - c s2 over the
allowed profiles. That bounds L from above:

| N | allowed profiles (by (T)) | c | largest s3 - c s2 (profile) | so L is at most |
|---|---|---|---|---|
| 20 | (4,7,9) (4,8,8) (5,6,9) (5,7,8) (6,6,8) (6,7,7) | 6.5 | -280.5 (6,7,7) | 30 |
| 21 | (3,9,9) (4,7,10) (4,8,9) (5,6,10) (5,7,9) (5,8,8) (6,6,9) (6,7,8) (7,7,7) | 7 | -336 (7,7,7) | 38 |
| 22 | (3,9,10) (4,7,11) (4,8,10) (4,9,9) (5,6,11) (5,7,10) (5,8,9) (6,6,10) (6,7,9) (6,8,8) (7,7,8) | 7 | -364 (7,7,8) | 33 |
| 23 | (3,9,11) (3,10,10) (4,7,12) (4,8,11) (4,9,10) (5,6,12) (5,7,11) (5,8,10) (5,9,9) (6,6,11) (6,7,10) (6,8,9) (7,7,9) (7,8,8) | 7.4 | -422.8 (7,8,8) | 38 |

Each bound is below 40, which contradicts L >= 40.

Worked line for N = 23 (all integers after multiplying by 5): the sum of 5 s3 - 37 s2 over the 40
classes equals 5(7084 + 9L) - 37 * 3289 = 45L - 86273. Every allowed profile has
5 s3 - 37 s2 <= -2114 (attained by (7,8,8): 5 * 147 - 37 * 77), so 45L - 86273 <= -84560, i.e.
L <= 38. The full list of profile values is in the test file and re-computed by the script.

Worked line for N = 22: with c = 7, the profile values s3 - 7 s2 are -383, -370, -373, -374,
-365, -367, -368, -365, -365, -365, -364 in the order of the table. So
6160 + 9L - 21021 <= 40 * (-364), i.e. L <= 33.

## What the method does not reach

At N = 24 the conditions can be met on paper: 28 classes of profile (8,8,8) and 12 of profile
(5,6,13) give the right sum of C(n_H, 2) (3588) and enough triples (8496 >= 8456). So closing
24, 25 and 26 needs more structure than these two moments. `N24.md` (2026-11-23) adds one such
piece: lines parallel to a class must fit inside its sections, which rules out (8,8,8) and caps a
24-point set at 41 lines, but does not exclude N = 24. The queued exact computation (`kmin43/`,
the k_4 = 27 test) remains the way to decide k_4; this bound shrinks what it must rule out.

## Sanity checks

- In F_3^2 the method gives exactly 7 = K(2,3) (Blokhuis and Mazzocca's formula).
- In F_3^3 it gives exactly 13, the known value of K(3,3) (first found by computer search; here
  N = 12 is excluded with c = 3.2, which leaves room for at most 12 lines where 13 are needed).
- In F_3^3 at N = 13 the conditions can be met, for example when every plane class has profile
  (3,4,6), so the method stops exactly at the true value.

## Prior art and honest limits

- The two moment identities are standard double counting (the "standard equations" of finite
  geometry). We did not find them applied to Kakeya sets in F_3^4 in the sources we checked:
  Lev's 2009 comment on Tao's blog, the Bukh-Chao paper (arXiv:2108.00074), and the
  2025-26 papers listed in `CORRECTIONS.md`.
- **The Polymath wiki's "Kakeya problem" page** (last modified 5 June 2009; read on 2026-11-22
  through the Internet Archive snapshot of 23 December 2020, because michaelnielsen.org now
  returns 404) states k_1 = 3, k_2 = 7, k_3 = 13 and k_4 <= 27 by computer, and "it seems likely
  that k_4 = 27 holds". It gives **no numerical lower bound for k_4 or k_5**. Its "slices argument"
  uses the same three slices as condition (T) and, with a Katz-Tao estimate, concludes
  k_n >= 3^(6(n-1)/11), about 6.0 for n = 4 and 11.0 for n = 5. So the bounds here are new as far
  as we can tell, and the page confirms the conjecture k_4 = 27 that this bound narrows towards.
- The work was carried out with AI agents working in MV2's name. No human expert has checked it.
