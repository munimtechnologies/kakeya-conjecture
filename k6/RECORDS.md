# Kakeya sets in F_3^6 with 105 points, in F_3^7 with 211, in F_3^8 with 422, in F_3^9 with 844 and in F_3^10 with 1702

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-26. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

A set K in F_3^n is a Kakeya set if it contains a full line {a, a+d, a+2d} in every
direction d. Let k_n be the least size of a Kakeya set in F_3^n.

* `k6/sets/kakeya6_105.json` is a Kakeya set in F_3^6 with **105 points**, so **k_6 <= 105**.
* `k6/sets/kakeya7_211.json` is a Kakeya set in F_3^7 with **211 points**, so **k_7 <= 211**.
* `k6/sets/kakeya8_422.json` is a Kakeya set in F_3^8 with **422 points**, so **k_8 <= 422**
  (classical construction: 511; Lev's pattern value: 427).
* `k6/sets/kakeya9_844.json` is a Kakeya set in F_3^9 with **844 points**, so **k_9 <= 844**
  (classical construction: 1023; Lev's pattern value: 853).
* `k6/sets/kakeya10_1702.json` is a Kakeya set in F_3^10 with **1702 points**, so
  **k_10 <= 1702** (classical construction: 2047; Lev's pattern value: 1707).
* Earlier sets, kept and still checked: `kakeya6_106.json` (106 points),
  `kakeya7_213.json` (213 points), both found on 2026-11-25, and `kakeya8_427.json`
  (427 points) and `kakeya9_846.json` (846 points), found on 2026-11-26.

All sets are checked at every commit by `k6/test_sets.py`, a stand-alone checker that
shares no code with the search: it decodes each integer i = sum_k x_k 3^(n-1-k) to the
point (x_0, ..., x_(n-1)) and checks all 364 (n = 6), 1093 (n = 7), 3280 (n = 8) or
9841 (n = 9) directions one by one. For n = 10 it collects the covered directions from
pairs of points instead (three distinct points of F_3^n are collinear iff they sum to 0)
and checks that all 29524 directions appear; the two methods are compared on small cases.
It also checks that no single point of the 105- or 211-point set can be deleted.

The 422-, 846-, 844- and 1702-point files are stored in two parts. The layer x_0 = 0 of the
422-point set is the file `kakeya7_211.json`, and its other 211 points are listed as
`free_points` p, standing for the integer 2187 + p (the file also keeps the compact code the
search printed, `points_code`, and the checker verifies that it decodes to the same list).
The layer x_0 = 0 of the 846- and 844-point sets is the file `kakeya8_422.json`, and their
other points are stored only as the search's compact code (`points_code`: sorted values,
gaps in base 32, see the file), standing for 6561 + p. The layer x_0 = 0 of the 1702-point
set is the 846-point set, and its other 856 points stand for 19683 + p. Both checkers decode
the codes with their own decoders.

The layered sets have a second, independent proof in `k6/test_layers.py`. The layer
x_0 = 0 of the 211-point set is exactly the 105-point set of F_3^6, which covers every
direction (0, d). A direction (1, v), v in F_3^6, is covered iff some x in that layer has
x + v in the layer x_0 = 1 and x + 2v in the layer x_0 = 2; the test checks all 729
vectors v. The other two layers have 56 and 50 points, together **106 = 2 x 53**, and no
point of either can be removed. In the same way:

| set | layer x_0 = 0 | layers x_0 = 1, 2 | vectors v checked |
|---|---|---|---|
| 211 in F_3^7 | 105-point set of F_3^6 | 56 + 50 = 106 | 729 |
| 422 in F_3^8 | 211-point set of F_3^7 | 96 + 115 = 211 | 2187 |
| 846 in F_3^9 | 422-point set of F_3^8 | 204 + 220 = 424 | 6561 |
| **844** in F_3^9 | 422-point set of F_3^8 | 207 + 215 = **422** | 6561 |
| **1702** in F_3^10 | 846-point set of F_3^9 | 443 + 413 = **856** | 19683 |

So the chains 105 -> 211 -> 422 -> 844 and 105 -> 211 -> 422 -> 846 -> 1702 have a proof
that uses only the direct check of the 105-point set in F_3^6. For the 1702-point set the
covered vectors are collected as {y - x : x in L_0, y in L_1, -(x + y) in L_2}.

## Why it matters

The best constructions we found in the literature for these dimensions are the classical
set {0,1}^n u {0,2}^n with 2^(n+1) - 1 points: 127 for n = 6, 255 for n = 7, 511 for
n = 8, 1023 for n = 9 and 2047 for n = 10 (Polymath wiki, Kakeya problem page; V. F. Lev's
2009 comment on T. Tao's blog). Products of smaller Kakeya sets are worse (3 x 53 = 159,
7 x 27 = 189, 13 x 13 = 169; in F_3^8, 27 x 27 = 729 and 3 x 211 = 633; in F_3^9, even
with our own records as factors, 3 x 422 = 1266, 7 x 211 = 1477, 13 x 105 = 1365 and
27 x 53 = 1431, all above the classical 1023). Kopparty, Lev, Saraf and Sudan
([arXiv:1003.3736](https://arxiv.org/abs/1003.3736), Theorem 8) have a random construction
of size (3 / 2^(2/3))^(n + O(sqrt n)), growth base about 1.8899, but no explicit sets for
small n. Chung, Du and Wesley ([arXiv:2608.23691](https://arxiv.org/abs/2608.23691)) found 53
points in F_3^5 and treated dimensions up to 5 only.

In 2009 Lev ([comment 994](https://terrytao.wordpress.com/2009/03/04/dhj3-900-999-density-hales-jewett-type-numbers/comment-page-3/#comment-36694))
computed k_3 = 13 and k_4 <= 27 and wrote "I am very curious to know whether k_5 = 53:
notice the pattern in 3, 7, 13, 27, 53, ...". Chung, Du and Wesley state this as the guess
that the recurrence k_n = k_(n-1) + 2 k_(n-2) continues and say that whether it continues
remains open. The pattern continues 107, 213, 427, 853, 1707.

**Consequence.** k_6 <= 105 < 107, k_7 <= 211 < 213, k_8 <= 422 < 427,
k_9 <= 844 < 853 and k_10 <= 1702 < 1707: our sets are below the pattern in all five
dimensions (by 2, 2, 5, 9 and 5 points). Since k_2 = 7 and k_3 = 13, the recurrence cannot
hold at all of n = 4, 5 and 6 (otherwise k_6 = 107). With k_4 = 27 (Lev's set and our
computer proof `PROOF26.md`, not human-checked) it fails at n = 5 or at n = 6: either
k_5 <= 52, or k_5 = 53 and k_6 <= 105 = k_5 + 2 k_4 - 2. We do not know which.

**The layer observation.** Each record is built on the previous one: a record set of
F_3^(n-1) in one layer plus a search over the two other layers. The extra cost was
54 on the 53-point layer in F_3^6 (total 107), 106 on the 105-point layer in F_3^7 (211),
211 on the 211-point layer in F_3^8 (422), 422 on the 422-point layer in F_3^9 (844; the
other runs gave 422, 423, 424, 424, 424) and 856 on the 846-point layer in F_3^10 (1702,
from searches of about five minutes that were still improving). Where the searches had
time to settle, the extra cost is the size of the layer itself plus 0 to 1 points, so the
sizes double from one dimension to the next: 211 = 2 x 105 + 1, 422 = 2 x 211,
844 = 2 x 422. In ratio to 2^n this is 844 / 512, about 1.648, against 5/3 for the pattern
and 2 for the classical construction. This is an observation from stalled searches, not a
theorem; if a proof of k_(n+1) <= 2 k_n existed it would give k_n <= 1.649 x 2^n for all
n >= 9, which is still base 2 and says nothing about large n, where the random construction
of Kopparty, Lev, Saraf and Sudan grows more slowly. The observation is also false in
F_3^3: a layered Kakeya set of F_3^3 needs at least 14 points (`k8/CENSUS.md`, proved by
hand and checked exhaustively), while k_3 = 13.

The free layers do not copy the fixed one in any simple way: for the 422-point set the
union of its free layers (96 + 115 points) is not a Kakeya set of F_3^7, the two free
layers meet the 211-point layer in only 7 and 9 points, and the best translate of the 211-point set covers 111 of its
211 points (`k8/probe9.py`, logged in `k6/results/records.json`, round 6).

Lower bounds: the general bound of Bukh and Chao, 3^n / (5/3)^(n-1), gives 57 (n = 6),
103 (n = 7), 184 (n = 8), 331 (n = 9) and 596 (n = 10). Nothing here is a lower bound,
and nothing says 105, 211, 422, 844 or 1702 is optimal.

## How the sets were found

`k6/fastwalk.py` (WalkSAT-style local search with incremental per-point costs) ran five
workers per dimension for about seven minutes each, in three rounds; rounds 4 to 7 used
`k8/layer.py`, which keeps a record set fixed in the layer x_0 = 0 and searches only the
two other layers (full log in `k6/results/records.json`):

| round | n | start | best | first reached at |
|---|---|-------|------|------------------|
| 1 | 6 | greedy from empty | 124 | 306 s |
| 1 | 6 | 53-point set frozen in the layer x_0 = 0 | 107 | 3 s |
| 1 | 6 | F_3 x (53-point set), 159 points | 106 | 415 s |
| 1 | 6 | classical set, 127 points | 106 | 46 s |
| 1 | 6 | (27-point set) x (7-point plane set), 189 points | 106 | 164 s |
| 2 | 6 | the 106-point set, two seeds | 106 | - |
| 2 | 6 | (13-point set) x (13-point set), 169 points | 123 | 366 s |
| 2 | 6 | classical set, 127 points | 106 | 25 s |
| 2 | 6 | F_3 x (53-point set), 159 points, noise 0.15 | **105** | 150 s |
| 3 | 6 | the 105-point set, two seeds | 105 | - |
| 3 | 6 | F_3 x (53-point set), three seeds | 106, 106, 106 | - |
| 1 | 7 | greedy from empty | 280 | 34 s |
| 1 | 7 | classical set, 255 points | 213 | 307 s |
| 1 | 7 | (53-point set) x (7-point plane set), 371 points | 214 | 268 s |
| 1 | 7 | (27-point set) x (13-point set), 351 points | 267 | 239 s |
| 2 | 7 | the 213-point set, two seeds | 213 | - |
| 2 | 7 | F_3 x (106-point set), 318 points, two seeds | 213, 216 | 313 s |
| 2 | 7 | classical set, 255 points | 215 | 316 s |
| 3 | 7 | 105-point set frozen in the layer x_0 = 0, then freed | **211** | 94 s |
| 3 | 7 | the same, another seed | 214 | - |
| 3 | 7 | 106-point set frozen in the layer x_0 = 0, then freed | 213 | - |
| 3 | 7 | F_3 x (105-point set), 315 points | 214 | 195 s |
| 3 | 7 | the 213-point set | 213 | - |
| 4 | 7 | 105-point layer, k8/layer.py | 211 | 156 s |
| 4 | 8 | 211-point layer, k8/layer.py, noise 0.10 | 427 | - |
| 4 | 8 | 211-point layer, k8/layer.py, noise 0.20 | 601 | 339 s |
| 5 | 8 | 211-point layer, from the 427-point set, noise 0.05 and 0.10 | 427, 427 | - |
| 5 | 8 | 211-point layer, fresh, noise 0.10 (seed 13) | 438 | 250 s |
| 5 | 8 | 211-point layer, fresh, noise 0.07 (seed 14) | 431 | 364 s |
| 5 | 8 | 211-point layer, fresh, noise 0.10 (seed 15) | **422** | 212 s |
| 5 | 9 | 427-point layer, fresh / classical start | 1321 / 1281 | - |
| 6 | 8 | 211-point layer, fresh, seeds 16, 17, 18 | 435, 423, 423 | - |
| 6 | 8 | 211-point layer, from the 422-point set | 422 | - |
| 6 | 9 | 422-point layer, lift start (1055 points), two seeds | 1055, 1054 | - |
| 6 | 9 | 422-point layer, start F_3 x (422-point set), 1266 points | 846 | 160 s |
| 7 | 9 | 422-point layer, start F_3 x K8, noise 0.10 / 0.13 | 846 / 844 | - / 413 s |
| 7 | 9 | 422-point layer, start K8 + t, K8 - t, noise 0.07 (seed 35) | **844** | 253 s |
| 7 | 9 | 422-point layer, start K8 + t, K8 - t, noise 0.10 | 845 | 253 s |
| 7 | 10 | 846-point layer, start F_3 x K9 (2538 points), noise 0.07 (seed 42) | **1702** | 224 s |
| 7 | 10 | 846-point layer, start F_3 x K9, noise 0.10 | 1756 | 277 s |
| 7 | 10 | 846-point layer, start K9 + t, K9 - t, noise 0.10 / 0.13 | 1747 / 1744 | 291 / 296 s |

The lift start puts F_3 x A and F_3 x B in the two free layers, where A and B are the free
layers of the 422-point set; it is a Kakeya set by construction, but the search could not
shrink it. The product start F_3 x K (also Kakeya: for every v the middle point of a line
of K in direction v works) and the shifted start K + t, K - t (for v put u = v - t; the
middle point a + u of a line {a, a + u, a + 2u} of K works) both shrink to about twice the
size of K.

The 105-point set has layers of 37, 31 and 37 points in x_0, 437 full lines, 305 of
its 364 directions covered by exactly one line, no point that can be deleted, and
|K cap -K| = 10. The 211-point set has 1316 full lines, 913 of its 1093 directions
covered by exactly one line, no point that can be deleted, and |K cap -K| = 12. The
427-point set has 4096 full lines, 2608 of its 3280 directions covered by exactly one
line, and no point that can be deleted. The 422-point set has 3888 full lines, 2854 of
its 3280 directions covered by exactly one line, and no point that can be deleted. The
846-point set has 11689 full lines, 8503 of its 9841 directions covered by exactly one
line, and no point that can be deleted.

## Next

`k8/probe11.py` asks whether the 844-point set has a deletable point, runs more layered
searches in F_3^9 on the 422-point layer, and gives the layered searches in F_3^10 more
time on top of the new 844-point layer (doubling would give 1688).
