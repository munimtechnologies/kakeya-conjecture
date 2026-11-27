# Kakeya sets in F_3^6 with 105 points, in F_3^7 with 211, in F_3^8 with 422, in F_3^9 with 844, in F_3^10 with 1689, and below Lev's pattern in every dimension

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-26, updated 2026-11-27. Made with AI agents
working in MV2's name. **No human mathematician has checked this yet.** Text CC BY 4.0, code
Apache-2.0.

## Update 2026-11-27: the doubling

The 105-point set has a "split" of excess 1 (`k9/splits/kakeya6_105_split.json`), and by the
doubling and lifting lemmas of `k9/SPLIT.md` this gives a Kakeya set of **1689 points in
F_3^10** (2 x 844 + 1; the layered search below had reached 1693) and explicit Kakeya sets
of **845 x 2^(n-9) - 1 points in F_3^n for every n >= 9**: 3379 in F_3^11, 6759 in F_3^12,
13519 in F_3^13, ... So **every dimension n >= 6 now has an explicit Kakeya set below Lev's
pattern**, by 2, 2, 5, 9, 18, 34, 68, 134, ... points, a gap that grows without bound. The
sets in F_3^10, F_3^11 and F_3^12 are checked directly by computer (`k9/test_splits.py`,
`k9/test_doubled_big.py`); from F_3^13 on they rest on the hand-proved Lemma 1. Numbers:
`k9/results/doubling.json`. The layered-search records below are unchanged and still checked
at every commit; `k6/results/records.json` keeps the layered-search log (its entry for n = 10
is the layered record 1693).

## Result

A set K in F_3^n is a Kakeya set if it contains a full line {a, a+d, a+2d} in every
direction d. Let k_n be the least size of a Kakeya set in F_3^n.

* `k6/sets/kakeya6_105.json` is a Kakeya set in F_3^6 with **105 points**, so **k_6 <= 105**.
* `k6/sets/kakeya7_211.json` is a Kakeya set in F_3^7 with **211 points**, so **k_7 <= 211**.
* `k6/sets/kakeya8_422.json` is a Kakeya set in F_3^8 with **422 points**, so **k_8 <= 422**
  (classical construction: 511; Lev's pattern value: 427).
* `k6/sets/kakeya9_844.json` is a Kakeya set in F_3^9 with **844 points**, so **k_9 <= 844**
  (classical construction: 1023; Lev's pattern value: 853).
* The doubled 844-point set (`k9/`, rebuilt by `k9/test_splits.py` from the stored split) is a
  Kakeya set in F_3^10 with **1689 points**, so **k_10 <= 1689** (classical construction:
  2047; Lev's pattern value: 1707). The best layered-search set,
  `k6/sets/kakeya10_1693.json`, has 1693 points.
* k_11 <= 3379, k_12 <= 6759 (checked directly), and k_n <= 845 x 2^(n-9) - 1 for all n >= 9.
* Earlier sets, kept and still checked: `kakeya6_106.json` (106 points),
  `kakeya7_213.json` (213 points), both found on 2026-11-25, and `kakeya8_427.json`
  (427 points), `kakeya9_846.json` (846 points), `kakeya10_1702.json` (1702 points) and
  `kakeya10_1695.json` (1695 points), found on 2026-11-26.

All sets in `k6/sets/` are checked at every commit by `k6/test_sets.py`, a stand-alone
checker that shares no code with the search: it decodes each integer
i = sum_k x_k 3^(n-1-k) to the point (x_0, ..., x_(n-1)) and checks all 364 (n = 6),
1093 (n = 7), 3280 (n = 8) or 9841 (n = 9) directions one by one. For n = 10 it collects the
covered directions from pairs of points instead (three distinct points of F_3^n are
collinear iff they sum to 0) and checks that all 29524 directions appear; the two methods
are compared on small cases. It also checks that no single point of the 105- or 211-point
set can be deleted.

The 422-, 846-, 844-, 1702-, 1695- and 1693-point files are stored in two parts. The layer
x_0 = 0 of the 422-point set is the file `kakeya7_211.json`, and its other 211 points are
listed as `free_points` p, standing for the integer 2187 + p (the file also keeps the
compact code the search printed, `points_code`, and the checker verifies that it decodes to
the same list). The layer x_0 = 0 of the 846- and 844-point sets is the file
`kakeya8_422.json`, and their other points are stored only as the search's compact code
(`points_code`: sorted values, gaps in base 32, see the file), standing for 6561 + p. The
layer x_0 = 0 of the 1693- and 1695-point sets is the 844-point set (of the earlier
1702-point set, the 846-point set), and their other points stand for 19683 + p. Both
checkers decode the codes with their own decoders.

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
| 1702 in F_3^10 | 846-point set of F_3^9 | 443 + 413 = 856 | 19683 |
| 1695 in F_3^10 | 844-point set of F_3^9 | 398 + 453 = 851 | 19683 |
| 1693 in F_3^10 | 844-point set of F_3^9 | 430 + 419 = 849 | 19683 |
| **1689** in F_3^10 (doubling, `k9/`) | 844-point set of F_3^9 | A + B of its split, 845 | all 29524 directions by pairs |

So the chain 105 -> 211 -> 422 -> 844 -> 1693 has a proof that uses only the direct check of
the 105-point set in F_3^6. For the sets in F_3^10 the covered vectors are collected as
{y - x : x in L_0, y in L_1, -(x + y) in L_2}.

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
remains open. The pattern continues 107, 213, 427, 853, 1707, 3413, 6827, 13653, ...

**Consequence.** k_6 <= 105 < 107, k_7 <= 211 < 213, k_8 <= 422 < 427,
k_9 <= 844 < 853, k_10 <= 1689 < 1707, and k_n <= 845 x 2^(n-9) - 1 below the pattern for
every n >= 9. Since k_2 = 7 and k_3 = 13, the recurrence cannot hold at all of n = 4, 5
and 6 (otherwise k_6 = 107). With k_4 = 27 (Lev's set and our computer proof `PROOF26.md`,
not human-checked) it fails at n = 5 or at n = 6: either k_5 <= 52, or k_5 = 53 and
k_6 <= 105 = k_5 + 2 k_4 - 2. We do not know which.

**The layer observation.** Each record is built on the previous one: a record set of
F_3^(n-1) in one layer plus a search over the two other layers. The extra cost was
54 on the 53-point layer in F_3^6 (total 107), 106 on the 105-point layer in F_3^7 (211),
211 on the 211-point layer in F_3^8 (422), 422 on the 422-point layer in F_3^9 (844; the
other logged runs gave 422, 423, 423, 423, 424, 424) and 849 on the 844-point layer in
F_3^10 (1693, twice; the other runs gave 850, 850, 851, 852, 853, 855, 855, 857 and 906).
Where the searches had time to settle, the extra cost is the size of the layer itself plus
0 to 5 points. The doubling lemma of `k9/SPLIT.md` now explains the "+1": any layer with a
split of excess 1 can be extended by exactly |layer| + 1 points. The searches sometimes did
one better (422 = 2 x 211, 844 = 2 x 422). The observation is false in F_3^3: a layered
Kakeya set of F_3^3 needs at least 14 points (`k8/CENSUS.md`, proved by hand and checked
exhaustively), while k_3 = 13.

**Split structure (search round 9).** Most points of the free layers lie in one translate
pair of the fixed layer L_0: for the 211-point set, all 56 points of L_1 lie in L_0 + t and
48 of the 50 points of L_2 in L_0 - t (one t); for the 422-point set 89 of 96 and 108 of
115; for the 844-point set 188 of 207 and 202 of 215 (for the sets found from the product
start, t = 0). A layered search that may only use points of L_0 itself in the free layers
reached 211 on the 105-point layer (the record), 423 on the 211-point layer (record 422),
845 on the 422-point layer (record 844) and 1696 on the 844-point layer (record 1693). Full
numbers: `k6/results/records.json`, `split_structure`. This led to the doubling lemma of
`k9/SPLIT.md`.

Lower bounds: the general bound of Bukh and Chao, 3^n / (5/3)^(n-1), gives 57 (n = 6),
103 (n = 7), 184 (n = 8), 331 (n = 9), 596 (n = 10), 1072 (n = 11), 1929 (n = 12) and
3471 (n = 13). Nothing here is a lower bound, and nothing says any of these sizes is optimal.

## How the sets were found

`k6/fastwalk.py` (WalkSAT-style local search with incremental per-point costs) ran five
workers per dimension for about seven minutes each, in three rounds; rounds 4 to 9 used
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
| 7 | 10 | 846-point layer, start F_3 x K9 (2538 points), noise 0.07 (seed 42) | 1702 | 224 s |
| 7 | 10 | 846-point layer, start F_3 x K9, noise 0.10 | 1756 | 277 s |
| 7 | 10 | 846-point layer, start K9 + t, K9 - t, noise 0.10 / 0.13 | 1747 / 1744 | 291 / 296 s |
| 8 | 9 | 422-point layer, start K8 + t, K8 - t / F_3 x K8 | 845 / 845 | 53 / 310 s |
| 8 | 10 | 844-point layer, start F_3 x K9 (2532 points), noise 0.07 (seed 61) | 1695 | 309 s |
| 8 | 10 | 844-point layer, start F_3 x K9, noise 0.05 / 0.07 / 0.10 | 1696 / 1699 / 1750 | 317 / 300 / 327 s |
| 8 | 10 | 844-point layer, start K9 + t, K9 - t, noise 0.07 / 0.05 | 1697 / 1699 | 341 / 314 s |
| 9 | 10 | 844-point layer, from the 1695-point set's free layers, noise 0.05 / 0.07 | 1694 / 1694 | 28 / 179 s |
| 9 | 10 | 844-point layer, start F_3 x K9, noise 0.07 (seed 73) | 1693 | 325 s |
| 9 | 10 | 844-point layer, start F_3 x K9, noise 0.06 (seed 74) | 1693 | 341 s |
| 9 | 10 | 844-point layer, start K9 + t, K9 - t, noise 0.07 | 1701 | 278 s |

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
line, and no point that can be deleted. The 844-point set has 11552 full lines, 8647 of
its 9841 directions covered by exactly one line, and no point that can be deleted.

## Next

`k9/probe14.py` runs layered searches in F_3^10 on the 844-point layer started from the
doubled set's own free layers (845 points). Any layered set on that layer with fewer than
845 free points would, by Lemma 2, carry a split of excess 1 again and lower the whole family
for n >= 10.
