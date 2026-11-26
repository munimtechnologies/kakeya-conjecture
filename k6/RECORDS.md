# Kakeya sets in F_3^6 with 105 points, in F_3^7 with 211, in F_3^8 with 422 and in F_3^9 with 846

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-26. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

A set K in F_3^n is a Kakeya set if it contains a full line {a, a+d, a+2d} in every
direction d. Let k_n be the least size of a Kakeya set in F_3^n.

* `k6/sets/kakeya6_105.json` is a Kakeya set in F_3^6 with **105 points**, so **k_6 <= 105**.
* `k6/sets/kakeya7_211.json` is a Kakeya set in F_3^7 with **211 points**, so **k_7 <= 211**.
* `k6/sets/kakeya8_422.json` is a Kakeya set in F_3^8 with **422 points**, so **k_8 <= 422**
  (classical construction: 511; Lev's pattern value: 427).
* `k6/sets/kakeya9_846.json` is a Kakeya set in F_3^9 with **846 points**, so **k_9 <= 846**
  (classical construction: 1023; Lev's pattern value: 853).
* Earlier sets, kept and still checked: `kakeya6_106.json` (106 points),
  `kakeya7_213.json` (213 points), both found on 2026-11-25, and `kakeya8_427.json`
  (427 points, found on 2026-11-26).

All sets are checked at every commit by `k6/test_sets.py`, a stand-alone checker that
shares no code with the search: it decodes each integer i = sum_k x_k 3^(n-1-k) to the
point (x_0, ..., x_(n-1)) and checks all 364 (n = 6), 1093 (n = 7), 3280 (n = 8) or
9841 (n = 9) directions. It also checks that no single point of the 105- or 211-point set
can be deleted.

The 422- and 846-point files are stored in two parts. The layer x_0 = 0 of the 422-point
set is the file `kakeya7_211.json`, and its other 211 points are listed as `free_points` p,
standing for the integer 2187 + p (the file also keeps the compact code the search printed,
`points_code`, and the checker verifies that it decodes to the same list). The layer
x_0 = 0 of the 846-point set is the file `kakeya8_422.json`, and its other 424 points are
stored only as the search's compact code (`points_code`: sorted values, gaps in base 32,
see the file), standing for 6561 + p; both checkers decode it with their own decoders.

The 211-, 422-, 427- and 846-point sets have a second, independent proof in
`k6/test_layers.py`. The layer x_0 = 0 of the 211-point set is exactly the 105-point set
of F_3^6, which covers every direction (0, d). A direction (1, v), v in F_3^6, is covered
iff some x in that layer has x + v in the layer x_0 = 1 and x + 2v in the layer x_0 = 2;
the test checks all 729 vectors v. The other two layers have 56 and 50 points, together
**106 = 2 x 53**, and no point of either can be removed. In the same way the layer x_0 = 0
of the 422-point set is exactly the 211-point set of F_3^7 (all 2187 vectors v checked;
the other two layers have 96 and 115 points, together **211**), and the layer x_0 = 0 of
the 846-point set is exactly the 422-point set of F_3^8 (all 6561 vectors v checked; the
other two layers have **424** points together). So the chain 105 -> 211 -> 422 -> 846 has
a proof that uses only the direct check of the 105-point set in F_3^6.

## Why it matters

The best constructions we found in the literature for these dimensions are the classical
set {0,1}^n u {0,2}^n with 2^(n+1) - 1 points: 127 for n = 6, 255 for n = 7, 511 for
n = 8 and 1023 for n = 9 (Polymath wiki, Kakeya problem page; V. F. Lev's 2009 comment on
T. Tao's blog). Products of smaller Kakeya sets are worse (3 x 53 = 159, 7 x 27 = 189,
13 x 13 = 169; in F_3^8, 27 x 27 = 729 and 3 x 211 = 633; in F_3^9, even with our own
records as factors, 3 x 422 = 1266, 7 x 211 = 1477, 13 x 105 = 1365 and 27 x 53 = 1431,
all above the classical 1023). Kopparty, Lev, Saraf and Sudan
([arXiv:1003.3736](https://arxiv.org/abs/1003.3736), Theorem 8) have a random construction
of size (3 / 2^(2/3))^(n + O(sqrt n)), growth base about 1.8899, but no explicit sets for
small n. Chung, Du and Wesley ([arXiv:2608.23691](https://arxiv.org/abs/2608.23691)) found 53
points in F_3^5 and treated dimensions up to 5 only.

In 2009 Lev ([comment 994](https://terrytao.wordpress.com/2009/03/04/dhj3-900-999-density-hales-jewett-type-numbers/comment-page-3/#comment-36694))
computed k_3 = 13 and k_4 <= 27 and wrote "I am very curious to know whether k_5 = 53:
notice the pattern in 3, 7, 13, 27, 53, ...". Chung, Du and Wesley state this as the guess
that the recurrence k_n = k_(n-1) + 2 k_(n-2) continues and say that whether it continues
remains open. The pattern continues 107, 213, 427, 853.

**Consequence.** k_6 <= 105 < 107, k_7 <= 211 < 213, k_8 <= 422 < 427 and
k_9 <= 846 < 853: our sets are below the pattern in all four dimensions, and the gap grows
(2, 2, 5, 7 points). Since k_2 = 7 and k_3 = 13, the recurrence cannot hold at all of
n = 4, 5 and 6 (otherwise k_6 = 107). With k_4 = 27 (Lev's set and our computer proof
`PROOF26.md`, not human-checked) it fails at n = 5 or at n = 6: either k_5 <= 52, or
k_5 = 53 and k_6 <= 105 = k_5 + 2 k_4 - 2. We do not know which.

**The layer observation.** Each record is built on the previous one: a record set of
F_3^(n-1) in one layer plus a search over the two other layers. The extra cost was
54 on the 53-point layer in F_3^6 (total 107), 106 on the 105-point layer in F_3^7 (211),
211 on the 211-point layer in F_3^8 (422; nine other runs of about seven minutes gave 212,
212, 216, 216, 216, 220, 224, 227 and 390) and 424 on the 422-point layer in F_3^9 (846).
So far the extra cost is the size of the layer itself plus 0 to 2 points, which makes the
sizes roughly double from one dimension to the next: 846 / 2^9 is about 1.65, against
5/3 for the pattern and 2 for the classical construction. This is an observation from
stalled searches, not a theorem, and it says nothing about large n, where the random
construction of Kopparty, Lev, Saraf and Sudan grows more slowly. It is also false in
F_3^3: a layered Kakeya set of F_3^3 needs at least 14 points (`k8/CENSUS.md`, proved by
hand and checked exhaustively), while k_3 = 13.

The free layers do not copy the fixed one in any simple way: for the 422-point set the
union of its free layers (96 + 115 points) is not a Kakeya set of F_3^7, the two free
layers meet the 211-point layer in only 7 and 9 points, and the best translate of the 211-point set covers 111 of its
211 points (`k8/probe9.py`, logged in `k6/results/records.json`, round 6).

Lower bounds: the general bound of Bukh and Chao, 3^n / (5/3)^(n-1), gives 57 (n = 6),
103 (n = 7), 184 (n = 8) and 331 (n = 9). Nothing here is a lower bound, and nothing says
105, 211, 422 or 846 is optimal.

## How the sets were found

`k6/fastwalk.py` (WalkSAT-style local search with incremental per-point costs) ran five
workers per dimension for about seven minutes each, in three rounds; rounds 4 to 6 used
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
| 6 | 9 | 422-point layer, start F_3 x (422-point set), 1266 points | **846** | 160 s |

The lift start puts F_3 x A and F_3 x B in the two free layers, where A and B are the free
layers of the 422-point set; it is a Kakeya set by construction, but the search could not
shrink it. The product start F_3 x K8 (also Kakeya: for every v the middle point of a line
of K8 in direction v works) shrank from 1266 to 846 points.

The 105-point set has layers of 37, 31 and 37 points in x_0, 437 full lines, 305 of
its 364 directions covered by exactly one line, no point that can be deleted, and
|K cap -K| = 10. The 211-point set has 1316 full lines, 913 of its 1093 directions
covered by exactly one line, no point that can be deleted, and |K cap -K| = 12. The
427-point set has 4096 full lines, 2608 of its 3280 directions covered by exactly one
line, and no point that can be deleted. The 422-point set has 3888 full lines, 2854 of
its 3280 directions covered by exactly one line, and no point that can be deleted.

## Next

`k8/probe10.py` asks whether the 846-point set has a deletable point, runs more layered
searches in F_3^9, and starts the first layered sets in F_3^10 on top of the 846-point
layer from F_3 x (846-point set) and its shifted variant K + t, K - t (classical 2047,
pattern 1707).
