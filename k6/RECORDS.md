# Kakeya sets in F_3^6 with 105 points and in F_3^7 with 211 points

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-26. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

A set K in F_3^n is a Kakeya set if it contains a full line {a, a+d, a+2d} in every
direction d. Let k_n be the least size of a Kakeya set in F_3^n.

* `k6/sets/kakeya6_105.json` is a Kakeya set in F_3^6 with **105 points**, so **k_6 <= 105**.
* `k6/sets/kakeya7_211.json` is a Kakeya set in F_3^7 with **211 points**, so **k_7 <= 211**.
* Earlier sets, kept and still checked: `kakeya6_106.json` (106 points) and
  `kakeya7_213.json` (213 points), both found on 2026-11-25.

All sets are checked at every commit by `k6/test_sets.py`, a stand-alone checker that
shares no code with the search: it decodes each integer i = sum_k x_k 3^(n-1-k) to the
point (x_0, ..., x_(n-1)) and checks all 364 (n = 6) or 1093 (n = 7) directions. It also
checks that no single point of a record set can be deleted.

The 211-point set has a second, independent proof in `k6/test_layers.py`. Its layer
x_0 = 0 is exactly the 105-point set of F_3^6, which covers every direction (0, d). A
direction (1, v), v in F_3^6, is covered iff some x in that layer has x + v in the layer
x_0 = 1 and x + 2v in the layer x_0 = 2; the test checks all 729 vectors v. The other two
layers have 56 and 50 points, together **106 = 2 x 53**, and no point of either can be
removed.

## Why it matters

The best constructions we found in the literature for these dimensions are the classical
set {0,1}^n u {0,2}^n with 2^(n+1) - 1 points: 127 for n = 6 and 255 for n = 7
(Polymath wiki, Kakeya problem page; V. F. Lev's 2009 comment on T. Tao's blog). Products
of smaller Kakeya sets are worse (3 x 53 = 159, 7 x 27 = 189, 13 x 13 = 169). Kopparty,
Lev, Saraf and Sudan ([arXiv:1003.3736](https://arxiv.org/abs/1003.3736), Theorem 8) have
a random construction of size (3 / 2^(2/3))^(n + O(sqrt n)), growth base about 1.8899, but
no explicit sets for small n. Chung, Du and Wesley
([arXiv:2608.23691](https://arxiv.org/abs/2608.23691)) found 53 points in F_3^5 and
treated dimensions up to 5 only.

In 2009 Lev ([comment 994](https://terrytao.wordpress.com/2009/03/04/dhj3-900-999-density-hales-jewett-type-numbers/comment-page-3/#comment-36694))
computed k_3 = 13 and k_4 <= 27 and wrote "I am very curious to know whether k_5 = 53:
notice the pattern in 3, 7, 13, 27, 53, ...". Chung, Du and Wesley state this as the guess
that the recurrence k_n = k_(n-1) + 2 k_(n-2) continues and say that whether it continues
remains open. The pattern continues 107, 213, 427.

**Consequence.** k_6 <= 105 < 107 and k_7 <= 211 < 213. Since k_2 = 7 and k_3 = 13, the
recurrence cannot hold at all of n = 4, 5 and 6 (otherwise k_6 = 107). With k_4 = 27 (Lev's
set and our computer proof `PROOF26.md`, not human-checked) it fails at n = 5 or at n = 6:
either k_5 <= 52, or k_5 = 53 and k_6 <= 105 = k_5 + 2 k_4 - 2. We do not know which.

**The layer observation.** 211 = 105 + 2 x 53 is exactly what the recurrence gives when it
is applied to the best values known now (k_6 <= 105, k_5 <= 53), and the set is built that
way: a record set of F_3^6 in one layer plus 2 x 53 points. Searches that froze a record
F_3^5 set (53 points) in one layer of F_3^6 stopped at 107 = 53 + 2 x 27. So, in every case
we tried, completing a record layer cost about 2 k_(n-2) points. If this held in general it
would give k_n <= k_(n-1) + 2 k_(n-2) as an inequality, and in F_3^8 a set of
211 + 2 x 105 = 421 points, below the pattern's 427 and the classical 511. This is an
observation from stalled searches, not a theorem; `k8/layer.py` tests it directly.

Lower bounds: the general bound of Bukh and Chao, 3^n / (5/3)^(n-1), gives 57 (n = 6) and
103 (n = 7). Nothing here is a lower bound, and nothing says 105 or 211 is optimal.

## How the sets were found

`k6/fastwalk.py` (WalkSAT-style local search with incremental per-point costs) ran five
workers per dimension for about seven minutes each, in three rounds (full log in
`k6/results/records.json`):

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

The 105-point set has layers of 37, 31 and 37 points in x_0, 437 full lines, 305 of
its 364 directions covered by exactly one line, no point that can be deleted, and
|K cap -K| = 10. The 211-point set has 1316 full lines, 913 of its 1093 directions
covered by exactly one line, no point that can be deleted, and |K cap -K| = 12.

## Next

`k8/layer.py` searches only the two free layers for a fixed record layer, a problem with
3^(n-1) conditions and |L_0| candidate lines each instead of the full space. It asks
whether about 2 x 105 = 210 extra points are enough around the 211-point layer in F_3^8
(total 421), and whether fewer than 2 k_(n-2) are enough in F_3^5, F_3^6 and F_3^7.
