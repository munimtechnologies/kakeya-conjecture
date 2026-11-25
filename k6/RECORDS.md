# Kakeya sets in F_3^6 with 106 points and in F_3^7 with 213 points

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

A set K in F_3^n is a Kakeya set if it contains a full line {a, a+d, a+2d} in every
direction d. Let k_n be the least size of a Kakeya set in F_3^n.

* `k6/sets/kakeya6_106.json` is a Kakeya set in F_3^6 with **106 points**, so **k_6 <= 106**.
* `k6/sets/kakeya7_213.json` is a Kakeya set in F_3^7 with **213 points**, so **k_7 <= 213**.

Both sets are checked at every commit by `k6/test_sets.py`, a stand-alone checker that
shares no code with the search: it decodes each integer i = sum_k x_k 3^(n-1-k) to the
point (x_0, ..., x_(n-1)) and checks all 364 (n = 6) or 1093 (n = 7) directions.

## Why it matters

The best constructions we found in the literature for these dimensions are the classical
set {0,1}^n u {0,2}^n with 2^(n+1) - 1 points: 127 for n = 6 and 255 for n = 7
(Polymath wiki, Kakeya problem page; V. F. Lev's 2009 comment on T. Tao's blog). Products
of smaller Kakeya sets are worse (3 x 53 = 159, 7 x 27 = 189, 13 x 13 = 169). Chung, Du and
Wesley ([arXiv:2608.23691](https://arxiv.org/abs/2608.23691)) found 53 points in F_3^5 and
treated dimensions up to 5 only.

In 2009 Lev ([comment 994](https://terrytao.wordpress.com/2009/03/04/dhj3-900-999-density-hales-jewett-type-numbers/comment-page-3/#comment-36694))
computed k_3 = 13 and k_4 <= 27 and wrote "I am very curious to know whether k_5 = 53:
notice the pattern in 3, 7, 13, 27, 53, ...". Chung, Du and Wesley state this as the guess
that the recurrence k_n = k_(n-1) + 2 k_(n-2) continues and say that whether it continues
remains open. The pattern continues 107, 213.

**Consequence.** k_6 <= 106 < 107. Since k_2 = 7 and k_3 = 13, the recurrence cannot hold
at all of n = 4, 5 and 6 (otherwise k_6 = 107). With k_4 = 27 (Lev's set and our computer
proof `PROOF26.md`, not human-checked) it fails at n = 5 or at n = 6: either
k_5 <= 52, or k_5 = 53 and k_6 <= 106 < k_5 + 2 k_4. We do not know which.

In F_3^7 the 213-point set only matches the pattern's value 213. Since k_6 <= 106, the
recurrence applied to the true values would now ask for k_7 <= k_6 + 2 k_5 <= 212.

Lower bounds: the general bound of Bukh and Chao, 3^n / (5/3)^(n-1), gives 57 (n = 6) and
103 (n = 7). Nothing here is a lower bound, and nothing says 106 or 213 is optimal.

## How the sets were found

`k6/fastwalk.py` (WalkSAT-style local search with incremental per-point costs) ran five
workers per dimension for about seven minutes each (log in `k6/results/records.json`):

| n | start | best | first reached at |
|---|-------|------|------------------|
| 6 | greedy from empty | 124 | 306 s |
| 6 | 53-point set frozen in the layer x_0 = 0 | 107 | 3 s |
| 6 | F_3 x (53-point set), 159 points | **106** | 415 s |
| 6 | classical set, 127 points | **106** | 46 s |
| 6 | (27-point set) x (7-point plane set), 189 points | **106** | 164 s |
| 7 | greedy from empty | 280 | 34 s |
| 7 | classical set, 255 points | **213** | 307 s |
| 7 | (53-point set) x (7-point plane set), 371 points | 214 | 268 s |
| 7 | (27-point set) x (13-point set), 351 points | 267 | 239 s |

Three different structured starts in F_3^6 all stopped at exactly 106. The second worker
kept the 53-point set fixed in the layer x_0 = 0 and reached 107 = 53 + 2 x 27 points
within seconds, but never fewer. The 106-point set is not of that shape: its layers
x_0 = 0, 1, 2 have 35, 34 and 37 points.

## Next

`k6/probe2.py` continues from both record sets with new seeds and other structured
starts, and reports structure numbers of the records (lines, directions covered once,
deletable points, hyperplane sections).
