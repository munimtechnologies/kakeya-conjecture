# The 53-point Kakeya set in F_3^5 is locally optimal (no 3-swap improvement)

*MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name;
not checked by a human mathematician. Text CC BY 4.0, code Apache-2.0.*

The smallest Kakeya set known in F_3^5 has 53 points. It was found by Chung, Du and Wesley
([arXiv:2608.23691](https://arxiv.org/abs/2608.23691)) and is published as
`kakeya_F3_d5_53.npy` in [dualverse-ai/station_data_v2](https://github.com/dualverse-ai/station_data_v2)
(`artifacts/finite_kakeya`). We download it at run time; it is not copied into this repository.
V. F. Lev's 2009 guess k_n = k_(n-1) + 2 k_(n-2) predicts exactly 53. The proven range is
45 <= k_5 <= 53 (lower bound: [QUAD45.md](QUAD45.md)).

## What we checked (k5/cdw.py)

* The set is Kakeya. It contains 131 lines: 111 directions have one line and 10 have two. No
  point can be removed, and the set has no centre of symmetry.
* Its 121 hyperplane profiles (sizes of the three parallel sections) range from (17,18,18) to
  (11,16,26). It has no section split of the shape (27,13,13) that a layer-by-layer version
  of Lev's recurrence would suggest.
* **Exact neighbourhood.** For r = 1, 2, 3 we tried every way to delete r of the 53 points and
  add at most r - 1 new ones (53, 1,378 and 23,426 deletions; for each deletion every way to
  restore the lost directions is enumerated). No Kakeya set came out. For r = 4 we covered
  36,750 of the 292,825 deletions (12.5%) before the time limit, also with no hit.
* **Controls.** With one or two extra points added, the same search recovers the 53-point set
  the expected number of times, and around the published 27-point set in F_3^4 it finds no
  26-point neighbour, as k_4 = 27 requires. The search is also checked against brute force
  in F_3^3 (`test_cdw_unit.py`).
* **Local search.** Five WalkSAT-style walkers (`walk.py`) started at the 53-point set with a
  52-point budget for 445 seconds each (about 160,000 steps each) never left size 53.

Numbers: [results/local53.json](results/local53.json). `test_local53.py` re-downloads the set,
re-checks its invariants and re-runs the complete r <= 2 neighbourhood.

## What it means

A 52-point Kakeya set, if one exists, is not a small modification of the known one. This is
evidence, not a proof: k_5 = 53 remains open.
