# k8: layered search for small Kakeya sets in F_3^n

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-26. Made with AI agents working in MV2's name;
not checked by a human mathematician. Code Apache-2.0, text CC BY 4.0.

Write F_3^n = F_3 x F_3^m (m = n - 1) and fix a Kakeya set L0 of F_3^m in the layer
x_0 = 0. Then {0} x L0 u {1} x L1 u {2} x L2 is a Kakeya set of F_3^n iff for every
v in F_3^m some x in L0 has x + v in L1 and x + 2v in L2 (`k8/layer.py`, docstring).

Why: the 211-point Kakeya set of F_3^7 (`k6/sets/kakeya7_211.json`) is of this form, with L0
the 105-point set of F_3^6 and |L1| + |L2| = 106 = 2 x 53, and a frozen 53-point layer in
F_3^6 led to 53 + 54 = 107. If completing a record layer always costs at most 2 k_(n-2)
points, then k_n <= k_(n-1) + 2 k_(n-2), and F_3^8 would have a Kakeya set of
211 + 210 = 421 points (Lev's 2009 pattern: 427; classical construction: 511).

* `layer.py`: incremental WalkSAT-style search over L1 u L2 for a fixed L0.
* `check.py`: stand-alone Kakeya check from the definition, used on every reported set.
* `test_layer.py`: tables against tuple arithmetic, incremental counters against a recount,
  the layer criterion against `check.py`, small searches in F_3^3 and F_3^4.
* `test_oracle_layer.py`: a report-only run (always fails on purpose, nothing committed):
  L0 = the 211-point set (F_3^8), the 105-point set (F_3^7), and the published 53- and
  27-point sets of Chung, Du and Wesley (arXiv:2608.23691; downloaded at run time from
  github.com/dualverse-ai/station_data_v2, not copied here) for F_3^6 and F_3^5.
