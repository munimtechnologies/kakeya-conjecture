# Field-symmetric Kakeya sets: first results (ORACLE-SYM, 2026-11-25)

MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; not checked
by a human mathematician. Text CC BY 4.0.

`sym/oracle.py` identifies F_3^n with the field F_(3^n) and, for each group
G = H_r x Frobenius_f (H_r the multiplicative subgroup of order r, f | n) with at most 60
point orbits, finds the least G-invariant Kakeya set with the branch-and-bound of `sym/bb.py`
(re-checked by `sym/check.py`). It ran in the publish gate of run 20261125T214710Z
(2026-11-25 22:07-22:13 UTC, 392 s). Only the end of its report was kept by the gate log, so
the table below lists the per-dimension best and the n = 8 lines that survived.

## Least invariant Kakeya set found, per dimension

| n | least invariant size | group (r, f) | status | best known Kakeya set (any) |
|---|---|---|---|---|
| 3 | 13 | (13, 1) | exact minimum in that class | 13 (= k_3) |
| 4 | 27 | (1, 2) | best found (search hit the time limit) | 27 (= k_4) |
| 5 | 55 | (11, 5) | exact minimum in that class | 53 (Chung, Du and Wesley) |
| 6 | 112 | (7, 6) | exact minimum in that class | 105 (k6/sets) |
| 7 | 1093 | (1093, 1) | exact; only very large groups have <= 60 orbits | 211 (k6/sets) |
| 8 | 615 | (205, 1) | exact minimum in that class | 511 (classical) |

Selected n = 8 classes (exact minima unless marked "best"): r = 205, f = 1: 615;
r = 41, f = 8: 657; r = 80, f = 8: 720; r = 40, f = 4: 760 (best); r = 160, f = 1: 800 (best);
r = 410, f = 8: 1231; r = 820, f = 1: 1640.

## Reading

No symmetry class with at most 60 orbits reaches the best known sizes for n = 5 to 8: the
gaps are +2 (n = 5), +7 (n = 6) and at least +104 (n = 8, against the classical 511). The
record sets are not invariant under any of these large field symmetries. These are minima
*within a symmetry class* only. The compete job ksym-v1 (CP-SAT, up to 400 orbits) has
been waiting for memory on every machine since 2026-11-25 22:07 UTC.
