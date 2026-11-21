# Corrections

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. Text CC BY 4.0.

## 2026-11-21: the 28-point Kakeya set in F_4^3 was not new
README versions 1.1-1.3 (section 7 and section 12) list "28 (new)" for F_4^3 and say the
28-point set is smaller than every explicit construction quoted there. That is wrong for F_4^3.
For even q, the construction of Kopparty, Saraf and Sudan has
2^-(n-1) q^n + (1 - 2^-(n-1)) q^(n-1) points (quoted in Bukh and Chao, "Sharp density bounds on
the finite field Kakeya problem", arXiv:2108.00074), which is 16 + 12 = 28 for n = 3, q = 4.
So 28 points in F_4^3 was already known; the README should have cited it.

What is new is the other direction: the computation kmin34-v1 (`results/kmin34-v1/`) found no
Kakeya set of 27 points, so 28 is the minimum by that computation (HiGHS in all 22 cases, CP-SAT
in 16; an independent DRAT-certified re-check is queued as kmin34cert-v1).

This note does not re-examine the other "new" entries (27 points in F_3^4, 56 in F_3^5, 64 in
F_4^4); for F_4^4 the formula above gives 88, so the 64-point set is below that construction.
A full re-check of those entries against the literature is still to be done.
