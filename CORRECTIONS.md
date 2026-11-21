# Corrections

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. Text CC BY 4.0.

## 2026-11-21 (second entry): K(3,3) = 13, the 27-point set in F_3^4 and the 56-point set in F_3^5 were not new

A literature re-check, promised in the entry below, found three over-claims in README versions 1.1-1.3
(sections 7 and 12 and the "New in version 1.1" line). We missed two sources:

1. **Vsevolod F. Lev, 2009.** In comment 994 on Terence Tao's Polymath thread "DHJ3: 900-999 (density
   Hales-Jewett type numbers)" (14 March 2009,
   https://terrytao.wordpress.com/2009/03/04/dhj3-900-999-density-hales-jewett-type-numbers/comment-page-3/#comment-36694),
   Lev, writing as "Seva", defined k_r as the smallest Kakeya set in F_3^r and wrote: "Using a computer, I also found
   k_3=13 and k_4≤27". He added "I suspect that, indeed, k_4=27 holds" and asked whether k_5 = 53, pointing to the
   pattern 3, 7, 13, 27, 53. The Polymath wiki page "Kakeya problem" recorded the same values.
2. **Chung, Du and Wesley, "Autonomous Mathematical Discovery in an Open-World Multi-Agent Environment",
   arXiv:2608.23691 (v1 24 Aug 2026, v2 14 Sep 2026).** Their Section 6.1 reports a 53-point Kakeya set in F_3^5,
   so k_5 <= 53. That matches Lev's pattern. They also cite Lev for k_3 = 13 and for k_4 <= 27 being believed sharp.

What this changes:

| README claim | Correct status |
|---|---|
| "Theorem (new). K(3,3) = 13", "the first exact finite-field Kakeya number in dimension three over a field with more than two elements" | Not first. Lev found k_3 = 13 by computer in 2009 (his comment gives no method details). Our exhaustive check of all 1,307,504 candidate 12-point sets, together with two solvers, is an independent, re-runnable **confirmation** of his value. |
| 27-point Kakeya set in F_3^4 listed as new | Lev had k_4 <= 27 in 2009. Ours is an independent example, not a new bound. |
| 56-point Kakeya set in F_3^5 listed as new | Superseded. Chung, Du and Wesley's 53 points is better. |
| 64-point Kakeya set in F_4^4 listed as new | Still not found anywhere else. The Kopparty-Saraf-Sudan even-q formula gives 88, and the Gold-function construction of Kyureghyan, Mueller and Wang (arXiv:1302.5591) is larger. We have not found a smaller published set, but we can no longer call our search exhaustive, so treat "new" here with caution. |

Consequences for the open computations:

- **kmin43-v2**, the queued exact computation of K(4,3) (the minimum in F_3^4, known range 18 to 27), is a test
  of Lev's 2009 conjecture k_4 = 27. If it proves that no Kakeya set of 26 points exists, it **confirms a
  conjecture open since 2009**. If it finds one, it refutes the conjecture. The pre-registered predictions and
  labels are unchanged; only the description of what the result would mean is corrected here.
- **K(3,4) = 28** (F_4^3, `results/kmin34-v1/`) is unaffected. Neither source, nor anything else we found,
  gives a lower bound above 21 for F_4^3. The 28-point upper bound is the Kopparty-Saraf-Sudan construction
  (entry below).

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
A full re-check of those entries against the literature is still to be done. (Done: see the entry above.)
