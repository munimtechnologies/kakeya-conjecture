# kmin34-v1 result: no Kakeya set of 27 points in F_4^3 (HiGHS, all 22 cases)

MV2 of Munim, Inc. (mv2@munimtech.com), graded 2026-11-21.
Made with AI agents working in MV2's name. **No human mathematician has checked this.**
Text CC BY 4.0. Pre-registration: `kmin34/prereg-kmin34.md`. Raw output: `kmin34.json` (this
folder), per-case table `cases.csv`, summary `summary.json`.

## What came out
- The job ran 7,074 s on 2 cores with no bug flags (AGL(2,4) has 2,880 distinct permutations,
  orbit sizes add up to C(16, k) for k = 1..7, the known 28-point set is accepted by both
  solvers at 28 and rejected at 27, its embedded copy equals the published certificate).
- At bound 27 the reduction leaves 22 cases (one per AGL(2,4)-orbit of a smallest plane
  section, sizes 1 to 6). **HiGHS proved all 22 infeasible.** CP-SAT proved 16 of them; it
  timed out (900 s) on the case n0 = 4, R = {0, 1, 4, 5} and on all five cases with n0 = 5.
  No case was open for both solvers and no set of 27 or fewer points was found.
- So, by this computation, **K(3, 4) = 28**: the smallest Kakeya set in F_4^3 has 28 points.

## Grading against the pre-registration
The pre-registered labels were K34-EQ-28 (both solvers prove every case, 45%), a smaller set
(20%), K34-EQ-28-CPSAT-ONLY (10%), UNDECIDED ("some case open for both solvers", 10%) and BUG
(15%). The outcome, every case proved by HiGHS and 16 of 22 by CP-SAT, matches none of them
exactly: the job's code falls back to the label **UNDECIDED**, while by the table's own wording
no case is open for both solvers. We report it as what it is: the both-solver bar of K34-EQ-28
was **not met**; the HiGHS-only mirror image of K34-EQ-28-CPSAT-ONLY was. The prediction "a
smaller set exists" (20%) failed.

## What it means
- The upper bound 28 is attained by the even-q construction of Kopparty, Saraf and Sudan, whose
  size 2^-(n-1) q^n + (1 - 2^-(n-1)) q^(n-1) is 16 + 12 = 28 for n = 3, q = 4 (as quoted by
  Bukh and Chao, arXiv:2108.00074). If the computation is right, that construction is exactly
  optimal in F_4^3. The Bukh-Chao density bound gives 21 here.
- As far as we could find (searches on 2026-11-21; Lund, Saraf and Wolf, arXiv:1609.01048,
  note that exact values were known only in the plane), this is the first exact value of a
  Kakeya number K(3, q) with q > 3, after K(3, 3) = 13 in this repository.

## Correction to earlier text in this repository
The README (v1.1-v1.3) calls the 28-point set in F_4^3 "new" and smaller than every explicit
construction it quotes. That was wrong for F_4^3: the Kopparty-Saraf-Sudan construction already
has 28 points there (our 28-point set matches it in size; we have not checked whether the two
are equivalent). See `CORRECTIONS.md`. The other entries (F_3^4: 27, F_3^5: 56, F_4^4: 64) are
not re-examined here; for F_4^4 the same formula gives 32 + 56 = 88.

## Limits and next step
The 6 cases proved only by HiGHS rest on one floating-point MIP solver, with no certificate.
An independent re-check is filed as `kmin34cert-v1` (code in `kmin34/cert/`): a new CNF
encoding written from scratch, solved by kissat with DRAT proofs, every proof checked by
drat-trim, with solver-free consistency checks against brute force. The value should be
treated as computer-assisted and unconfirmed until that job reports.
