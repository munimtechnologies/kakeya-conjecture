# Addendum v2 to the kmin43 pre-registration (K(4,3): 26 or 27?)

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-21, written before kmin43-v2 runs.
Made with AI agents working in MV2's name. **No human mathematician has checked this note.**
Text CC BY 4.0; code Apache-2.0.

## What happened to kmin43-v1

kmin43-v1 ran on 2026-11-21 (16:43-16:49 UTC) and stopped in its unit tests, before any search.
That was a pre-registered BUG outcome, and nothing is claimed from it. Two tests failed, both
for the same reason: the helper `classical(n)`, which builds the classical Kakeya set used only as
a sanity check, returned 35 points in F_3^4 instead of 31.

**Cause.** The classical set is built by induction: K_{m+1} = {(t, d^2 + t d)} ∪ ({0} x K_m'),
where K_m' must be an affine copy of K_m that contains the slice {0,1}^m of the first part. v1
used K_m itself with no affine map, so the overlap was too small and the sizes went 3, 7, 16, 35
instead of 3, 7, 15, 31. v2 uses K_m' = phi(K_m) with phi(t, y) = (t - 1, 2y). This affine map
sends {1,2} x {0,2}^(m-1), which lies inside K_m, onto {0,1}^m. The docstring of `classical()`
gives the proof of |K_n| = 2^(n+1) - 1.

The other ten v1 tests passed: the group order, orbit counts, the slice/brute-force agreement, the
known 27-point set (accepted at bound 27, rejected when fixed at bound 26 by both solvers), and
infeasibility for a 1-point smallest section. **The reduction and the 0/1 model are unchanged.**

## Changes in v2 (all fixed before the run)

1. `classical()` fixed as above. New tests check sizes 3, 7, 15, 31 against a second, independent
   brute-force Kakeya test, and check the embedded cube that the induction step needs.
2. **Deadline.** The main search stops starting new cases after 260 minutes (job timeout 360,
   leaving room for host pauses). Cases it does not reach are reported as `NOT-RUN`, and the
   verdict is then UNDECIDED (never EQ-27).
3. **Sanity at bound 27** now runs *before* the main search. It goes from smallest-section size
   9 down to 1 and stops at the first set that CP-SAT finds and the brute-force test confirms.
   The pre-registered requirement ("CP-SAT finds at least one set at bound 27") is unchanged.
4. HiGHS/CP-SAT time limits in one unit test were raised (60 s to 120/180 s) so that a slow
   machine does not cause a false failure.

Verdicts, solver limits for the main search (CP-SAT 900 s, HiGHS 600 s per case) and the
predictions of the original pre-registration stand. One prediction is updated for v2 only:
P(BUG) 15% -> 8%, since the failing piece was found and the rest passed in v1.
