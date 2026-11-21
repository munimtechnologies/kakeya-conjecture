# Pre-registration: kmin34cert-v1 (a certificate for K(3, 4) >= 28)
MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-21, written before the job runs.
Made with AI agents working in MV2's name. **No human mathematician has checked this note.**
Text CC BY 4.0; code `kmin34cert.py` and tests `test_kmin34cert.py` in this folder, Apache-2.0.

## Why
kmin34-v1 (`results/kmin34-v1/`) reduced "is there a Kakeya set of <= 27 points in F_4^3?" to
22 small 0/1 models. HiGHS proved all 22 infeasible; CP-SAT proved 16 and timed out on 6. The
claim K(3, 4) = 28 therefore rests on one floating-point MIP solver for 6 cases. This job
re-checks every case by a different route that leaves a checkable proof.

## Method
- Same reduction as `prereg-kmin34.md` (smallest plane section moved to slice 0, then to an
  AGL(2, 4)-orbit representative). The 22 representatives are copied from kmin34-v1, but the
  group is rebuilt by new code and the job re-checks that their orbits cover all C(16, k)
  subsets for k = 1..6.
- New CNF encoding written from scratch (nothing imported from kmin34.py): 48 point variables,
  AND-auxiliaries for the 16 transversal and 5 horizontal directions, Sinz sequential counters
  for "every plane section >= |R|" and "|K| <= 27".
- Each case: kissat (built from source at the commit recorded in the output) writes a DRAT
  proof; drat-trim checks it. Limits 1,800 s for kissat and 2,400 s for drat-trim per case,
  2 cases in parallel. The 22 CNF files are kept in the output for anyone to re-solve.

## Self-checks (any failure means BUG)
- F_4 axioms; AGL(2, 4) = 2,880 distinct line-preserving permutations; orbit coverage.
- Kakeya brute force: the known 28-point set and F_4^3 are Kakeya; 200 random 20-point sets
  are not (they cannot be, by the Bukh-Chao bound 21).
- Pipeline: pigeonhole 5 -> 4 is UNSAT with a VERIFIED proof; a satisfiable formula with an
  empty proof is not verified.
- Counters: exhaustive check against popcount for 5 literals (kissat, and solver-free in tests).
- The normalised known 28-point set is SAT at 28 and UNSAT (verified) at 27 with its slices
  fixed; the free formula for its representative at 28 is SAT with a brute-force-valid witness.
- 90 perturbed/random sets: kissat, a solver-free canonical assignment and unit propagation
  must all agree with brute force (Kakeya and size <= bound), with at least 3 of each answer.

## Outcomes and predictions (written before any result)
| verdict | meaning | probability |
|---|---|---|
| CERT-EQ-28 | all 22 cases UNSAT with drat-trim-VERIFIED proofs: K(3, 4) = 28 with a checkable certificate | 60% |
| PARTIAL | some cases time out or proofs are too big to check in the limits | 22% |
| FOUND-27 | kissat finds a brute-force-valid Kakeya set of <= 27 points (HiGHS was wrong) | 2% |
| BUG | a self-check fails (setup, tool output parsing or encoding) | 16% |

On CERT-EQ-28 we update the article "kakeya-sets-in-f4-cubed-need-28-points" and the results
with the certificate; on PARTIAL we file -v2 for the open cases only (`--only`), with longer
limits or symmetry breaking; on FOUND-27 we retract the K(3, 4) = 28 claim publicly.
