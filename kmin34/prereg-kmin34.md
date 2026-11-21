# Pre-registration: kmin34-v1 (how small can a Kakeya set in F_4^3 be?)
MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-21, written before the job runs.
Made with AI agents working in MV2's name. **No human mathematician has checked this note.**
Text CC BY 4.0; code `kmin34.py` and tests `test_kmin34.py` in this folder, Apache-2.0.

## Question
A Kakeya set in F_q^n contains a full line in every direction; K(n, q) is the least size. Exact
values are classical in the plane, and this repository proved K(3, 3) = 13 (release v1.1.0).
For F_4^3 the known range is **21 <= K(3, 4) <= 28**: 21 is the Bukh-Chao bound recorded in
`data/kakeya-min-v1/results.json`, 28 is the certificate found there (CP-SAT stopped with
lower bound 24 after 1,200 s on the full 64-point model). The job asks whether 28 is optimal
and, if it is not, finds the true minimum. A sister job (`kmin43/`) asks the same for F_3^4.

## Method (exact, not a heuristic)
Write points as (y, s) with y in F_4^2, s in F_4, and K_s = {y : (y, s) in K}.
1. A transversal direction (w, 1) has a line in K iff some a in K_0 has a + t w in K_t for
   t = 1, 2, 3 (16 directions). A horizontal direction (w, 0) has a line iff some slice K_s
   contains a line of F_4^2 in direction w (5 directions).
2. Let |K| <= b and let n0 be the least size of K n H over the 84 affine planes H. Four
   parallel planes partition the space, so n0 <= b // 4 (6 for b = 27), and n0 >= 1.
3. An affine map of F_4^3 moves a smallest plane to {s = 0}; an affine map of F_4^2 applied to
   y alone (which is affine on F_4^3) moves K_0 to the canonical representative R of its
   AGL(2, 4)-orbit (2,880 elements; orbits enumerated by canonical augmentation with the
   check that orbit sizes add up to C(16, k)).
4. The image solves MODEL(R, b): 48 0/1 point variables for K_1, K_2, K_3, the transversal and
   horizontal conditions, all 84 plane sections >= |R|, and |K| <= b. If MODEL(R, 27) is
   infeasible for every representative R with |R| <= 6, then **K(3, 4) = 28**.
5. Each model is solved by OR-Tools CP-SAT (900 s) and HiGHS (300 s). Every solution is
   re-checked by a brute-force Kakeya test and by a solver-free row check. If a set of size
   m <= 27 is found, the search restarts at bound m - 1 (representatives already proved
   infeasible by both solvers stay infeasible, since the model only tightens).

## Self-checks (any failure means BUG)
- F_4 field axioms with the encoding of `kakeya.py` (2 = w, 2*2 = 3); 5 and 21 directions.
- AGL(2, 4) is 2,880 distinct permutations preserving the 20 lines; orbit sizes add to C(16, k).
- The slice test agrees with the brute-force test on random sets; a normalised Kakeya set
  passes the model's rows at its own size and fails one below; a non-Kakeya set fails.
- The known 28-point set (embedded copy checked by hand; the downloaded certificate is
  preferred) is accepted by both solvers at 28 and rejected at 27 with its slices fixed; the
  free model for its representative is feasible at 28; at bound 28 at least one
  representative is feasible.

## Outcomes and predictions (written before any result)
| verdict | meaning | probability |
|---|---|---|
| K34-EQ-28 | both solvers: no Kakeya set of <= 27 points; K(3, 4) = 28 | 45% |
| K34-EQ-m or K34-LE-m (m <= 27) | a smaller Kakeya set exists (new record), exact if the next pass closes | 20% |
| K34-EQ-28-CPSAT-ONLY | CP-SAT proves every case, HiGHS times out on some | 10% |
| UNDECIDED | some case open for both solvers within the limits | 10% |
| BUG | a self-check fails | 15% |

## What we will publish either way
`results/kmin34-v1/kmin34.json` with every case, solver times and certificates. A decided value
would, as far as we know, be the first exact value of K(3, q) for q = 4 (a 2016 paper,
arXiv:1609.01048, notes that exact minima are known only in the plane).
Limits: the infeasibility proofs come from solvers (no DRAT certificate yet); a solver-free
replication is the planned follow-up.
