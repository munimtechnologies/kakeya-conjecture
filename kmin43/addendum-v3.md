# kmin43 addendum v3 (2026-11-24): decide 25 first, then 26

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
No human mathematician has checked this. Text CC BY 4.0, code Apache-2.0.

The exact search of `prereg-kmin43.md` (v2, `addendum-v2.md`) waited three days in the compute
queue because it asked for 2 GB and 4 cores. Since then k_4 >= 25 has been proved
(`moments/SD.md`). Version 3 keeps the reduction and the model unchanged and makes four changes.

1. **Order.** Bound 25 is decided first (no solution for any case proves k_4 >= 26), then
   bound 26 (no solution proves k_4 = 27, Lev's conjecture, and would need an independent
   check before any claim).
2. **Smaller smallest section.** For any N points of F_3^4,
   sum over the 120 hyperplanes of |K n H|^2 = 40 N + 13 N (N - 1) (`moments/min_section.py`,
   tested). If every parallel class had all three sections >= m, the sum would be at most
   40 (2 m^2 + (N - 2m)^2). Hence a set of 25 points has a section of at most 6 points and a set
   of 26 points one of at most 7. The search runs n0 = 1..6 for bound 25 and n0 = 1..7 for
   bound 26 (v2 ran up to 8).
3. **Symmetry breaking** |K_1| <= |K_2|: the linear map (y, s) -> (y, -s) fixes slice 0 and swaps
   the other two slices, so this loses no solution. The known 27-point set is still found by the
   free search with and without this constraint (checked in the job).
4. **Solvers.** CP-SAT first on every case, then HiGHS on the same cases as an independent check;
   2 workers and 1 GB, so the job can start sooner.

Verdicts: `K4-GE-26` (both solvers, bound 25), `K4-GE-26-CPSAT-ONLY`, `K4-EQ-27` (both solvers,
bound 26), `K4-EQ-27-CPSAT-ONLY`, `FOUND-<size>` (a verified set, re-checked by brute force),
`UNDECIDED`, `BUG`. A verified set of 24 or fewer points would contradict `moments/SD.md` and is
reported as a bug. Solver INFEASIBLE answers are not certificates; a DRAT-certified run follows
any positive verdict before an article claims it as a theorem.
