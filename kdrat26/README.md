# kdrat26: an independent SAT + DRAT re-check of the 26-point search

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. No human
mathematician has checked this. Code licence: Apache-2.0; text CC BY 4.0.

`kprof26/` reduces Lev's conjecture k_4 = 27 (the smallest Kakeya set in F_3^4 has 27 points)
to a finite list of cases: P1-P4 (some hyperplane cut has a positive weight) and Z1-Z5 (every
cut has an 8-point section and there is exactly one line in each of the 40 directions; Z3-Z5 are
impossible by arithmetic). There, each case is decided by CP-SAT and re-solved by HiGHS.

This folder decides the same cases a third way, so that a "no" answer comes with a certificate
that a small independent program can check:

1. Each case and AGL(3,3)-orbit representative of slice 0 (the same model spec as kprof26,
   built by `kprof25.make_spec`) is written as a CNF formula: exact totalizer counters for the
   slice sizes and for the sections of all 40 hyperplane classes, one clause forbidding each
   pair of section sizes that is not in the allowed-cut table, a variable for each line (true
   exactly when its three points are), at least one line per direction, and in Z cases at most
   one. The moment identities are not encoded; they hold for every point set, so leaving them
   out keeps the formula a correct model.
2. CaDiCaL solves the formula and writes a DRAT proof when the answer is UNSAT.
3. drat-trim checks every proof. A representative counts as certified only when drat-trim
   prints `s VERIFIED`.
4. Any satisfying assignment is turned back into a point set and re-checked by brute force (the
   Kakeya property, all 40 cuts, and the case conditions).

Sanity checks in every run: the known 27-point set (fixed slices) is satisfiable and is the
set returned; forbidding its own cut makes it unsatisfiable; for each of the 27 sets "known
set minus one point" the formula agrees with brute force in both the all-cuts and the Z model;
a pigeonhole formula (8 into 7) is refuted with a verified proof; and a free search for 27
points finds a set that passes the brute-force check. The unit tests evaluate the encoding by
unit propagation against brute force, without any solver.

Verdicts: `CERT-K4-EQ-27` (every representative of P1-P4, Z1, Z2 is UNSAT with a verified
proof), `CERT-PART-CLOSED`, `UNSAT-PROOF-INCOMPLETE`, `UNDECIDED`, `FOUND-26` (a 26-point
Kakeya set, which would refute Lev's conjecture) or `BUG`.

Even a full `CERT-K4-EQ-27` certifies only the search. The reduction to the cases (the moment
identities, the sum-difference lemma and the arithmetic remark) is the hand argument in
`kprof26/README.md` and `moments/`, and the orbit representatives come from `kprof25`.
