# kprof25: is there a 25-point Kakeya set in F_3^4?

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
No human mathematician has checked this. Code Apache-2.0, text CC BY 4.0.

Known: 25 <= k_4 <= 27 (`moments/SD.md`). `moments/BAD25.md` proves that every 25-point Kakeya set
has at least two hyperplane classes whose cut is (6,6,13), (5,7,13) or (6,7,12), because these are
the only allowed cuts with positive weight in the exact identity

    sum over the 40 classes of (7 s3 - 51 s2 + 3302) = 63 L - 2420,   L >= 40.

So, after an affine map, a 25-point Kakeya set would be in one of three positions:

| case | cut of the last class | excluded everywhere | slice 0 |
|---|---|---|---|
| A | (6,6,13) | nothing | canonical 6-point set |
| B | (5,7,13) | (6,6,13) | canonical 5-point set |
| C | (6,7,12) | (6,6,13), (5,7,13) | canonical 6-point set |

`kprof25.py` decides every case and every AGL(3,3)-orbit of slice 0 with CP-SAT. The model has:

- one Boolean per point of slices 1 and 2;
- one Boolean per line of F_3^4 (the AND of its points), with a line in every one of the 40
  directions;
- a table of allowed cuts (with their s2 and s3) in each of the 40 classes;
- the two moment identities with the exact line count.

The HiGHS MILP solver then re-solves every case that CP-SAT finds infeasible, using a different
0/1 encoding of the table. Before the search, the model must accept the known 27-point set, reject
it once one of its own cuts is forbidden, and re-find a 27-point set by free search. Any set it
finds is re-checked by brute force.

Verdicts: `K4-GE-26` (both solvers: no 25-point Kakeya set, so k_4 >= 26), `K4-GE-26-CPSAT-ONLY`,
`FOUND-25` (a verified 25-point Kakeya set), `UNDECIDED` (time ran out) or `BUG`.

Even a K4-GE-26 verdict is a computer search. It is not a checked proof until a DRAT certificate
or an independent re-implementation confirms it.

Results go to `results/kprof25-v1/` once the job has run. The job (`kprof25-v1`, 1 GB, 2 CPUs)
runs `python -m unittest -v test_kprof25` first.
