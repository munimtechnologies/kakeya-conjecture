# The pencil LP does not reach k_4 = 27 (negative result)

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's
name; not checked by a human mathematician. Text CC BY 4.0, code Apache-2.0.

## Question

The class-level counting of [QUAD45.md](QUAD45.md) and [RESIDUE.md](RESIDUE.md) cannot prove
k_4 = 27, k_5 >= 46 or k_6 >= 79: RESIDUE.md gives integer points satisfying every
class-level constraint at 26, 45 and 78 points. [PENCIL.md](PENCIL.md) proposed the next level:
a linear program over *pencils* (the 4 hyperplane classes through a codimension-2 subspace,
which cut a set into a 3 x 3 array of cells), with the exact cubic pencil identity, a Kakeya
filter per pencil and pair-incidence rows. Before using it on F_3^5 we tested it where the
answer is known.

## What we ran

* `pencil_gate.py` builds the LP of `pencil.py` for every quadric type of Q_K, first at type
  level, then for each u = 0..H (u = number of classes with an unbalanced profile), and solves
  it with the standard-library simplex `lpsolve.py`. A branch counts as excluded only when
  integer Farkas multipliers pass the exact check `pencil.check_certificate`.
* Calibration: at 12 points in F_3^3 every branch is excluded with checked certificates
  (`test_pencil_n3.py`), which reproves Lev's k_3 >= 13 without a search.
* Test case: 26 points in F_3^4, all 9 quadric types (`test_pencil26_results.py`).

## Result

None of the 9 types is excluded ([results/pencil26.json](results/pencil26.json)). Every
type-level LP is feasible, and between 5 and 22 values of u per type stay feasible. Feasibility
is the floating-point verdict of the solver; only exclusions are certified, so the precise
statement is "not excluded".

So the pencil LP (rule sd2, 17 profiles, 1902 pencil features after pruning) is too weak to
reprove k_4 = 27. Two exhaustive searches ([PROOF26.md](../PROOF26.md)) are still the only
proof, and we do not expect the same LP to give k_5 >= 46. Any further lower-bound progress
needs constraints that see more than one codimension-2 subspace at a time (for example
projections to F_3^3, 27 cells), or a different idea.

## Reproduce

`python -m unittest k5.test_pencil26_results` (about 80 s, Python 3.11 standard library).
