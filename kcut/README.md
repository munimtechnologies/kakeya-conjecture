# kcut: a solver-free exhaustive search for 25- and 26-point Kakeya sets in F_3^4

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Code Apache-2.0, text CC BY 4.0.
Python standard library only (`kcut.py`, `kcutc.py`, tests in `test_kcut.py`, `test_kcutc.py`).

## Why

k_4, the smallest size of a Kakeya set in F_3^4, is 25, 26 or 27 (`moments/SD.md` proves 25; a
27-point set has been known since 2009, and Lev conjectured 27 is optimal). The exact searches
`kprof25`, `kprof26` and `kdrat26` use CP-SAT and SAT solvers. This folder is a third,
independent route with **no solver at all**: plain enumeration plus an exact set-cover search,
small enough to read and re-implement.

Adding a point to a Kakeya set keeps it Kakeya, so **the 26-point search alone decides the
question**: if no 26-point Kakeya set exists, k_4 = 27 (`moments/PRODUCT.md`, section 4).

## The reduction

Cut F_3^4 = F_3^3 x F_3 by the class {s = 0, 1, 2}, with sections A, B, C of F_3^3. K is Kakeya iff

1. for every d in F_3^3 some a in A has a + d in B and a + 2d in C (the 27 crossing directions), and
2. each of the 13 directions of F_3^3 has a line inside A, B or C.

Affine maps of F_3^4 let us put the smallest section at s = 0 as a fixed AGL(3,3)-orbit
representative, the middle one at s = 1 containing the origin (a shear (x, s) -> (x - s v, s)
translates B and fixes A), and the largest at s = 2. For each such A and every B, the program
searches for C exactly. C must meet each set W[d] = {-(a + b) : b - a = d}; the search branches on
the smallest unmet W[d], excludes earlier choices in later branches, prunes with a disjoint-packing
bound, and then completes the parallel directions line by line.

Which cuts must be searched comes from the moment identities (`moments/BAD25.md`,
`moments/min_section.py`, and the factored form in `moments/PRODUCT.md`):

- **N = 25:** some class has cut (6,6,13), (5,7,13) or (6,7,12).
- **N = 26:** some class has cut (6,6,14), (5,7,14), (6,7,13) or (7,7,12). Otherwise every class
  has an 8-point section and some hyperplane has at most 7 points, which leaves a class with cut
  (5,8,13), (6,8,12) or (7,8,11). `moments/test_product.py` re-derives this list from (T) and (SD)
  over all cuts of 26 points.

## kcutc: the same search with a symmetry reduction (added 2026-11-24)

With A fixed, every affine map g(x) = M x + t fixing A gives a bijection
(x, s) -> (M x + t + s v, s) of F_3^4 that fixes A, sends B to M B + (t + v) and C to
M C + (t + 2v). Since v is free, B and every M (B - b) (M a linear part of the stabiliser of A,
b in B) have the same answer. `kcutc.py` solves for C only when B is the smallest of these images,
which keeps exactly one B per orbit and cuts the number of C-searches by about |B| times the
stabiliser order. Its tests compare it with the unreduced search and with the k_3 = 13 control, and
check the symmetry claim directly. Job `kcut26c-v1` runs it on the seven N = 26 cuts.

## Controls (in the jobs and the tests)

- In F_3^3, the same engine searches every cut of 12 and of 13 points. It must find nothing at 12
  and a set at 13, because k_3 = 13.
- In F_3^2 sections, `solve` and `run` are compared with brute force over all C (and all B) on
  random instances. The B count must match exactly.
- The known 27-point set: its cut (6,8,13) must be found by the full enumeration, and no 26-point
  subset of it is Kakeya.
- The cut (6,6,12) must give nothing, since k_4 >= 25 is already proved.
- Every reported set is re-checked against the plain definition of a Kakeya set.
- kcutc: orbit-stabiliser check |Stab(A)| * |orbit| = |AGL(3,3)| for every representative.

## Note on testing (2026-11-24)

MV2's publishing pipeline runs `unittest discover` (Python 3.11) at commit time, and that only
enters folders that are Python packages. This folder has no `__init__.py`, so `test_kcut.py` and
`test_kcutc.py` were **not** run when they were committed; they run inside the compute jobs. The
folder `gate26/` is a package that loads `kcut.py` by file path and runs the controls and the
seven 26-point cuts at commit time, one commit per cut (see `gate26/README.md`).

## Status

Jobs `kcut-v1` and `kcut26c-v1` are queued (on 2026-11-24 every compute host was short of free
memory, so no job was running). Their results, and those of the `gate26/` commits, will be added
here; no claim is made until then.
