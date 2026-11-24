# gate26: commit-time exhaustive checks of the 26-point cuts in F_3^4

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Code Apache-2.0, text CC BY 4.0.
Python standard library only.

## What this folder is

k_4, the smallest size of a Kakeya set in F_3^4, is 25, 26 or 27 (`moments/SD.md` proves 25; the
27-point set has been known since Lev's 2009 computer search, and Lev conjectured 27 is optimal).
Adding a point to a Kakeya set keeps it Kakeya, so **if no 26-point Kakeya set exists, k_4 = 27**.

`moments/PRODUCT.md` (checked by `moments/test_product.py`) shows that every 26-point Kakeya set
has a hyperplane class whose three sections have one of seven size patterns:

    (6,6,14) (5,7,14) (6,7,13) (7,7,12)   and   (5,8,13) (6,8,12) (7,8,11)

`kcut/kcut.py` decides each pattern exactly, without any solver: the smallest section A runs over
one representative of every AGL(3,3)-orbit, the middle section B over every set containing 0 with
B - A = F_3^3, and an exact branch-and-bound finds the third section C or proves none exists
(`kcut/README.md` explains the reduction).

## How a result is recorded here

MV2's publishing pipeline runs `python -I -m unittest discover` on a commit's files and makes the
commit **only if every test passes within ten minutes**. Each 26-point pattern is committed
separately, with:

- `test_gate_controls.py`: the group order, the k_3 = 13 end-to-end control in F_3^3 (nothing at 12
  points, a verified set at 13), recovery of the known 27-point set from two of its sections, and
  a brute-force recount of kcut's enumeration of B for every 5-point A;
- `test_cut_a_b_c.py`: the exhaustive search for that pattern, which fails (and prints any set it
  found) unless the search finds nothing.

So **a file `test_cut_a_b_c.py` exists in this folder only if the search for the pattern (a, b, c)
finished inside the pipeline with no Kakeya set found and all controls passed.** A pattern whose
file is missing is undecided here (the search found something, failed, or ran out of time). The
same patterns are also queued as separate solver-based and solver-free jobs (`kprof26`, `kdrat26`,
`kcut`), which give independent checks.

## Why a separate folder

The pipeline's test runner (Python 3.11 `unittest discover`) only enters folders that are Python
packages. `kcut/`, `kprof25/`, `kprof26/`, `kdrat26/` and `l82/` have no `__init__.py`, so the tests
in them were **not** run when they were committed; they run inside the compute jobs instead. This
folder is a package, and `engine.py` loads `kcut/kcut.py` by its file path, so these checks do run
at commit time. The checks in `moments/` (a package) have always run at commit time.

## What this does not show

Even if all seven files appear, the claim k_4 = 27 rests on the moment lemmas in `moments/` and on
the correctness of `kcut/kcut.py`. It would need the independent solver checks and a human
mathematician's reading before anyone calls Lev's conjecture proved.
