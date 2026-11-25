# The pencil LP for Kakeya sets in F_3^n

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.
Code: `k5/pencil.py` (exact, standard library), `k5/pencil_lp.py` (HiGHS via scipy, with exact
Farkas certificates); checks: `k5/test_pencil.py`.

## Why

`RESIDUE.md` showed that counting over single hyperplane classes (moment identities, the pair
condition (T), SD2 and all mod-3 information) has integer solutions at 26 points in F_3^4 and
45 points in F_3^5. A proof of k_4 = 27 without an exhaustive search, or of k_5 >= 46, has to
use how classes fit together. The smallest such structure is the **pencil**.

## Set-up

Let K be a set of N points in F_3^n, n >= 3, and V a subspace of codimension 2, with
annihilator spanned by linear forms xi1, xi2. The 9 cosets of V ("cells") carry the counts
m(a, b) = #{x in K : xi1.x = a, xi2.x = b}, a 3 x 3 array indexed by the plane F_3^2. The four
hyperplane classes containing V are xi1, xi2, xi1 + xi2, xi1 + 2 xi2; the sections of each are
the three parallel lines of the array in one of its four directions, so the four class profiles
are the line sums of the array. Conversely the array is determined by its line sums:
m(x) = (sum of the four line sums through x - N) / 3.

Constants: H = (3^n - 1)/2 classes, T = 3^(n-1), a1 = (3^(n-1) - 1)/2, a2 = (3^(n-2) - 1)/2,
P = [n choose 2]_3 pencils, B1 = [n-1 choose 2]_3, B2 = [n-2 choose 2]_3, D = 3^(n-2). For
n = 4: P = 130, a1 = 13, B1 = 13, B2 = 1, D = 9. For n = 5: P = 1210, a1 = 40, B1 = 130,
B2 = 13, D = 27.

## Facts (with proofs)

1. **Incidences.** A class (hyperplane W) lies in a1 pencils (the hyperplanes of W). Two
   distinct classes lie in exactly one common pencil (the intersection of their hyperplanes).
2. **Cell triples.** sum over pencils of sum over cells C(m, 3) = B2 C(N, 3) + (B1 - B2) L,
   with L the number of lines in K. Proof: three distinct points lie in one coset of V iff V
   contains the span of their differences, which is 1-dimensional for the L collinear triples
   (B1 subspaces V contain it) and 2-dimensional otherwise (B2 contain it).
3. **Pencil Kakeya filter.** Fix a class through V. The D directions inside its hyperplane but
   not inside V each need a line of K, and such a line lies in one section of the class,
   meeting its three cells. In a section whose cells have a, b, c points, lines in distinct
   directions use distinct pairs between any two cells, so at most min(ab, bc, ca) directions
   are served. Hence, summed over the 3 sections, sum of min(ab, bc, ca, D) >= D, for each of
   the 4 classes. The (D - 1)/2 directions inside V need lines inside single cells; a cell
   with m points is a copy of AG(n-2, 3) and serves at most max_dirs(m) directions (exact
   table for n = 4, checked by brute force; for n = 5 at most m(m-1)/6, and at most 12 when
   m < 13 since k_3 = 13).
4. **Line counts.** A line of K lies inside a cell for B1 pencils and crosses three cells (one
   section of one class through V) for P - B1 pencils; inside a section for a1 classes and
   crossing for T classes. A crossing line is fixed by its pair in any two of its cells, and a
   set of m points holds at most m(m-1)/6 lines (exact table in AG(2,3)).
5. **Monotonicity.** Adding points to a Kakeya set keeps it Kakeya, so excluding N excludes
   every smaller size.

## The LP

Variables (all >= 0): x_p, the number of classes with sorted profile p (p allowed by (T) and
SD2); y_f, the number of pencils whose array has feature f = (multiset of its four sorted line
sums, sum C(m,3), lines-inside bound, lines-crossing bound), over all arrays that pass the
filter of fact 3; L' = L - H. Rows: the two class moment identities; sum y = P;
sum_f (number of directions of f with profile p) y_f = a1 x_p for every p (fact 1); the cell
triple identity (fact 2); the four line-count inequalities (fact 4); the residue groups of a
quadric type z (`QUAD45.md`); and, for u = number of classes with an unbalanced profile,
sum_f C(k_U(f), 2) y_f = C(u, 2) (fact 1: every pair of classes shares one pencil). Loop over
all 2n + 1 types z and all u = 0..H. If every branch is infeasible, no Kakeya set of N points
exists.

Arrays are enumerated up to the affine group AGL(2, 3) (it permutes the 4 directions as S4 and
rows and columns freely), and dominated features (same profiles and cell triples, fewer
allowed lines) are dropped; both reductions keep infeasibility valid, and the first is checked
against brute force in `test_pencil.py`.

**Certificates.** For each infeasible branch the job extracts a Farkas vector (lam, mu >= 0)
and checks exactly, in rational arithmetic, that
b.lam + h.mu - sum over columns with r_j < 0 of r_j U_j < 0, where r = A^T lam + G^T mu and U_j
are valid column bounds (x_p <= H, y_f <= P, L' <= C(N,2)/3). This proves the branch has no
solution without trusting the floating-point solver.

## Checks (`test_pencil.py`)

- incidence counts of fact 1 for n = 3, 4 (every pair of classes in exactly one pencil);
- the AG(2,3) tables by brute force over all 512 subsets;
- enumeration up to symmetry equals brute force over all arrays (n = 3, N = 12, 13, 15);
- on Lev's 27-point set in F_3^4, on random minimal Kakeya sets in F_3^3 and F_3^4 and on a
  random Kakeya set in F_3^5: every pencil array passes the filters, its feature is in the
  enumeration, the identities hold, and the true vector (x, y, L') satisfies every LP row
  exactly in its own branch (z, u);
- the exact Farkas checker on toy systems (accepts a valid certificate, rejects bad ones, and
  never certifies a feasible system).

## Status

The LP runs as compute job `k5pencil-v1` on (n, N) = (3, 12), (3, 13), (4, 25), (4, 26),
(4, 27), (5, 45), then upward from 46 if 45 falls, with N = 13, 27 and 53 as controls that must
stay feasible (sets of those sizes exist). Results will be recorded here. A certified
infeasibility at (4, 26) would give a proof of k_4 = 27 that is independent of the two
exhaustive searches in `PROOF26.md` and checkable from a finite certificate; one at (5, 45)
would give k_5 >= 46.

## Credits

V. Lev (k_3 = 13 and the 27-point set in F_3^4, 2009, recorded on the Polymath wiki);
B. Bukh and T.-W. Chao (arXiv:2108.00074); S. Chung, W. Du and W. J. Wesley (53 points in F_3^5,
arXiv:2608.23691).
