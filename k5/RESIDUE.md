# Section sizes mod 3, and where the profile method stops (26 in F_3^4, 45 in F_3^5, 78 in F_3^6)

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.
Code: `k5/residue.py`; checks: `k5/test_residue.py`.

## Summary

1. **Residue lemma.** For any N points in F_3^n, the residues mod 3 of the three section sizes of
   *every* hyperplane class are fixed by two objects: a linear form l (from the sum of the points)
   and a quadratic form R (from the sum of their squares). The quadratic form Q_K used in
   `QUAD45.md` (which proved k_5 >= 45) is Q_K = N R - l^2. So **when 3 divides N, Q_K has rank at
   most 1**, and the section residues split the classes into four groups, not three.
2. **Negative result.** The profile method (the two moment identities, the pair condition (T),
   the refined sum-difference bound SD2 and *all* mod-3 information above) has explicit integer
   solutions at **26 points in F_3^4, 45 points in F_3^5 and 78 points in F_3^6**. So that method,
   on its own, cannot prove k_4 = 27, k_5 >= 46 or k_6 >= 79. This explains why the two
   pass-only tests filed on 2026-11-25 (`test_quad_k5_46.py`, `test_quad_k4_27.py`) did not land,
   and it shows that the exhaustive search behind k_4 = 27 (`PROOF26.md`) was really needed.
   The next gains have to use geometry beyond single hyperplane classes.

## 1. The residue lemma

Let K be a multiset of N points in F_3^n and xi a nonzero linear form. Put
n_j = #{x in K : xi.x = j}, l(xi) = sum over x of xi.x and R(xi) = sum over x of (xi.x)^2, all
read in F_3. Over F_3 the indicators of t = 0, 1, 2 are 1 - t^2, 2t^2 + 2t and 2t^2 + t. Summing
over K:

    n_0 = N - R,    n_1 = -R - l,    n_2 = -R + l      (mod 3).

l is linear in xi and R is a quadratic form in xi. Summing (xi.(x - y))^2 over ordered pairs gives
2N R - 2 l^2, so Q_K(xi) = sum over unordered pairs of (xi.(x - y))^2 = N R - l^2.

**Consequences.**

- If 3 does not divide N, translate K so that its sum is 0 (possible since N is invertible mod 3).
  Then l = 0 and every class has sizes {N - r, -r, -r} mod 3, where r = R(xi). These three
  residue patterns are exactly the three groups that `quad.py` already uses, so for these N the
  LP in `QUAD45.md` already has all mod-3 information.
- If 3 divides N, then Q_K = -l^2 has rank 0 or 1: only the types (H, 0, 0) and
  ((3^(n-1) - 1)/2, 0, 3^(n-1)) occur, not all 2n + 1. Classes with l(xi) = 0 have all three sizes
  congruent to -R(xi); classes with l(xi) != 0 have three different residues. On the hyperplane
  l = 0 of the dual space, R is translation invariant, so its type there (a quadratic form in
  n - 1 variables) is a further constraint.

For a hypothetical 45-point Kakeya set in F_3^5 this says: either every hyperplane section has
size congruent mod 3 to the other two sections of its class, or there is one special direction
s such that this holds for the 40 classes through s and fails (three different residues) for the
other 81. Any future search over 45-point sets can use this to prune cuts.

`test_residue.py` checks the lemma on every subset of F_3^2, on random sets in F_3^3 to F_3^5
(including N divisible by 3), and on Lev's 27-point Kakeya set; and it checks the rank-1
statement on random sets with N = 9, 12, 27, 45.

## 2. Integer points the profile method cannot exclude

Notation: a class has profile (n_1, n_2, n_3); H = (3^n - 1)/2 classes; T = 3^(n-1);
a1 = (3^(n-1) - 1)/2; a2 = (3^(n-2) - 1)/2. The method asks for profile counts with

- total H classes, each profile allowed by (T) and SD2 (`k5/SD2.md`);
- sum of s2 = a1 C(N, 2) and sum of s3 = a2 C(N, 3) + (a1 - a2) L for an integer L >= H
  (L = number of lines in K);
- residues: some linear l and quadratic R whose predicted residues match the profiles class by
  class (here checked as equality of the two multisets of residue patterns, which is equivalent).

| space | N | profile counts | l, R | L |
|---|---|---|---|---|
| F_3^4 | 26 | 10 x (6,6,14), 1 x (6,8,12), 1 x (6,9,11), 28 x (8,9,9) | 0, 0 | 53 |
| F_3^5 | 45 | 23 x (9,12,24), 1 x (12,12,21), 1 x (12,15,18), 96 x (15,15,15) | 0, 0 | 125 |
| F_3^6 | 78 | 39 x (18,18,42), 1 x (18,21,39), 3 x (17,19,42), 234 x (25,26,27), 3 x (24,26,28), 3 x (24,25,29), 81 x (26,26,26) | xi_1, xi_2^2 | 364 |

`residue.check_point` verifies every condition exactly, and the test also checks that the LP in
`quad.py` does not exclude the quadric type of each point (weak duality says it cannot).

**Reading the table.**

- F_3^4, 26 points: the profile LP is satisfied by ten classes cut as (6,6,14). So k_4 = 27 cannot
  come from counting profiles alone. The two exhaustive searches of `PROOF26.md` are needed, and
  they exclude these patterns by geometry.
- F_3^5, 45 points: with every section size divisible by 3, about 17 to 23 classes cut as
  (9,12,24) satisfy all counting. A short calculation (the hand version of the LP for this group)
  shows any solution of this residue type needs between 17 and 23 such classes. To prove
  k_5 >= 46 one has to show that a Kakeya set cannot have that many sections of this shape,
  or rule out the other residue types in the same way.
- F_3^6, 78 points: the old bound k_6 >= 78 (moment LP with SD2) does not improve with the residue
  information either.

## 3. What would go further

The natural next level is the **pencil**: the four hyperplane classes containing a fixed
codimension-2 subspace V. They cut K into 9 cells (the cosets of V), and the 4 profiles are the
line sums of a 3 x 3 array. Over all pencils there is an exact cubic identity
sum over pencils of sum over cells C(m, 3) = a C(N, 3) + b L (for n = 5: 13 C(N, 3) + 117 L),
which the class-level LP does not see, and each pencil carries its own Kakeya condition: for each
of its 4 classes, the 3^(n-2) directions in that class but not in V need lines inside the
class's sections, each running through three cells of one row of the array. A pencil-level LP
(variables = numbers of pencils of each cell type) is the planned next step.

## What is new and what is not

The indicator identities over F_3 are elementary, and summing polynomial identities over a set
is standard (it is the degree-2 part of the polynomial method). We found no prior use of the
section-residue structure for Kakeya sets (searches on 2026-11-25). New here: the explicit
form Q_K = N R - l^2, the rank-1 consequence for N divisible by 3, and the three explicit
counting solutions that mark where single-class counting stops. The solutions are not Kakeya
sets. They only show that a proof needs more than these constraints.

## Credits

V. Lev (k_3 = 13 and the 27-point set in F_3^4, 2009, recorded on the Polymath wiki);
B. Bukh and T.-W. Chao (arXiv:2108.00074); S. Chung, W. Du and W. J. Wesley (53 points in F_3^5,
arXiv:2608.23691).
