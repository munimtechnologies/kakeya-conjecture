# No Kakeya set in F_3^5 has 44 points (k_5 >= 45), by a quadric constraint on sections

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

A Kakeya set in F_3^5 is a set of points containing a full line in each of the 121 directions.
Let k_5 be the least size of such a set. Known before: k_5 <= 53 (Chung, Du and Wesley,
arXiv:2608.23691) and, from our own moment and sum-difference method, k_5 >= 44
(`moments/SD.md`, `k5/SD2.md`). **Theorem.** No Kakeya set in F_3^5 has 44 points. Since adding
points keeps a set Kakeya, k_5 >= 45.

The proof below is by hand, apart from finite checks of the list of allowed section profiles,
which `k5/test_sd2.py` and `k5/test_quad.py` do exactly. It does **not** depend on the
exhaustive search `k5cut.py` (still queued); that search is now an independent second route to
the same statement.

## The new ingredient: section sizes mod 3 come from a quadratic form

For a hyperplane class xi (a nonzero linear form on F_3^5, up to sign) the three parallel
hyperplanes xi = 0, 1, 2 cut K into sections of sizes (n1, n2, n3), the *profile* of xi. Let
s2(xi) = C(n1,2) + C(n2,2) + C(n3,2): it counts the pairs {x, y} of K with xi(x - y) = 0. Over
F_3, the indicator of t = 0 is 1 - t^2. Hence

    s2(xi) = C(N, 2) - Q_K(xi)   (mod 3),    Q_K(xi) = sum over pairs {x, y} of K of xi(x - y)^2,

and **Q_K is a quadratic form in xi**. Quadratic forms over F_3 are classified by rank and
discriminant, so the numbers (z0, z1, z2) of classes where Q_K takes the values 0, 1, 2 form one
of eleven vectors for 5 variables:

| form | (z0, z1, z2) |
|---|---|
| 0 | (121, 0, 0) |
| rank 1 | (40, 81, 0), (40, 0, 81) |
| rank 2 | (67, 27, 27) hyperbolic, (13, 54, 54) elliptic |
| rank 3 | (40, 54, 27), (40, 27, 54) |
| rank 4 | (49, 36, 36) hyperbolic, (31, 45, 45) elliptic |
| rank 5 | (40, 45, 36), (40, 36, 45) |

`test_quad.py` recomputes this table from representatives, checks the classification
exhaustively for 3 and 4 variables and on random forms in 5, and checks the congruence on random
sets and on Lev's 27-point set.

## Proof of the theorem

Let K be a Kakeya set with N = 44 points and L >= 121 lines.

**Step 1 (two identities).** Summing over all 121 classes,
sum s2 = 40 C(44,2) = 37840 and sum s3 = 13 C(44,3) + 27 L, where s3 = sum C(n_i, 3) (a pair lies
in a common section for 40 classes; a non-collinear triple for 13, a line for 40).

**Step 2 (allowed profiles).** Every profile satisfies (T) (any two sections multiply to at least
81) and SD2 (`k5/SD2.md`). Put Delta2 = s2 - 301, Delta3 = s3 - 1274 (the values for the profile
(14, 15, 15), called *base*) and f = Delta3 - 15 Delta2. Then

    sum over classes of Delta2 = 1419,    sum over classes of f = 27 (L - 121) >= 0.

Among allowed profiles (exact enumeration in `test_sd2.py` and `test_quad.py`): the base has
Delta2 = 0 and f = 0; the profile (10, 11, 23), called *heavy*, has Delta2 = 52 and f = 2; every
other profile (a *defect*) has Delta2 >= 1 and f <= -1. Writing t_i = n_i - 16, one checks
f = 1 - 2 Delta2 + t1 t2 t3 / 2, and from this **f = 2 (mod 3) if Delta2 = 1 (mod 3), and
f = 0 (mod 3) otherwise**. So a defect with Delta2 not = 1 (mod 3) has f <= -3.

**Step 3 (budget).** Let k be the number of heavy classes and d the number of defects. From
sum f >= 0: the defects have total |f| <= 2k. From sum Delta2 = 1419: 52 k <= 1419, so k <= 27,
and the defects have total Delta2 = 1419 - 52 k. So d <= 2k and d <= 1419 - 52k.

**Step 4 (the quadric).** Since 946 - 301 = 645 is divisible by 3, Delta2(xi) = -Q_K(xi) (mod 3).
So base classes have Q_K = 0, heavy classes Q_K = 2, and every class with Q_K = 1 is a defect with
Delta2 = 2 (mod 3), hence f <= -3. With Step 3: 3 z1 <= 2k <= 54, so **z1 <= 18**. By the table
only two forms qualify:

- Q_K = 0, (121, 0, 0): no class can be heavy (heavy needs Q_K = 2), so k = 0, then d <= 2k = 0,
  every class is base and sum Delta2 = 0, not 1419. Impossible.
- Q_K = 2 l^2, (40, 0, 81): base classes lie among the 40 zeros, so d >= 121 - k - 40 = 81 - k.
  With d <= 2k this needs k >= 27, so k = 27 and d >= 54; but d <= 1419 - 52 * 27 = 15.
  Impossible.

So no 44-point Kakeya set exists. QED.

## The same argument as a linear program (code: `quad.py`)

For a quadric type z, a profile p can only sit on a class with Q_K = (C(N,2) - s2(p)) mod 3.
The profile counts x_p >= 0 must fill exactly z_c classes of each residue c, meet the pair
identity exactly and the triple identity with L >= 121. Weak duality gives, for every real mu,

    sum x_p s3(p) <= mu * S2 + sum_c z_c * max over residue-c profiles of (s3 - mu s2),

so a mu that pushes this below the needed triple total excludes the type. `quad.excluded(5, 44)`
finds such a mu for each of the eleven types, and `check_certificate` re-checks each one directly
over all profiles. The proof above is the hand version (mu = 15 for most types). As controls, the
same code must **not** exclude 13 points in F_3^3 or 27 points in F_3^4, where Lev's Kakeya sets
exist, and Lev's 27-point set must be a feasible point for its own quadric type.

The job `k5quad-v1` runs the same LP for F_3^3 to F_3^7 and reports whether it gives more
(k_5 >= 46, or new bounds for k_6 >= 78 and k_7 >= 137). Two extra tests,
`test_quad_k5_46.py` and `test_quad_k4_27.py`, are committed only if they pass: they assert that
the LP alone excludes 45 points in F_3^5, and 26 points in F_3^4 (which would be a search-free
proof of k_4 = 27).

## What is new and what is not

The moment identities are standard counting (they go back to the Polymath project's slice ideas
for Kakeya sets over F_3, 2009). The congruence s2 = C(N,2) - Q_K (mod 3) is an elementary
observation; we found no use of it for Kakeya sets in the literature we searched (2026-11-25).
What is new here is using it to restrict which hyperplane profiles can occur together, and the
bound k_5 >= 45 that follows.

## Limits

Carried out with AI agents working in MV2's name; no human mathematician has checked it. The
proof relies on the classification of quadratic forms over finite fields (standard; e.g. Lidl and
Niederreiter, *Finite Fields*, Chapter 6), on SD2 (`k5/SD2.md`, also unchecked by a human), and
on exact enumerations of section profiles done by the tests. The gap for k_5 is now 45 to 53.

## Credits

V. Lev found k_3 = 13 and the 27-point set in F_3^4 (2009, recorded on the Polymath wiki).
B. Bukh and T.-W. Chao proved the best general lower bound (arXiv:2108.00074).
S. Chung, W. Du and W. J. Wesley built the 53-point set in F_3^5 (arXiv:2608.23691).
