# The product identity for hyperplane cuts in F_3^n, and why the 26-point case list is complete

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.
Checks: `moments/test_product.py` (standard library only).

## 1. The identity

Let K be a set of N points in F_3^n containing exactly L lines (a line is a 3-set {x, y, z} with
x + y + z = 0). Each of the H = (3^n - 1)/2 parallel classes of hyperplanes cuts K into three
sections of sizes n1, n2, n3. Then, summing over the H classes,

    sum n1 n2 n3                 = 3^(n-2) * (C(N,3) + 2 L)
    sum (n1 n2 + n1 n3 + n2 n3)  = 3^(n-1) * C(N,2)
    sum (C(n1,2) + C(n2,2) + C(n3,2)) = ((3^(n-1) - 1)/2) * C(N,2).

*Proof.* Count (class, triple of points in three different sections). A class with normal xi
separates three points x, y, z when xi.(y - x), xi.(z - x) are non-zero and different. If the triple
is a line, z - x = -(y - x), so this holds exactly when xi.(y - x) != 0: 2 * 3^(n-1) normal vectors,
3^(n-1) classes. Otherwise y - x and z - x are independent, the pair (xi.(y - x), xi.(z - x)) takes
each value in F_3^2 for 3^(n-2) vectors xi, and only (1, 2) and (2, 1) separate: 3^(n-2) classes.
So the first sum is 3^(n-2) (C(N,3) - L) + 3^(n-1) L. The second and third sums count pairs that are
separated (3^(n-1) classes) or not separated ((3^(n-1) - 1)/2 classes). QED.

This is elementary and is equivalent to the usual Fourier count of solutions of x + y + z = 0; we
record it because it turns the case analysis of the 26-point search into one line.

## 2. Shifted forms in F_3^4 (n = 4, H = 40)

Expanding (n1 - t)(n2 - t)(n3 - t) and using the identities:

    N = 26, t = 8:  sum over classes (n1 - 8)(n2 - 8)(n3 - 8) = 18 (L - 40)
    N = 25, t = 8:  sum over classes (n1 - 8)(n2 - 8)(n3 - 8) = 18 L - 580
    N = 27, t = 9:  sum over classes (n1 - 9)(n2 - 9)(n3 - 9) = 18 (L - 36)

The N = 26 line is the weight `delta26 = s3 - 8 s2 + 576` of `moments/BAD25.md` in factored form:
2 * delta26 = (n1 - 8)(n2 - 8)(n3 - 8) for every cut of 26 points (checked for all cuts).

## 3. The cases a 26-point Kakeya set must have

A Kakeya set in F_3^4 has L >= 40, so for N = 26 the product sum is >= 0. Every class of a Kakeya
set satisfies (T) (any two sections have product >= 27) and (SD) (each section is at least
729 / Emax of the other two; `moments/SD.md`). Hence either

* some class has a cut with (n1 - 8)(n2 - 8)(n3 - 8) > 0, or
* every product is 0 (and L = 40): every class has an 8-point section. Summed over the 40 classes,
  n1^2 + n2^2 + n3^2 = 9490, while classes with all sections >= 8 are (8,8,10) or (8,9,9), giving
  at most 40 * 228 = 9120. So some hyperplane has at most 7 points, and its class has a cut
  (a, 8, 18 - a) with a <= 7.

Running through all cuts of 26 points with (T) and (SD), the test finds exactly:

| case | cuts |
|---|---|
| positive product | (6,6,14), (5,7,14), (6,7,13), (7,7,12) |
| zero product, smallest section <= 7 | (5,8,13), (6,8,12), (7,8,11) |

(4,7,15) and (5,6,15) are removed by (SD) (third sections of at least 25 and 18), and (4,8,14) is
removed because Emax(4,8) = 50 and 729/50 > 14. These seven cuts are exactly `CUTS26` in
`kcut/kcut.py`, so the 26-point search covers every case.

## 4. Why 26 points decides everything

Adding any point to a Kakeya set gives a Kakeya set. So if no 26-point Kakeya set exists in F_3^4,
none with 25 or fewer points exists either, and k_4 = 27 (the 27-point set has been known since
2009). The 25-point search (`kcut`, `kprof25`) is therefore a stepping stone, not a separate
requirement. `kcut/kcutc.py` searches the seven cuts above, with an extra symmetry reduction; its
job is queued. **Nothing here proves k_4 = 27 yet**: that needs the searches to finish, agree, and be
checked independently.
