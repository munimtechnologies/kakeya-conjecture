# Collinear plane sections of Kakeya sets: a bound for every q (hand proof)

MV2 of Munim, Inc. Text CC BY 4.0. Written 2026-11-22 by AI agents working in MV2's name;
no human mathematician has checked it. It extends `collinear-section-lemma.md` (same folder)
from q = 4 to all q. Its use is to prune cases in exact Kakeya-number searches (F_4^3 now,
F_5^3 later). It is a small lemma, not a research result on its own.

## Setting (as in collinear-section-lemma.md)

K is a Kakeya set in F_q^n. H is a hyperplane, and B = K ∩ H has k points lying on one line
of H. With the collinear-section lemma,

    |K| >= k + q^(n-2) · m(q, k),

where m(q,k) is the least value of Σ_{t ≠ 0} |R_t|. Here R_t = ∪_b (β_b + t A_b) ⊂ F_q,
(A_b) is a partition of F_q into parts indexed by the k base points (empty parts allowed),
and β_b are the k distinct positions of the base points on their line.

## Lemma (every q)

    m(q, k) >= C(q,2) + s(q,k),

where C(q,2) = q(q-1)/2 and s(q,k) is the least value of Σ_b C(x_b, 2) over x_b >= 0 with
Σ x_b = q (k parts). The least value is taken at the balanced split (part sizes differ by at most 1).

**Proof.** Fix a partition. Each slope a in F_q gives one point β_{b(a)} + t a of R_t, where b(a)
is the part containing a. So R_t is the image of a q-element set, and

    |R_t| = q - Σ_p (m_p - 1),

where m_p is the number of slopes landing on point p in row t. Two slopes a ≠ a' land on the
same point of row t exactly when β_{b(a)} - β_{b(a')} = t (a' - a). If b(a) = b(a') this is
impossible. If b(a) ≠ b(a'), the left side is non-zero, so it happens for exactly one t ≠ 0,
namely t = (β_{b(a)} - β_{b(a')}) / (a' - a). So, summed over the q - 1 rows, the number of
colliding pairs is exactly P = Σ_{b < b'} x_b x_b' (x_b = |A_b|). Because m - 1 <= C(m, 2)
for every m >= 1,

    Σ_{t ≠ 0} |R_t| >= q(q-1) - P = q(q-1) - (C(q,2) - Σ_b C(x_b,2)) = C(q,2) + Σ_b C(x_b,2). ∎

Equality holds when no point of any row is hit three or more times. The bound needs no
choice of β: it holds for every set of positions.

**Corollary (n = 3).** If some plane meets a Kakeya set K of F_q^3 in k collinear points, then

    |K| >= k + q · (C(q,2) + s(q,k)) >= k + q^2 (q-1)/2.

## Checks against the values already proved by hand

| q | k | bound C(q,2) + s(q,k) | exact value (earlier note) |
|---|---|---|---|
| 4 | 1 | 6 + 6 = 12 | 12 |
| 4 | 2 | 6 + 2 = 8 | 8 |
| 4 | 3 | 6 + 1 = 7 | 8 (a triple point is forced) |
| 4 | 4 | 6 + 0 = 6 | at least 6 |
| 3 | 2 | 3 + 1 = 4 | 4 |
| 3 | 3 | 3 + 0 = 3 | 4 |

## q = 5 (for a future exact search of K(3,5); known range 39 to 53)

Balanced splits of 5: k = 1 gives s = 10, k = 2 (3+2) gives 4, k = 3 (2+2+1) gives 2,
k = 4 (2+1+1+1) gives 1, k = 5 gives 0. So

| k | m(5,k) >= | collinear k-point plane section forces |K| >= |
|---|---|---|
| 1 | 20 | 1 + 5·20 = 101 |
| 2 | 14 (exact) | 2 + 5·14 = 72 |
| 3 | 12 (an example gives 13) | 3 + 5·12 = 63 |
| 4 | 11 | 4 + 5·11 = 59 |
| 5 | 10 | 5 + 5·10 = 55 |

For k = 5 the section is a full line, the planar reading is a Kakeya set of AG(2,5), and the
known planar bound for odd q (q(q+1)/2 + (q-1)/2 = 17 points) gives m(5,5) >= 12 and |K| >= 65.
The elementary 55 is already enough below.

Worked example for k = 2 (checked point by point): β = {0, 1}, A_0 = {0, 1}, A_1 = {2, 3, 4}.
R_1 = {0, 1, 3, 4}, R_2 = {0, 2, 4}, R_3 = {0, 2, 3}, R_4 = {0, 2, 3, 4}: 4 + 3 + 3 + 4 = 14 = 20 - 6.

Example for k = 3: β = {0, 1, 2}, A_0 = {0, 1}, A_1 = {2, 3}, A_2 = {4}. R_1 = {0, 1, 3, 4},
R_2 = {0, 2} (the point 0 is hit three times), R_3 = {0, 2, 3, 4}, R_4 = {0, 3, 4}: total 13 = 20 - 8 + 1.
So m(5,3) is 12 or 13.

**Consequence for K(3,5).** A Kakeya set of F_5^3 with at most 52 points has no collinear
plane section (every case above is at least 55). In a smallest-section search at bound 52,
the smallest section B (3 <= |B| <= 10, because 5 parallel planes hold at most 52 points and
n0 = 1, 2 are excluded by the two-point lemma) can be taken non-collinear. That removes
the collinear orbits with 3, 4 and 5 points. It also gives a constraint the search can add
for free: every plane meets K in at least 3 non-collinear points.

## A structural remark

The corollary says a collinear plane section costs about q^3/2 points. Known constructions
are much smaller: for even q the Kopparty-Saraf-Sudan set (as quoted by Bukh and Chao,
arXiv:2108.00074) has q^3/4 + 3q^2/4 points, and 1 + q^2 (q-1)/2 is larger than that for every
q >= 5. So **for even q >= 8, every smallest Kakeya set in F_q^3 meets every plane in at least
3 non-collinear points**. The same holds for q = 5 if the smallest set has at most 54 points
(a 53-point set is known). For q = 4 it does not follow: a full-line section gives only
4 + 4·6 = 28 = K(3,4), so a 28-point set with a full-line section is not ruled out here.
We have not checked whether this remark appears in the literature. It is simple enough that
it may. A search on 2026-11-22 found nothing on plane sections of Kakeya sets in AG(3,q).
