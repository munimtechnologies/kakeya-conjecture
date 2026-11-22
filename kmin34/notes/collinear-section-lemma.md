# Collinear-section lemma for Kakeya sets in F_q^n (hand proof)

MV2 of Munim, Inc. Text CC BY 4.0. Written 2026-11-22 by AI agents working in MV2's name;
no human mathematician has checked it. It is a hand sanity check for the computer-assisted
bound K(3,4) >= 28 (results/kmin34-v1), not a replacement for that computation or for the
queued DRAT certificate (kmin34cert-v1).

## Setting

K is a Kakeya set in F_q^n. H is the hyperplane "last coordinate 0" (layer 0), and B = K ∩ H
has k points. Every one of the q^(n-1) slanted directions (d, 1), d in F_q^(n-1), has a line
of K in it. That line meets H in exactly one point, which must lie in B. Pick one such line
for each d. That gives a partition of the direction space F_q^(n-1) into parts D_b (b in B):
the line in direction (d,1) is b + t(d,1), and in layer t it gives the point b + t d.
So for t ≠ 0, layer t of K contains S_t = ∪_b (b + t D_b), and

    |K| >= k + Σ_{t ≠ 0} |S_t|.                                            (1)

## Lemma (collinear section)

Suppose B lies on one line of H, with direction v. Write b = b_1 + β_b v (the β_b are
distinct elements of F_q). Then

    |K| >= k + q^(n-2) · m(q, k),

where m(q,k) is the least value of Σ_{t ≠ 0} |∪_b (β_b + t A_b)| over all partitions
(A_b) of F_q into parts indexed by the k base points.

**Proof.** Split the direction space into the q^(n-2) cosets C = c + <v>. A direction
c + s v in D_b gives the point b + t(c + s v) = b_1 + t c + (β_b + t s) v in layer t. For
fixed t ≠ 0 all of these lie in the single coset b_1 + t c + <v>, and different C give
different cosets. So S_t splits into q^(n-2) disjoint pieces, one per C, each a copy of
∪_b (β_b + t A_b), where A_b = {s : c + s v in D_b} partitions F_q. Sum over t and C and use (1). ∎

**Planar reading.** For one coset, the points (β_b + t a, t) for a in A_b, t in F_q form, in an
abstract plane F_q^2, a union of lines in all q non-horizontal directions, which meets row 0
exactly in the k points β_b. So m(q,k) + k is the least size of such a planar set. If k = q,
row 0 is a full line, the planar set is a planar Kakeya set, and the elementary bound
(q+1 lines in different directions, any two meeting at most once) gives
m(q,q) >= q(q+1)/2 - q = q(q-1)/2.

## Values for q = 4 (F_4 = {0, 1, w, w^2}, w^2 = w + 1)

- **k = 1:** m = 12 (a full fan).
- **k = 2:** per coset 3(x_1 + x_2) - x_1 x_2 with x_1 + x_2 = 4, so m = 8. This agrees with the
  two-point base lemma (|K| >= 34 in F_4^3).
- **k = 3:** take β_1 = 0. The part sizes are (4,0,0), (3,1,0), (2,2,0) or (2,1,1).
  - (4,0,0) gives 12, and (3,1,0) gives at least 9, since every layer has at least 3 points.
  - (2,2,0): t A_1 ∪ (β_2 + t A_2) = t (A_1 ∪ (A_2 + β_2 u)) with u = 1/t. A_1 and A_2 are the two
    cosets of an additive subgroup {0, e} of F_4. Two of the three translations swap them (size 2)
    and one fixes them (size 4). Total 8.
  - (2,1,1) with A_1 = {s_1, s_2}, A_2 = {a}, A_3 = {a'}: point 2 lies at t(a + β_2 u). As u runs over
    F_4*, a + β_2 u runs over F_4 \ {a} = A_1 ∪ {a'}, so it leaves t A_1 at exactly one t, namely
    β_2 u = a' - a. Point 3 leaves t A_1 exactly when β_3 u = a - a' (= a' - a, because the
    characteristic is 2). Those two layers differ because β_2 ≠ β_3. At each of them the escaping
    point is a new point (the other one is inside t A_1). Total 6 + 2 = 8.
  - So **m(4,3) = 8**.
- **k = 4 (a full line):** m(4,4) >= 4·3/2 = 6.

## Consequence for K(3,4) (n = 3, q = 4, at most 27 points)

The 22 cases of kmin34 are the AGL(2,4)-orbits of the smallest plane section B, with
|B| = n0 = 1..6. By the lemma:

- n0 = 1: |K| >= 1 + 4·12 = 49. n0 = 2: |K| >= 2 + 4·8 = 34.
- **n0 = 3, B collinear:** |K| >= 3 + 4·8 = **35**.
- **n0 = 4, B a full line:** |K| >= 4 + 4·6 = **28**.

So 4 of the 22 cases are impossible at 27 by hand: both n0 <= 2 orbits, the collinear
n0 = 3 orbit, and the full-line n0 = 4 orbit. (For n0 = 3 there are exactly two orbits,
collinear and triangle, because AGL(2,4) is transitive on non-collinear triples and on
collinear triples.)

## The triangle case is not closed by fan counting alone

Take B = {0, e_1, e_2} and u = 1/t. Then S_t = t·W_u, where
W_u = D_1 ∪ (D_2 + u e_1) ∪ (D_3 + u e_2). Inclusion-exclusion summed over u gives
Σ_u |W_u| = 48 - P + T. Here P counts pairs (d, d') in different parts D_i, D_j whose
difference is parallel to b_j - b_i. T counts homothetic triangles {d, d + u e_1, d + u e_2}
with one vertex in each of D_1, D_2, D_3. The partition

- D_1 = the diagonal {(a, a)},
- D_2 = {(0,1),(1,0),(0,w),(w,0),(0,w^2),(w^2,0)},
- D_3 = {(w,w^2),(w^2,w),(1,w^2),(w^2,1),(1,w),(w,1)}

has P = 12 + 12 = 24 and T = 0, so |W_1| = |W_w| = |W_{w^2}| = 8 (checked point by point).
That makes (1) equal to exactly 3 + 24 = 27. This particular K fails anyway. The only full line
inside any W_u is the diagonal, so the horizontal directions e_1, e_2, (1,w) and (1,w^2) have no
line of K, and adding points goes past 27. But a hand exclusion of the triangle case would have
to use the horizontal directions (and the other partitions with P - T = 24). We leave it to the
solvers. kmin34-v1 closed both n0 = 3 cases with both solvers in at most 4.1 s.

## Small extra for F_3^4 (Lev's k_4 = 27, bound 26)

For q = 3 a collinear 3-point section is a full line. The planar bound gives m(3,3) >= 3.
A direct check gives m(3,3) = 4. With parts (1,1,1), A_b = {σ(β_b)} for a permutation
σ of F_3. Every permutation of F_3 is affine, σ(β) = cβ + e with c ≠ 0, so the layer-t
image of β ↦ (1 + tc)β + te has 1 point when t = -1/c and 3 points at the other t, giving
4 in total. Parts (2,1,0) give at least 2 + 2 = 4, and (3,0,0) gives 6. So a full-line
section forces |K| >= 3 + 9·4 = 39 > 26 in F_3^4.
