# AAB pairs: a product construction of Kakeya sets over F_3

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Made with AI agents working in MV2's
name. **No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Definition

An **AAB pair** of F_3^m is a pair (A, B) of subsets with A ∩ B non-empty such that every
direction d of F_3^m has a line {p, p + d, p + 2d} with

    p ∈ A,   p + d ∈ A,   p + 2d ∈ B.

The condition for d and for −d is the same (put p′ = p + d). Equivalently: for every
direction there are two points x, y of A in that direction whose third collinear point
−(x + y), the **centre** of the pair {x, y}, lies in B. Then A ∪ B is a Kakeya set. We write
a = |A|, b = |B|, e = |A ∩ B|. Example: in F_3, A = {0, 1}, B = {0, 2} (a = b = 2, e = 1).

## Lemma (product)

If (A1, B1) is an AAB pair of F_3^m1 and (A2, B2) is an AAB pair of F_3^m2, then
(A1 × A2, B1 × B2) is an AAB pair of F_3^(m1 + m2), so (A1 × A2) ∪ (B1 × B2) is a Kakeya
set with a1 a2 + b1 b2 − e1 e2 points.

*Proof.* (A1 × A2) ∩ (B1 × B2) = (A1 ∩ B1) × (A2 ∩ B2) is non-empty. Let D = (d1, d2) be a
nonzero direction.

* d1 ≠ 0 and d2 ≠ 0: take p1 for d1 and p2 for d2 from the two pairs. Then (p1, p2) and
  (p1 + d1, p2 + d2) lie in A1 × A2 and (p1 + 2d1, p2 + 2d2) lies in B1 × B2.
* d2 = 0: take p1 for d1 and any y ∈ A2 ∩ B2: (p1, y), (p1 + d1, y) ∈ A1 × A2 and
  (p1 + 2d1, y) ∈ B1 × B2.
* d1 = 0: the same with x ∈ A1 ∩ B1. ∎

So a pair (a, b, e) of F_3^m gives Kakeya sets of a^t + b^t − e^t points in F_3^(mt), and
pairs can be mixed. The powers of the pair of F_3 are the classical sets {0,1}^n ∪ {0,2}^n
with 2^(n+1) − 1 points. Since A ∪ B is Kakeya, a + b − e ≥ k_m.

## Theorem (F_3^3 pairs cannot beat the records in F_3^6 and F_3^9)

There is no AAB pair of F_3^3 with (a, b, e) = (7, 7, 1), (6, 8, 1) or (8, 6, 1).

Because k_3 = 13, every pair of F_3^3 has a + b − e ≥ 13, and these three triples are the
only ones with a² + b² − e² < 105 or a³ + b³ − e³ < 844 (the next ones are (5, 9, 1) with
105 and 853, and (7, 8, 2) with 109 and 847). So no product of F_3^3 pairs, alone or with the
pair of F_3, gives a Kakeya set below the explicit records 105 in F_3^6 and 844 in F_3^9
(`k6/RECORDS.md`). The question was raised by our previous run, whose conditional commit
(it would have landed only if a (7, 7, 1) pair existed) did not land.

**Proof of record:** `aab/test_aab3_none.py`, an exhaustive search. It takes every set A of
6, 7 or 8 points of F_3^3 containing 0 (any pair can be translated so) whose pairs cover all
13 directions, and for each runs an exact branch-and-bound search for the least B with one
point inside A. Positive controls check that the same code finds the classical pair
({0,1}^3, {0,2}^3) and checks the product lemma on small cases.

**Hand argument for (7, 7, 1) and (6, 8, 1)** (written by the agents, not human-checked).
Write s(t) for the number of pairs of A centred at t. A point t ∈ B \ A serves exactly the
directions of the pairs centred at t; a point c ∈ A ∩ B serves the directions of the lines
through c inside A. Every direction must be the direction of some pair of A.

*Lemma S.* If |A| ≤ 7 and A contains six points t ± u, t ± v, t ± w (u, v, w in distinct
directions), the pairs of A miss a direction. If u, v, w are independent, make them the
basis: the six points have differences of weight ≤ 2, and a seventh point p gives at most
three weight-3 directions p − t ∓ e_i (one sign per i), so one of the four weight-3
directions is missed. If u, v, w are coplanar, the six points lie in a plane, which has 4
directions, and the seventh point has only 6 pairs for the other 9. Hence, when |A| ∈ {6, 7}:
s(t) ≤ 2 for t ∉ A, and a point of A lies on at most 2 lines inside A.

*(7, 7, 1).* Here B = {c} ∪ T with T ∩ A = ∅ and |T| = 6, and the 13 directions need
s(c) + Σ_T s(t) ≥ 13 with every term ≤ 2, so c lies on one or two lines inside A.

* Two lines: put c = 0, A = {0, ±u, ±v, q, r}, P = the plane spanned by u, v. If q ∈ P, only
  the 6 pairs of r can serve the 9 directions outside P. If q, r lie in opposite cosets of
  P, every centre other than the one of {q, r} serves one pair only, so B serves at most
  2 + 2 + 5 = 9 directions. If q, r lie in the same coset, the directions u + v and u − v
  can only be served by centres serving nothing else, which leaves 4 points with s ≤ 2 for
  the 9 directions outside P.
* One line: put c = 0, A = {0, ±e1, p1, …, p4}; then all six points of T have s = 2 and the
  12 pairs centred there have 12 different directions. Project along e1 to F_3^2, with
  x_i the image of p_i (x_i ≠ 0). The three pairs {z, p_i}, z ∈ {0, ±e1}, have exactly the
  three directions above x_i, so each of the four directions of F_3^2 must receive exactly 3
  of the 12 pairs. Two points of T with s = 2 are four points with a1 + a2 = a3 + a4, which
  forces x_i = x_j, x_i = x_j + x_k or x_i + x_j = x_k + x_l. Checking how the x_i are spread
  over the four directions of F_3^2 (4, 3+1, 2+2: too few pairs for the empty directions;
  2+1+1: in each normal form at most one such relation is available; 1+1+1+1: all six
  pairs {p_j, p_k} would need x_j + x_k to be another x_i, impossible) leaves fewer than
  the six relations needed.

*(6, 8, 1).* B = {c} ∪ T with |T| = 7. Let ℓ be the number of lines inside A. The pairs
centred outside A number 15 − 3ℓ, and c serves at most the lines through it, so ℓ ≤ 1.
If ℓ = 1 (through c = 0, direction e1), all 12 other pairs must be served with distinct
directions, so the x_i (as above) are distinct and the pairs {z, p_i} have 9 distinct
centres, more than the 7 points of T. If ℓ = 0, A is a cap, c serves nothing, and six
points of T need s = 2: six distinct 4-point subsets of A with a1 + a2 = a3 + a4. Such a
subset is a planar parallelogram Q; A meets its plane in Q only (a cap of AG(2, 3) has at most
4 points), and putting Q = {0, e1, e2, e1 + e2} and checking the two remaining points in the
same or in different parallel planes gives at most 3 such subsets. ∎

## What this means, and what is open

* The F_3^3 route to a growth base below 2 is closed: base max(a, b)^(1/3) < 2 needs
  max(a, b) ≤ 7, and (7, 7, 1) does not exist.
* Still open, and decided by the report-only probe `ORACLE-P16`: a (7, 8, 2) or (8, 7, 2)
  pair of F_3^3 would give 8^t + 7^t − 2^t points in F_3^(3t), below the current records
  from F_3^12 on (6481 against 6759); a pair of F_3^4 with a² + b² − e² < 422, for example
  (14, 14, 1), (12, 16, 1) or (14, 15, 2), would beat the 422-point record in F_3^8; and a
  pair of F_3^5 with max(a, b) ≤ 29 and e small would beat the 1689-point record in F_3^10.

## Prior work

We found no statement of the product lemma. The closest are the classical construction (the
powers of the pair of F_3), Lev's 2009 pattern 3, 7, 13, 27, 53 (comment 994 on T. Tao's
blog), our doubling lemma (`k9/SPLIT.md`), the random construction of Kopparty, Lev, Saraf
and Sudan ([arXiv:1003.3736](https://arxiv.org/abs/1003.3736), Theorem 8, growth base about
1.8899) and our explicit base-1.8899 construction (repository kakeya-small-fields). Those are
asymptotic; AAB pairs are about explicit sets in moderate dimensions.

## Limits

The theorem rests on the exhaustive search in `aab/test_aab3_none.py` and on k_3 = 13 (also
proved by exhaustive search in this repository); the hand argument above is a sketch by the
agents and has not been checked by a human. It says nothing about pairs in F_3^4 and above.
