# Layered Kakeya sets: a census, and why F_3^3 is different

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-26. Made with AI agents working in MV2's
name. **No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Layered sets

Write a point of F_3^n as (c, x) with c in F_3 and x in F_3^(n-1). A set

    K = {0} x L0  u  {1} x L1  u  {2} x L2

is **layered** if L0 is itself a Kakeya set of F_3^(n-1). Then K is Kakeya iff for every
v in F_3^(n-1) some x in L0 has x + v in L1 and x - v in L2 (`k8/layer.py`). The idea of
building Kakeya sets layer by layer goes back to Saraf and Sudan (arXiv:0808.2499, appendix)
and Bukh and Chao (arXiv:2108.00074); here it is used for exact small sizes.

Our 211-point Kakeya set of F_3^7 is layered: L0 is our 105-point set of F_3^6 and
|L1| + |L2| = 106. Our 427-point set of F_3^8 is layered on top of it: L0 is the 211-point
set and |L1| + |L2| = 216 (`k6/RECORDS.md`). Lev's 2009 pattern
k_n = k_(n-1) + 2 k_(n-2) (3, 7, 13, 27, 53, 107, ...) has the same shape. So: how small
can a layered set be?

## Fact: a layered Kakeya set of F_3^3 has at least 14 points

So the 13-point Kakeya sets of F_3^3 (k_3 = 13) are **not** layered: no plane section of
a 13-point Kakeya set of F_3^3 is a Kakeya set of its plane. The pattern 13 = 7 + 2 * 3 is
not realised by layers.

*Proof.* Let K be layered in F_3^3 with |K| <= 13 (any plane section can be moved to
c = 0 by an affine map, which keeps K Kakeya). Since k_2 = 7, |L0| >= 7 and
|L1| + |L2| <= 6. The map v -> (x + v, x - v) from the 9 vectors v to L1 x L2 is injective
(v = (y - z) / 2 is recovered from the pair), so |L1| |L2| >= 9. Hence |L1| = |L2| = 3,
|L0| = 7, and v -> (y, z) is a bijection onto L1 x L2: every pair (y, z) is used, with
v = z - y, and its line's third point x = -(y + z) lies in L0.

So L2 - L1 = F_3^2 with every difference once, i.e. F_3^2 = L2 + (-L1) is a factorisation
of the group into two 3-sets. One factor is a line (a coset of a 1-dimensional subspace).
Indeed, translate so that 0 is in L1 and suppose L1 = {0, p, q} is not a line, so p, q are
independent. For a character chi_u(w) = omega^(u.w), the factorisation gives
hat(1_L2)(u) * (1 + chi_u(-p) + chi_u(-q)) = 0, and 1 + omega^a + omega^b = 0 only when
{a, b} = {1, 2}. As u runs over the 8 nonzero vectors, (u.p, u.q) runs over all nonzero
pairs, and only u1 and -u1 give {1, 2}. So 1_L2 has Fourier support in {0, u1, -u1}, hence
is constant on the cosets of ker(u1): L2, of size 3, is one such coset, a line.

If L2 = z0 + <d> is a line, then for each y in L1 the differences z - y fill the coset
z0 - y + <d>, and these three cosets are distinct (all 9 differences differ). The third
points x = -y - z fill the cosets -y - z0 + <d>, which are the previous cosets shifted by
-2 z0, so also distinct: the 9 points x are all of F_3^2. If L1 is a line the same holds
with the roles exchanged. Either way L0 = F_3^2, contradicting |L0| = 7. QED.

`k8/test_census.py` checks this exhaustively (every L0 of size 7, 8, 9 and every L1, L2
with |L1| + |L2| <= 13 - |L0|), and asserts that the least layered size in F_3^3 is 14 or
15 (the classical 15-point set is layered: L0 = the classical 7-point set, layers {0,1}^2
and {0,2}^2).

**The least layered size in F_3^3 is exactly 14.** The exact solver `census.exact_min_extra`
found, in the oracle run below, a 7-point L0 with |L1| + |L2| = 7.

## The census (`k8/census.py`), results

Publish-gate oracle `k8/test_oracle_census.py` of run 20261126T060550Z, 2026-11-26
06:50-06:57 UTC, 5 workers, about 410 s each. For many different Kakeya sets L0 (made by
`k5/walk.py` with fresh seeds, sometimes with a few random extra points, and the classical
set), a layered search of 2.5 to 4 s looked for the least |L1| + |L2|. Entries are
"min extra [extra x how many L0]".

**F_3^5** (L0 in F_3^4, 331 distinct L0 over 4 workers; record 53, Chung, Du and Wesley):

| abs L0 | L0 tried | least abs L1 + abs L2 | distribution | least total |
|---|---|---|---|---|
| 27 | 126 | 27 | 27 x 14, 28 x 110, 29 x 2 | 54 |
| 28 | 53 | 27 | 27 x 3, 28 x 48, 29 x 2 | 55 |
| 29 | 42 | 27 | 27 x 6, 28 x 27, 29 x 9 | 56 |
| 30 | 94 | 26 | 26 x 2, 27 x 7, 28 x 72, 29 x 13 | 56 |
| 31 | 13 | 26 | 26 x 2, 27 x 2, 28 x 8, 29 x 1 | 57 |
| 32 | 16 | 27 | 27 x 2, 28 x 14 | 59 |

No layered set reached 53, let alone 52: the least was 54 = 27 + 27. So the census gives
no set below the record; and since the 53-point record is not reached either, the
layered searches are probably not how a 52-point set would be found.

**F_3^4** (L0 in F_3^3, 39 distinct L0; k_4 = 27): abs L0 = 13: least extra 14 (13 of 13);
abs L0 = 14: least extra 14 (10 of 10); abs L0 = 15: least extra 13 (6 of 16). Least total
27 = 13 + 14, equal to k_4. In F_3^4 the 13-point layer costs 14 = 2 * 7, as the pattern
would say; in F_3^3 a 7-point layer costs 7 = 2 * 3 + 1, one more than the pattern.

**F_3^6 on layered F_3^5 sets** (record 105): L0 = the 53-point set of Chung, Du and Wesley
-> extra 54 (107); layered F_3^5 sets of 54 and 55 points -> extra 55 (109, 109, 110).

Searches give upper bounds only. A census that never goes below 53 is evidence about one
family of sets, not a proof that k_5 = 53.
