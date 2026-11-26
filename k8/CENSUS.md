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
|L1| + |L2| = 106. Lev's 2009 pattern k_n = k_(n-1) + 2 k_(n-2) (3, 7, 13, 27, 53, 107, ...)
has the same shape. So: how small can a layered set be?

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
with |L1| + |L2| <= 13 - |L0|), and computes the least layered size in F_3^3; the classical
15-point set is layered (L0 = the classical 7-point set, layers {0,1}^2 and {0,2}^2), so the
least size is 14 or 15.

## The census (`k8/census.py`)

For many different Kakeya sets L0 of F_3^4 (from `k5/walk.py` with fresh seeds, sometimes
with a few random extra points, and the classical 31-point set), a layered search looks for
the least |L1| + |L2|. A verified layered set of at most 52 points in F_3^5 would beat the
53-point record of Chung, Du and Wesley (arXiv:2608.23691). The same is done for F_3^4
(L0 in F_3^3, record 27) and recursively for F_3^6: layered sets of F_3^5 found this way,
and the published 53-point set, are used as L0 (record 105). Results go to `CENSUS.md`
once the gate oracle has run.

Searches give upper bounds only. A census that never goes below 53 is evidence about one
family of sets, not a proof that k_5 = 53.
