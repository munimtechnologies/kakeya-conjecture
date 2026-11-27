# Rainbow triples of F_3^3 cannot beat the Kopparty–Lev–Saraf–Sudan base

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Made with AI agents working in MV2's
name. **No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Theorem

**Every rainbow triple (A, M, C) of F_3^3 has |A|·|M|·|C| ≥ 308.**

Since (3/2^(2/3))^9 = (27/4)^3 = 307.55…, a triple of F_3^3 would give, through the tripling
lemma of `TRIPLES.md`, a Kakeya growth base over F_3 below the 1.88988… of Kopparty, Lev,
Saraf and Sudan ([arXiv:1003.3736](https://arxiv.org/abs/1003.3736)) exactly when its product
is at most 307. The theorem says there is no such triple: the base of dimension 3 is at
least 308^(1/9) = 1.8902…, and products with triples of dimensions 1 and 2 (bases 2 and 2)
cannot go lower. In October 2026 (`TRIPLES.md`) an energy argument had left the size pattern
(6, 7, 7), product 294, as the only one below 308, and CP-SAT returned UNKNOWN on it after
600 s; the theorem rules it out.

## Proof

Definitions, the normal form N1–N3, the pair bound and the hyperplane lemma H are in
`SCAN.md`. Take a triple with product P ≤ 307 in normal form, so |A| ≤ |M| ≤ |C|.

1. |A|³ ≤ P ≤ 307 < 343, so |A| ≤ 6.
2. If A lies in a plane, lemma H gives P ≥ 9·τ(2) = 576. So A spans F_3^3, which needs at
   least 4 points: |A| ∈ {4, 5, 6}.
3. |M| ≥ max(|A|, ⌈27/|A|⌉) and |A|·|M|² ≤ P leave |M| ∈ {7, 8} for |A| = 4 and
   |M| ∈ {6, 7} for |A| = 5 and 6.
4. `rainbow/test_tau3_308.py` scans every full-dimensional affine class of A of these sizes,
   both sizes of M, and every M in normal form N3 with M − A = F_3^3, and computes the least
   C exactly. No (A, M) gives |A|·|M|·|C| < 308. ∎

The scan does not use the energy lemma, so it re-checks every size pattern below 308. Its
controls are in `rainbow/test_scan_controls.py` (see `SCAN.md`).

## What is open

* The exact value of τ(3): `TAU3.md` and `test_tau3_448.py` (bound 448, attained by the
  (7, 8, 8) triple of October 2026).
* Dimension 4: a triple with |A||M||C| ≤ 2075 (< (27/4)^4 = 2075.9) would beat the
  Kopparty–Lev–Saraf–Sudan base. By lemma H, all three of its sets would span F_3^4, since
  9·τ(3) ≥ 9·308 = 2772. The best known product in F_3^4 is 448·8 = 3584.

## Limits

The theorem rests on an exhaustive computer search and on short reductions that are proved
in `SCAN.md` but have not been checked by a human. It says nothing about dimensions 4 and
higher, where a rainbow triple could still beat the Kopparty–Lev–Saraf–Sudan base.
