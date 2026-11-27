# An exhaustive scan for rainbow triples of F_3^m

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Made with AI agents working in MV2's
name. **No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Rainbow triples and the Kakeya base over F_3

A triple (A, M, C) of subsets of F_3^m is **rainbow** if every d ∈ F_3^m, d = 0 included,
has a witness a ∈ A with a + d ∈ M and a + 2d ∈ C. Write τ(m) for the least |A|·|M|·|C|.
By the tripling lemma of `TRIPLES.md` (October 2026), a rainbow triple of F_3^m gives a
Kakeya set of 3·|A||M||C| points in F_3^(3m), and products of rainbow triples are rainbow,
so the Kakeya growth base over F_3 is at most τ(m)^(1/(3m)) for every m. The best known
upper bound for that base is 3/2^(2/3) = 1.88988… (Kopparty, Lev, Saraf and Sudan,
[arXiv:1003.3736](https://arxiv.org/abs/1003.3736)). A triple beats it exactly when
|A||M||C| < (27/4)^m: at most 307 for m = 3 and at most 2075 for m = 4.

## Normal form

Every rainbow triple can be moved to one with the same multiset of sizes, still rainbow,
that satisfies N1 to N3.

* **N1 (role permutations).** If a witnesses d for (A, M, C), then a + d witnesses d for
  (M, C, A), and a + 2d witnesses −d for (C, M, A). A 3-cycle and a reversal generate all
  six orders, so we may take |A| ≤ |M| ≤ |C|.
* **N2 (affine maps).** An affine bijection φ applied to all three sets keeps the triple
  rainbow (φ(a) witnesses the image direction), so A may be the canonical representative
  of its affine class. `rscan.canon` takes, over every ordered affine frame of A, the affine
  map sending it to (0, e_1, …, e_k) and keeps the least sorted image; `rscan.classes`
  builds every class of x-sets from the classes of (x − 1)-sets. The class lists are
  checked by orbit–stabiliser counts (`test_scan_controls.py`, `test_tau3_448.py`).
* **N3 (translations).** (A, M + t, C + 2t) is rainbow when (A, M, C) is (the witness of d
  serves d + t). So, for any m ∈ M, (A, M − m, C − 2m) is rainbow with 0 ∈ M − m; we take
  the m that makes the bit mask of M − m least.

**Size bounds.** d ↦ (a_d, a_d + d) is injective from F_3^m to A × M, so |A||M| ≥ 3^m.
With N1, |M| ≥ max(|A|, ⌈3^m/|A|⌉), and a product P forces |A|·|M|² ≤ P.

**Hyperplane lemma H.** *If one of A, M, C lies in an affine hyperplane of F_3^m, then
|A||M||C| ≥ 9 τ(m − 1).* Proof: by N1 take it to be A; translate so that A ⊆ H, a linear
hyperplane, and pick v ∉ H. Write points as h + s v (h ∈ H, s ∈ F_3) and let
M_s = {h : h + s v ∈ M}, C_s = {h : h + s v ∈ C}. For d = d_H + s v and a ∈ A ⊆ H,
a + d ∈ M iff a + d_H ∈ M_s, and a + 2d ∈ C iff a + 2d_H ∈ C_{2s}. So for each s the
triple (A, M_s, C_{2s}) is rainbow in H ≅ F_3^(m−1), and |M_s||C_{2s}| ≥ τ(m − 1)/|A|.
As s ↦ 2s permutes F_3, Cauchy–Schwarz gives
|M||C| = (Σ_s |M_s|)(Σ_s |C_{2s}|) ≥ (Σ_s (|M_s||C_{2s}|)^(1/2))² ≥ 9 τ(m − 1)/|A|. ∎

For m = 3, τ(2) = 64 (`TRIPLES.md`, re-found here) gives 576. So below 576 all three sets
of a triple of F_3^3 span F_3^3, and in particular have at least 4 points.

## The scan

For a fixed A and size s, `rscan.scan_task` walks through the sorted sets M = {0 < p_1 < …}
of size s with M − A = F_3^m (pruned when an uncovered direction can no longer be covered,
or when more directions are uncovered than the remaining points can serve), keeps the M
that are least among their translates (N3), and computes the least C exactly. (A, M, C) is
rainbow iff C meets every set

    P_d = {a + 2d : a ∈ A, a + d ∈ M},   d ∈ F_3^m,

so the least C is a minimum hitting set of the P_d (`rscan.min_hit`: branching on the
smallest unmet set, pruned by a greedy packing of pairwise disjoint sets, and checked
against brute force on 300 random instances). A triple with product below the bound is
reported if the least C gives |A|·s·|C| < bound. This is a different method from the
October CP-SAT models and from the witness-by-witness branch and bound tried earlier: it
enumerates M directly, so its running time does not depend on how the solver branches.

## Controls (`test_scan_controls.py`)

* The checker `check_triple` (straight from the definition) accepts ({0,1}, {0,1}, {0,2})
  and its products, and rejects near misses.
* τ(1) = 8 and τ(2) = 64 are found by the scan, every triple it reports is checked, and
  nothing is found below them.
* A inside a line of the plane gives nothing below 9·τ(1) = 72 (lemma H one dimension down).
* The (7, 8, 8) triple of the October CP-SAT run is rainbow; put in normal form N3, the
  scan re-finds a product ≤ 448 from it, both from the whole M and from its first three
  points.

## First case of the F_3^3 bound

The same test file proves: **no rainbow triple of F_3^3 whose smallest set has 4 points has
product below 308.** By lemma H such a set spans F_3^3, so it is an affine frame (one class);
the normal form leaves |M| ∈ {7, 8}, and the scan finds no M that works.

The full bound (smallest sets of 4, 5 and 6 points) is `test_tau3_308.py` and
`RAINBOW3.md`; the exact value τ(3) = 448 is `test_tau3_448.py` and `TAU3.md`. Each of them
is committed only after its test passes.

## Limits

Everything rests on the reductions above, proved here but not checked by a human, and on
the correctness of the scan code, which the controls test but cannot prove.
