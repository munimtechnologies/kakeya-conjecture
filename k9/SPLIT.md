# Splits of Kakeya sets in F_3^m and a doubling lemma

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-26 (lemmas), 2026-11-27 (result).
Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Result

The 105-point Kakeya set of F_3^6 (`k6/sets/kakeya6_105.json`) has a split with
**excess 1** (`k9/splits/kakeya6_105_split.json`). By Lemma 2 it lifts to splits with
excess 1 of the layered records 211, 422, 844 and 1693, and Lemma 1 on the 844-point set
of F_3^9 then gives explicit Kakeya sets of

    845 · 2^(n−9) − 1   points in F_3^n,  for every n ≥ 9:

1689 in F_3^10 (previous record 1693), 3379 in F_3^11, 6759 in F_3^12, 13519 in F_3^13, ...
Together with the records 105, 211, 422 and 844 in dimensions 6 to 9, **every dimension
n ≥ 6 has an explicit Kakeya set below Lev's 2009 pattern** (5 · 2^n + (−1)^n)/3, and the gap

    (5 · 2^n + (−1)^n)/3 − (845 · 2^(n−9) − 1) = 25 · 2^n / 1536 + 1 + (−1)^n / 3

grows without bound (18 points in F_3^10, 34 in F_3^11, 68 in F_3^12, 134 in F_3^13).
In ratio to 2^n the family tends to 845/512 ≈ 1.6504, against 5/3 for the pattern and 2
for the classical construction. Table: `k9/results/doubling.json`.

What is checked by computer, without the lemmas: the 1689-point set (all 29524 directions,
`k9/test_splits.py`), the 3379-point set (all 88573 directions) and the 6759-point set
(layer criterion over all 177147 vectors on the 3379-point layer, `k9/test_doubled_big.py`).
From F_3^13 on the sets rest on Lemma 1, whose proof is below; the split it needs is
verified for the 6759-point set.

## Definitions

Let K be a Kakeya set of F_3^m (a full line {a, a+d, a+2d} in every direction d). A
**split** of K is a pair (A, B) of subsets with A ∪ B = K and

    C(A, B) := { b − a : a ∈ A, b ∈ B, −(a + b) ∈ K } = F_3^m.

Its **excess** is e = |A ∩ B|. Always e ≥ 1, because the zero vector lies in C(A, B) only
through a pair a = b (then −(a + b) = a). For u ≠ 0 the condition asks for two points
a ∈ A, b ∈ B on a line of K with b − a = u (three distinct points of F_3^m are collinear
iff they sum to 0). For example, 2-colour K so that every direction has a line of K with
both colours, then put one point into both A and B: that is a split with excess 1.

## Lemma 1 (doubling)

If (A, B) is a split of K, then

    K⁺ = {0} × K ∪ {1} × A ∪ {2} × B

is a Kakeya set of F_3^(m+1) with |K| + |A| + |B| = 2|K| + e points.

*Proof.* Directions (0, d) are covered by {0} × K. A direction (1, v) is covered iff some
x ∈ K has x + v ∈ A and x + 2v ∈ B. Given a ∈ A, b ∈ B with x := −(a + b) ∈ K and
v := b − a, we get x + v = −2a = a and x + 2v = b − 3a = b. So the covered v are exactly
C(A, B) = F_3^m. Every direction of F_3^(m+1) is a multiple of some (0, d) or (1, v). ∎

## Lemma 2 (lifting)

Let K′ = {0} × K ∪ {1} × L₁ ∪ {2} × L₂ be **any** Kakeya set of F_3^(m+1) whose layer
x₀ = 0 is K, and let (A, B) be a split of K. Then

    A′ = {0} × A ∪ {1} × L₁,   B′ = {0} × B ∪ {2} × L₂

is a split of K′ with the same excess.

*Proof.* A′ ∪ B′ = K′ and A′ ∩ B′ = {0} × (A ∩ B). For U = (0, u), the pair (a, b) of K for
u gives ((0, a), (0, b)), whose third point (0, −(a + b)) is in {0} × K. Every other U is
±D with D = (1, v). K′ is Kakeya, so it holds a line p₀ = (0, x), p₁ = (1, x + v),
p₂ = (2, x + 2v), and p₁ ∈ A′, p₂ ∈ B′. Then p₂ − p₁ = D, with third point p₀ ∈ K′. For −D:
if x ∈ A take (p₀, p₂) (p₂ − p₀ = −D, third point p₁); if x ∈ B take (p₁, p₀)
(p₀ − p₁ = −D, third point p₂). One of the two holds because A ∪ B = K. ∎

## Corollary

K⁺ is layered with layer K, so Lemma 2 gives it the split (A ∪ {1} × A, B ∪ {2} × B), again
with excess e, and Lemma 1 applies again. Hence a split of a Kakeya set K ⊂ F_3^m with
excess e gives Kakeya sets of

    2^j (|K| + e) − e   points in F_3^(m+j),  for every j ≥ 0.

By Lemma 2 a split of one layer is enough for a whole tower of layered sets: our records
105 → 211 → 422 → 844 → 1693 (`k6/RECORDS.md`) each have the previous record as their layer
x₀ = 0. So a split of the 105-point set of F_3^6 with excess e gives splits of all of
them with excess e. Moreover **every** Kakeya set of F_3^(m+1) that a layered search builds
on one of these layers inherits a split of excess 1, so any future layered record on this
chain doubles in the same way: a layered set of s points in F_3^n gives 2^j (s + 1) − 1 in
F_3^(n+j).

Starting from K = F_3 with A = {0, 1}, B = {0, 2} (excess 1) the lemmas give 2^(n+1) − 1
points in F_3^n, exactly the classical construction {0,1}^n ∪ {0,2}^n. Our earlier sticky
family 31 · 2^(n−4) − 1 (repository kakeya-small-fields) has the same shape, from a 30-point
sticky set of F_3^4.

## Code and checks

* `k9/split.py`: the definitions, `verify_split` (checks A ∪ B = K and computes C(A, B) from
  the definition, with coordinate tuples), `double` (Lemma 1 and the split of the
  Corollary), `lift` (Lemma 2), and `find_split`, a WalkSAT-style local search over the
  colourings "A only / B only / both" of K that lowers the excess step by step.
* `k9/splits/kakeya6_105_split.json`, `k9/splits/kakeya7_211_split.json`: the stored splits
  (one state per point, encoding in the files).
* `k9/test_split.py`: rebuilds the classical sizes 7, 15, 31, 63 with the lemmas and checks
  each set by directions; searches a split of the 105-point set with its own seeds and
  requires excess at most 4.
* `k9/test_splits.py`: decodes the stored splits twice, verifies them from the definition,
  checks the doubled sets 211 (F_3^7) and 423 (F_3^8), lifts the 105-point split through
  211 → 422 → 844 → 1693 with excess 1 at every step, and checks the 1689-point set of F_3^10
  over all 29524 directions with its own pair checker (no code shared with `k9/split.py`).
* `k9/test_doubled_big.py`: checks the 3379-point set of F_3^11 over all 88573 directions,
  the 6759-point set of F_3^12 by the layer criterion over all 177147 vectors, the split of
  the 6759-point set over all 531441 vectors, and every number in `k9/results/doubling.json`.

## Prior work

We found no statement of these lemmas. Lev's 2009 pattern 3, 7, 13, 27, 53 (comment 994 on
T. Tao's blog), the classical construction {0,1}^n ∪ {0,2}^n and the layered criterion of
Saraf and Sudan ([arXiv:0808.2499](https://arxiv.org/abs/0808.2499)) are the closest. Kopparty,
Lev, Saraf and Sudan ([arXiv:1003.3736](https://arxiv.org/abs/1003.3736)) prove growth base
about 1.8899 by a random construction, so for large n sets of size c · 2^n are far from
optimal; the doubling is about explicit sets in moderate dimensions.
