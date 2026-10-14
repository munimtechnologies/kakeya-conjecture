# The Kakeya conjecture in October 2026

**What is proved, what changed in 2025-2026, and what is open.** This is an exposition of the
Kakeya set conjecture in Euclidean space and over finite fields, with complete proofs of the
elementary results (the planar Córdoba-type bound, Dvir's finite-field theorem, and several small
finite-field facts), a dated record of the dimension bounds, and a tested, dependency-free Python
companion (`kakeya.py`) that checks Kakeya sets in AG(*n*, *q*), builds the two explicit
constructions below, and computes every finite-field bound quoted here.

**New in version 1.2:** explicit Kakeya sets in **F**₃^*n* with base below 2. Cut a point of
**F**₃^3*m* into three blocks of length *m*; the set of points in which one block is a 0/1 vector
with at most *m*/3 ones and each of the other two blocks misses one of the three symbols is a
Kakeya set. This gives *K*(*n*, 3) ≤ 189 · 432^(*n*/9) ≈ 189 · 1.9626^*n*, where the explicit
constructions quoted in the literature have base 2 (section 7, last part; proofs in
[`TRIPLES.md`](TRIPLES.md), code and tests in `triples.py` and `test_triples.py`).

**New in version 1.1:** the first exact finite-field Kakeya number in dimension three over a field
with more than two elements. The smallest Kakeya set in **F**₃³ has exactly **13** points
(section 7), proved by an exhaustive search that runs in seconds with no solver. The release also
gives smaller Kakeya sets than the known constructions in **F**₃⁴, **F**₃⁵, **F**₄³ and **F**₄⁴,
with every set published and re-checked.

By MV2 of Munim, Inc. This article and the code were researched and written with AI agents working in
MV2's name. Every statement below is either proved here or cited to a source; section 10 separates
what is proved from what is open, and section 12 lists what is new here. Corrections are welcome by
issue or pull request.

Version 1.2, 2026-10-14. Text and data CC BY 4.0, code Apache-2.0.

Related problem record: Justin Sun Prize JSP-000015, "Kakeya dimension conjecture (general
dimension)" (MathDB #315983).

---

## 1. The question

A **Kakeya set** (or Besicovitch set) in **R**^*n* is a compact set containing a unit line segment
in every direction.

**Kakeya set conjecture.** Every Kakeya set in **R**^*n* has Hausdorff dimension *n* (and hence
also Minkowski dimension *n*).

The Justin Sun Prize catalogue states it as: "Must a set containing a unit line segment in every
direction have full Hausdorff dimension in its Euclidean ambient space?" Its title says "general
dimension", so the open content in October 2026 is *n* ≥ 4 (section 10).

History in three lines. Kakeya (1917) asked for the smallest area of a planar region in which a
unit needle can be turned around. Besicovitch (1919, published again in 1928) showed that a planar
set containing a unit segment in every direction can have measure zero, which answers the needle
problem in a surprising way: the area can be made arbitrarily small. Since such sets can be
null, the right question is how *large* they must be in the sense of dimension, and that is the
conjecture above. It sits underneath a chain of central problems in harmonic analysis: the
Kakeya *maximal function* conjecture, the Fourier restriction conjecture (JSP-001041 in the same
pilot) and the Bochner-Riesz conjecture. Bochner-Riesz implies restriction (Tao 1999), restriction
implies the Kakeya maximal conjecture, and that implies the set conjecture. Fefferman's 1971 disc multiplier theorem used Besicovitch's
construction to show that the ball multiplier is unbounded on *L*^*p*(**R**^*n*) for *p* ≠ 2,
*n* ≥ 2; that is where Kakeya sets entered Fourier analysis.

## 2. Two notions of size

For a bounded set *E* ⊂ **R**^*n* let *N*(*E*, δ) be the least number of δ-balls covering *E*.

- **Minkowski (box) dimension**: the lower and upper limits of log *N*(*E*, δ) / log(1/δ) as
  δ → 0. Equivalently, the δ-neighbourhood *E*_δ has volume about δ^(*n* − dim).
- **Hausdorff dimension**: the infimum of *s* for which *E* can be covered by balls of radii *r_i*
  with Σ *r_i*^*s* arbitrarily small. It is always at most the lower Minkowski dimension, so a
  Hausdorff lower bound is the stronger statement.

Both are computed through **δ-tubes**: a δ-tube is a cylinder of length 1 and radius δ around a
unit segment (in the plane, a 1 × δ rectangle). If *K*
is a Kakeya set, its δ-neighbourhood contains a δ-tube in every direction, so a lower bound on
|∪ *T*| for δ-tubes pointing in δ-separated directions gives a lower bound on |*K*_δ| and hence
on the Minkowski dimension. Hausdorff bounds need the same estimate when each tube is replaced by
a fixed fraction of itself (Zahl's survey, Conjecture 1.2(B)).

## 3. The plane, with a proof

**Theorem (Davies 1971).** Every Kakeya set in **R**² has Hausdorff dimension 2.

Here is the complete argument for the Minkowski version, which shows the mechanism that every later
proof tries to extend.

**Lemma 1.** Two δ-tubes in **R**² whose directions make an angle θ ∈ (0, π/2] intersect in a set
of area at most δ²/sin θ ≤ (π/2) δ²/θ.

*Proof.* Each tube lies in an infinite strip of width δ. Two strips of width δ at angle θ intersect
in a parallelogram whose sides have length δ/sin θ and whose height is δ, so its area is
δ²/sin θ. Finally sin θ ≥ 2θ/π on [0, π/2]. ∎

**Theorem 2.** Let *T*_1, ..., *T_N* be δ-tubes in **R**² with directions θ_*i* = *i*δ,
*N* = ⌊π/δ⌋. Then |∪ *T_i*| ≥ *N*δ / (1 + π(1 + ln *N*)), which is about 1/log(1/δ).

*Proof.* Let *f* = Σ 1_{*T_i*}. By Cauchy-Schwarz, (∫ *f*)² ≤ |supp *f*| · ∫ *f*², so
|∪ *T_i*| ≥ (Σ |*T_i*|)² / Σ_{*i*,*j*} |*T_i* ∩ *T_j*|. Here Σ |*T_i*| = *N*δ, since each
tube is a 1 × δ rectangle. For fixed *i* and *j* ≠ *i*, put *m* = min(|*i* − *j*|, *N* − |*i* − *j*|);
since *N*δ ≤ π, the angle between the two tubes (taken in (0, π/2]) is at least *m*δ, and each
value *m* = 1, ..., ⌊*N*/2⌋ occurs for at most two values of *j*. By Lemma 1,
Σ_*j* |*T_i* ∩ *T_j*| ≤ δ + 2 Σ_{*m* ≤ *N*/2} (π/2)δ²/(*m*δ) ≤ δ(1 + π(1 + ln *N*)).
Summing over *i* and dividing gives the bound. ∎

Since the δ-neighbourhood of a planar Kakeya set contains such a family (it contains a 1 × δ
rectangle around each of its unit segments), |*K*_δ| ≳ 1/log(1/δ), and
*N*(*K*, δ) ≳ δ^(−2)/log(1/δ): the Minkowski dimension is 2. The logarithm is a real loss, not an
artefact: Besicovitch's sets show that |*K*_δ| can tend to 0. Córdoba (1977) turned this
computation into the sharp *L*² bound for the planar Kakeya maximal function.

**Why the argument stops at *n* = 2.** In **R**^*n* the same computation only controls pairs of
tubes. Tubes can now meet in complicated configurations (many tubes through one point lying near a
plane, or arranged along a ruled surface), and pair counting alone gives only the dimension
bound (*n* + 1)/2 (Drury 1983; Christ, Duoandikoetxea and Rubio de Francia 1986). Every later
record comes from a new way to rule out such configurations.

## 4. The record in R^n, dated

| Year | Who | Result |
| --- | --- | --- |
| 1971 | Davies | dimension 2 in **R**² |
| 1977 | Córdoba | sharp Kakeya maximal estimate in **R**² |
| 1983, 1986 | Drury; Christ-Duoandikoetxea-Rubio de Francia | (*n* + 1)/2 in **R**^*n* |
| 1991 | Bourgain | first improvement on (*n* + 1)/2 ("bush" argument) |
| 1995 | Wolff | (*n* + 2)/2 in **R**^*n* (hairbrush; 5/2 in **R**³, 3 in **R**⁴) |
| 1999 | Bourgain | (13*n* + 12)/25 via additive combinatorics ("arithmetic" method) |
| 1999 | Katz-Tao | Minkowski (4*n* + 3)/7, Hausdorff (6*n* + 5)/11 |
| 2000 | Katz-Łaba-Tao | upper Minkowski 5/2 + ε in **R**³ |
| 2001 | Łaba-Tao | upper Minkowski (*n* + 2)/2 + ε_*n* for *n* ≥ 4 |
| 2002 | Katz-Tao | Minkowski *n*/α + (α − 1)/α with α³ − 4α + 2 = 0 (α = 1.675...); Hausdorff (2 − √2)(*n* − 4) + 3 |
| 2018 | Guth-Zahl, with Zahl and Katz-Rogers (polynomial Wolff axioms) | 3 + 1/40 in **R**⁴ |
| 2019 | Katz-Zahl | Hausdorff 5/2 + ε in **R**³ |
| 2021 | Katz-Zahl (planebrush) | 3.059 in **R**⁴ |
| 2019-2021 | Zahl; Hickman-Rogers-Zhang | polynomial-partitioning ("broad") estimates; new bounds in some dimensions *n* ≥ 7 |
| 2025 | **Wang-Zahl** | **dimension 3 in R³: the conjecture holds for *n* = 3** |
| 2026 | Rai Choudhuri (preprint, not refereed) | 13/4 − (12951 − 54√10593)/41972 = 3.0738... in **R**⁴ |

The 2002 Katz-Tao bounds are quoted from Theorem 1.1 of their paper; the paper notes that its
Hausdorff bound is better than its Minkowski bound for *n* ≤ 23. The **R**⁴ value 3.0738 is
Theorem 1 of arXiv 2610.12480 (submitted 30 August 2026), which combines a slightly improved
Guth-Zahl trilinear estimate, obtained by applying Wang-Zahl inside thickened neighbourhoods of
varieties, with the Katz-Zahl planebrush. The author says some AI assistance was used. We give its
number exactly as stated and have not checked the proof.

For *n* ≥ 5 we do not attempt a per-dimension table: the best bounds come from a mix of the
Katz-Tao arithmetic bounds and the polynomial-method estimates, and Zahl's 2025 survey (ICM 2026)
does not tabulate them. All known bounds for large *n* grow like *cn* with *c* well below 1 (for
example *n*/1.675 for the Minkowski dimension), so the gap to *n* widens with the dimension.

## 5. What Wang and Zahl proved, and how

**Theorem (Wang-Zahl 2025).** Every Kakeya set in **R**³ has Hausdorff dimension 3 and Minkowski
dimension 3.

The proof is spread over three papers: the *sticky* case ("Sticky Kakeya sets and the sticky
Kakeya conjecture", J. Amer. Math. Soc. 39 (2026) 515-585), an Assouad-dimension paper
(Invent. Math. 241 (2025)), and the general case, arXiv 2502.17655 (February 2025), which reduces
the general case to the sticky one by induction on scales. Guth, Wang and Zahl posted a streamlined
proof in January 2026 (arXiv 2601.14411); Guth's Séminaire Bourbaki survey appeared in April 2026
(arXiv 2604.03416), and Alex Cohen circulated detailed checking notes the same month. Hong Wang
received a 2026 Fields Medal at the International Congress in Philadelphia on 23 July 2026, with
this work cited. As of the sources we read, the main 2025 paper is still listed as a preprint.

The shape of the argument, following Guth's survey:

1. **Stickiness.** A family of δ-tubes is *sticky* if, at every coarser scale ρ, the tubes group
   into about ρ^(1−*n*) fat ρ-tubes, the minimum possible. Sticky families are the hardest case.
2. **Grains.** In a worst-case sticky family the union locally looks like a product of a small
   square with a fractal set *A*: the tubes stack into "grains" (Katz-Łaba-Tao).
3. **The obstruction is complex conjugation.** Following grains around a closed cycle forces an
   approximate identity *A* + *s*(*x*)*y* + *s̃*(*y*)*x* ≈ *A*. The model solution is
   **R** + *x̄y* + *ȳx* = **R** inside **C**, and it really happens: the Heisenberg group in **C**³
   gives a tube arrangement over **C** that is sticky and grainy with the 5/2 behaviour. This is
   why every method that ignores the order structure of **R** stalls at 5/2.
4. **Sum-product kills it over R.** Bourgain's discretised sum-product theorem says **R** has no
   such approximate subrings, which rules out the conjugation structure; repeating a small gain
   across scales gives the sharp result.
5. **Reduction to sticky.** In a non-sticky family, the densest convex sets containing many tubes
   are handled by induction on scales (the "convex Wolff axioms" theorem, Guth's Theorem 1.2).

What it does **not** prove: the Kakeya *maximal function* conjecture in **R**³ (Wang-Zahl say so
explicitly, and Zahl's survey lists obstructions to adapting the argument; a special case for SL₂
lines was proved by Green, Harris and Ou, Analysis & PDE 19 (2026)), restriction in **R**³, or
anything in **R**⁴ beyond the improved bounds of section 4.

## 6. The finite-field problem, with proofs

Wolff (1999) proposed a model problem. A **Kakeya set** in **F**_*q*^*n* is a set *K* containing a
full line {*b* + *td* : *t* ∈ **F**_*q*} for every direction *d* ≠ 0. Write *K*(*n*, *q*) for the
minimum size. Wolff asked whether |*K*| ≥ *c_n q^n*.

**Theorem (Dvir 2008).** Every Kakeya set in **F**_*q*^*n* has at least C(*q* + *n* − 1, *n*)
points. In particular |*K*| ≥ *q^n*/*n*!.

*Proof.* The polynomials in *n* variables of degree at most *q* − 1 form a space of dimension
C(*q* − 1 + *n*, *n*). If |*K*| were smaller, a nonzero such *P* would vanish on *K* (more unknowns
than linear conditions). Let *d* = deg *P*, so 1 ≤ *d* ≤ *q* − 1 (*d* = 0 would make *P* a nonzero
constant vanishing on the nonempty set *K*), and let *P_d* be its top-degree homogeneous part. For
a direction *v*, *K* contains a line *b* + *tv*, so the one-variable polynomial *t* ↦ *P*(*b* + *tv*)
has degree at most *d* < *q* and *q* roots; it is identically zero. Its coefficient of *t^d* is
*P_d*(*v*), so *P_d*(*v*) = 0 for every *v* ≠ 0, and *P_d*(0) = 0 because *d* ≥ 1. So *P_d* is a
nonzero polynomial of degree less than *q* in each variable that vanishes on all of **F**_*q*^*n*,
which is impossible (induct on *n*: a nonzero one-variable polynomial of degree < *q* has fewer
than *q* roots). ∎

The proof is ten lines. Nothing like it is known in **R**^*n*, where δ-tubes are not algebraic
objects, so the finite-field model turned out to be much softer than the real problem. The sharp
constant came later:

- Dvir, Kopparty, Saraf and Sudan (2009, method of multiplicities): |*K*| ≥ *q^n*/2^*n*.
- **Bukh and Chao (Discrete Analysis 2021):** |*K*| ≥ *q^n* / (2 − 1/*q*)^(*n*−1), so the density
  is at least 2^(1−*n*) asymptotically. This is sharp: constructions with density
  2^(1−*n*) + *O*(1/*q*) exist (Dvir-Kopparty-Saraf-Sudan, Saraf-Sudan; the parabola set below).

So for fixed *n* the first-order behaviour is settled: *K*(*n*, *q*) = (2^(1−*n*) + *o*(1)) *q^n*
as *q* → ∞ (through odd *q* by the parabola set; Kyureghyan, Müller and Wang give constructions of
the same density for even *q*).

**The plane is solved exactly.** Blokhuis and Mazzocca (2008) proved
*K*(2, *q*) = *q*(*q* + 1)/2 + (*q* − 1)/2 for odd *q*; for even *q* the minimum is *q*(*q* + 1)/2
(the smallest examples are classified: Blokhuis, De Boeck, Mazzocca and Storme 2014). Dover and
Mellinger found by exhaustive search every attainable size of a Kakeya set in the Desarguesian
planes of order *q* ≤ 9. For even *q*, Dvir's bound C(*q* + 1, 2) = *q*(*q* + 1)/2 is therefore
exact in the plane.

Four further elementary facts, each checked by the companion code:

**(a) The counting bound.** *K*(*n*, *q*) ≥ the least *s* with *s*(*s* − 1) ≥ *q*(*q^n* − 1).
*Proof.* *K* contains (*q^n* − 1)/(*q* − 1) lines in distinct directions; two distinct lines share
at most one point, so the C(*q*, 2) pairs of points on each line are disjoint sets of pairs inside
the C(*s*, 2) pairs of *K*. ∎ It is weak for large *q* but exact for *q* = 2 and *n* ≤ 4.

**(b) Products.** If *A* ⊂ **F**_*q*^*n* and *B* ⊂ **F**_*q*^*m* are Kakeya sets, so is
*A* × *B*, hence *K*(*n* + *m*, *q*) ≤ *K*(*n*, *q*) *K*(*m*, *q*). *Proof.* For a direction
(*d*, *e*) with *d*, *e* ≠ 0, take lines *a* + *td* ⊂ *A* and *b* + *te* ⊂ *B*; then
(*a* + *td*, *b* + *te*) ⊂ *A* × *B*. If *e* = 0 use (*a* + *td*, *b*) for any *b* ∈ *B*, and
symmetrically if *d* = 0. ∎

**(c) Over F₂ a Kakeya set is a difference basis.** A line in **F**₂^*n* is a pair {*x*, *x* + *d*},
so *K* is Kakeya exactly when every nonzero vector is a sum *a* + *b* with *a*, *b* ∈ *K*: *K* is a
**difference basis** of the group *C*₂^*n*. Banakh and Gavrylkiv computed the difference sizes of all
abelian groups of order below 96, which gives

| *n* | 1 | 2 | 3 | 4 | 5 | 6 |
| --- | --- | --- | --- | --- | --- | --- |
| *K*(*n*, 2) | 2 | 3 | 5 | 6 | 10 | 14 |
| counting bound (a) | 2 | 3 | 5 | 6 | 9 | 12 |

(Acta Sci. Math. (Szeged) 85 (2019), Tables 1 and 3, where the group is written (*C*₂)^*n*). For
*n* = 4 the optimum is a perfect Sidon set: {0, *e*₁, *e*₂, *e*₃, *e*₄, *e*₁ + *e*₂ + *e*₃ + *e*₄}
has 15 pairwise sums, and they are exactly the 15 nonzero vectors. The companion code reproves
*n* ≤ 4 by exhaustion. Their Proposition 7.3 gives, for *n* ≥ 2, the two-sided bound
√2 < *K*(*n*, 2)/2^(*n*/2) < 3/√2 (and < 2 for even *n*).

**(d) Two explicit constructions.**

- *The parabola set* (*q* odd). Write points as (*x*, *t*) with *x* ∈ **F**^(*n*−1). For each
  *b* ∈ **F**^(*n*−1) take the line {((*tb_i* − *b_i*²)_*i*, *t*) : *t* ∈ **F**}, of direction (*b*, 1),
  and add the hyperplane *t* = 0 for the directions (*b*, 0). For *t* ≠ 0 the slice at height *t* is
  {*x* : *t*²/4 − *x_i* is a square for every *i*}, because *t*²/4 − (*tb* − *b*²) = (*b* − *t*/2)².
  So the set has (*q* − 1)((*q* + 1)/2)^(*n*−1) + *q*^(*n*−1) points, density 2^(1−*n*) + *O*(1/*q*),
  and for *n* = 2 it has (*q*² + 2*q* − 1)/2 points: exactly the Blokhuis-Mazzocca minimum.
- *The Smyth set.* *K* = (**F**^*)^*n* ∪ {0, 1}^*n*, with (*q* − 1)^*n* + 2^*n* − 1 points.
  This size was reported by Clifford Smyth in a 2008 comment on Tao's blog; here is a proof that
  this set works. Given a direction *d*, let *Z* be the set of coordinates where *d* is zero (not all
  of them) and 1_*Z* the 0/1 vector supported on *Z*. On the line 1_*Z* + *td* the point *t* = 0
  is 1_*Z* ∈ {0, 1}^*n*, and for *t* ≠ 0 every coordinate is nonzero: 1 on *Z*, *td_i* ≠ 0 off *Z*. ∎
  For *q* = 3 it gives *K*(*n*, 3) ≤ 2^(*n*+1) − 1, which beats the parabola set's
  2^*n* + 3^(*n*−1) for every *n* ≥ 3, and no 0/1 vector other than (1, ..., 1) can be dropped
  from it (checked by the tests for *n* = 3).

## 7. Small cases: what is known exactly

| (*n*, *q*) | points | best lower bound (which) | best verified upper bound (which) |
| --- | --- | --- | --- |
| (2, *q*) | *q*² | exact: Blokhuis-Mazzocca | exact |
| (*n*, 2), *n* ≤ 6 | 2^*n* | exact: Banakh-Gavrylkiv | exact |
| (3, 3) | 27 | **13 (new, this repository)** | **13 (new)**; Smyth 15 |
| (4, 3) | 81 | 18 (Bukh-Chao) | **27 (new)**; Smyth 31 |
| (5, 3) | 243 | 32 (Bukh-Chao) | **56 (new)**; Smyth 63 |
| (3, 4) | 64 | 21 (Bukh-Chao) | **28 (new)**; Smyth 34 |
| (4, 4) | 256 | 48 (Bukh-Chao) | **64 (new)**; Smyth 96 |
| (3, 5) | 125 | 39 (Bukh-Chao) | 53 (AlphaEvolve formula; found again here); parabola 61 |
| (3, 7) | 343 | 100 (Bukh-Chao) | 128 (found here); parabola 145 |

Every number in the table is recomputed by `python kakeya.py bounds N Q`, every upper bound is a
set that the tests build or load and check line by line, and the lower bound 13 is re-proved by the
tests (below). The 2025 AlphaEvolve paper (Georgiev, Gómez-Serrano, Tao and Wagner,
arXiv 2511.02864) gives constructions in **F**_*p*³ for primes *p* ≡ 1 (mod 4) of size
*p*³/4 + 7*p*²/8 − 1/8, which is 53 at *p* = 5; our solver found a 53-point set in **F**₅³
independently, and did not beat it in 25 minutes.

**Theorem (new). *K*(3, 3) = 13: the smallest Kakeya set in F₃³ has exactly 13 points.**

*Upper bound.* The set

    (0,0,0) (1,0,0) (2,0,0) (1,1,0) (0,0,1) (2,0,1) (0,1,1)
    (1,1,1) (0,2,1) (2,2,1) (0,0,2) (2,1,2) (0,2,2)

contains a line in each of the 13 directions of **F**₃³ (over **F**₃ a line is a 3-term
progression *b*, *b* + *d*, *b* + 2*d*):

| direction | line in the set | direction | line in the set |
| --- | --- | --- | --- |
| (1,0,0) | (0,0,0) (1,0,0) (2,0,0) | (1,1,0) | (0,0,1) (1,1,1) (2,2,1) |
| (0,1,0) | (0,0,1) (0,1,1) (0,2,1) | (1,2,0) | (0,2,1) (1,1,1) (2,0,1) |
| (0,0,1) | (0,0,0) (0,0,1) (0,0,2) | (1,1,1) | (0,0,2) (1,1,0) (2,2,1) |
| (0,1,1) | (0,0,0) (0,1,1) (0,2,2) | (1,1,2) | (0,2,1) (1,0,0) (2,1,2) |
| (0,1,2) | (2,0,0) (2,1,2) (2,2,1) | (1,2,1) | (0,2,2) (1,1,0) (2,0,1) |
| (1,0,1) | (0,0,2) (1,0,0) (2,0,1) | (1,2,2) | (0,2,2) (1,1,1) (2,0,0) |
| (1,0,2) | (0,1,1) (1,1,0) (2,1,2) | | |

*Lower bound.* Suppose *K* is a Kakeya set with at most 12 points. Translating *K* so that its line
in direction *e*₁ becomes the axis {(0,0,0), (1,0,0), (2,0,0)} keeps it a Kakeya set of the same
size, and adding points keeps it Kakeya, so there would be a 12-point Kakeya set containing the
axis. There are C(24, 9) = 1,307,504 such 12-point sets, and none of them contains a line in all
13 directions. `kakeya_sets_through_axis` checks every one with bitmasks, with no solver, in a few
seconds; `test_k33_is_thirteen` runs it on every test run. Separately, two independent
integer-programming solvers (OR-Tools CP-SAT, with and without symmetry fixing, and HiGHS through
SciPy) each proved 13 optimal in MV2's run `kakeya-min-v1`. ∎

So the sequence *K*(*n*, 3) begins 3, 7, 13. The minimum set is not a "lift" of a planar one in
the obvious way: its slices by the last coordinate have 4, 6 and 3 points, and none of them is a
planar Kakeya set (that needs 7 points).

**Why this is new.** Lund, Saraf and Wolf (arXiv 1609.01048) wrote that no exact value is known in
any dimension other than 2. The *q* = 2 values follow from Banakh and Gavrylkiv's difference bases
(section 6(c)); for *q* ≥ 3 and *n* ≥ 3 we found no exact value in the papers cited here, in the
AlphaEvolve paper (which treats asymptotic constructions), on OEIS, or in answer to the 2008
question on Tao's blog about the smallest subset of **F**₃^*n* with a 3-term progression in every
direction. If an earlier computation exists, we would like to cite it.

**New upper bounds.** The same run found Kakeya sets of 27 points in **F**₃⁴, 56 in **F**₃⁵,
28 in **F**₄³ and 64 in **F**₄⁴. Each is smaller than every explicit construction quoted in this
article: the Smyth set (31, 63, 34, 96), the parabola set (43 and 113 for *q* = 3), products of
smaller optimal sets (for example 3 · 13 = 39 and 7 · 7 = 49 in **F**₃⁴, 10 · 10 = 100 in **F**₄⁴),
and the Kyureghyan-Müller-Wang construction for *q* = 4 (46 in **F**₄³ by their formula (2)).
They are not claimed optimal: the solvers stopped at their time limits. CP-SAT reported a lower
bound of 24 for **F**₄³; we do not use it, because it was not certified independently. The
27-point set in **F**₃⁴ has slices of 13, 6 and 8 points by its last coordinate, which suggests a
construction lifting **F**₃³ sets to **F**₃⁴; MV2 is testing that next.

All sets, solver logs and timings are in `data/kakeya-min-v1/` (`results.json`: for every case the
runs of each solver and a `certificate`, the smallest set found as a list of coordinates;
`results.txt`: one line per case). `solve.py` reproduces the run
(`pip install -r requirements-solve.txt`, then `python solve.py --out-dir out`; about 3.3 hours on
4 cores).

**Open, as far as we can tell:** every other exact value with *n* ≥ 3 and *q* ≥ 3. The smallest
open cases are *K*(4, 3) ∈ [18, 27] and *K*(3, 4) ∈ [21, 28].

**Growth over F₃: explicit sets below base 2 (new in version 1.2).** For fixed *q* = 3 and
growing *n*, the limit β = lim *K*(*n*, 3)^(1/*n*) exists, and Kopparty, Lev, Saraf and Sudan
(arXiv 1003.3736) showed 1.8 ≤ β ≤ 3/2^(2/3) ≈ 1.890. Their upper bound uses randomly chosen
linear maps; the explicit constructions they and later papers give for *q* = 3 have base 2 (the
Smyth set above, 2^(*n*+1) − 1, and 3 · 2^*n*). Here is an explicit one with base below 2.

Let *W* ⊂ **F**₃^*m* be the vectors that use at most two of the three symbols, and *B* the 0/1
vectors with at most *m*/3 ones.

**Theorem (new).** The set *K_m* = (*W* × *W* × *B*) ∪ (*W* × *B* × *W*) ∪ (*B* × *W* × *W*) ⊂ **F**₃^3*m*
is a Kakeya set. Consequently *K*(*n*, 3) ≤ 189 · 432^(*n*/9) for every *n*, and
432^(1/9) = 2^(4/9) · 3^(1/3) ≈ 1.96262.

*Proof.* First, for every *d* ∈ **F**₃^*m* there is *a* with *a* ∈ *W*, *a* + *d* ∈ *W* and
*a* + 2*d* ∈ *B*: let *j* be a least frequent symbol of *d* (it occurs at most *m*/3 times), let
*e* be the 0/1 vector with ones where *d_i* = *j*, and put *a* = *d* + *e*. Then *a* uses only
*j* + 1 and *j* + 2; *a* + *d* = 2*d* + *e* uses only 2*j* + 1 and 2*j* + 2; and
*a* + 2*d* = *e* ∈ *B*. Now take a direction *d* = (*d*₁, *d*₂, *d*₃) of **F**₃^3*m*, let *a_i*
be as above for *d_i*, and put *x* = (*a*₁, *a*₂ + *d*₂, *a*₃ + 2*d*₃). Then *x* ∈ *W* × *W* × *B*,
*x* + *d* = (*a*₁ + *d*₁, *a*₂ + 2*d*₂, *a*₃) ∈ *W* × *B* × *W* and
*x* + 2*d* = (*a*₁ + 2*d*₁, *a*₂, *a*₃ + *d*₃) ∈ *B* × *W* × *W*. For the size, |*W*| = 3 · 2^*m* − 3
and |*B*| ≤ 2^(*m* h(1/3)) = (27/4)^(*m*/3) (the binomial tail bound, h the binary entropy), so
|*K_m*| ≤ 3 · 9 · 4^*m* · (27/4)^(*m*/3) = 27 · 432^(*m*/3). Other *n* follow by products with
*K*(1, 3) = 3 and *K*(2, 3) = 7. ∎

The exact size is |*K_m*| = 3*w*²*b* − 3*wb*² + *b*³ with *w* = |*W*|, *b* = |*B*|, because *B* lies
inside *W*. It first beats the Smyth set at *m* = 17 (*n* = 51), and the clean bound beats it for
every *n* = 3*m* ≥ 138. The construction does **not** improve the probabilistic 1.890; it is, as far
as we found, the first explicit one below 2. [`TRIPLES.md`](TRIPLES.md) puts it in a general
framework ("rainbow triples": the growth base equals the best base of a single small triple), with
the tripling lemma, an energy lemma, the exact plane values and the open targets in **F**₃³.

## 8. Variants that are open or recently settled

- **Z/NZ and p-adic rings.** Over the rings **Z**/*N***Z** the analogue of Dvir's theorem was
  proved by Dhar (Advances in Combinatorics 2024, general *N*), and Arsovski proved the *p*-adic
  Kakeya conjecture (J. Amer. Math. Soc. 37 (2024)).
- **Nikodym and Furstenberg sets.** The Furstenberg set conjecture in the plane was settled by Ren
  and Wang (2024), building on Orponen-Shmerkin; Guth cites it as a hint that sticky cases are key.
- **Plany sets in F_q⁴.** Łaba, Rai Choudhuri and Zahl (Res. Math. Sci. 13 (2026)) give a
  finite-field exposition of the Katz-Zahl planebrush.

## 9. Claimed proofs

Short claimed proofs of the general conjecture appear from time to time. A claim deserves attention
only if it handles the Heisenberg obstruction of section 5 (any argument that works verbatim over
**C** proves something false) and recovers the Wang-Zahl theorem as a special case. We do not list
individual unrefereed claims.

## 10. Proved versus open, as of 2026-10-14

**Proved.**
- *n* = 2: every Kakeya set has dimension 2 (Davies 1971); the maximal conjecture holds (Córdoba 1977).
- *n* = 3: every Kakeya set has Hausdorff and Minkowski dimension 3 (Wang-Zahl 2025; three papers,
  one published in JAMS 2026, the main one still a preprint; a streamlined proof exists).
- All *n*: Hausdorff dimension ≥ max((*n* + 2)/2, (2 − √2)(*n* − 4) + 3) and the Katz-Tao Minkowski
  bound; in **R**⁴, at least 3.059 (Katz-Zahl, published) and 3.0738 (2026 preprint).
- Finite fields: |*K*| ≥ *q^n*/(2 − 1/*q*)^(*n*−1), sharp to first order; the plane exactly; *q* = 2
  exactly for *n* ≤ 6; *K*(3, 3) = 13 (new here, by exhaustive computation that anyone can rerun);
  for *q* = 3, 1.8 ≤ β ≤ 1.890 (KLSS, the upper bound probabilistic) and explicit sets with
  *K*(*n*, 3) ≤ 189 · 1.9627^*n* (new here, section 7).

**Open.**
- The Kakeya set conjecture for every *n* ≥ 4 (this is the open content of JSP-000015).
- The Kakeya maximal function conjecture for every *n* ≥ 3.
- The restriction and Bochner-Riesz conjectures for *n* ≥ 3.
- Exact finite-field values *K*(*n*, *q*) for *n* ≥ 3, *q* ≥ 3 other than *K*(3, 3); the
  second-order term in *K*(*n*, *q*).
- The growth base β for *q* = 3, and an explicit construction matching the random 1.890.

## 11. Why this matters beyond the problem

Kakeya estimates control how much wave packets pointing in different directions can overlap, which
is the core difficulty in restriction theory, Bochner-Riesz summation, local smoothing for the wave
equation and decoupling. Dvir's finite-field theorem became a standard tool in theoretical computer
science (randomness mergers and extractors) and set off the modern wave of the polynomial method in
combinatorics, which went on to the joints theorem and the cap-set bound.

## 12. What is original here

- **Explicit Kakeya sets in F₃^*n* with base below 2** (section 7, last part): *K*(*n*, 3) ≤
  189 · 432^(*n*/9), base 1.9626, from a two-symbol "rainbow triple" and a tripling lemma
  (`TRIPLES.md`). As far as we found, the explicit constructions in the literature for *q* = 3
  have base 2, and the best bound, 1.890, is probabilistic.
- **The exact value *K*(3, 3) = 13** (section 7), with a solver-free proof that reruns in seconds,
  an explicit 13-point set and a table of its 13 lines. As far as we found, it is the first exact
  finite-field Kakeya number in dimension at least 3 over a field with more than two elements.
- **Smaller Kakeya sets** in **F**₃⁴ (27), **F**₃⁵ (56), **F**₄³ (28) and **F**₄⁴ (64) than the
  constructions we found in the literature, all published as coordinate lists and re-checked.
- A dated, sourced record of the dimension bounds through October 2026, including the August 2026
  **R**⁴ preprint and the 2026 publication and recognition of the Wang-Zahl work.
- Self-contained proofs of the planar bound, Dvir's theorem, the counting bound, the product lemma,
  the parabola construction and its exact size, and the Smyth construction (whose proof we did not
  find written out).
- The identification *K*(*n*, 2) = Δ[*C*₂^*n*] with difference bases, which imports exact values
  for *n* ≤ 6 from group theory. We did not find this identification stated in the Kakeya
  literature we read; it is elementary.
- The small-case table of section 7 with every bound recomputed and every upper bound verified by
  machine, and a precise statement of which exact values are open.

## 13. Using the companion code

`kakeya.py` needs only Python 3.8+ and no packages.

    python kakeya.py bounds 3 3       # counting, Dvir, Bukh-Chao, parabola and Smyth sizes
    python kakeya.py smyth 3 3        # print the 15-point Kakeya set in F_3^3
    python kakeya.py parabola 3 5     # print the 61-point Kakeya set in F_5^3
    python kakeya.py check 3 3 set.txt  # check a set: one point per line, coordinates separated by spaces
    python kakeya.py brute 2 5        # exact minimum by brute force (small planes only)
    python kakeya.py exhaust 3 3 12   # no Kakeya set of 12 points in F_3^3 (a few seconds)
    python kakeya.py exhaust 3 3 13   # print a 13-point one
    python kakeya.py records          # the exact values and bounds of section 7
    python -m unittest -v test_kakeya

`triples.py` (Python 3.8+, no packages) holds the rainbow-triple tools of `TRIPLES.md`:

    python triples.py wwb 40          # sizes and bases of the two-symbol construction, m = 1..40
    python triples.py plane           # the exhaustive F_3^2 table behind tau(2) = 64
    python triples.py energy 3 512    # size patterns in F_3^3 that the energy lemma leaves open
    python -m unittest -v test_triples

Fields supported: **F**_*p* for primes *p*, and **F**₄, **F**₈, **F**₉, **F**₁₆. The tests check the
field axioms, the line structure of AG(*n*, *q*), the plane minima for *q* ≤ 5 by brute force and
by the exhaustive search, *K*(3, 2) = 5 and *K*(4, 2) = 6 by exhaustion, the equivalence (c) on all
256 subsets of **F**₂³, both constructions in every case used in section 7 (line by line, with
their exact sizes), the product lemma, the exhaustive search against a slow independent check,
*K*(3, 3) = 13, and every one of the 18 solver certificates in `data/kakeya-min-v1/results.json`.
`test_triples` checks the two-symbol theorem exhaustively for *m* ≤ 5, the Kakeya property of
*K_m* line by line for *n* = 3*m* ≤ 9 and on random directions up to *n* = 60, and every size
quoted above.

`solve.py` (needs `ortools`, `scipy` and `numpy`, see `requirements-solve.txt`) is the
integer-programming search that produced the data: one 0/1 variable per point and per line,
exactly one chosen line per direction, every point of a chosen line in the set, minimise the size;
solved with CP-SAT and HiGHS, and every set returned is re-checked with `is_kakeya`.
`search_s3.py` and `anneal.py` (CP-SAT and simulated annealing) search for small rainbow triples in
**F**₃³ and beyond.

## Sources

- Z. Dvir, "On the size of Kakeya sets in finite fields", J. Amer. Math. Soc. 22 (2009) 1093-1097; arXiv 0803.2336.
- Z. Dvir, S. Kopparty, S. Saraf, M. Sudan, "Extensions to the method of multiplicities, with applications to Kakeya sets and mergers", SIAM J. Comput. 42 (2013).
- S. Kopparty, V. F. Lev, S. Saraf, M. Sudan, "Kakeya-type sets in finite vector spaces", J. Algebraic Combin. 34 (2011) 337-355; arXiv 1003.3736.
- G. Kyureghyan, P. Müller, Q. Wang, "On the size of Kakeya sets in finite vector spaces", Electron. J. Combin. 20(3) (2013) #P36; arXiv 1302.5591.
- B. Bukh, T.-W. Chao, "Sharp density bounds on the finite field Kakeya problem", Discrete Analysis 2021:26; arXiv 2108.00074.
- A. Blokhuis, F. Mazzocca, "The finite field Kakeya problem", in Building Bridges, Bolyai Soc. Math. Stud. 19 (2008); arXiv 0911.4370.
- A. Blokhuis, M. De Boeck, F. Mazzocca, L. Storme, "The Kakeya problem: a gap in the spectrum and classification of the smallest examples", Des. Codes Cryptogr. 72 (2014) 21-31.
- B. Lund, S. Saraf, C. Wolf, "Finite field Kakeya and Nikodym sets in three dimensions", arXiv 1609.01048.
- T. Banakh, V. Gavrylkiv, "Difference bases in finite Abelian groups", Acta Sci. Math. (Szeged) 85 (2019) 119-137; arXiv 1704.02471.
- B. Georgiev, J. Gómez-Serrano, T. Tao, A. Z. Wagner, "Mathematical exploration and discovery at scale", arXiv 2511.02864 (2025).
- T. Tao, "Dvir's proof of the finite field Kakeya conjecture", blog post, 2008-03-24, https://terrytao.wordpress.com/2008/03/24/dvirs-proof-of-the-finite-field-kakeya-conjecture/ (comments include the Smyth construction and the F_3 question).
- L. Perron, F. Didier, "CP-SAT", Google OR-Tools, https://developers.google.com/optimization; Q. Huangfu, J. A. J. Hall, "Parallelizing the dual revised simplex method", Math. Prog. Comp. 10 (2018) (HiGHS), used through scipy.optimize.milp.
- R. O. Davies, "Some remarks on the Kakeya problem", Proc. Cambridge Philos. Soc. 69 (1971) 417-421.
- A. Córdoba, "The Kakeya maximal function and the spherical summation multipliers", Amer. J. Math. 99 (1977).
- S. Drury, "L^p estimates for the X-ray transform", Illinois J. Math. 27 (1983).
- M. Christ, J. Duoandikoetxea, J. L. Rubio de Francia, Duke Math. J. 53 (1986).
- J. Bourgain, "Besicovitch type maximal operators and applications to Fourier analysis", GAFA 1 (1991); "On the dimension of Kakeya sets and related maximal inequalities", GAFA 9 (1999).
- T. Wolff, "An improved bound for Kakeya type maximal functions", Rev. Mat. Iberoam. 11 (1995) 651-674.
- N. Katz, I. Łaba, T. Tao, "An improved bound on the Minkowski dimension of Besicovitch sets in R^3", Ann. of Math. 152 (2000) 383-446.
- I. Łaba, T. Tao, "An improved bound for the Minkowski dimension of Besicovitch sets in medium dimension", GAFA 11 (2001); "An X-ray transform estimate in R^n", Rev. Mat. Iberoam. 17 (2001).
- N. Katz, T. Tao, "Bounds on arithmetic projections, and applications to the Kakeya conjecture", Math. Res. Lett. 6 (1999).
- T. Tao, "The Bochner-Riesz conjecture implies the restriction conjecture", Duke Math. J. 96 (1999).
- C. Fefferman, "The multiplier problem for the ball", Ann. of Math. 94 (1971).
- N. Katz, T. Tao, "New bounds for Kakeya problems", J. Anal. Math. 87 (2002) 231-263; arXiv math/0102135.
- L. Guth, J. Zahl, "Polynomial Wolff axioms and Kakeya-type estimates in R^4", Proc. London Math. Soc. 117 (2018) 192-220; N. Katz, K. Rogers, "On the polynomial Wolff axioms", GAFA 28 (2018); J. Zahl, GAFA 28 (2018).
- N. Katz, J. Zahl, "An improved bound on the Hausdorff dimension of Besicovitch sets in R^3", J. Amer. Math. Soc. 32 (2019); arXiv 1704.07210. "A Kakeya maximal function estimate in four dimensions using planebrushes", Rev. Mat. Iberoam. 37 (2021) 317-359.
- J. Hickman, K. Rogers, R. Zhang, "Improved bounds for the Kakeya maximal conjecture in higher dimensions", Amer. J. Math. (2022).
- H. Wang, J. Zahl, "Sticky Kakeya sets and the sticky Kakeya conjecture", J. Amer. Math. Soc. 39 (2026) 515-585; "The Assouad dimension of Kakeya sets in R^3", Invent. Math. 241 (2025) 153-206; "Volume estimates for unions of convex sets, and the Kakeya set conjecture in three dimensions", arXiv 2502.17655.
- L. Guth, H. Wang, J. Zahl, "A streamlined proof of the Kakeya set conjecture in R^3", arXiv 2601.14411 (2026).
- L. Guth, "The Kakeya conjecture, after Wang and Zahl", Séminaire Bourbaki, arXiv 2604.03416 (2026).
- J. Zahl, "A survey of the Kakeya conjecture, 2000-2025", Proc. ICM 2026, vol. 4, 270-287; arXiv 2512.09397.
- M. Rai Choudhuri, "New bound on the Hausdorff dimension of Kakeya sets in R^4", arXiv 2610.12480 (2026, preprint).
- I. Łaba, M. Rai Choudhuri, J. Zahl, "A bound for plany Kakeya sets in F_q^4 using the planebrush method", Res. Math. Sci. 13 (2026); arXiv 2507.09605.
- B. Green, T. Harris, Y. Ou, "An L^{3/2} SL2 Kakeya maximal inequality", Analysis & PDE 19 (2026) 1295-1338.
- M. Dhar, "The Kakeya set conjecture for Z/NZ for general N", Advances in Combinatorics 2024:2; B. Arsovski, "The p-adic Kakeya conjecture", J. Amer. Math. Soc. 37 (2024).
- Fields Medal 2026: NYU news, July 2026, https://www.nyu.edu/about/news-publications/news/2026/july/nyu-professor-hong-wang-wins-fields-medal.html
- Justin Sun Prize catalogue, JSP-000015: https://github.com/TheJustinSunPrize/awards (problems/catalog-0001-0100.md, problems/mathdb-pilot.md).

## Licence

Text and data (`data/`): CC BY 4.0. Code: Apache-2.0 (see LICENSE). Cite as: MV2 (Munim, Inc.),
"The Kakeya conjecture in October 2026", version 1.2,
https://github.com/munimtechnologies/kakeya-conjecture.
