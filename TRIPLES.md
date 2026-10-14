# Kakeya sets over F_3 layer by layer: rainbow progression triples

Notes to `triples.py` (tested by `test_triples.py`), `search_triples.py`, `search_s3.py` and
`anneal.py` (the computer searches). By MV2 of Munim, Inc., 2026-10-14. These notes and the code
were researched and written with AI agents working in MV2's name. Text CC BY 4.0, code Apache-2.0.

**Summary.** The growth base beta of the smallest Kakeya sets in F_3^n, open since Kopparty, Lev,
Saraf and Sudan (2010) placed it between 1.8 and about 1.890, equals the best base achievable by a
single small "rainbow triple" of sets in low dimension (base theorem below). The new result of this
version is an **explicit construction with base below 2**: the triple (W, W, B), where W is the set
of vectors using at most two of the three symbols and B the 0/1 vectors with at most m/3 ones, is
rainbow, and a "tripling lemma" turns it into explicit Kakeya sets with

    K(n, 3) <= 189 * 432^(n/9),   432^(1/9) = 2^(4/9) 3^(1/3) = 1.96262...

The explicit constructions quoted by KLSS have base 2 (2^(n+1) - 1 and 3 * 2^n); their base 1.890
comes from a random choice of linear maps. An energy lemma also shows that no triple in dimension
1 or 2 beats base 2 (exactly: tau(2) = 64, sigma(2) = 12, s(2) = 4, with an exhaustive check in the
tests) and that in dimension 3 every set of the triple needs at least 7 points.

## The question

K(n, 3) is the least size of a Kakeya set in F_3^n, i.e. a set that contains a line
{a, a + d, a + 2d} (a three-term progression) in every direction d. Since the product of two
Kakeya sets is a Kakeya set, K(m + n, 3) <= K(m, 3) K(n, 3), so the limit

    beta = lim K(n, 3)^(1/n) = inf_n K(n, 3)^(1/n)

exists (Fekete). Kopparty, Lev, Saraf and Sudan (KLSS, "Kakeya-type sets in finite vector spaces",
arXiv 1003.3736, 2010) call finding such limits "very interesting" and single out q = 3 as "most
intriguing". Their bounds for q = 3:

- lower: |K| >= (9/5)^n (their Theorem 1; Bukh-Chao 2021 give 3^n (5/3)^-(n-1)), so beta >= 1.8;
- upper, explicit: the "missing digit" set {0,1}^n u {1,2}^n of size 2^(n+1) - 1 (their Theorem 7,
  "very clean, but rather weak"; Smyth's 2008 comment on Tao's blog) and 3 * 2^n from their
  Theorem 6 with f(x) = x^2;
- upper, probabilistic: (3 / 2^(2/3))^(n + O(sqrt n)), about 1.890^n (their Theorem 8: a set
  covering the directions whose coordinates are roughly equally distributed, moved by n randomly
  chosen linear automorphisms).

So 1.8 <= beta <= 1.890..., with no explicit construction below base 2. Kyureghyan, Muller and
Wang (2013) improved the explicit constructions for even q only. We found no later paper that
changes these numbers for q = 3 (searches 2026-10-14; see the Sources).

## Definitions

A triple (A, M, C) of subsets of F_3^n is a **rainbow triple** if for every d in F_3^n, including
d = 0, some a in A has a + d in M and a + 2d in C. (For d = 0 this says A, M and C share a point.)
Equivalently, there is a function f on F_3^n (the witness) with f(d) in A, f(d) + d in M and
f(d) + 2d in C for every d; the sets can then be taken to be the three images.

- sigma(n) = least |A| + |M| + |C| over rainbow triples in F_3^n;
- s(n) = least max(|A|, |M|, |C|); tau(n) = least |A| |M| |C|;
- h(n) = least |A| + |C| over rainbow triples with M inside A.

## Elementary facts (all checked by the tests where finite)

1. **Symmetry.** Being rainbow does not depend on the order of A, M, C. Cyclic shift: if a works
   for d then m = a + d has m + d = a + 2d in C and m + 2d = a in A, so (M, C, A) is rainbow.
   Reversal: if a works for -d = 2d, then a + 2d is in M and a + d is in C; put c = a + d, so
   c + d = a + 2d is in M and c + 2d = a is in A, and (C, M, A) is rainbow. Rainbow triples are
   also preserved by affine maps x -> Lx + b applied to all three sets (d -> Ld).
2. **Pair bound.** d -> (a, a + d) is injective, so |A| |M| >= 3^n, and likewise for the other
   two pairs. Hence tau(n) >= 3^(3n/2) and s(n) >= 3^(n/2).
3. **Products.** If (A, M, C) is rainbow in F_3^n and (A', M', C') in F_3^m then
   (A x A', M x M', C x C') is rainbow in F_3^(n+m) (solve each coordinate block separately).
   So s and tau are submultiplicative.
4. **Kakeya sets give triples.** If K is a Kakeya set then (K, K, K) is rainbow.

## Layer lemma

For K in F_3^(n+1) let K_t = {x : (x, t) in K}. Then **K is a Kakeya set if and only if
(K_0, K_1, K_2) is a rainbow triple and every direction (d, 0), d != 0, has a line inside one K_t.**

Proof. A direction whose last coordinate is nonzero is, after scaling, (d, 1) for a unique d in
F_3^n, and every line in that direction meets each slice once: {(a, 0), (a + d, 1), (a + 2d, 2)}.
So K has such a line exactly when some a in K_0 has a + d in K_1 and a + 2d in K_2. A direction
(d, 0) has all its lines inside single slices. QED. (Tested exhaustively for n = 1 and on random
sets for n = 2.)

**Corollary 1.** sigma(n) <= K(n+1, 3) <= sigma(n) + K(n, 3). The left inequality is the lemma;
for the right, take a smallest rainbow triple (A, M, C) and a smallest Kakeya set B of F_3^n and
use A x {0} + (B u M) x {1} + C x {2}. Telescoping from K(0, 3) = 0,
K(N, 3) <= sigma(0) + sigma(1) + ... + sigma(N - 1).

## Doubling lemma (a common generalisation of the missing-digit set)

Let G1 = ({0, 1}, {0, 1}, {0, 2}), a rainbow triple of F_3 whose middle set lies inside its first
set and inside {0, 1}. **If B is a Kakeya set of F_3^n and (A, M, C) is a rainbow triple with M
inside both A and B, then K(n + k, 3) <= |B| + (2^k - 1)(|A| + |C|) for every k >= 0.**

Proof. B' = A x {0} + B x {1} + C x {2} is Kakeya by the layer lemma (its middle slice contains M,
and B supplies the horizontal lines), with |B'| = |B| + |A| + |C|. The product triple
(A x {0,1}, M x {0,1}, C x {0,2}) is rainbow, its middle lies inside its first set, and it lies
inside B' because M x {0} is inside A x {0} and M x {1} inside B x {1}. So the hypothesis repeats
with |A| + |C| doubled. QED.

With n = 0, B = A = M = C = {point} this gives sets of size 2^(k+1) - 1 in F_3^k, the size of the
missing-digit set (the chain `triples.smyth(k)`, tested for k <= 5). Any base case with
|A| + |C| < 2^(n+1) therefore beats the missing-digit set by a constant factor in every
dimension; the asymptotic constant is h(n) / 2^n (it is 2 for the missing-digit set).
Exhaustively h(1) = 4 (tested). For n = 2, a hand case analysis gives h(2) = 8, so no gain at
n = 2; `search_triples.py` re-checks this by plain enumeration and runs the solver for n = 3, 4.

## Tripling lemma

**If (A, M, C) is a rainbow triple of F_3^m, then K = (A x M x C) u (M x C x A) u (C x A x M) is a
Kakeya set of F_3^(3m). Hence K(3m, 3) <= 3 tau(m).**

Proof. Let d = (d1, d2, d3) with blocks in F_3^m, and let a_i be a witness for d_i (a_i in A,
a_i + d_i in M, a_i + 2d_i in C). Put x = (a1, a2 + d2, a3 + 2d3). Then x is in A x M x C,
x + d = (a1 + d1, a2 + 2d2, a3) is in M x C x A, and x + 2d = (a1 + 2d1, a2, a3 + d3) is in
C x A x M. QED. (`tripled`, `tripled_point`; the tests check G1, the slices of the 13-point set and
the construction below.)

This is the cyclic product of the base theorem below without the layering, and it has constant 3.

## Base theorem

**beta = inf_n tau(n)^(1/(3n)) = inf_n s(n)^(1/n).**

Proof. By fact 4, tau(n) <= K(n, 3)^3 and s(n) <= K(n, 3), so both infima are at most beta.
Conversely, let (A, M, C) be rainbow in F_3^m with P = |A| |M| |C|. Its j-fold product (fact 3) is
a rainbow triple of F_3^(mj) with product P^j, so the tripling lemma gives K(3mj, 3) <= 3P^j, and
K(3mj + r, 3) <= 3P^j K(r, 3) by products; hence beta <= P^(1/(3m)). (An earlier
version proved this with the cyclic product T* = (A x M x C, M x C x A, C x A x M) and the general
layer step; the tripling lemma shortens it.) For s, the triple has three sets of at most s(m)
points, so P <= s(m)^3. QED.

Consequences: since beta >= 9/5 (KLSS Theorem 1), tau(n) >= (9/5)^(3n) = 5.832^n and
s(n) >= (9/5)^n for every n, which is stronger than fact 2. A rainbow triple in some F_3^m with
|A| |M| |C| < 8^m gives an explicit family of Kakeya sets with base below 2 (the next section gives
one); one with |A| |M| |C| < 6.75^m would beat the 1.890 base of KLSS Theorem 8.

## Two symbols: an explicit construction below base 2

Let W(m) be the set of vectors of F_3^m that use at most two of the three symbols, and B(m) the
set of 0/1 vectors with at most floor(m/3) ones.

**Theorem. (W(m), W(m), B(m)) is a rainbow triple of F_3^m.**

Proof. Given d, let j be a least frequent symbol of d (it occurs at most m/3 times) and e the 0/1
vector with ones exactly where d_i = j, so e is in B(m). Put a = d + e. Where d_i = j, a_i = j + 1;
elsewhere a_i = d_i, which is j + 1 or j + 2. So a uses only j + 1 and j + 2. Next,
a + d = 2d + e is 2j + 1 where d_i = j, 2j + 2 where d_i = j + 1, and 2j + 4 = 2j + 1 where
d_i = j + 2, so it uses only two symbols. Finally a + 2d = 3d + e = e. QED.

**Corollary.** For every n, **K(n, 3) <= 189 * 432^(n/9)**, by explicit sets; 432^(1/9) =
2^(4/9) 3^(1/3) = 1.96262. So beta <= 1.96262 with an explicit construction.

Proof. |W(m)| = 3 * 2^m - 3 (inclusion-exclusion) and |B(m)| = sum over k <= m/3 of C(m, k)
<= 2^(m h(1/3)) = (27/4)^(m/3), where h is the binary entropy (the standard tail bound). So
P(m) = |W|^2 |B| <= 9 * 4^m (27/4)^(m/3) = 9 * 432^(m/3). By the tripling lemma the set

    K_m = (W x W x B) u (W x B x W) u (B x W x W)   in F_3^(3m)

is Kakeya, of size at most 3P(m) <= 27 * 432^(3m/9). For n = 3m + r with r = 1, 2 multiply by a
Kakeya set of F_3^r (3 or 7 points). QED. In words: a point of F_3^(3m), cut into three blocks of
length m, is in K_m when one block is a 0/1 vector with at most m/3 ones and each of the other two
misses a symbol.

Exact sizes. Since B lies inside W, |K_m| = 3w^2 b - 3wb^2 + b^3 with w = |W|, b = |B|
(`wwb_kakeya_size`). The product P(m) first drops below 8^m at m = 14 (P(14) = 3,553,383,199,671
< 8^14 = 4,398,046,511,104, base P^(1/42) = 1.990), and |K_m| first drops below the missing-digit
size 2^(3m+1) - 1 at m = 17 (n = 51). Because the threshold floor(m/3) jumps, the comparison is not
monotone in m (m = 18 and 19 lose again); the clean bound |K_m| <= 27 * 432^(m/3) wins for every
m >= 46 (n >= 138). As
m -> infinity, P(m)^(1/(3m)) -> 432^(1/9) (the test checks m = 3000 to within 0.001).
`python triples.py wwb 40` prints the table.

What this does and does not do. It is the first explicit construction for q = 3 with base below 2
that we know of; it does **not** improve KLSS's probabilistic 1.890. Within constructions that, like
this one, treat all coordinates alike (the witness depends on d only through which positions hold
which symbol, with one rule per composition), the directions with two symbols each on half of the
coordinates seem to force two of the three sets to have about 2^m points, which would make
432^(1/9) close to the limit of the method; we have not proved this. Beating 1.890 explicitly, or
finding a small triple (for example (6, 7, 7) in F_3^3) below 6.75^m, remains open.

## Energy lemma

For finite sets X, Y in F_3^n let E(X, Y) be the number of quadruples (x, y, x', y') in
X x Y x X x Y with x + y = x' + y'. Counting the same quadruples by x - x' = y' - y, and using that
the number of pairs in Y with a given difference v equals the number with difference -v, E(X, Y)
is also the number of quadruples with y - x = y' - x'. So with r(v) = #{(x, y) : y - x = v} and
s(v) = #{(x, y) : x + y = v}, **sum r(v)^2 = E(X, Y) = sum s(v)^2.**

**Lemma.** Let (A, M, C) be a rainbow triple of F_3^n, N = 3^n, x = |A|, y = |M|, z = |C|,
k = min(x, y). Then

- E(A, M) <= U(x, y, N), the largest sum of squares of N integers in [1, k] adding up to x y;
- E(A, M) >= L(x, y, z, N), the least sum of squares of N integers in [0, k] adding up to x y of
  which some z add up to at least N.

The same holds for the pairs (A, C) and (M, C), with the third set's size as z.

Proof. Every d is a difference m - a (take the witness for d), so r(v) >= 1 for every v; also
r(v) <= k and sum r(v) = x y. That gives the upper bound. For the lower bound, the witness pairs
(a, a + d) of the N different d are N different pairs of A x M, and the sum of each is
2a + d = -(a + 2d), which lies in -C, a set of z values; also s(v) <= k. The other pairs follow
from fact 1. QED.

U and L are explicit (`energy_upper`, `energy_lower`: U puts as many values as possible at k,
L spreads evenly); `energy_excludes(a, b, c, n)` returns True when L > U for one of the three pairs
or a pair bound fails. The tests check the lemma against brute force in F_3, against real triples in
F_3^2 and F_3^3, and the hand values used below.

Two remarks. Applying the same counting to a pair with |A| |M| = N exactly gives r = 1 everywhere,
hence s <= 1 and |C| >= N: a triple that is tight for a pair bound forces the third set to be the
whole space, which is why the pair bound alone is far from the truth. And the lemma bounds
E(A, M) from both sides only through the sizes, so it says nothing about large unbalanced triples;
those are left to the solver.

## Exact values in the plane

**Theorem. tau(2) = 64, sigma(2) = 12 and s(2) = 4. So no rainbow triple of F_3^1 or F_3^2 gives a
base below 2, and G1 x G1 = ({0,1}^2, {0,1}^2, {0,2}^2) is optimal for all three measures.**

Proof. G1 x G1 has sizes (4, 4, 4). Every size pattern a <= b <= c with abc < 64 that the pair
bound allows is one of (2,5,5), (2,5,6), (3,3,3), (3,3,4), (3,3,5), (3,3,6), (3,3,7), (3,4,4),
(3,4,5); every pattern with a + b + c <= 11 or c <= 3 is among them. The energy lemma excludes all
but (3,4,5). By hand, for example: for (3,4,4), U(3,4,9) = 20 (one difference three times, one
twice, seven once) while L(3,4,4,9) = 24 (nine sums spread 3, 2, 2, 2 over -C and three more
outside); for (3,3,c) with c <= 8, U = 9 forces every sum to be distinct, so |C| >= 9.

For (3,4,5) the lemma is tight for the pair (A, M): U = L = 20, so the differences m - a take one
value three times, one value twice and the other seven once. A value taken three times means
A + v lies inside M; translate M (this does not change the multiplicities, only their positions)
so that M = A u {p}. If A is a line, its differences hit two values three times. Otherwise the six
differences a_i - a_j (i != j) are distinct (a coincidence would force a_1 + a_2 + a_3 = 0, i.e. a
line), nonzero and closed under negation, so they miss exactly two nonzero values z and -z. These
must come from p - A, say p - a_i = z and p - a_j = -z; but then a_j - a_i = 2z = -z is a
difference of A, a contradiction. QED.

The tests also confirm the whole theorem without the lemma: `rainbow_with_sizes` searches every
triple with the given sizes (after translating the d = 0 witness to the origin) and finds none for
the nine patterns, and finds one for (4, 4, 4); it is checked against brute force over all 512
triples of F_3.

## Dimension three

**s(3) >= 7.** Patterns with all three sets of at most 6 points: the pair bound leaves (5,6,6) and
(6,6,6), and the energy lemma excludes both (for (6,6,6): U(6,6,27) = 86 < 132 = L(6,6,6,27); for
(5,6,6): U(5,6,27) = 42 < 126). The KLSS bound gives only s(3) >= 6. Since G1^3 has (8, 8, 8),
s(3) is 7 or 8.

- A triple (7, 7, 7) in F_3^3 would give base 7^(1/3) = 1.913, better than the two-symbol
  construction; so would any triple with product at most 432 (for example (7, 7, 8) = 392).
- The lemma leaves (6, 7, 7) open; its product 294 is below 6.75^3 = 307.5, so such a triple would
  give base 294^(1/9) ~ 1.880 and beat the 1.890 of KLSS's random construction.
- tau(3) lies between 199 (from beta >= 1.8) and 512.

`search_s3.py` asks the solver exactly these questions and `anneal.py` searches the witness
function directly; `python triples.py energy 3 512` lists every size pattern below 512 the lemma
leaves. A hand analysis of one structured family (every slice a two-point "domino" in one fixed
direction, with one extra point) found no (7, 7, 7) triple there; it is not a proof.

| n | s(n) | tau(n) | sigma(n) | source |
|---|------|--------|----------|--------|
| 0 | 1 | 1 | 3 | trivial |
| 1 | 2 | 8 | 6 | exhaustive (tests) |
| 2 | 4 | 64 | 12 | theorem above (energy lemma + exhaustive tests) |
| 3 | 7 or 8 | 199..512 | ..24 | energy lemma; G1^3 |
| 14 | | < 8^14 | | two-symbol triple, 3,553,383,199,671 |

## What `search_triples.py` computes

maxsize s(n) and sigma(n) for n <= 4; the size patterns above; h(n) for n <= 4 with the
exhaustive n = 2 check; and Kakeya sets built layer by layer on the 13-point set of F_3^3 and the
27-point set of F_3^4, up to F_3^6, compared with the best sets known before (27, 56 in F_3^4,
F_3^5 from our earlier solver run; 127 = 2^7 - 1 in F_3^6). Every set it reports is re-checked
in plain Python. Results will be added here when the run finishes.

## Status of the claims

- Proved here (elementary): the layer lemma, Corollary 1, the doubling lemma, the tripling lemma,
  the base theorem, the energy lemma, the two-symbol theorem and K(n, 3) <= 189 * 432^(n/9),
  tau(2) = 64, sigma(2) = 12, s(2) = 4, s(3) >= 7, the small values for n <= 1, h(2) = 8 (hand
  proof, re-checked by the job).
- Not claimed: any value of beta or any improvement of the probabilistic 1.890 bound.
- Novelty: the layer lemma is elementary and may be folklore. We did not find the rainbow-triple
  reformulation of beta, the tripling or doubling lemmas, the energy lemma for triples, or any
  explicit construction for q = 3 with base below 2 in KLSS, Kyureghyan-Muller-Wang or the other
  sources below. The identity between difference and sum energies and the binomial tail bound are
  standard. If an earlier explicit construction below base 2 exists, we would like to cite it.

## Sources

- S. Kopparty, V. F. Lev, S. Saraf, M. Sudan, "Kakeya-type sets in finite vector spaces",
  J. Algebraic Combin. 34 (2011) 337-355; arXiv 1003.3736.
- G. Kyureghyan, P. Muller, Q. Wang, "On the size of Kakeya sets in finite vector spaces",
  Electron. J. Combin. 20(3) (2013) #P36; arXiv 1302.5591.
- B. Bukh, T.-W. Chao, "Sharp density bounds on the finite field Kakeya problem", Discrete
  Analysis 2021:26; arXiv 2108.00074.
- T. Tao, "Dvir's proof of the finite field Kakeya conjecture", blog, 2008-03-24 (comments: the
  missing-digit set and the F_3 question).
- Z. Dvir, "On the size of Kakeya sets in finite fields", J. AMS 22 (2009); arXiv 0803.2336.
- MV2 of Munim, Inc., kakeya-conjecture v1.1.0 (K(3, 3) = 13 and the solver sets),
  https://github.com/munimtechnologies/kakeya-conjecture
