"""Rainbow progression triples in F_3^n and the layer construction of Kakeya sets over F_3.

Over F_3 a line is a three-term progression {a, a + d, a + 2d}. A triple (A, M, C) of subsets of
F_3^n is a *rainbow triple* if for every d in F_3^n (d = 0 included) some a in A has a + d in M and
a + 2d in C. Layer lemma: a set K in F_3^(n+1) is a Kakeya set exactly when its three slices
K_t = {x : (x, t) in K}, t = 0, 1, 2, form a rainbow triple (K_0, K_1, K_2) and every direction
(d, 0) with d != 0 has a line inside one slice. (A line in direction (d, 1) through (a, 0) is
{(a, 0), (a + d, 1), (a + 2d, 2)}; every direction with last coordinate nonzero is some (d, 1).)

Consequences proved in TRIPLES.md and checked by test_triples.py:
- sigma(n) <= K(n+1, 3) <= sigma(n) + K(n, 3), where sigma(n) is the least |A| + |M| + |C| of a
  rainbow triple in F_3^n; so K(N, 3) <= sigma(0) + ... + sigma(N - 1).
- Doubling: if B in F_3^n is Kakeya and (A, M, C) is rainbow with M inside A and B, then
  K(n + k, 3) <= |B| + (2^k - 1)(|A| + |C|). Smyth's set, of size 2^(k+1) - 1, is the case n = 0.
- Base: lim K(n, 3)^(1/n) = inf over n of tau(n)^(1/(3n)), tau(n) the least |A| |M| |C|.
- Energy lemma: the additive energy of (A, M) is at most energy_upper(...) because every d is a
  difference m - a, and at least energy_lower(...) because 3^n of the sums a + m lie in -C.
  energy_excludes(a, b, c, n) applies it to all three pairs.
- In the plane: tau(2) = 64, sigma(2) = 12, s(2) = 4 (rainbow_with_sizes, exhaustive), so no
  triple in F_3^1 or F_3^2 beats base 2; in F_3^3, s(3) >= 7 (energy lemma).
- Tripling lemma: for a rainbow triple (A, M, C) of F_3^m, the union of A x M x C, M x C x A and
  C x A x M is a Kakeya set of F_3^(3m) (tripled).
- Two symbols: (W, W, B) is a rainbow triple of F_3^m, where W is the set of vectors using at most
  two of the three symbols and B the 0/1 vectors with at most m // 3 ones. Tripled, it gives
  explicit Kakeya sets with K(n, 3) <= 189 * 432^(n/9), about 189 * 1.9626^n: base below 2.

    python triples.py smyth 4        # the doubling chain from F_3^0, checked to be Kakeya
    python triples.py check          # rainbow / Kakeya checks of the built-in examples
    python triples.py plane          # the exhaustive F_3^2 table behind tau(2) = 64
    python triples.py energy 3 512   # size patterns with abc < 512 in F_3^3 the energy lemma leaves
    python triples.py wwb 40         # sizes of the two-symbol construction for m = 1..40

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0.
"""

import math
import sys
from itertools import combinations, product

# The 13-point Kakeya set of F_3^3 (kakeya.py K33; K(3, 3) = 13).
K33 = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (1, 1, 0), (0, 0, 1), (2, 0, 1), (0, 1, 1), (1, 1, 1),
       (0, 2, 1), (2, 2, 1), (0, 0, 2), (2, 1, 2), (0, 2, 2)]

# The one-dimensional rainbow triple used by the doubling chain: sizes 2, 2, 2, middle inside
# the first and inside {0, 1}.
G1 = ({(0,), (1,)}, {(0,), (1,)}, {(0,), (2,)})


def vectors(n):
    return list(product(range(3), repeat=n))


def add(u, v):
    return tuple((x + y) % 3 for x, y in zip(u, v))


def scale(t, u):
    return tuple((t * x) % 3 for x in u)


def directions(n):
    """One vector per direction of F_3^n: nonzero, first nonzero coordinate 1."""
    out = []
    for d in vectors(n):
        nz = [x for x in d if x]
        if nz and nz[0] == 1:
            out.append(d)
    return out


def is_kakeya(K, n):
    """True if K (tuples in F_3^n) contains a line {a, a + d, a + 2d} in every direction d."""
    K = set(K)
    return all(any(add(a, d) in K and add(a, scale(2, d)) in K for a in K) for d in directions(n))


def missing_differences(A, M, C, n):
    """The d in F_3^n with no a in A such that a + d in M and a + 2d in C."""
    M, C = set(M), set(C)
    return [d for d in vectors(n)
            if not any(add(a, d) in M and add(a, scale(2, d)) in C for a in A)]


def is_rainbow(A, M, C, n):
    return not missing_differences(A, M, C, n)


def lift(A, M, C):
    """A x {0} + M x {1} + C x {2}, a subset of F_3^(n+1)."""
    return ({tuple(a) + (0,) for a in A} | {tuple(m) + (1,) for m in M}
            | {tuple(c) + (2,) for c in C})


def slices(K):
    """The slices (K_0, K_1, K_2) of K in F_3^(n+1) by the last coordinate."""
    return tuple({p[:-1] for p in K if p[-1] == t} for t in range(3))


def is_kakeya_by_layers(K, n):
    """The layer lemma's test for K in F_3^(n+1); must agree with is_kakeya(K, n + 1)."""
    S = slices(K)
    if not is_rainbow(S[0], S[1], S[2], n):
        return False
    return all(any(is_line_in(S[t], d) for t in range(3)) for d in directions(n))


def is_line_in(S, d):
    S = set(S)
    return any(add(a, d) in S and add(a, scale(2, d)) in S for a in S)


def product_triple(T, U):
    """Coordinatewise product of a rainbow triple of F_3^n and one of F_3^m."""
    return tuple({tuple(a) + tuple(b) for a in S for b in R} for S, R in zip(T, U))


def doubling_chain(B, T, n, k):
    """Apply the doubling step k times to a Kakeya set B of F_3^n and a rainbow triple T = (A, M, C)
    with M inside A and inside B. Returns the Kakeya set of F_3^(n+k), of size
    |B| + (2^k - 1)(|A| + |C|), and the last triple."""
    A, M, C = (set(S) for S in T)
    B = set(B)
    if not (M <= A and M <= B):
        raise ValueError("the middle set must lie inside the first set and inside B")
    for _ in range(k):
        B = lift(A, B, C)
        A, M, C = product_triple((A, M, C), G1)
    return B, (A, M, C)


def layered_step(B, T):
    """General step: a Kakeya set of F_3^(n+1) from a Kakeya set B of F_3^n and any rainbow triple
    (A, M, C) of F_3^n, of size |A| + |B u M| + |C|."""
    A, M, C = T
    return lift(A, set(B) | set(M), C)


def smyth(k):
    """The doubling chain from F_3^0: a Kakeya set of F_3^k with 2^(k+1) - 1 points."""
    return doubling_chain({()}, ({()}, {()}, {()}), 0, k)[0]


def sizes(T):
    return tuple(len(S) for S in T)


# Energy lemma. For a rainbow triple (A, M, C) of F_3^n put N = 3^n, x = |A|, y = |M|, z = |C|,
# k = min(x, y) and E = #{(a, m, a', m') : a + m = a' + m'} = #{(a, m, a', m') : m - a = m' - a'}.
# Upper bound: r(v) = #{(a, m) : m - a = v} is at least 1 for every v (each d is a difference) and
# at most k, with total x y. Lower bound: the witness pairs (a, a + d) of the N differences d are
# distinct, and each has sum 2a + d = -(a + 2d) in -C, a set of z values; s(v) <= k as well.


def spread(total, slots, cap):
    """Least sum of squares of `slots` integers in [0, cap] adding up to `total` (None if none)."""
    if total < 0 or total > slots * cap:
        return None
    if slots == 0:
        return 0
    q, r = divmod(total, slots)
    return r * (q + 1) ** 2 + (slots - r) * q * q


def energy_upper(x, y, N):
    """Largest sum of r(v)^2 over N values with 1 <= r(v) <= min(x, y) and total x y (None if the
    pair bound x y >= N fails)."""
    k = min(x, y)
    e = x * y - N
    if e < 0:
        return None
    if k == 1:
        return N if e == 0 else None
    q, rem = divmod(e, k - 1)
    if q + (1 if rem else 0) > N:
        return None
    return N + q * (k * k - 1) + (1 + rem) ** 2 - 1


def energy_lower(x, y, z, N):
    """Least sum of s(v)^2 over N values with 0 <= s(v) <= min(x, y), total x y, and at least N of
    the total on some z values (None if impossible)."""
    k = min(x, y)
    best = None
    for inside in range(N, x * y + 1):
        a = spread(inside, z, k)
        b = spread(x * y - inside, N - z, k)
        if a is not None and b is not None and (best is None or a + b < best):
            best = a + b
    return best


def energy_excludes(a, b, c, n):
    """True if the energy lemma (with the pair bound) rules out a rainbow triple of F_3^n with set
    sizes a, b, c (in any order: being rainbow does not depend on the order)."""
    N = 3 ** n
    for x, y, z in ((a, b, c), (a, c, b), (b, c, a)):
        hi = energy_upper(x, y, N)
        lo = energy_lower(x, y, z, N)
        if hi is None or lo is None or lo > hi:
            return True
    return False


def surviving_patterns(n, below):
    """Size patterns a <= b <= c with abc < below that the energy lemma does not exclude."""
    N = 3 ** n
    out = []
    for a in range(1, N + 1):
        for b in range(a, N + 1):
            for c in range(b, N + 1):
                if a * b * c >= below:
                    break
                if not energy_excludes(a, b, c, n):
                    out.append((a, b, c))
    return sorted(out, key=lambda p: (p[0] * p[1] * p[2], p))


def rainbow_with_sizes(a, b, c, n):
    """A rainbow triple of F_3^n whose sets have exactly a, b, c points, or None. Exhaustive: the
    d = 0 witness lies in all three sets, and translating all three by the same vector keeps a
    triple rainbow, so 0 may be assumed to lie in each set."""
    pts = vectors(n)
    zero, rest = pts[0], pts[1:]
    for A0 in combinations(rest, a - 1):
        A = (zero,) + A0
        for M0 in combinations(rest, b - 1):
            M = set(M0)
            M.add(zero)
            needs = []
            for d in pts:
                need = {add(p, scale(2, d)) for p in A if add(p, d) in M}
                if not need:
                    break
                needs.append(need)
            else:
                for C0 in combinations(rest, c - 1):
                    C = set(C0)
                    C.add(zero)
                    if all(not C.isdisjoint(need) for need in needs):
                        return set(A), M, C
    return None


# Every size pattern a <= b <= c in F_3^2 with abc < 64 = 8^2 that the pair bound ab >= 9 allows.
PLANE_PATTERNS = [(2, 5, 5), (2, 5, 6), (3, 3, 3), (3, 3, 4), (3, 3, 5), (3, 3, 6), (3, 3, 7),
                  (3, 4, 4), (3, 4, 5)]


# Tripling lemma. If (A, M, C) is rainbow in F_3^m, then K = A x M x C u M x C x A u C x A x M is a
# Kakeya set of F_3^(3m): for d = (d1, d2, d3) take witnesses a_i of d_i; the point
# x = (a1, a2 + d2, a3 + 2 d3) has x in A x M x C, x + d in M x C x A and x + 2d in C x A x M.


def tripled(T):
    """The tripling lemma's set A x M x C u M x C x A u C x A x M for a triple T = (A, M, C)."""
    A, M, C = T
    out = set()
    for X, Y, Z in ((A, M, C), (M, C, A), (C, A, M)):
        out |= {tuple(x) + tuple(y) + tuple(z) for x in X for y in Y for z in Z}
    return out


def tripled_point(d, m, witness):
    """The point x of the tripling lemma for the direction d of F_3^(3m), given a witness function
    of the triple (a = witness(d') has a in A, a + d' in M, a + 2d' in C)."""
    d1, d2, d3 = tuple(d[:m]), tuple(d[m:2 * m]), tuple(d[2 * m:])
    return witness(d1) + add(witness(d2), d2) + add(witness(d3), scale(2, d3))


# Two symbols. W(m): the vectors of F_3^m that use at most two of the three symbols. B(m): the 0/1
# vectors with at most m // 3 ones. For d let j be a least frequent symbol of d and e the 0/1
# indicator of {i : d_i = j}, so e is in B(m). Then a = d + e uses only j + 1 and j + 2, a + d = 2d + e
# uses only 2j + 1 and 2j + 2, and a + 2d = e. So (W(m), W(m), B(m)) is a rainbow triple.

TWO_SYMBOL_BASE = 432 ** (1.0 / 9)   # 2^(4/9) 3^(1/3) = 2^((2 + h(1/3)) / 3), about 1.96262


def in_w(v):
    return len(set(v)) <= 2


def in_b(v):
    return all(x in (0, 1) for x in v) and sum(v) <= len(v) // 3


def wwb_witness(d):
    """The point a with a in W, a + d in W and a + 2d in B."""
    counts = [sum(1 for x in d if x == s) for s in range(3)]
    j = counts.index(min(counts))
    return tuple((x + 1) % 3 if x == j else x for x in d)


def wwb_triple(m):
    """(W(m), W(m), B(m)) as sets of tuples (small m only)."""
    W = {v for v in vectors(m) if in_w(v)}
    B = {v for v in vectors(m) if in_b(v)}
    return W, set(W), B


def wwb_sizes(m):
    """(|W(m)|, |W(m)|, |B(m)|) for m >= 1, by counting."""
    w = 3 * 2 ** m - 3
    return w, w, sum(math.comb(m, k) for k in range(m // 3 + 1))


def wwb_product(m):
    w, _, b = wwb_sizes(m)
    return w * w * b


def wwb_kakeya_size(m):
    """Exact size of the tripled set W x W x B u W x B x W u B x W x W in F_3^(3m) (B lies in W)."""
    w, _, b = wwb_sizes(m)
    return 3 * w * w * b - 3 * w * b * b + b ** 3


def wwb_base(m):
    """(|W|^2 |B|)^(1/(3m)): by the base theorem, lim K(n, 3)^(1/n) is at most this."""
    return 2 ** (math.log2(wwb_product(m)) / (3 * m))


def main(argv):
    if len(argv) >= 2 and argv[0] == "smyth":
        k = int(argv[1])
        K = sorted(smyth(k))
        print("doubling chain from F_3^0: %d points in F_3^%d, Kakeya: %s"
              % (len(K), k, is_kakeya(K, k)))
        for p in K:
            print(" ".join(map(str, p)))
        return 0
    if argv and argv[0] == "check":
        S = slices(K33)
        print("K33 slices:", sizes(S), "rainbow:", is_rainbow(*S, 2))
        print("G1 rainbow:", is_rainbow(*G1, 1))
        for k in range(6):
            K = smyth(k)
            print("chain k=%d: %d points, Kakeya %s" % (k, len(K), is_kakeya(K, k)))
        return 0
    if argv and argv[0] == "plane":
        for p in PLANE_PATTERNS + [(4, 4, 4)]:
            T = rainbow_with_sizes(*p, 2)
            print("sizes %s (abc %d): energy lemma excludes %s; exhaustive: %s"
                  % (p, p[0] * p[1] * p[2], energy_excludes(*p, 2),
                     "none" if T is None else [sorted(S) for S in T]))
        return 0
    if len(argv) >= 3 and argv[0] == "energy":
        n, below = int(argv[1]), int(argv[2])
        for p in surviving_patterns(n, below):
            P = p[0] * p[1] * p[2]
            print("%s abc %d base %.4f" % (p, P, P ** (1.0 / (3 * n))))
        return 0
    if len(argv) >= 2 and argv[0] == "wwb":
        print("m |W| |B| |W|^2|B|/8^m base n=3m Kakeya/(2^(n+1)-1)")
        for m in range(1, int(argv[1]) + 1):
            w, _, b = wwb_sizes(m)
            print("%d %d %d %.4f %.5f %d %.4f"
                  % (m, w, b, wwb_product(m) / 8 ** m, wwb_base(m), 3 * m,
                     wwb_kakeya_size(m) / (2 ** (3 * m + 1) - 1)))
        print("limit 432^(1/9) = %.5f" % TWO_SYMBOL_BASE)
        return 0
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
