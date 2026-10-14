"""Finite-field Kakeya sets: exact arithmetic, line enumeration, checking, constructions and bounds.

A Kakeya set in F_q^n is a set of points containing a full line {b + t*d : t in F_q} in every
direction d. This module builds the finite field F_q (q prime or q in {4, 8, 9, 16}), the points,
directions and lines of AG(n, q), checks the Kakeya property, builds the parabola construction and
product sets, computes the classical lower bounds, finds the minimum size of a Kakeya set by
brute force when the search space is small, and decides by exhaustive search whether a Kakeya set
of a given size exists (this proves K(3, 3) = 13).

Points are encoded as integers: the point (c_0, ..., c_{n-1}) is sum c_i * q**i.

    python kakeya.py check N Q FILE     check a Kakeya set (one point per line, coordinates c_0 ... c_{n-1})
    python kakeya.py bounds N Q         lower bounds and the parabola upper bound for AG(N, Q)
    python kakeya.py parabola N Q       print the parabola Kakeya set in AG(N, Q), Q odd
    python kakeya.py smyth N Q          print the set (F*)^N plus the 0/1 vectors in AG(N, Q)
    python kakeya.py brute N Q          exact minimum by brute force (small cases only)
    python kakeya.py exhaust N Q K      is there a Kakeya set of at most K points? (small cases only)
    python kakeya.py records            the exact values and bounds computed by MV2 (section 7)

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. Apache-2.0 licence.
"""

import json
import sys
from fractions import Fraction
from itertools import combinations, product
from math import comb

# Monic irreducible polynomials for the non-prime fields, lowest coefficient first.
IRREDUCIBLE = {
    4: (2, 2, [1, 1, 1]),        # x^2 + x + 1 over F_2
    8: (2, 3, [1, 1, 0, 1]),     # x^3 + x + 1 over F_2
    9: (3, 2, [1, 0, 1]),        # x^2 + 1 over F_3
    16: (2, 4, [1, 1, 0, 0, 1]),  # x^4 + x + 1 over F_2
}

# Minimum Kakeya set sizes in F_2^n. Over F_2 a Kakeya set is exactly a difference basis of the
# group C_2^n (see is_difference_cover_f2); these are the difference sizes computed by Banakh and
# Gavrylkiv, Acta Sci. Math. (Szeged) 85 (2019), Tables 1 and 3.
KNOWN_F2 = {1: 2, 2: 3, 3: 5, 4: 6, 5: 10, 6: 14}

# MV2's computation kakeya-min-v1 (October 2026), raw output in data/kakeya-min-v1/results.json:
# (n, q) -> (best lower bound, smallest Kakeya set found). The lower bounds are Bukh-Chao except
# (3, 3), where kakeya_sets_through_axis proves 13 with no solver; every upper bound is a set in
# that file, re-checked line by line by the tests.
COMPUTED = {(3, 3): (13, 13), (4, 3): (18, 27), (5, 3): (32, 56), (3, 4): (21, 28),
            (4, 4): (48, 64), (3, 5): (39, 53), (3, 7): (100, 128)}

# A 13-point Kakeya set in F_3^3, coordinates (c_0, c_1, c_2); it is the smallest possible.
K33 = [(0, 0, 0), (1, 0, 0), (2, 0, 0), (1, 1, 0), (0, 0, 1), (2, 0, 1), (0, 1, 1), (1, 1, 1),
       (0, 2, 1), (2, 2, 1), (0, 0, 2), (2, 1, 2), (0, 2, 2)]


def is_prime(m):
    if m < 2:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return True


def _digits(a, p, k):
    out = []
    for _ in range(k):
        out.append(a % p)
        a //= p
    return out


def _undigits(ds, p):
    a = 0
    for c in reversed(ds):
        a = a * p + c
    return a


def _poly_mulmod(a, b, p, k, mod):
    prod = [0] * (2 * k - 1)
    for i in range(k):
        if a[i]:
            for j in range(k):
                prod[i + j] = (prod[i + j] + a[i] * b[j]) % p
    for d in range(2 * k - 2, k - 1, -1):
        c = prod[d]
        if c:
            for i in range(k + 1):
                prod[d - k + i] = (prod[d - k + i] - c * mod[i]) % p
    return prod[:k]


class GF:
    """The field with q elements, elements 0..q-1, with full addition and multiplication tables."""

    def __init__(self, q):
        self.q = q
        if is_prime(q):
            p, k, mod = q, 1, None
        elif q in IRREDUCIBLE:
            p, k, mod = IRREDUCIBLE[q]
        else:
            raise ValueError("unsupported field size %r" % q)
        self.p, self.k = p, k
        if k == 1:
            self.add = [[(a + b) % p for b in range(q)] for a in range(q)]
            self.mul = [[(a * b) % p for b in range(q)] for a in range(q)]
        else:
            ds = [_digits(a, p, k) for a in range(q)]
            self.add = [[_undigits([(x + y) % p for x, y in zip(ds[a], ds[b])], p)
                         for b in range(q)] for a in range(q)]
            self.mul = [[_undigits(_poly_mulmod(ds[a], ds[b], p, k, mod), p)
                         for b in range(q)] for a in range(q)]
        self.neg = [self.add[a].index(0) for a in range(q)]


class AffineSpace:
    """AG(n, q): points, normalised directions and all lines."""

    def __init__(self, n, q):
        self.n, self.q = n, q
        self.F = GF(q)
        self.size = q ** n
        self.coords = [tuple((i // q ** j) % q for j in range(n)) for i in range(self.size)]
        self.directions = [c for c in self.coords if any(c) and c[_first_nonzero(c)] == 1]
        self.lines = [self._lines_in_direction(d) for d in self.directions]

    def encode(self, c):
        return sum(x * self.q ** j for j, x in enumerate(c))

    def point_on_line(self, b, t, d):
        add, mul = self.F.add, self.F.mul
        return tuple(add[b[j]][mul[t][d[j]]] for j in range(self.n))

    def line(self, b, d):
        return frozenset(self.encode(self.point_on_line(b, t, d)) for t in range(self.q))

    def _lines_in_direction(self, d):
        seen = set()
        out = []
        for b in self.coords:
            L = self.line(b, d)
            if L not in seen:
                seen.add(L)
                out.append(L)
        return sorted(out, key=lambda L: sorted(L))

    def unit(self, j):
        return tuple(1 if i == j else 0 for i in range(self.n))


def _first_nonzero(c):
    for i, x in enumerate(c):
        if x:
            return i
    return -1


def is_kakeya(points, space):
    """True if the set of encoded points contains a full line in every direction."""
    S = set(points)
    return all(any(L <= S for L in lines) for lines in space.lines)


def missing_directions(points, space):
    """Directions (as coordinate tuples) in which the set contains no full line."""
    S = set(points)
    return [d for d, lines in zip(space.directions, space.lines) if not any(L <= S for L in lines)]


def is_difference_cover_f2(points, n):
    """True if every nonzero vector of F_2^n is a sum a + b of two points of the set.

    Over F_2 the line through x in direction d is the pair {x, x + d}, so this is exactly the
    Kakeya property in F_2^n. With the integer encoding, addition in F_2^n is XOR.
    """
    S = set(points)
    sums = {a ^ b for a in S for b in S}
    return all(v in sums for v in range(1, 2 ** n))


def symmetry_lines(space):
    """Lines that may be fixed without loss of generality (n >= 2).

    Translate a Kakeya set K so that its line in direction e1 is L0 = span(e1). Translating
    along e1 keeps L0 and moves K's line in direction e2 to {b + s*e2} with b_1 = b_2 = 0.
    If b = 0 that line is span(e2); otherwise a linear map fixing e1 and e2 sends b to e3,
    so the line becomes e3 + span(e2). Translations and linear maps preserve the size and the
    Kakeya property, so some minimum Kakeya set contains L0 and one of these two lines.
    Returns {direction index: allowed lines}.
    """
    zero = (0,) * space.n
    e1, e2 = space.unit(0), space.unit(1)
    allowed2 = [space.line(zero, e2)]
    if space.n >= 3:
        allowed2.append(space.line(space.unit(2), e2))
    return {space.directions.index(e1): [space.line(zero, e1)],
            space.directions.index(e2): allowed2}


def brute_force_minimum(space, allowed=None):
    """Minimum Kakeya set size by trying every choice of one line per direction.

    `allowed` optionally restricts the lines tried in some directions (see symmetry_lines).
    """
    allowed = allowed or {}
    best = space.size
    best_choice = None
    for choice in product(*[allowed.get(i, lines) for i, lines in enumerate(space.lines)]):
        u = frozenset().union(*choice)
        if best_choice is None or len(u) < best:
            best = len(u)
            best_choice = u
    return best, sorted(best_choice)


def line_masks(space):
    """Every direction's lines as bitmasks over the encoded points (bit p stands for point p)."""
    return [[sum(1 << p for p in L) for L in lines] for lines in space.lines]


def mask_points(mask):
    """The encoded points of a bitmask, in increasing order."""
    return [p for p in range(mask.bit_length()) if mask >> p & 1]


def kakeya_sets_through_axis(space, size):
    """Yield, as bitmasks, every Kakeya set of exactly `size` points that contains span(e1).

    Exhaustive and solver-free. Every Kakeya set has a translate containing span(e1) (move its
    line in direction e1 there), and adding points to a Kakeya set keeps it Kakeya, so AG(n, q)
    has a Kakeya set of at most `size` points if and only if this yields at least one set. It
    tries binom(q^n - q, size - q) sets, so it is practical only for small cases: proving
    K(3, 3) > 12 takes binom(24, 9) = 1,307,504 sets.
    """
    if size < space.q or size > space.size:
        return
    axis = space.line((0,) * space.n, space.unit(0))
    base = sum(1 << p for p in axis)
    rest = [1 << p for p in range(space.size) if p not in axis]
    e1 = space.directions.index(space.unit(0))
    checks = [lines for i, lines in enumerate(line_masks(space)) if i != e1]
    for combo in combinations(rest, size - space.q):
        s = base | sum(combo)
        for lines in checks:
            for L in lines:
                if (s & L) == L:
                    break
            else:
                break
        else:
            yield s


def has_kakeya_set_of_size(space, size):
    """True if AG(n, q) has a Kakeya set of at most `size` points (exhaustive; small cases only)."""
    return next(kakeya_sets_through_axis(space, size), None) is not None


def parabola_kakeya(n, q):
    """The parabola Kakeya set in AG(n, q), q odd. Returns (space, sorted encoded points).

    Write points as (x, t) with x in F^(n-1) and t the last coordinate. For each b in F^(n-1)
    take the line {(t*b_i - b_i^2)_i, t) : t in F}, which has direction (b, 1); together these
    lines cover every direction whose last coordinate is nonzero. Add the hyperplane t = 0 for
    the remaining directions (b, 0). Size: parabola_size(n, q).
    """
    if q % 2 == 0:
        raise ValueError("parabola_kakeya needs odd q")
    S = AffineSpace(n, q)
    add, mul, neg = S.F.add, S.F.mul, S.F.neg
    pts = {S.encode(c) for c in S.coords if c[-1] == 0}
    for b in product(range(q), repeat=n - 1):
        for t in range(q):
            x = tuple(add[mul[t][bi]][neg[mul[bi][bi]]] for bi in b)
            pts.add(S.encode(x + (t,)))
    return S, sorted(pts)


def parabola_size(n, q):
    """(q - 1) ((q + 1)/2)^(n - 1) + q^(n - 1): for t != 0 the slice of the parabola set at height t
    is the set of x with t^2/4 - x_i a square for every i, which has ((q + 1)/2)^(n - 1) points."""
    return (q - 1) * ((q + 1) // 2) ** (n - 1) + q ** (n - 1)


def smyth_kakeya(n, q):
    """The set (F^*)^n together with the 0/1 vectors, in AG(n, q). Returns (space, sorted points).

    For a direction d with zero coordinates Z, the line through the 0/1 vector 1_Z in direction
    d has nonzero coordinates everywhere except at that one point. Size: smyth_size(n, q).
    """
    S = AffineSpace(n, q)
    pts = {S.encode(c) for c in S.coords if all(c) or all(x in (0, 1) for x in c)}
    return S, sorted(pts)


def smyth_size(n, q):
    """(q - 1)^n + 2^n - 1: the all-ones vector is counted in both parts."""
    return (q - 1) ** n + 2 ** n - 1


def product_set(points_a, n, points_b, q):
    """Encoded points of A x B in AG(n + m, q), where A lies in AG(n, q) and B in AG(m, q).

    If A and B are Kakeya sets, so is A x B (section 6 of the README).
    """
    shift = q ** n
    return sorted(a + shift * b for a in points_a for b in points_b)


def plane_minimum(q):
    """Known minimum size of a Kakeya set in AG(2, q) (Blokhuis-Mazzocca; even q: q(q+1)/2)."""
    if q % 2 == 0:
        return q * (q + 1) // 2
    return q * (q + 1) // 2 + (q - 1) // 2


def counting_bound(n, q):
    """Least s with s(s - 1) >= q(q^n - 1).

    A Kakeya set contains (q^n - 1)/(q - 1) lines in distinct directions; two distinct lines share
    at most one point, so their q(q - 1)/2 pairs of points are disjoint sets of pairs.
    """
    need = q * (q ** n - 1)
    s = 1
    while s * (s - 1) < need:
        s += 1
    return s


def dvir_bound(n, q):
    """Dvir's polynomial-method lower bound: binom(q + n - 1, n)."""
    return comb(q + n - 1, n)


def bukh_chao_bound(n, q):
    """Bukh-Chao lower bound q^n / (2 - 1/q)^(n-1), rounded up to an integer."""
    v = Fraction(q ** (2 * n - 1), (2 * q - 1) ** (n - 1))
    return -((-v.numerator) // v.denominator)


def load_certificates(path):
    """Read a solver results file (data/kakeya-min-v1/results.json).

    Returns a list of (n, q, best_found, encoded points of the smallest set found).
    """
    with open(path) as f:
        cases = json.load(f)
    out = []
    for case in cases:
        n, q = case["n"], case["q"]
        pts = [sum(x * q ** j for j, x in enumerate(c)) for c in case["certificate"]]
        out.append((n, q, case["best_found"], pts))
    return out


def read_points(path, space):
    pts = []
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].split()
            if not line:
                continue
            c = tuple(int(x) for x in line)
            if len(c) != space.n or not all(0 <= x < space.q for x in c):
                raise ValueError("bad point %r" % (c,))
            pts.append(space.encode(c))
    return pts


def main(argv):
    if argv[:1] == ["records"]:
        print("(n, q)   lower  upper  (MV2 computation kakeya-min-v1; see README section 7)")
        for (n, q), (lo, hi) in sorted(COMPUTED.items()):
            print("(%d, %d)   %5d  %5d%s" % (n, q, lo, hi, "  exact" if lo == hi else ""))
        return 0
    if len(argv) < 3 or argv[0] not in ("check", "bounds", "parabola", "smyth", "brute", "exhaust"):
        print(__doc__)
        return 2
    cmd, n, q = argv[0], int(argv[1]), int(argv[2])
    if cmd == "bounds":
        print("AG(%d,%d): %d points, %d directions" % (n, q, q ** n, (q ** n - 1) // (q - 1)))
        print("counting bound   %d" % counting_bound(n, q))
        print("Dvir bound       %d" % dvir_bound(n, q))
        print("Bukh-Chao bound  %d" % bukh_chao_bound(n, q))
        if q % 2:
            print("parabola set     %d" % parabola_size(n, q))
        print("Smyth set        %d" % smyth_size(n, q))
        if n == 2:
            print("known minimum    %d" % plane_minimum(q))
        if q == 2 and n in KNOWN_F2:
            print("known minimum    %d" % KNOWN_F2[n])
        if (n, q) in COMPUTED:
            lo, hi = COMPUTED[(n, q)]
            print("MV2 computation  %d..%d" % (lo, hi))
        return 0
    if cmd in ("parabola", "smyth"):
        S, pts = parabola_kakeya(n, q) if cmd == "parabola" else smyth_kakeya(n, q)
        for p in pts:
            print(" ".join(str(x) for x in S.coords[p]))
        return 0
    S = AffineSpace(n, q)
    if cmd == "brute":
        best, pts = brute_force_minimum(S, symmetry_lines(S) if n >= 2 else None)
        print("minimum %d" % best)
        for p in pts:
            print(" ".join(str(x) for x in S.coords[p]))
        return 0
    if cmd == "exhaust":
        k = int(argv[3])
        s = next(kakeya_sets_through_axis(S, k), None)
        if s is None:
            print("no Kakeya set of at most %d points in AG(%d,%d)" % (k, n, q))
            return 1
        print("a Kakeya set of %d points in AG(%d,%d):" % (k, n, q))
        for p in mask_points(s):
            print(" ".join(str(x) for x in S.coords[p]))
        return 0
    pts = read_points(argv[3], S)
    miss = missing_directions(pts, S)
    print("%d distinct points; %s" % (len(set(pts)), "Kakeya set" if not miss else
                                       "NOT a Kakeya set, missing directions %s" % miss))
    return 0 if not miss else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
