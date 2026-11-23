# sd_bound.py -- hyperplane moments plus a sum-difference lemma: lower bounds for the smallest
# Kakeya set k_n in F_3^n.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.
#
# Setting. K is a Kakeya set in F_3^n with N points: it contains a full line (3 points) in each of
# the H = (3^n - 1)/2 directions. L >= H is the number of lines contained in K. The hyperplanes of
# F_3^n fall into H parallel classes of 3; for a class, its profile (n1, n2, n3) lists |K cap H_i|.
#
# Moment identities (each pair of points lies on one line, each line lies in a1 hyperplane classes;
# a non-collinear triple spans a plane, which lies in a2 classes; collinear triples are the L lines):
#   sum over classes of C(n1,2)+C(n2,2)+C(n3,2) = a1 C(N,2)
#   sum over classes of C(n1,3)+C(n2,3)+C(n3,3) = a2 C(N,3) + (a1 - a2) L
# with a1 = (3^(n-1) - 1)/2 and a2 = (3^(n-2) - 1)/2.
#
# Per-class conditions (T = 3^(n-1) = number of directions NOT parallel to the class):
#  (T)  n_i n_j >= T for every two sections: each non-parallel direction has a line of K, which meets
#       H_i and H_j in one point each, and distinct directions give distinct pairs.
#  (SD) n_k >= ceil(T^2 / Emax(n_i, n_j)) for every choice of two sections i, j and the third k.
#       Proof (see SD.md): pick one line of K per non-parallel direction; its points p in H_i and q in
#       H_j give a set R of T pairs whose differences q - p are distinct and cover all T possible
#       values. The third point of the line is -(p + q), so the sums p + q of R take at most n_k
#       values, and Cauchy-Schwarz gives T^2 <= n_k * #{(r, r') in R^2 : equal sums}. That count is at
#       most the additive energy E(S_i, S_j) = #{p + q = p' + q'}, which equals the difference energy
#       #{q - p = q' - p'} = sum_d r(d)^2, where r(d) = #{(p, q) : q - p = d} satisfies
#       1 <= r(d) <= min(n_i, n_j) and sum_d r(d) = n_i n_j. Emax is the largest such sum of squares.
#  (HZ) only for n = 4: the 13 directions parallel to a class are each carried by a line inside one
#       section, a subset of a copy of AG(3,3); a section of m <= 8 points carries at most D3(m)
#       directions (computed exhaustively below), so D3(n1) + D3(n2) + D3(n3) >= 13.
#
# Linear relaxation. Averaging over the H classes, a Kakeya set of size N needs a probability
# distribution on allowed profiles with mean s2 = a1 C(N,2) / H and mean s3 >= (a2 C(N,3) +
# (a1 - a2) H) / H. That exists iff the upper concave envelope of the points (s2, s3) reaches the
# required mean s3 at the required mean s2 (at most two profiles are needed). feasible() decides this
# in exact integer arithmetic. If it is infeasible, no Kakeya set of size N exists.

from itertools import product
from math import comb


def flat_counts(n):
    """(H, a1, a2): directions (= hyperplane classes) of F_3^n, classes containing a given direction,
    classes containing a given 2-dimensional direction space."""
    H = (3 ** n - 1) // 2
    a1 = (3 ** (n - 1) - 1) // 2
    a2 = (3 ** (n - 2) - 1) // 2 if n >= 2 else 0
    return H, a1, a2


def energy_cap(x, y, T):
    """Largest sum of r(d)^2 over T values d with 1 <= r(d) <= min(x, y) and sum r(d) = x*y.
    None if no such r exists."""
    m = min(x, y)
    total = x * y
    if total < T or total > T * m:
        return None
    e = total - T
    if e == 0:
        return T
    cap = m - 1  # >= 1 here, because T < total <= T*m forces m >= 2
    q, rem = divmod(e, cap)
    ones = T - q - (1 if rem else 0)
    return ones + q * m * m + ((1 + rem) ** 2 if rem else 0)


def sd_floor(x, y, T):
    """Least possible size of the third section when two sections have x and y points (condition SD).
    None if x, y are impossible (x*y < T)."""
    E = energy_cap(x, y, T)
    if E is None:
        return None
    return -(-(T * T) // E)


# ---------------------------------------------------------------- AG(3,3): directions in small sets

def ag33_families():
    """fam[d] = sorted list of the 9 lines (bitmasks over 27 points) of AG(3,3) in direction d."""
    pts = list(product(range(3), repeat=3))
    idx = {p: i for i, p in enumerate(pts)}
    dirs = [v for v in pts if any(v) and next(c for c in v if c) == 1]
    fam = []
    for v in dirs:
        s = set()
        for p in pts:
            m = 0
            for t in range(3):
                m |= 1 << idx[tuple((p[i] + t * v[i]) % 3 for i in range(3))]
            s.add(m)
        fam.append(sorted(s))
    return fam


def _pc(x):
    return bin(x).count("1")


def min_union(k, fam=None):
    """f(k): least size of a union of k lines of AG(3,3) in pairwise distinct directions, 1 <= k <= 13.
    Exhaustive. The configuration may be assumed to contain the line through the origin in direction 0
    (the affine group of F_3^3 is transitive on lines and permutes directions); the other lines are
    taken in increasing direction order."""
    fam = fam or ag33_families()
    nd = len(fam)
    start = next(m for m in fam[0] if m & 1)
    best = [3 * k + 1]

    def rec(mask, count, d):
        if count == k:
            u = _pc(mask)
            if u < best[0]:
                best[0] = u
            return
        if count + (nd - d) < k:
            return
        for dd in range(d, nd):
            for m in fam[dd]:
                nm = mask | m
                if _pc(nm) < best[0]:
                    rec(nm, count + 1, dd + 1)

    rec(start, 1, 1)
    return best[0]


def d3_table(kmax=5):
    """f(1..kmax) and the map m -> D3(m) for m <= f(kmax) - 1 (largest number of directions carried
    by an m-point subset of AG(3,3)); f is nondecreasing, so D3(m) = max{k : f(k) <= m} there."""
    fam = ag33_families()
    f = [0] + [min_union(k, fam) for k in range(1, kmax + 1)]
    limit = f[kmax] - 1
    table = {}
    for m in range(0, limit + 1):
        table[m] = max(k for k in range(0, kmax + 1) if f[k] <= m)
    return f[1:], table


_D3 = None


def D3(m):
    """Upper bound for the number of directions carried by m points of AG(3,3): exact for m <= 8,
    the trivial 13 above."""
    global _D3
    if _D3 is None:
        _D3 = d3_table(5)[1]
    return _D3[m] if m in _D3 else 13


# ---------------------------------------------------------------- the linear relaxation

def profile_ok(p, n, sd=True, hz=False):
    T = 3 ** (n - 1)
    if any(c < 0 or c > T for c in p):
        return False
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        x, y, z = p[i], p[j], p[k]
        if x * y < T:
            return False
        if sd:
            f = sd_floor(x, y, T)
            if f is None or z < f:
                return False
    if hz:
        if n != 4:
            raise ValueError("condition HZ is implemented for n = 4 only")
        if D3(p[0]) + D3(p[1]) + D3(p[2]) < 13:
            return False
    return True


def profiles(n, N, sd=True, hz=False):
    out = []
    for a in range(0, N // 3 + 1):
        for b in range(a, (N - a) // 2 + 1):
            c = N - a - b
            if c >= b and profile_ok((a, b, c), n, sd, hz):
                out.append((a, b, c))
    return out


def moments_of(p):
    return sum(comb(c, 2) for c in p), sum(comb(c, 3) for c in p)


def feasible(n, N, sd=True, hz=False):
    """True iff the linear relaxation admits a Kakeya set of size N in F_3^n (exact)."""
    H, a1, a2 = flat_counts(n)
    X = a1 * comb(N, 2)                       # H * required mean of s2
    Y = a2 * comb(N, 3) + (a1 - a2) * H       # H * least required mean of s3 (L = H)
    best = {}
    for p in profiles(n, N, sd, hz):
        s2, s3 = moments_of(p)
        if s3 > best.get(s2, -1):
            best[s2] = s3
    pts = sorted(best.items())
    for x1, y1 in pts:
        if H * x1 > X:
            break
        if H * x1 == X and H * y1 >= Y:
            return True
        for x2, y2 in pts:
            if H * x2 <= X:
                continue
            if (H * y1 - Y) * (x2 - x1) + (y2 - y1) * (X - H * x1) >= 0:
                return True
    return False


def threshold(n, sd=True, hz=False, stop=10000):
    """Least N for which the relaxation is feasible; every Kakeya set in F_3^n has at least this many
    points."""
    for N in range(1, stop):
        if feasible(n, N, sd, hz):
            return N
    raise RuntimeError("no feasible N below stop")


def certificate_max(n, N, alpha, beta, sd=True, hz=False):
    """max over allowed profiles of alpha*s3 - beta*s2 (an explicit dual certificate check)."""
    return max(alpha * moments_of(p)[1] - beta * moments_of(p)[0] for p in profiles(n, N, sd, hz))


def certificate_refutes(n, N, alpha, beta, sd=True, hz=False):
    """True if summing alpha*s3 - beta*s2 over the H classes, with L = H lines (the fewest), already
    exceeds H * (max over allowed profiles); then no Kakeya set of size N exists (alpha > 0)."""
    H, a1, a2 = flat_counts(n)
    lhs = alpha * (a2 * comb(N, 3) + (a1 - a2) * H) - beta * a1 * comb(N, 2)
    return alpha > 0 and lhs > H * certificate_max(n, N, alpha, beta, sd, hz)


if __name__ == "__main__":
    print("f(1..5) in AG(3,3):", d3_table(5)[0])
    for n in (2, 3, 4, 5, 6):
        t_only = threshold(n, sd=False)
        with_sd = threshold(n, sd=True)
        line = "n=%d: moments+(T) give k_n >= %d; with (SD) k_n >= %d" % (n, t_only, with_sd)
        if n == 4:
            line += "; with (SD)+(HZ) k_4 >= %d" % threshold(4, sd=True, hz=True)
        print(line)
