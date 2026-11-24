# plane_lp.py -- a plane-class relaxation for Kakeya sets in F_3^4 (3x3 arrays of plane sections).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.
#
# Setting (see PLANE.md). K is a Kakeya set in F_3^4 with N points and L >= 40 lines. For each of the
# 130 two-dimensional direction spaces W, the 9 planes parallel to W cut K into a 3x3 array M of
# section sizes (cells = cosets of W, indexed by (u.x, v.x) for two normals u, v spanning W-perp).
# The 12 lines of this AG(2,3) are the 12 hyperplanes of the 4 hyperplane classes containing W, so
# every line sum of M is a section size of a hyperplane class.
#
# Exact identities (sums over the 130 classes W; C_k = C(N,k)):
#   (I2) sum c2 = 13 C_2                       c_k(M) = sum over cells of C(m,k)
#   (I3) sum c3 = C_3 + 12 L                   (a non-collinear triple spans one 2-flat direction)
#   (H3) sum h3 = 13 (4 C_3 + 9 L)             h_k(M) = sum over the 12 lines of C(line sum, k)
#   (I4) sum (h4 - 39 c4) = 13 C_4             (4 points span a 3-flat, in one hyperplane, or lie in a
#                                               2-flat, in 4 hyperplanes; each hyperplane class lies
#                                               in 13 of the classes W)
#   (Ik) sum (hk - 39 ck) <= 13 C_k, k >= 5    (k points may also span all of F_3^4)
#   (LW) sum linecap >= 13 L                   lines of K with direction in W lie inside cells; a cell
#                                               of m points holds at most ML2[m] lines
#   (LX) sum crosscap >= 117 L                 a line with direction not in W crosses three collinear
#                                               cells and is fixed by its points in two of them
# Per-array conditions (each is necessary for the array of a real Kakeya set):
#   every line-sum triple is an allowed hyperplane profile (sd_bound: (T), (SD), (HZ));
#   cells have 0..9 points; for each of the 4 line directions some line has 3 non-empty cells
#   (the image of K in F_3^4 / W is a Kakeya set of AG(2,3));
#   sum over cells of D2[m] >= 4 (the 4 directions of W are carried inside cells);
#   for each line direction, sum over its 3 lines of cross_cap >= 9 (the 9 directions of U \ W,
#   U = the hyperplane class of that direction, are carried by lines crossing 3 collinear cells).
#
# Relaxation: a Kakeya set of size N gives 130 arrays satisfying all of the above. plane_feasible
# decides, with an exact rational simplex (Bland's rule), whether nonnegative weights on the valid
# arrays (total 130) and some L >= Lmin satisfy the identities. If not, no such Kakeya set exists.

from fractions import Fraction
from itertools import permutations, product
from math import comb, gcd
from operator import mul

from moments import sd_bound as S

N_CLASSES = 130
KMAX = 13
D2 = (0, 0, 0, 1, 1, 2, 3, 4, 4, 4)     # most directions carried by m points of AG(2,3)
ML2 = (0, 0, 0, 1, 1, 2, 3, 5, 8, 12)   # most lines contained in m points of AG(2,3)


def _ag2_lines():
    out = []
    for v in ((0, 1), (1, 0), (1, 1), (1, 2)):
        fam = set()
        for i in range(3):
            for j in range(3):
                fam.add(tuple(sorted(3 * ((i + t * v[0]) % 3) + (j + t * v[1]) % 3
                                     for t in range(3))))
        out.append(tuple(sorted(fam)))
    return tuple(out)


LINES2 = _ag2_lines()   # 4 directions x 3 parallel lines, cells numbered 3*i + j

_PROFILE = {}


def profile_allowed(p):
    q = tuple(sorted(p))
    r = _PROFILE.get(q)
    if r is None:
        r = S.profile_ok(q, 4, sd=True, hz=True)
        _PROFILE[q] = r
    return r


def cross_cap(a, b, c):
    """Upper bound for the number of the 9 crossing directions that three parallel plane sections
    of sizes a, b, c inside one AG(3,3) can carry: distinct directions need distinct pairs in every
    two sections, and carrying all 9 requires the sum-difference condition with T = 9."""
    t = min(9, a * b, a * c, b * c)
    if t == 9:
        for x, y, z in ((a, b, c), (a, c, b), (b, c, a)):
            f = S.sd_floor(x, y, 9)
            if f is None or z < f:
                return 8
    return t


def array_ok(M):
    if any(m < 0 or m > 9 for m in M):
        return False
    if sum(D2[m] for m in M) < 4:
        return False
    for fam in LINES2:
        sums = [M[a] + M[b] + M[c] for a, b, c in fam]
        if not profile_allowed(sums):
            return False
        if not any(M[a] and M[b] and M[c] for a, b, c in fam):
            return False
        if sum(cross_cap(M[a], M[b], M[c]) for a, b, c in fam) < 9:
            return False
    return True


_CB = [[comb(s, k) for k in range(KMAX + 1)] for s in range(28)]


def raw_moments(M):
    """(c, h, linecap, crosscap) with c[k], h[k] for k = 0..KMAX."""
    c = [0] * (KMAX + 1)
    for m in M:
        row = _CB[m]
        for k in range(KMAX + 1):
            c[k] += row[k]
    h = [0] * (KMAX + 1)
    cross = 0
    for fam in LINES2:
        for a, b, d in fam:
            row = _CB[M[a] + M[b] + M[d]]
            for k in range(KMAX + 1):
                h[k] += row[k]
            cross += min(M[a] * M[b], M[a] * M[d], M[b] * M[d])
    return c, h, sum(ML2[m] for m in M), cross


def features(M):
    c, h, lc, xc = raw_moments(M)
    return ((c[2], c[3], h[3], h[4] - 39 * c[4])
            + tuple(h[k] - 39 * c[k] for k in range(5, KMAX + 1)) + (lc, xc))


def _comps(s):
    return [(x, y, s - x - y) for x in range(0, min(9, s) + 1)
            for y in range(0, min(9, s - x) + 1) if s - x - y <= 9]


_VALID = {}


def valid_arrays(N):
    """All valid arrays whose row sums are nondecreasing (every valid array is a row permutation,
    an affine map of AG(2,3), of one of these; features are invariant under affine maps)."""
    if N not in _VALID:
        _VALID[N] = _valid_arrays(N)
    return _VALID[N]


def _valid_arrays(N):
    rows = S.profiles(4, N, sd=True, hz=True)
    ordered = sorted({q for p in rows for q in permutations(p)})
    out = []
    for r0, r1, _ in rows:
        c0, c1 = _comps(r0), _comps(r1)
        for a in c0:
            for b in c1:
                for col in ordered:
                    z0 = col[0] - a[0] - b[0]
                    if z0 < 0 or z0 > 9:
                        continue
                    z1 = col[1] - a[1] - b[1]
                    if z1 < 0 or z1 > 9:
                        continue
                    z2 = col[2] - a[2] - b[2]
                    if z2 < 0 or z2 > 9:
                        continue
                    M = a + b + (z0, z1, z2)
                    if array_ok(M):
                        out.append(M)
    return out


_FEAT = {}


def feature_set(N):
    if N not in _FEAT:
        _FEAT[N] = sorted({features(M) for M in valid_arrays(N)})
    return _FEAT[N]


# ---------------------------------------------------------------- exact phase-I simplex

def phase1_feasible(cols, b):
    """Is there x >= 0 with sum_j x_j cols[j] = b?  cols: integer columns, b: integers.
    Exact revised simplex, Bland's rule (terminates)."""
    m = len(b)
    cols = [tuple(c) for c in cols]
    b = list(b)
    for i in range(m):
        if b[i] < 0:
            b[i] = -b[i]
            cols = [c[:i] + (-c[i],) + c[i + 1:] for c in cols]
    n = len(cols)
    allcols = cols + [tuple(int(i == k) for i in range(m)) for k in range(m)]
    cost = [0] * n + [1] * m
    basis = [n + i for i in range(m)]
    Binv = [[Fraction(int(i == j)) for j in range(m)] for i in range(m)]
    xB = [Fraction(v) for v in b]
    while True:
        cB = [cost[j] for j in basis]
        y = [sum((cB[i] * Binv[i][k] for i in range(m)), Fraction(0)) for k in range(m)]
        D = 1
        for v in y:
            D = D * v.denominator // gcd(D, v.denominator)
        Y = [int(v * D) for v in y]
        inb = set(basis)
        entering = None
        for j in range(n + m):
            if j in inb:
                continue
            if cost[j] * D - sum(map(mul, Y, allcols[j])) < 0:
                entering = j
                break
        if entering is None:
            return sum(xB[i] for i in range(m) if basis[i] >= n) == 0
        a = allcols[entering]
        u = [sum((Binv[i][k] * a[k] for k in range(m) if a[k]), Fraction(0)) for i in range(m)]
        r = None
        best = None
        for i in range(m):
            if u[i] > 0:
                ratio = xB[i] / u[i]
                if r is None or ratio < best or (ratio == best and basis[i] < basis[r]):
                    r, best = i, ratio
        if r is None:
            raise RuntimeError("phase I unbounded (impossible)")
        piv = u[r]
        Binv[r] = [v / piv for v in Binv[r]]
        xB[r] = xB[r] / piv
        for i in range(m):
            if i != r and u[i] != 0:
                f = u[i]
                Binv[i] = [Binv[i][k] - f * Binv[r][k] for k in range(m)]
                xB[i] = xB[i] - f * xB[r]
        basis[r] = entering


def lp_system(N, Lmin=40, F=None):
    """Columns and right-hand side of the relaxation (weights on feature vectors, s = L - Lmin,
    slacks)."""
    if F is None:
        F = feature_set(N)
    C = [comb(N, k) for k in range(KMAX + 1)]
    nk = KMAX - 4                                  # inequalities k = 5..KMAX
    m = 5 + nk + 2
    cols = []
    for f in F:
        c2, c3, h3, w4 = f[0], f[1], f[2], f[3]
        wk = f[4:4 + nk]
        lc, xc = f[4 + nk], f[5 + nk]
        cols.append((1, c2, c3, h3, w4) + tuple(wk) + (lc, xc))
    cols.append((0, 0, -12, -117, 0) + (0,) * nk + (-13, -117))           # s = L - Lmin
    for k in range(nk):                                                   # slacks of (Ik)
        e = [0] * m
        e[5 + k] = 1
        cols.append(tuple(e))
    for r in (5 + nk, 6 + nk):                                            # surplus of (LW), (LX)
        e = [0] * m
        e[r] = -1
        cols.append(tuple(e))
    b = ([N_CLASSES, 13 * C[2], C[3] + 12 * Lmin, 13 * (4 * C[3] + 9 * Lmin), 13 * C[4]]
         + [13 * C[k] for k in range(5, KMAX + 1)] + [13 * Lmin, 117 * Lmin])
    return cols, b


def plane_feasible(N, Lmin=40):
    """False means: no Kakeya set in F_3^4 with N points and at least Lmin lines exists."""
    F = feature_set(N)
    if not F:
        return False
    cols, b = lp_system(N, Lmin, F)
    return phase1_feasible(cols, b)


# ---------------------------------------------------------------- tools for real point sets

def _normals():
    return [v for v in product(range(3), repeat=4) if any(v) and next(c for c in v if c) == 1]


def _norm(w):
    lead = next(c for c in w if c)
    inv = 1 if lead == 1 else 2
    return tuple((c * inv) % 3 for c in w)


def plane_classes():
    """One pair (u, v) of normals per 2-dimensional space of normals (130 of them)."""
    seen, out = set(), []
    vs = _normals()
    for i, u in enumerate(vs):
        for v in vs[i + 1:]:
            key = frozenset([u, v, _norm(tuple((a + c) % 3 for a, c in zip(u, v))),
                             _norm(tuple((a + 2 * c) % 3 for a, c in zip(u, v)))])
            if key not in seen:
                seen.add(key)
                out.append((u, v))
    return out


def array_of(K, u, v):
    M = [0] * 9
    for x in K:
        i = sum(a * b for a, b in zip(u, x)) % 3
        j = sum(a * b for a, b in zip(v, x)) % 3
        M[3 * i + j] += 1
    return tuple(M)


def count_lines(K):
    Ks = set(K)
    pts = sorted(Ks)
    n = 0
    for i, p in enumerate(pts):
        for q in pts[i + 1:]:
            r = tuple((-a - b) % 3 for a, b in zip(p, q))
            if r in Ks and r > q:
                n += 1
    return n


def is_kakeya(K):
    Ks = set(K)
    for d in _normals():
        if not any(tuple((a + b) % 3 for a, b in zip(x, d)) in Ks
                   and tuple((a + 2 * b) % 3 for a, b in zip(x, d)) in Ks for x in Ks):
            return False
    return True


def sort_rows(M):
    rows = sorted([M[0:3], M[3:6], M[6:9]], key=sum)
    return tuple(rows[0]) + tuple(rows[1]) + tuple(rows[2])


# A 27-point Kakeya set of F_3^4 found by the kakeya-min-v1 integer program
# (data/kakeya-min-v1/results.json); slices by the last coordinate have 13, 6 and 8 points.
K27 = [(0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
       (0, 1, 1, 0), (1, 1, 1, 0), (0, 2, 1, 0), (2, 2, 1, 0), (0, 0, 2, 0), (0, 1, 2, 0),
       (1, 2, 2, 0), (1, 2, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (0, 0, 2, 1),
       (2, 1, 2, 1), (2, 0, 0, 2), (0, 1, 0, 2), (0, 2, 0, 2), (1, 0, 1, 2), (2, 0, 1, 2),
       (1, 2, 1, 2), (0, 1, 2, 2), (2, 1, 2, 2)]


if __name__ == "__main__":
    for N in (24, 25, 26, 27):
        print(N, len(feature_set(N)), plane_feasible(N))
