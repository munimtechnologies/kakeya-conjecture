# bad25.py -- the hyperplane cuts a 25-point Kakeya set in F_3^4 must have.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.
#
# K is a Kakeya set in F_3^4 with N = 25 points and L >= 40 lines. For each of the 40 parallel
# classes of hyperplanes, its cut (profile) (n1, n2, n3) lists the sizes of K's three sections;
# s2 = sum C(n_i, 2) and s3 = sum C(n_i, 3). For EVERY 25-point set (Kakeya or not):
#     sum over classes of s2 = 13 C(25, 2) = 3900,
#     sum over classes of s3 = 4 C(25, 3) + 9 L = 9200 + 9 L,
# so with delta(cut) = 7 s3 - 51 s2 + 3302,
#     sum over the 40 classes of delta = 63 L - 2420.                                   (*)
# A Kakeya set has L >= 40, so the sum is at least 100. Conditions (T) and (SD) of SD.md leave 14
# possible cuts (allowed_profiles(25)); delta is positive on only three of them:
#     (6,7,12): 25,   (5,7,13): 60,   (6,6,13): 76,
# zero on (7,7,11) and (5,8,12), and negative on the other nine. Hence every 25-point Kakeya set has
# at least two classes with one of these three cuts, and if no section has 13 points, at least four
# classes with cut (6,7,12).

from itertools import product
from math import comb

T = 27
HAND25 = [(4, 9, 12), (4, 10, 11), (5, 7, 13), (5, 8, 12), (5, 9, 11), (5, 10, 10), (6, 6, 13),
          (6, 7, 12), (6, 8, 11), (6, 9, 10), (7, 7, 11), (7, 8, 10), (7, 9, 9), (8, 8, 9)]
BAD25 = [(6, 7, 12), (5, 7, 13), (6, 6, 13)]

# N = 26: sum s2 = 13 C(26,2) = 4225, sum s3 = 4 C(26,3) + 9 L = 10400 + 9 L, so with
# delta26(cut) = s3 - 8 s2 + 576 the sum over the 40 classes is 9 (L - 40) >= 0.
HAND26 = [(4, 9, 13), (4, 10, 12), (4, 11, 11), (5, 7, 14), (5, 8, 13), (5, 9, 12), (5, 10, 11),
          (6, 6, 14), (6, 7, 13), (6, 8, 12), (6, 9, 11), (6, 10, 10), (7, 7, 12), (7, 8, 11),
          (7, 9, 10), (8, 8, 10), (8, 9, 9)]
BAD26 = [(7, 7, 12), (6, 7, 13), (5, 7, 14), (6, 6, 14)]


def energy_cap(x, y, T=T):
    """Largest sum of r(d)^2 over T values with 1 <= r(d) <= min(x, y) and sum r(d) = x y
    (None if impossible). Same definition as sd_bound.energy_cap, re-implemented here."""
    m = min(x, y)
    total = x * y
    if total < T or total > T * m:
        return None
    e = total - T
    if e == 0:
        return T
    q, rem = divmod(e, m - 1)
    ones = T - q - (1 if rem else 0)
    return ones + q * m * m + ((1 + rem) ** 2 if rem else 0)


def sd_floor(x, y, T=T):
    E = energy_cap(x, y, T)
    if E is None:
        return None
    return -(-(T * T) // E)


def allowed_profiles(N):
    """Sorted cuts (a <= b <= c, a + b + c = N) that pass (T) and (SD) for every pair."""
    out = []
    for a in range(0, N // 3 + 1):
        for b in range(a, (N - a) // 2 + 1):
            c = N - a - b
            if c < b:
                continue
            p = (a, b, c)
            ok = True
            for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
                f = sd_floor(p[i], p[j])
                if f is None or p[k] < f:
                    ok = False
                    break
            if ok:
                out.append(p)
    return out


def moments_of(p):
    return sum(comb(c, 2) for c in p), sum(comb(c, 3) for c in p)


def delta(p):
    s2, s3 = moments_of(p)
    return 7 * s3 - 51 * s2 + 3302


def delta26(p):
    s2, s3 = moments_of(p)
    return s3 - 8 * s2 + 576


# ---------------------------------------------------------------- brute force on actual sets

def points():
    return list(product(range(3), repeat=4))


def directions():
    return [v for v in points() if any(v) and next(c for c in v if c) == 1]


def cuts(K):
    """The 40 cuts of K (unsorted), one per direction u: sizes of {x : u.x = c}, c = 0, 1, 2."""
    out = []
    for u in directions():
        cnt = [0, 0, 0]
        for x in K:
            cnt[sum(a * b for a, b in zip(u, x)) % 3] += 1
        out.append(tuple(cnt))
    return out


def count_lines(K):
    Ks = set(K)
    n = 0
    for v in directions():
        for x in Ks:
            if tuple((a + b) % 3 for a, b in zip(x, v)) in Ks and \
                    tuple((a + 2 * b) % 3 for a, b in zip(x, v)) in Ks:
                n += 1
    return n // 3  # each line is found once from each of its 3 points


def is_kakeya(K):
    Ks = set(K)
    for v in directions():
        if not any(tuple((a + b) % 3 for a, b in zip(x, v)) in Ks and
                   tuple((a + 2 * b) % 3 for a, b in zip(x, v)) in Ks for x in Ks):
            return False
    return True


# The 27-point Kakeya set of kakeya-min-v1 (data/kakeya-min-v1/results.json in this repository).
KNOWN27 = [(0, 0, 0, 0), (1, 0, 0, 0), (2, 0, 0, 0), (2, 2, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0),
           (0, 1, 1, 0), (1, 1, 1, 0), (0, 2, 1, 0), (2, 2, 1, 0), (0, 0, 2, 0), (0, 1, 2, 0),
           (1, 2, 2, 0), (1, 2, 0, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 2, 1, 1), (0, 0, 2, 1),
           (2, 1, 2, 1), (2, 0, 0, 2), (0, 1, 0, 2), (0, 2, 0, 2), (1, 0, 1, 2), (2, 0, 1, 2),
           (1, 2, 1, 2), (0, 1, 2, 2), (2, 1, 2, 2)]


if __name__ == "__main__":
    for p in allowed_profiles(25):
        print(p, moments_of(p), delta(p), "BAD" if p in BAD25 else "")
