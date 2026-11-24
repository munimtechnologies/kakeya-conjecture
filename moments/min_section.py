# min_section.py -- every set of N points in F_3^4 has a small hyperplane section.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.
#
# Identity. For any set K of N points in F_3^4,
#     sum over the 120 affine hyperplanes H of |K n H|^2 = 40 N + 13 N (N - 1),
# because each point lies in 40 hyperplanes (one per parallel class) and each pair of distinct
# points lies in exactly 13 (the hyperplanes containing the line through them: one in each class
# whose direction space contains the line's direction, and there are 13 such 3-dimensional
# subspaces of F_3^4 through a 1-dimensional one).
#
# Lemma. If every one of the 40 parallel classes had all three sections >= m, each class would
# contribute at most 2 m^2 + (N - 2m)^2 to the sum (a convex function of the three sizes is
# largest at a vertex of {a + b + c = N, a, b, c >= m}). So if 40 (2 m^2 + (N - 2m)^2) is smaller
# than the identity's value, some section has fewer than m points.
#   N = 25: sum = 8800; with m = 7 the cap is 40 * (98 + 121) = 8760 < 8800, so some section <= 6.
#   N = 26: sum = 9490; with m = 8 the cap is 40 * (128 + 100) = 9120 < 9490, so some section <= 7.
# This bounds the smallest section n0 in the exact search kmin43 (v3).


def second_moment(N):
    return 40 * N + 13 * N * (N - 1)


def nmax_for_size(N):
    """Largest m for which the lemma does NOT exclude 'all sections >= m' (0 if N < 3)."""
    S = second_moment(N)
    best = 0
    for m in range(1, N // 3 + 1):
        if 40 * (2 * m * m + (N - 2 * m) ** 2) >= S:
            best = m
    return best


def points():
    return [(a, b, c, d) for a in range(3) for b in range(3) for c in range(3) for d in range(3)]


def directions():
    out = []
    for v in points():
        nz = [x for x in v if x]
        if nz and nz[0] == 1:
            out.append(v)
    return out


def section_sizes(K):
    """The 120 section sizes |K n H|, grouped by class: list of 40 triples."""
    out = []
    for u in directions():
        cnt = [0, 0, 0]
        for x in K:
            cnt[sum(a * b for a, b in zip(u, x)) % 3] += 1
        out.append(tuple(cnt))
    return out
