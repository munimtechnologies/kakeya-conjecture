# n24.py -- the case N = 24 for Kakeya sets in F_3^4 (see moments/N24.md).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only,
# exact rational arithmetic.

from fractions import Fraction
from math import comb

import kakeya_moments as km

# Upper bounds on D(m), the most directions of lines inside an m-point subset of AG(3,3)
# (hand lemma in N24.md, checked by ag33_lines.fits in test_n24.py). For m >= 9 we use the
# trivial bound 13, so the "covers" condition below only ever uses the proved values.
D_UPPER = {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 3, 7: 4, 8: 4}


def d_upper(m):
    return D_UPPER.get(m, 13)


def covers(profile):
    """Horizontal condition (H): the 13 directions parallel to a hyperplane class each carry a
    line of K, and such a line lies inside one of the three hyperplanes of the class. So the
    three sections together must carry lines in all 13 directions."""
    return sum(d_upper(x) for x in profile) >= 13


def best_avg_s3(profs, C, S2):
    """Largest average of s3 over C classes whose average s2 is exactly S2/C, using real mixing
    weights over the given profiles (an LP with two equality rows: singles and pairs suffice)."""
    x = Fraction(S2, C)
    pts = {}
    for P in profs:
        s2, s3 = km.moments(P)
        pts[s2] = max(pts.get(s2, -1), s3)
    best = None
    for a2, a3 in pts.items():
        if a2 == x:
            best = Fraction(a3) if best is None else max(best, Fraction(a3))
    for a2, a3 in pts.items():
        for b2, b3 in pts.items():
            if a2 < x < b2:
                v = a3 + Fraction(b3 - a3) * (x - a2) / (b2 - a2)
                best = v if best is None else max(best, v)
    return best


def max_lines(N, profs):
    """Most lines L a Kakeya set of N points in F_3^4 can contain if every hyperplane class has
    a profile from profs: sum C(n_H,3) = 4 C(N,3) + 9 L is at most C * best_avg_s3."""
    C, S2, _ = km.targets(3, 4, N)
    b = best_avg_s3(profs, C, S2)
    if b is None:
        return None
    return (C * b - 4 * comb(N, 3)) // 9


def deficit(profile):
    """18 s3 - 141 s2 measured from the line through (7,8,9) and (5,6,13); 0 on that line."""
    s2, s3 = km.moments(profile)
    return -8835 - (18 * s3 - 141 * s2)


def main():
    all24 = km.profiles(3, 4, 24)
    ok24 = [P for P in all24 if covers(P)]
    print("profiles for N = 24:", len(all24), "; failing (H):", [P for P in all24 if not covers(P)])
    print("lines at most, (T) only:", max_lines(24, all24), "; (T) and (H):", max_lines(24, ok24))
    for P in sorted(ok24, key=deficit):
        print(P, km.moments(P), "deficit", deficit(P))


if __name__ == "__main__":
    main()
