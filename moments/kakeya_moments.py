# kakeya_moments.py -- hyperplane-moment lower bounds for Kakeya sets in F_q^n.
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0.
#
# Standard library only, exact integer arithmetic only.
#
# The argument (moments/README.md has the full proof):
#   Let K be a Kakeya set in F_q^n with N points, and let n_H = |K n H| for every affine
#   hyperplane H. Hyperplanes fall into C = (q^n - 1)/(q - 1) parallel classes of q.
#   (T) Transversal condition. In each class, a line whose direction is not parallel to the
#       class meets each of its q hyperplanes once, and two of those points fix the line.
#       K has a line in each of the q^(n-1) such directions, so for any two hyperplanes
#       H, H' of one class: n_H * n_H' >= q^(n-1).
#   (M) Moment identities, from counting pairs and triples of points of K by the hyperplanes
#       that contain them:
#         sum_H C(n_H, 2) = mu1 * C(N, 2)
#         sum_H C(n_H, 3) = mu2 * C(N, 3) + (mu1 - mu2) * T
#       where mu1 = (q^(n-1) - 1)/(q - 1) hyperplanes contain a given line,
#             mu2 = (q^(n-2) - 1)/(q - 1) hyperplanes contain a given affine plane,
#       and T >= C * C(q, 3) is the number of collinear triples of K (one line per direction,
#       q points on a line, and two lines share at most one point).
#   Write each class as its "profile" (the sorted list of its q section sizes). The
#   profiles of the C classes are C points (s2, s3) = (sum C(n_i,2), sum C(n_i,3)) whose
#   average must have s2 exactly mu1*C(N,2)/C and s3 at least (mu2*C(N,3) + (mu1-mu2)*Tmin)/C.
#   If no convex combination of allowed profile points does that, no Kakeya set of N points
#   exists. (That is a linear programme with two equality rows, so an optimal mix uses at
#   most two profiles; feasible() checks single profiles and pairs.)

import argparse
import json
from math import comb


def gauss1(q, k):
    """Number of 1-dimensional subspaces of F_q^k (0 for k = 0)."""
    return (q ** k - 1) // (q - 1)


def profiles(q, n, N):
    """All sorted tuples of q section sizes, each in [1, q^(n-1)], summing to N, with every
    pairwise product >= q^(n-1) (condition (T))."""
    cap = q ** (n - 1)
    out = []

    def rec(prefix, remaining, lo, parts_left):
        if parts_left == 1:
            v = remaining
            if lo <= v <= cap and all(p * v >= cap for p in prefix):
                out.append(tuple(prefix + [v]))
            return
        v = max(lo, 1)
        while v * parts_left <= remaining and v <= cap:
            if all(p * v >= cap for p in prefix):
                prefix.append(v)
                rec(prefix, remaining - v, v, parts_left - 1)
                prefix.pop()
            v += 1

    rec([], N, 1, q)
    return out


def moments(profile):
    return (sum(comb(x, 2) for x in profile), sum(comb(x, 3) for x in profile))


def targets(q, n, N):
    """(C, S2, S3min): number of classes, exact sum of C(n_H,2), lower bound on sum of C(n_H,3)."""
    C = gauss1(q, n)
    mu1 = gauss1(q, n - 1)
    mu2 = gauss1(q, n - 2)
    tmin = C * comb(q, 3)
    return C, mu1 * comb(N, 2), mu2 * comb(N, 3) + (mu1 - mu2) * tmin


def feasible(q, n, N):
    """True if the moment conditions can be met by C profiles (real mixing weights).
    False is a proof that no Kakeya set of exactly N points exists in F_q^n."""
    C, S2, S3 = targets(q, n, N)
    best = {}
    for P in profiles(q, n, N):
        s2, s3 = moments(P)
        if s3 > best.get(s2, -1):
            best[s2] = s3
    if not best:
        return False
    lows, highs = [], []
    for s2, s3 in sorted(best.items()):
        if C * s2 == S2:
            if C * s3 >= S3:
                return True
        elif C * s2 < S2:
            lows.append((s2, s3))
        else:
            highs.append((s2, s3))
    for a2, a3 in lows:
        for b2, b3 in highs:
            d = b2 - a2
            # value of the segment at x = S2/C is a3 + (b3-a3)(x-a2)/d; compare with S3/C,
            # multiplied through by C*d > 0
            if C * a3 * d + (b3 - a3) * (S2 - C * a2) >= S3 * d:
                return True
    return False


def lower_bound(q, n):
    """Smallest N for which the moment conditions are feasible: every Kakeya set in F_q^n
    has at least this many points. (A Kakeya set of N points gives one of N+1, so excluding
    each N below the bound separately is exactly what is needed.)"""
    for N in range(1, q ** n + 1):
        if feasible(q, n, N):
            return N
    return None


def certificate(q, n, N, c_num, c_den):
    """Hand-checkable dual certificate for one excluded N: with c = c_num/c_den, return
    (max over profiles of c_den*s3 - c_num*s2, total needed, lines bound) where the moment
    identities force C * max >= c_den*S3(T) - c_num*S2. Returned 'lmax' is the largest
    number of lines (q = 3) or collinear triples consistent with the inequality."""
    C, S2, S3 = targets(q, n, N)
    mu1, mu2 = gauss1(q, n - 1), gauss1(q, n - 2)
    vals = [c_den * moments(P)[1] - c_num * moments(P)[0] for P in profiles(q, n, N)]
    top = max(vals)
    # c_den * (mu2*C(N,3) + (mu1-mu2)*T) - c_num*S2 <= C*top  =>  T <= (C*top + c_num*S2 - c_den*mu2*C(N,3)) / (c_den*(mu1-mu2))
    num = C * top + c_num * S2 - c_den * mu2 * comb(N, 3)
    den = c_den * (mu1 - mu2)
    return {"N": N, "c": f"{c_num}/{c_den}", "max_profile_value": top,
            "collinear_triples_at_most": num // den, "collinear_triples_needed": C * comb(q, 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    res = {}
    for q, n in [(3, 2), (3, 3), (3, 4), (3, 5), (4, 3), (5, 3)]:
        res[f"q={q},n={n}"] = lower_bound(q, n)
        print(f"q={q} n={n}: every Kakeya set has at least {res[f'q={q},n={n}']} points")
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
