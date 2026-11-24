# plane_lp2.py -- a sharper plane-class relaxation for Kakeya sets in F_3^4 (see PLANE2.md).
#
# MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name.
# No human mathematician has checked this. Licence: Apache-2.0. Standard library only.
#
# It keeps everything of plane_lp.py (PLANE.md) and adds two exact facts about three parallel
# planes P1, P2, P3 inside one hyperplane (a copy of AG(3,3)). Identify each plane with F_3^2 and
# put P1, P2, P3 at heights 0, 1, 2. Points x in P1, y in P2, z in P3 are collinear exactly when
# x + y + z = 0 in F_3^2, and the line then has direction (y - x, 1). So for sections A, B, C:
#   * the lines of K that meet all three planes are triples (x, y, z) in A x B x C with x + y + z = 0;
#     their number is at most T3[|A|][|B|][|C|], the exact maximum over all subsets of those sizes;
#   * the directions they carry are the distinct differences y - x; their number is at most
#     G[|A|][|B|][|C|], again the exact maximum.
# Both tables are computed exhaustively (triple_tables). They replace min(ab, ac, bc) (the old LX
# capacity) and cross_cap (the old condition 4) of plane_lp.py.
#
# New per-array condition (lines per class). For a plane class W with array M, every line of K
# either has its direction in W (then it lies in a cell, at most ML2[m] lines per cell) or crosses
# three collinear cells (at most T3 lines per quotient line). Hence, for EVERY class,
#   linecap(M) + crosscap3(M) >= L >= Lmin,
# so arrays with too little room for lines cannot occur at all.
#
# plane2_feasible(N, Lmin) is the relaxation of plane_lp.plane_feasible with these changes. If it is
# infeasible, no Kakeya set in F_3^4 with N points and at least Lmin lines exists.

from math import comb

from moments import plane_lp as P


def _idx(a, b):
    return 3 * (a % 3) + (b % 3)


# F_3^2 on indices 0..8 (index 3*a + b is the vector (a, b)); ADD[i][j] = i + j, SUB[i][j] = j - i
ADD = [[_idx(i // 3 + j // 3, i % 3 + j % 3) for j in range(9)] for i in range(9)]
SUB = [[_idx(j // 3 - i // 3, j % 3 - i % 3) for j in range(9)] for i in range(9)]
POP = [bin(s).count("1") for s in range(512)]
BITS = [[i for i in range(9) if (s >> i) & 1] for s in range(512)]


def triple_tables():
    """(G, T3), each indexed [a][b][c] with 0 <= a, b, c <= 9.
    T3[a][b][c] = max #{(x, y, z) in A x B x C : x + y + z = 0} over A, B, C in F_3^2 of sizes a, b, c.
    G[a][b][c]  = max #{y - x : (x, y, z) such a triple} over the same sets.
    Exhaustive. Translating A by t, B by s and C by -(s + t) keeps the triples and shifts every
    difference by s - t, so A and B may be assumed to contain 0. Swapping A and B negates the
    differences, so pairs with B < A are covered by (B, A). For fixed A, B the best C of size k
    is the set of negatives of k sums, chosen among all 2^9 subsets of sums."""
    G = [[[0] * 10 for _ in range(10)] for _ in range(10)]
    T3 = [[[0] * 10 for _ in range(10)] for _ in range(10)]
    reps = [s for s in range(1, 512) if s & 1]
    U = [0] * 512
    for A in reps:
        a = POP[A]
        bitsA = BITS[A]
        for B in reps:
            if B < A:
                continue
            b = POP[B]
            bitsB = BITS[B]
            dm = [0] * 9
            cnt = [0] * 9
            for x in bitsA:
                addx = ADD[x]
                subx = SUB[x]
                for y in bitsB:
                    s = addx[y]
                    dm[s] |= 1 << subx[y]
                    cnt[s] += 1
            best = [0] * 10
            for S in range(1, 512):
                low = S & -S
                u = U[S ^ low] | dm[low.bit_length() - 1]
                U[S] = u
                k = POP[S]
                p = POP[u]
                if p > best[k]:
                    best[k] = p
            cs = sorted(cnt, reverse=True)
            tb = [0] * 10
            run = 0
            for k in range(1, 10):
                run += cs[k - 1]
                tb[k] = run
            for p, q in ((a, b), (b, a)):
                gpq = G[p][q]
                tpq = T3[p][q]
                for k in range(1, 10):
                    if best[k] > gpq[k]:
                        gpq[k] = best[k]
                    if tb[k] > tpq[k]:
                        tpq[k] = tb[k]
    return G, T3


_TAB = []


def tables():
    if not _TAB:
        _TAB.append(triple_tables())
    return _TAB[0]


def caps2(M):
    """(linecap, crosscap3): most lines of K inside cells, most lines of K crossing three cells."""
    T3 = tables()[1]
    lc = sum(P.ML2[m] for m in M)
    xc = 0
    for fam in P.LINES2:
        for a, b, c in fam:
            xc += T3[M[a]][M[b]][M[c]]
    return lc, xc


def extra_ok(M, Lmin):
    """The new per-array conditions (for arrays that already pass plane_lp.array_ok)."""
    G = tables()[0]
    for fam in P.LINES2:
        if sum(G[M[a]][M[b]][M[c]] for a, b, c in fam) < 9:
            return False
    lc, xc = caps2(M)
    return lc + xc >= Lmin


def array_ok2(M, Lmin=40):
    return P.array_ok(M) and extra_ok(M, Lmin)


def features2(M):
    c, h, lc, _ = P.raw_moments(M)
    lc2, xc = caps2(M)
    return ((c[2], c[3], h[3], h[4] - 39 * c[4])
            + tuple(h[k] - 39 * c[k] for k in range(5, P.KMAX + 1)) + (lc2, xc))


def valid_arrays2(N, Lmin=40):
    return [M for M in P.valid_arrays(N) if extra_ok(M, Lmin)]


def feature_set2(N, Lmin=40):
    return sorted({features2(M) for M in valid_arrays2(N, Lmin)})


def plane2_feasible(N, Lmin=40):
    """False means: no Kakeya set in F_3^4 with N points and at least Lmin lines exists."""
    F = feature_set2(N, Lmin)
    if not F:
        return False
    cols, b = P.lp_system(N, Lmin, F)
    return P.phase1_feasible(cols, b)


def real_set_guard():
    """Cheap guard run by every readout test: all 130 arrays of the real 27-point Kakeya set K27 pass
    every condition (old and new) with its true number of lines, and its line totals meet (LW), (LX)."""
    K = P.K27
    if not P.is_kakeya(K):
        return False
    L = P.count_lines(K)
    tot_lc = tot_xc = 0
    for u, v in P.plane_classes():
        M = P.array_of(K, u, v)
        if not array_ok2(M, L):
            return False
        lc, xc = caps2(M)
        tot_lc += lc
        tot_xc += xc
    return L >= 40 and tot_lc >= 13 * L and tot_xc >= 117 * L


def hand_values():
    """Small values checked by hand in PLANE2.md."""
    return {"T3(1,1,1)": 1, "T3(2,2,2)": 3, "G(2,2,2)": 3, "T3(3,3,3)": 9, "G(1,9,9)": 9,
            "T3(9,9,9)": 81, "G(9,9,9)": 9}


if __name__ == "__main__":
    G, T3 = tables()
    print("G(3,3,3) =", G[3][3][3], " T3(3,3,3) =", T3[3][3][3])
    for N in (24, 25, 26, 27):
        print(N, len(feature_set2(N)), plane2_feasible(N), comb(N, 2))
