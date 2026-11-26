# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 MV2 of Munim, Inc. (mv2@munimtech.com)
# Made with AI agents working in MV2's name. Not yet checked by a human mathematician.
"""Split extensions of Kakeya sets in F_3^m and a doubling lemma (standard library only).

Definitions.  Let K be a Kakeya set of F_3^m.  A *split* of K is a pair (A, B) with
A u B = K and

    C(A, B) := { b - a : a in A, b in B, -(a + b) in K }  =  F_3^m.

Its *excess* is e = |A n B|; e >= 1, because 0 is in C(A, B) only through a = b.

Lemma 1 (doubling).  If (A, B) is a split of K, then
    K+ = {0} x K  u  {1} x A  u  {2} x B
is a Kakeya set of F_3^(m+1) with |K| + |A| + |B| = 2|K| + e points.
Proof.  Directions (0, d) are covered by {0} x K.  A direction (1, v) is covered iff some
x in K has x + v in A and x + 2v in B.  Given a in A, b in B with x := -(a + b) in K and
v := b - a, we get x + v = -2a = a and x + 2v = b - 3a = b.  So the covered v are exactly
C(A, B) = F_3^m, and every direction of F_3^(m+1) is a multiple of some (0, d) or (1, v).

Lemma 2 (lifting).  Let K' = {0} x K u {1} x L1 u {2} x L2 be any Kakeya set of F_3^(m+1)
whose layer x_0 = 0 is K, and let (A, B) be a split of K.  Then
    A' = {0} x A u {1} x L1,   B' = {0} x B u {2} x L2
is a split of K' with the same excess.
Proof.  A' u B' = K' and A' n B' = {0} x (A n B).  For U = (0, u): take the a, b of K for u;
(0, a), (0, b) work, with third point (0, -(a + b)) in {0} x K.  For U = +-D, D = (1, v):
K' is Kakeya, so it holds a line p0 = (0, x), p1 = (1, x + v), p2 = (2, x + 2v), with
p1 in A' and p2 in B'.  Then p2 - p1 = D with third point p0.  For -D: if x is in A, take
(p0, p2) (p2 - p0 = -D, third point p1); if x is in B, take (p1, p0) (p0 - p1 = -D, third
point p2).  x is in A or B because A u B = K.

Corollary.  K+ is layered with layer K, so Lemma 2 gives it the split
(A u {1} x A, B u {2} x B) with excess e again, and Lemma 1 applies again: a split of K
with excess e gives Kakeya sets of 2^j (|K| + e) - e points in F_3^(m+j) for all j >= 0.

This module finds splits by local search (WalkSAT-style, over the colourings A only / B only
/ both of the points of K, constraints = the nonzero vectors u, each needing an ordered pair
(a, b) on a line of K with b - a = u), lifts them through stored layered sets (Lemma 2),
and builds the doubled sets (Lemma 1).  Every claim it makes is checked by verify_split(),
which uses the definition of C(A, B) directly (coordinate tuples, no search code).

Points are integers i = sum_k x_k 3^(m-1-k) (first coordinate most significant), as in
k6/sets/.  A split of K is stored as the code of its states over sorted(K): state 1 = A
only, 2 = B only, 3 = both, three states per character of ALPHA (16 s0 + 4 s1 + s2)."""
import json
import os
import random
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETS = os.path.join(_ROOT, "k6", "sets")
ALPHA = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"
_HI = "ABCDEFGHIJKLMNOPQRSTUVWXYZ6789+/"
_LO = "abcdefghijklmnopqrstuvwxyz012345"


# ---- sets --------------------------------------------------------------------------------
def _unpack(code):
    vals, cur, acc = [], -1, 0
    for ch in code:
        k = _LO.find(ch)
        if k < 0:
            acc = 32 * acc + _HI.index(ch)
        else:
            cur += 32 * acc + k
            vals.append(cur)
            acc = 0
    return vals


def load_set(name):
    """(n, sorted indices) of a stored set; two-part files are rebuilt recursively."""
    with open(os.path.join(SETS, os.path.basename(name)), encoding="utf-8") as f:
        d = json.load(f)
    if "points" in d:
        return d["n"], sorted(d["points"])
    _, base = load_set(d["layer0_file"])
    free = d["free_points"] if "free_points" in d else _unpack(d["points_code"])
    return d["n"], sorted(base + [d["offset"] + p for p in free])


def digits(i, m):
    out = []
    for _ in range(m):
        i, r = divmod(i, 3)
        out.append(r)
    return tuple(reversed(out))


def index_of(x):
    i = 0
    for c in x:
        i = 3 * i + c
    return i


# ---- split codes -------------------------------------------------------------------------
def encode_split(states):
    s = list(states) + [0] * (-len(states) % 3)
    return "".join(ALPHA[16 * s[i] + 4 * s[i + 1] + s[i + 2]] for i in range(0, len(s), 3))


def decode_split(code, count):
    out = []
    for ch in code:
        v = ALPHA.index(ch)
        out += [v // 16, (v // 4) % 4, v % 4]
    assert all(s == 0 for s in out[count:]) and len(out) - count < 3, "code length"
    return out[:count]


def split_from_states(K, states):
    """K sorted indices, states 1/2/3 -> (A, B) as sets of indices."""
    A = set(k for k, s in zip(K, states) if s & 1)
    B = set(k for k, s in zip(K, states) if s & 2)
    return A, B


def states_of(K, A, B):
    return [(1 if k in A else 0) + (2 if k in B else 0) for k in K]


# ---- independent verification (definition of C(A, B), coordinate tuples) -----------------
def verify_split(m, K, A, B):
    """(ok, excess, number of vectors of F_3^m missing from C(A, B))."""
    Ks, As, Bs = set(K), set(A), set(B)
    if not (As | Bs == Ks and As <= Ks and Bs <= Ks):
        return False, len(As & Bs), -1
    T = set(digits(k, m) for k in Ks)
    At = [digits(a, m) for a in sorted(As)]
    Bt = [digits(b, m) for b in sorted(Bs)]
    cov = set()
    for a in At:
        for b in Bt:
            if tuple((-x - y) % 3 for x, y in zip(a, b)) in T:
                cov.add(tuple((y - x) % 3 for x, y in zip(a, b)))
    missing = 3 ** m - len(cov)
    return missing == 0, len(As & Bs), missing


def double(m, K, A, B):
    """Lemma 1 and the split of Corollary: K+ in F_3^(m+1) and its split."""
    N = 3 ** m
    Kp = sorted(list(K) + [N + a for a in A] + [2 * N + b for b in B])
    Ap = set(A) | set(N + a for a in A)
    Bp = set(B) | set(2 * N + b for b in B)
    return Kp, Ap, Bp


def lift(m, Kprime, A, B):
    """Lemma 2: a split of the layered set Kprime of F_3^(m+1) from a split of its layer 0."""
    N = 3 ** m
    L0 = set(i for i in Kprime if i < N)
    if L0 != set(A) | set(B):
        raise ValueError("layer 0 of the set is not A u B")
    Ap = set(A) | set(i for i in Kprime if N <= i < 2 * N)
    Bp = set(B) | set(i for i in Kprime if i >= 2 * N)
    return Ap, Bp


# ---- lines of K --------------------------------------------------------------------------
def lines_of(m, K):
    """All full lines of K: (q0, q1, q2, dir) with q1 = q0 + D, q2 = q0 + 2D (positions in K),
    D the direction with first nonzero entry 1; and the number of directions met."""
    pts = [digits(k, m) for k in K]
    pos = {p: i for i, p in enumerate(pts)}
    dir_id = {}
    lines = []
    n = len(pts)
    for i in range(n):
        x = pts[i]
        for j in range(i + 1, n):
            y = pts[j]
            k = pos.get(tuple((-a - b) % 3 for a, b in zip(x, y)))
            if k is None or k <= j:
                continue
            d = tuple((b - a) % 3 for a, b in zip(x, y))
            lead = 0
            for c in d:
                if c:
                    lead = c
                    break
            if lead == 1:
                D, q = d, (i, j, k)
            else:
                D, q = tuple((2 * c) % 3 for c in d), (i, k, j)
            di = dir_id.get(D)
            if di is None:
                di = dir_id[D] = len(dir_id)
            lines.append((q[0], q[1], q[2], di))
    return lines, len(dir_id)


# colours: 0 = A only, 1 = B only, 2 = both.  TAB[9 c0 + 3 c1 + c2] = (number of ordered
# pairs (a, b) on the line with b - a = +D, the same for -D), a in A, b in B.
def _ina(c):
    return 0 if c == 1 else 1


def _inb(c):
    return 0 if c == 0 else 1


TAB = []
for _c0 in range(3):
    for _c1 in range(3):
        for _c2 in range(3):
            _q = (_c0, _c1, _c2)
            TAB.append((sum(_ina(_q[i]) * _inb(_q[(i + 1) % 3]) for i in range(3)),
                        sum(_ina(_q[(i + 1) % 3]) * _inb(_q[i]) for i in range(3))))


class SplitSearch:
    def __init__(self, n, lines, ndirs):
        self.n = n
        self.lines = lines
        self.ndirs = ndirs
        self.inc = [[] for _ in range(n)]       # (a, b, c, constraint, weight of q)
        self.by_dir = [[] for _ in range(ndirs)]
        for (a, b, c, d) in lines:
            self.by_dir[d].append((a, b, c))
            self.inc[a].append((a, b, c, 2 * d, 9))
            self.inc[b].append((a, b, c, 2 * d, 3))
            self.inc[c].append((a, b, c, 2 * d, 1))

    def reset(self, col):
        self.col = list(col)
        self.sat = [0] * (2 * self.ndirs)
        col = self.col
        for (a, b, c, d) in self.lines:
            p, q = TAB[9 * col[a] + 3 * col[b] + col[c]]
            self.sat[2 * d] += p
            self.sat[2 * d + 1] += q
        self.uns = []
        self.pos = [-1] * (2 * self.ndirs)
        for c in range(2 * self.ndirs):
            if self.sat[c] == 0:
                self.pos[c] = len(self.uns)
                self.uns.append(c)

    def _flag(self, c, bad):
        pos, uns = self.pos, self.uns
        if bad and pos[c] < 0:
            pos[c] = len(uns)
            uns.append(c)
        elif not bad and pos[c] >= 0:
            i = pos[c]
            last = uns.pop()
            if last != c:
                uns[i] = last
                pos[last] = i
            pos[c] = -1

    def delta(self, q, new):
        col, sat = self.col, self.sat
        old = col[q]
        brk = mk = 0
        for (a, b, c, cp, w) in self.inc[q]:
            i0 = 9 * col[a] + 3 * col[b] + col[c]
            p0, m0 = TAB[i0]
            p1, m1 = TAB[i0 + (new - old) * w]
            if p0 != p1:
                s = sat[cp]
                t = s - p0 + p1
                if s > 0 and t == 0:
                    brk += 1
                elif s == 0 and t > 0:
                    mk += 1
            if m0 != m1:
                s = sat[cp + 1]
                t = s - m0 + m1
                if s > 0 and t == 0:
                    brk += 1
                elif s == 0 and t > 0:
                    mk += 1
        return brk, mk

    def apply(self, q, new):
        col, sat = self.col, self.sat
        old = col[q]
        for (a, b, c, cp, w) in self.inc[q]:
            i0 = 9 * col[a] + 3 * col[b] + col[c]
            p0, m0 = TAB[i0]
            p1, m1 = TAB[i0 + (new - old) * w]
            if p0 != p1:
                sat[cp] += p1 - p0
                self._flag(cp, sat[cp] == 0)
            if m0 != m1:
                sat[cp + 1] += m1 - m0
                self._flag(cp + 1, sat[cp + 1] == 0)
        col[q] = new

    def candidates(self, c):
        out = []
        for (a, b, cc) in self.by_dir[c >> 1]:
            for q in (a, b, cc):
                if self.col[q] != 2:
                    out.append(q)
        return out

    def run(self, seed, max_steps, deadline=None, noise=0.15, stall=4000, tenure=3,
            reloc=0.02):
        """Search for a split with excess as small as possible (stops at excess 1).
        Returns (best colouring, its excess, steps)."""
        rng = random.Random(seed)
        n = self.n
        col = [rng.randrange(2) for _ in range(n)]
        col[rng.randrange(n)] = 2
        self.reset(col)
        while self.uns:                      # greedy: promote points to 'both'
            c = self.uns[rng.randrange(len(self.uns))]
            cand = sorted(set(self.candidates(c)))
            best, bl = -1, []
            for q in cand:
                mk = self.delta(q, 2)[1]
                if mk > best:
                    best, bl = mk, [q]
                elif mk == best:
                    bl.append(q)
            self.apply(rng.choice(bl), 2)
        best_col = list(self.col)
        best_e = sum(1 for x in best_col if x == 2)
        steps = 0
        last_flip = [-10 ** 9] * n
        while best_e > 1 and steps < max_steps:
            if deadline is not None and time.time() > deadline:
                break
            wild = [q for q in range(n) if self.col[q] == 2]
            if rng.random() < 0.5:
                w = rng.choice(wild)
                b0, b1 = self.delta(w, 0)[0], self.delta(w, 1)[0]
                nc = 0 if b0 < b1 or (b0 == b1 and rng.random() < 0.5) else 1
            else:
                opts, bb = [], None
                for w in wild:
                    for nc in (0, 1):
                        b = self.delta(w, nc)[0]
                        if bb is None or b < bb:
                            bb, opts = b, [(w, nc)]
                        elif b == bb:
                            opts.append((w, nc))
                w, nc = rng.choice(opts)
            self.apply(w, nc)
            e_now = best_e - 1
            last, bestu = steps, len(self.uns)
            while self.uns and steps - last < stall and steps < max_steps:
                steps += 1
                if deadline is not None and (steps & 1023) == 0 and time.time() > deadline:
                    break
                c = self.uns[rng.randrange(len(self.uns))]
                cand = self.candidates(c)
                if rng.random() < noise:
                    q = rng.choice(cand)
                else:
                    bq, key = [], None
                    for q in set(cand):
                        brk, mk = self.delta(q, 1 - self.col[q])
                        if brk > 0 and steps - last_flip[q] <= tenure:
                            continue
                        k = (brk, -mk)
                        if key is None or k < key:
                            key, bq = k, [q]
                        elif k == key:
                            bq.append(q)
                    q = rng.choice(bq) if bq else rng.choice(cand)
                self.apply(q, 1 - self.col[q])
                last_flip[q] = steps
                if reloc and rng.random() < reloc:
                    wild = [x for x in range(n) if self.col[x] == 2]
                    q2 = rng.choice(cand)
                    if wild and self.col[q2] != 2:
                        w = rng.choice(wild)
                        b0, b1 = self.delta(w, 0)[0], self.delta(w, 1)[0]
                        self.apply(w, 0 if b0 <= b1 else 1)
                        self.apply(q2, 2)
                if len(self.uns) < bestu:
                    bestu, last = len(self.uns), steps
            if not self.uns:
                best_col, best_e = list(self.col), e_now
            else:
                self.reset(best_col)
            if deadline is not None and time.time() > deadline:
                break
        return best_col, best_e, steps


def find_split(m, K, seed=1, max_steps=300000, deadline=None, **kw):
    """Search a split of the Kakeya set K (sorted indices in F_3^m).  Returns a dict with the
    states (over sorted K), the excess, the code and search statistics."""
    t0 = time.time()
    K = sorted(K)
    lines, ndirs = lines_of(m, K)
    if ndirs != (3 ** m - 1) // 2:
        raise ValueError("K is not a Kakeya set: %d of %d directions" % (ndirs, (3 ** m - 1) // 2))
    S = SplitSearch(len(K), lines, ndirs)
    t1 = time.time()
    col, e, steps = S.run(seed, max_steps, deadline=deadline, **kw)
    states = [c + 1 for c in col]
    single = sum(1 for d in S.by_dir if len(d) == 1)
    return {"m": m, "size": len(K), "excess": e, "states": states, "code": encode_split(states),
            "steps": steps, "lines": len(lines), "single": single,
            "setup_s": round(t1 - t0, 1), "search_s": round(time.time() - t1, 1)}
