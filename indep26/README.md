# indep26: an independent check that no Kakeya set in F_3^4 has 26 points

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0 (Python standard
library only).

## Why this folder exists

The seven commit-time checks in [`gate26/`](../gate26) (search code [`kcut/kcut.py`](../kcut/kcut.py))
all landed on 2026-11-24: none of the seven hyperplane cuts listed in `moments/PRODUCT.md` holds a
26-point Kakeya set. If that is right, then k_4 = 27, which would confirm the value Lev's 2009
computer search suggested and the Polymath wiki called likely. Before anyone relies on it, the same
conclusion has to be reached a second, independent way. This folder does that:

* **New code.** `kk.py` was written from scratch without reusing kcut/ or gate26/ code. Only the
  coordinates of Lev's 27-point set were copied, as test data, and they are re-verified here.
* **Fewer lemmas.** The case list uses only two exact counting identities and the pair-product
  bound (T). It does **not** use the sum-difference lemma of `moments/SD.md`, so it has to search
  ten cuts instead of seven: (4,7,15), (5,6,15) and (4,8,14) are searched here, not ruled out.
* **Different search.** Orbit representatives come from a closure of AGL(3,3) generators (union by
  breadth-first search), the middle section is enumerated by plain lexicographic search with a
  covering bound, and the third section by a branch and bound with exclusion and packing bounds.

## The argument

Write K for a Kakeya set in F_3^4: it contains a full line {x, x+d, x+2d} for each of the 40
directions d. Adding a point to a Kakeya set keeps it Kakeya, and a 27-point Kakeya set exists
(Lev, 2009). So **k_4 = 27 if and only if no Kakeya set has exactly 26 points.**

Each of the 40 classes of parallel hyperplanes {xi.x = t} (t = 0, 1, 2) cuts K into three sections.
After an affine change of coordinates, the class is t = the last coordinate, and the sections A, B, C
(at t = 0, 1, 2) are subsets of F_3^3. K is Kakeya exactly when

* **(X)** for every v in F_3^3 there is x in A with x + v in B and x + 2v in C (the 27 directions
  (v, 1)), and
* **(D)** each of the 13 directions of F_3^3 has a line inside A, B or C.

**(T)** From (X), B - A, C - A and C - B are each all of F_3^3, so every two sections have sizes
whose product is at least 27.

**Identities** (`moments/PRODUCT.md`, proved there in five lines; re-checked numerically in
`test_controls.py` on random sets and on Lev's set). For N points containing L lines, summed over
the 40 classes of F_3^4:

    sum n1 n2 n3 = 9 (C(N,3) + 2 L)        sum (n1 n2 + n1 n3 + n2 n3) = 27 C(N,2).

For N = 26 these give sum (n1 - 8)(n2 - 8)(n3 - 8) = 18 (L - 40), and a Kakeya set has L >= 40
(one line per direction). So:

1. **Some class has (n1-8)(n2-8)(n3-8) > 0.** All three factors cannot be positive (that needs
   27 points), so two sections have at most 7 points and one at least 9. Sorted, with (T):
   (4,7,15), (5,6,15), (5,7,14), (6,6,14), (6,7,13), (7,7,12).
2. **Or every product is 0** (and L = 40), so every class has an 8-point section. Summed over the
   classes, n1^2 + n2^2 + n3^2 = 40 * 676 - 2 * 8775 = 9490, while a class whose sections are all
   at least 8 is (8,8,10) or (8,9,9), with at most 228; 40 * 228 = 9120 < 9490. So some class has
   an 8-point section and a section of at most 7 points: (a, 8, 18 - a) with 4 <= a <= 7 by (T):
   (4,8,14), (5,8,13), (6,8,12), (7,8,11).

`cuts26_from_identity()` derives exactly these ten cuts and `test_cut_list` checks them.

**Normalising one class.** Translating t by 1 and replacing t by -t permute the three sections in
every possible way, so the smallest section may be put at t = 0 and the next smallest at t = 1. The
maps (x, t) -> (g x, t) with g in AGL(3,3) replace A by any set in its orbit, so A may be taken to be
an orbit representative. The shear (x, t) -> (x + t u, t) fixes A and translates B, so B may be
taken to contain the origin. All of these maps are affine bijections of F_3^4 that preserve the
class, so they map Kakeya sets to Kakeya sets of the same size.

**The search, for a cut (a, b, c).** For every orbit representative A of a-subsets of F_3^3, every
b-set B through the origin with B - A = F_3^3, decide exactly whether some C with |C| <= c satisfies
(X) and (D). Any solution is assembled into a set of F_3^4 and re-checked by brute force
(`is_kakeya`). A test file `test_cut_a_b_c*.py` passes only if no solution exists for its cut (or
for its part of the orbit representatives, index = part mod parts). The publishing pipeline commits
a test file only after it passes, inside a package folder, so it really runs at commit time (see
`gate26/README.md` for why that matters).

**Conclusion, if all test files below land:** no 26-point Kakeya set exists in F_3^4, so k_4 = 27.

## Controls (run before every search, and in `test_controls.py`)

* The generated linear groups have orders 48 = |GL(2,3)| and 11232 = |GL(3,3)|; with translations
  432 = |AGL(2,3)|. Orbit sizes of k-subsets add up to C(27, k) and divide 303264 = |AGL(3,3)|.
* The middle-section enumeration equals brute force in F_3^2 (all cases) and in F_3^3 (one case,
  230230 candidates).
* The third-section search returns exactly the brute-force minimum for every case in F_3^2 with
  |A| in {2, 3, 4} and |B| <= 5.
* **Lev's k_3 = 13, one dimension down, by the same pipeline:** no 12-point Kakeya set in F_3^3
  under any cut allowed by (T), and a 13-point one is found.
* **Lev's 27-point set in F_3^4:** it is Kakeya; for each of the 4 coordinate classes and each of the
  6 orders of its sections, the search completes the first two sections within the size of the
  third; and the full pipeline (orbit representative, B through the origin, C search) finds a
  Kakeya set of at most 27 points in the cut (6, 8, 13) from the orbit of its 6-point section.

## Files

| file | what it checks |
|---|---|
| `kk.py` | the search (standard library only) |
| `test_controls.py` | everything in "Controls" above, plus the identities and the cut list |
| `test_cut_4_7_15.py`, `test_cut_5_6_15.py`, `test_cut_4_8_14.py`, `test_cut_5_7_14.py`, `test_cut_6_6_14.py`, `test_cut_5_8_13.py`, `test_cut_6_7_13.py` | one cut each |
| `test_cut_7_7_12_p{0,1}of2.py`, `test_cut_6_8_12_p{0,1}of2.py`, `test_cut_7_8_11_p{0..3}of4.py` | one part of a cut each |

A missing test file means that cut (or part) is still undecided by this check, for example because
it did not finish within the pipeline's ten-minute limit.

`test_controls.py` lands with the first commit; each cut file lands in its own commit, which also
runs `kk.quick_controls()` (the k_3 = 13 and 27-point controls) before its search.

Run locally from the repository root: `python -m unittest discover -s indep26 -t . -p "test_*.py"`
runs this folder; `python -m unittest indep26.test_cut_6_6_14` runs one cut.

## Limits

This is a computer proof. It relies on the correctness of `kk.py`, which the controls test but
cannot prove, and on the hand proofs of the identities above. An agreement between `kcut` (gate26)
and `kk` (indep26), written independently, makes a shared bug unlikely, not impossible. A formal
certificate (for example DRAT proofs from a SAT encoding, queued as `kdrat26`) and a human reading
are still needed before k_4 = 27 should be called proved. Credit: Lev (2009) for k_3 = 13 and the
27-point set; the Polymath project wiki for recording them.
