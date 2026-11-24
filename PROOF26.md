# k_4 = 27: no Kakeya set in F_3^4 has 26 points (computer-assisted)

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-24. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

This file supersedes the bounds table in the October survey in `README.md` (which lists the F_3^4
case as open, 18 to 27) and should be read together with `CORRECTIONS.md`.

## Statement

A Kakeya set in F_3^4 is a set of points containing a full line {x, x+d, x+2d} in each of the 40
directions d. Write k_4 for the smallest size of such a set. Lev found a 27-point Kakeya set by
computer in 2009, and the Polymath wiki's Kakeya page (last edited June 2009) recorded "k_4 <= 27"
and that "it seems likely that k_4 = 27 holds". Before November 2026 the best proved lower bound
we could find was 18.

**Result.** Two independently written exhaustive searches find no 26-point Kakeya set in F_3^4.
Since adding a point to a Kakeya set keeps it Kakeya, no Kakeya set has 26 or fewer points, and
k_4 = 27.

## The argument in one page

1. **Identities** (`moments/PRODUCT.md`, gate-tested by `moments/test_product.py`; re-checked in
   `indep26/test_controls.py`). For N points with L lines, over the 40 hyperplane classes:
   sum n1 n2 n3 = 9 (C(N,3) + 2L) and sum (n1 n2 + n1 n3 + n2 n3) = 27 C(N,2). For N = 26 this
   gives sum (n1-8)(n2-8)(n3-8) = 18 (L - 40) >= 0, since a Kakeya set has L >= 40.
2. **Cases.** Either some class has (n1-8)(n2-8)(n3-8) > 0, so two sections have at most 7 points;
   or every product is 0, and a sum-of-squares count (9490 > 40 * 228) gives a class (a, 8, 18-a)
   with a <= 7. With the pair-product bound (T) (any two sections have product >= 27) this leaves
   ten cuts: (4,7,15), (5,6,15), (5,7,14), (6,6,14), (6,7,13), (7,7,12), (4,8,14), (5,8,13),
   (6,8,12), (7,8,11).
3. **Normal form.** For a class with sections A, B, C (smallest first), affine maps of F_3^4
   preserving the class make A an AGL(3,3)-orbit representative and put the origin in B.
4. **Search.** For every such (A, B), decide exactly whether a third section C of the allowed size
   gives a Kakeya set (27 transversal conditions plus 13 in-plane directions).

## The two checks

| check | code | cuts | lemmas used | where |
|---|---|---|---|---|
| gate26 | `kcut/kcut.py` | the 7 cuts left after the sum-difference lemma (`moments/SD.md`) | identities, (T), SD | `gate26/` (7 commits) |
| indep26 | `indep26/kk.py`, written without kcut code | all 10 cuts | identities, (T) only | `indep26/` (11 search commits) |

So 7 cuts are excluded by two independent searches, and the other 3 by one search plus the
sum-difference lemma. Every check is a unit test inside a Python package, run by the publishing
pipeline before the commit; a test file exists in the repository only if its search finished
within ten minutes and found nothing. Before every search the same code must reproduce Lev's
k_3 = 13 one dimension down and recover Lev's 27-point set. The per-cut record, with the commit
of every check, is `indep26/k4_result.json` (checked by `indep26/test_result_file.py`).

Pipeline check: `indep26/test_gate_must_reject.py` is a deliberately failing test, filed through
the same pipeline. It must never appear in this repository; if it does, the commit-time checks did
not run and none of the files above may be read as passed checks.

## Limits

This is a computer proof: it depends on the correctness of `kcut/kcut.py` and `indep26/kk.py`,
which their controls test but cannot prove, and on the hand proofs of the identities and of the
sum-difference lemma. Two independent programs agreeing makes a shared bug unlikely, not
impossible. A certificate-level re-check (SAT encodings with DRAT proofs, `kdrat26/`) is queued.
A human reading is still needed.

Credit: V. Lev (2009) for k_3 = 13 and the 27-point set; the Polymath project wiki for recording
them; Dvir (2008) and Dvir-Kopparty-Saraf-Sudan for the general finite-field Kakeya bounds.
