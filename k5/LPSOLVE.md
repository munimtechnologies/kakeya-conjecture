# The pencil LP without scipy: a standard-library solver with exact certificates

MV2 of Munim, Inc. (mv2@munimtech.com), 2026-11-25. Made with AI agents working in MV2's name.
**No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

`PENCIL.md` sets up the pencil LP for Kakeya sets in F_3^n. Its HiGHS script `pencil_lp.py`
needs scipy. This note adds a route that needs only Python's standard library, so that the
LP can run inside commit-time tests:

- `lpsolve.py`: a phase-1 revised simplex in floating point (Dantzig pricing with a fallback
  to Bland's rule, periodic refactorisation). It never proves anything by itself. When it
  ends with a positive phase-1 objective, the dual vector of the final basis is recomputed
  **exactly** (B^T y = c_B solved with Fractions) and turned into integer Farkas multipliers.
  The docstring proves why such a dual is a Farkas certificate.
- `pencil_gate.py`: builds each branch with `pencil.build_lp`, solves it, and counts a
  branch as excluded only if `pencil.check_certificate` (exact, rational, with the column
  bounds of `build_lp`) accepts the multipliers. A wrong floating-point decision can only
  lose a certificate, never create one.

## Tests

- `test_lpsolve.py`: toy systems; 150 random systems with a planted solution (must be
  feasible); 150 random systems made infeasible on purpose (must be certified); 200 random
  systems where every "infeasible" answer must carry a valid certificate. Controls: the
  pencil LP must stay feasible at 13 points in F_3^3 and on the (z, u) branch of Lev's
  27-point set in F_3^4, where Kakeya sets exist.
- `test_pencil_n3.py`: calibration; asserts that the pencil LP excludes 12 points in
  F_3^3 (Lev's k_3 = 13, without a search). It is committed only if it passes.
- `test_pencil26_z0.py` ... `test_pencil26_z8.py`: one file per quadric type of F_3^4
  (sorted order of `pencil.quadric_types(4)`). Each asserts that every branch of its type is
  excluded at 26 points with an exactly checked certificate, and is committed only if it
  passes. **If all nine files are in the repository, the pencil LP alone proves that no
  Kakeya set in F_3^4 has 26 points (k_4 = 27), independently of the exhaustive searches of
  `PROOF26.md`.** A missing file means its type was not excluded (or its test ran out of
  time); nothing is claimed for it.

The class-level LP (`quad.py`) cannot do this: at 26 points in F_3^4 it excludes none of the
nine quadric types (`RESIDUE.md` gives explicit class-level solutions). The pencil level adds
how the four classes through a codimension-2 subspace fit together.
