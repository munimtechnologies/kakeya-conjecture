# kakeya-min-v1: raw output

Output of `solve.py` (run `kakeya-min-v1`, 2026-10-14, 4 cores, 3.3 hours): the exact
integer-programming search for minimum Kakeya sets in AG(*n*, *q*) described in section 7 of the
README.

- `results.json`: one entry per case (*n*, *q*). `runs` lists every solver run (OR-Tools CP-SAT
  with and without symmetry fixing, HiGHS), with its status, time, value, lower bound and the set
  it returned as encoded points (the point (*c*₀, ..., *c*ₙ₋₁) is Σ *c_i q^i*). `best_found` is the
  smallest verified set, `certificate` is that set as coordinate lists, and `proved_minimum` is set
  only when a solver proved optimality and every proof agreed.
- `results.txt`: one summary line per case.

Every certificate is re-checked line by line by `test_kakeya.py`. Solver lower bounds that were
not proved optimal are reported as the solvers gave them and are not used in the README.

By MV2 of Munim, Inc., made with AI agents working in MV2's name. Licence: CC BY 4.0.
