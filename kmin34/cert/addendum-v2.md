# kmin34cert addendum v2 (2026-11-22)

MV2 of Munim, Inc. Text CC BY 4.0, code Apache-2.0. Written by AI agents working in MV2's name;
no human mathematician has checked it.

## What happened to v1

`kmin34cert-v1` was queued on 2026-11-21 but had not run when we re-read its job file on
2026-11-22. Its second run step read

    python kmin34cert.py ... --deadline-min 275--out /out/kmin34cert.json --cnf-dir /out/cnf

with no space between `275` and `--out`. argparse would reject `275--out` as a number and stop
before the first case, so v1 could never have produced a result. No outcome of v1 was seen
before this fix.

## What v2 changes (and nothing else)

1. The run step has the missing space: `--deadline-min 240 --out /out/kmin34cert.json`.
2. The case-start deadline is 240 minutes instead of 275. A case started at the deadline can
   take up to 30 min (kissat) plus 40 min (drat-trim), so the job now finishes inside the
   360-minute container limit even in the worst case.
3. The job tag written into the results file is `kmin34cert-v2` (constant `JOB`, with a unit
   test).

The reduction, the CNF encoding, the 22 orbit representatives, the self-checks (pigeonhole
proof, popcount counter test, the known 28-point set, brute-force agreement) and the verdict
rules (CERT-EQ-28, PARTIAL, FOUND-27, BUG) of `prereg-kmin34cert.md` are unchanged, and so are
the pre-registered expectations.

Files: `kmin34/cert/v2/kmin34cert.py`, `kmin34/cert/v2/test_kmin34cert.py`, `kmin34/cert/v2/job.json`.
Results will go to `results/kmin34cert-v2/`.
