# τ(3) = 448: the least rainbow triple of F_3^3

MV2, Munim, Inc. (mv2@munimtech.com), 2026-11-27. Made with AI agents working in MV2's
name. **No human mathematician has checked this yet.** Text CC BY 4.0, code Apache-2.0.

## Theorem

**Every rainbow triple (A, M, C) of F_3^3 has |A|·|M|·|C| ≥ 448, and 448 is attained.**

The bound is attained by the (7, 8, 8) triple found by CP-SAT on 2026-10-15
(`compete-results/kakeya-s3-v2`), listed in `rainbow/test_tau3_448.py` and checked there
from the definition:

    A = {000, 001, 021, 100, 121, 201, 222}
    M = {000, 012, 020, 102, 112, 200, 202, 210}
    C = {000, 002, 010, 012, 020, 102, 200, 212}

So τ(3) = 448, the dimension-3 base is 448^(1/9) = 1.9705…, and there is no rainbow triple
of F_3^3 whose three sets all have at most 7 points (7³ = 343). Through the tripling lemma of
`TRIPLES.md` this gives Kakeya sets of 3·448^k points in F_3^(9k); among triples of
dimension at most 3 none does better.

## Proof

As in `RAINBOW3.md`, with the bound 448: a triple in normal form with product P ≤ 447 has
|A| ≤ 7 (8³ = 512); lemma H (9·τ(2) = 576 > 447) makes A span F_3^3, so |A| ∈ {4, …, 7};
the size bounds leave |M| ∈ {7, …, 10}, {6, …, 9}, {6, 7, 8} and {7} for |A| = 4, 5, 6, 7.
`rainbow/test_tau3_448.py` scans all of them, with the class list of 7-sets checked by
orbit–stabiliser counts, and finds no triple with product below 448. ∎

## Limits

An exhaustive computer search on top of short hand-proved reductions (`SCAN.md`); not
checked by a human. It says nothing about τ(4) beyond lemma H: with τ(3) = 448, a triple of
F_3^4 with a set inside a hyperplane has product at least 9·448 = 4032, so every triple
below 4032 (in particular any with product ≤ 2075, which would beat the
Kopparty–Lev–Saraf–Sudan base) has all three sets spanning F_3^4. The best known triple of
F_3^4 is the product (7, 8, 8) × ({0,1}, {0,1}, {0,2}), with product 3584.
