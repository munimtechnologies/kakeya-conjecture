# Corrections

MV2 of Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name. Text CC BY 4.0.

## 2026-11-24: two priority claims in the October 2026 survey (README, version 1.1) were wrong

1. **K(3, 3) = 13 is not new.** Vsevolod Lev found by computer in 2009 that the smallest Kakeya set
   in F_3^3 has 13 points, and he also found a 27-point Kakeya set in F_3^4. Both are recorded on
   the Polymath project's "Kakeya problem" wiki page (last modified 5 June 2009; archived copy:
   https://web.archive.org/web/20201223142807/https://asone.ai/polymath/index.php?title=Kakeya_problem).
   Our exhaustive computation is an independent confirmation of k_3 = 13, and our 27-point set in
   F_3^4 matches the size of Lev's. The survey's words "new", "first exact finite-field Kakeya
   number in dimension three" and "27 (new)" should be read with this correction.
2. **56 points in F_3^5 was not smaller than the literature.** Chung, Du and Wesley
   (https://arxiv.org/abs/2608.23691) had already given a 53-point Kakeya set in F_3^5.

We have not found earlier constructions with 28 points in F_4^3 or 64 in F_4^4, but we have not
searched for those as carefully. Results after the survey are dated and list their prior work:
`moments/SD.md` (k_4 >= 25), `moments/BAD25.md`, `moments/PRODUCT.md`, and the searches in `kcut/`,
`kprof25/`, `kprof26/` and `kdrat26/`.
