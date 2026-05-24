---
fontsize: 10pt
output:
  pdf_document:
    keep_tex: true
    latex_engine: lualatex
    extra_dependencies: ["mathtools", "amsmath", "amssymb"]
monofont: "Menlo"
---

# WEeDS Background — Reorganization & Verification Report

**Subject document:** `WEeDS_Background_reorganized.md` (reorganized from `052326_WE_BidirectionalSeedingBackground.md`)
**Scope of this report:** (1) what was reorganized and where content moved; (2) independent verification of every worked calculation and the principal derivations; (3) assessment of statistical validity; (4) a complete catalog of internal-consistency and cross-reference issues, with an old $\to$ new section map. The first reorganization round altered only structure, ordering, headings, the Table of Contents, and an orientation note. A subsequent round — documented in §6 — additionally expanded the weighted-ensemble methodology section, added inline citations and a reference list, standardized load-bearing notation, split long paragraphs, converted one extensive inline expression to display math, and removed in-text section-number cross-references at the author's request. Substantive errors found during verification are reported here rather than silently patched, per the agreed workflow.

---

## 1. What was reorganized

The source file was two concatenated documents that collided on section numbering: a "background methods" document (§§0–4, with §2 placed *after* §3, an out-of-order §4, an orphaned "§5.5", and a corrupted Table of Contents wedged mid-text) and a self-contained brief, *"Bidirectional Seeding, Mode Excitation, and Pooling,"* that **restarted its own numbering at §1–§10**.

The reorganization unifies these into one continuous document with the following moves.

**Global ordering.** Restored the logical order Abstract $\to$ Program Statement $\to$ **Theoretical Foundations (§2)** $\to$ **Pipeline Architecture (§3)** $\to$ **WE Sampling (§4)**. In the source, §2 followed §3. The pipeline architecture figure was relocated to the head of §3 (it is the "Pipeline architecture" figure). The Table of Contents was rebuilt from scratch to match the actual contents (the original ToC referenced §§5–30 and mislabeled §7 as "Pipeline Integration" with §20.x subsection anchors).

**The brief folded into Chapter 4.** Because the downstream pipeline owns chapter slots §5–§30 (tICA, RiteWeight, density-ratio, DMD, kEDMD, PCCA+, GP-FES, VAMPnet, validation — none of which are present in this file; they live in the companion Pipeline document), the seeding/pooling brief was folded into Chapter 4, its true thematic home (WE experimental design). This is collision-free and matches the brief's own integration intent.

**Duplicate derivations merged to a canonical home.** The brief repeated several derivations already in §2. These were consolidated:

| Brief content | Canonical home in reorganized doc | Disposition |
|---|---|---|
| §1.1 data structure restated | §3.1 / §4.3 | dropped (duplicate) |
| §1.2 snapshot view | §2.6 | dropped (duplicate) |
| §1.3 the -1 shift / orthogonality | §2.6a | dropped (duplicate; see issue I5) |
| §1.4 density-relaxation expansion | §2.6 | dropped (duplicate) |
| §1.5 trajectory-pair view | §2.5 | dropped (duplicate) |
| §1.6 VAMP-2 | §2.5 | dropped (duplicate) |
| §1.7 estimator equivalence/failure modes | **§2.5.1** (relocated, retitled) | kept |
| §1.8 noise structure of density ratio | **§2.6b** (relocated, new home) | kept |
| §1.9 amplitude/temporal/spatial decomposition | **§2.6.1** (relocated) | kept |

**Editorial sections dissolved.** Brief §9 ("Suggested integration into the methods document") was consumed — its guidance is exactly the reorganization performed here — and removed. Brief §10 ("Loose ends and open questions") was removed *as a section*, but its five substantive open questions were preserved as clearly-marked `Caveat` notes folded into the relevant technical subsections (heterogeneous-coverage and bootstrap-interaction caveats $\to$ §4.14; bandwidth–pooling caveat $\to$ §4.12; higher-mode-cancellation and asymmetric-equilibrium caveats $\to$ §4.11). No information was lost; only the "loose ends" framing was removed.

**Subsection reordering inside §5 of the brief.** The brief's identifiability subsections appeared in file order 5.1, 5.2, **5.4, 5.3** (mis-numbered, though the *Condition 1–7* content flowed correctly). They are renumbered in reading order as §4.12.1–§4.12.4.

### 1.1 Old $\to$ new section map

Use this table to update any in-text cross-references (see issue I6).

| Source | Reorganized |
|---|---|
| §2 (was after §3) | §2 (now before §3) |
| — | §2.5.1 (was brief §1.7) |
| — | §2.6.1 (was brief §1.9) |
| — | §2.6b (was brief §1.8) |
| §4.7 Equilibrium WE / WEED | §4.5 |
| §4.5 PC choice & bidirectional seeding | §4.6 |
| §4.9 PC independence & limits | §4.7 |
| §5.5 What WESTPA produces | §4.8 |
| brief §2 Mode excitation | §4.9 (subs §2.1–2.6 $\to$ §4.9.1–4.9.6) |
| §4.6 Coverage vs crossings | §4.10 |
| brief §3.3 asymmetric advance | §4.10.1 |
| brief §4 Coverage–excitation tradeoff | §4.11 (subs §4.1–4.6 $\to$ §4.11.1–4.11.6) |
| brief §5 Identifiability conditions | §4.12 (§5.1 $\to$ 4.12.1, §5.2 $\to$ 4.12.2, §5.4 $\to$ 4.12.3, §5.3 $\to$ 4.12.4) |
| §4.5b Pooled-corpus architecture | §4.13 |
| brief §6.1 / §6.2 | §4.13.1 / §4.13.2 (brief §6.3 dropped, duplicate) |
| brief §7 Destructive interference | §4.14 (§7.1–7.6 $\to$ §4.14.1–4.14.6; §7.3.1–7.3.6 $\to$ §4.14.3.1–4.14.3.6) |
| brief §8 Operational recommendations | §4.15 (§8.1–8.3 $\to$ §4.15.1–4.15.3) |
| brief §9 Integration notes | removed (executed) |
| brief §10 Loose ends | removed; caveats folded into §4.11, §4.12, §4.14 |

Line count went from 2,398 to 1,987; the difference is the merged duplicate theory (brief §1.1–1.6), the dropped redundant brief §3.1/§3.2/§6.3, and the removed editorial framing. **No unique technical content was discarded.**

---

## 2. Mathematical verification

Every worked numerical example and the principal derivations were re-derived independently and checked numerically (symbolic/`numpy`). **All arithmetic is correct.**

### 2.1 Worked calculations (all verified exact)

**§4.11.1 — asymmetric eigenfunction normalization** ($\pi_F=0.9,\ \pi_U=0.1$).
From $\mathbb{E}_\pi[\psi_2]=\pi_F\psi_2^F+\pi_U\psi_2^U=0$, $\psi_2^U=-(\pi_F/\pi_U)\psi_2^F=-9\psi_2^F$ $\checkmark$. Normalization $\pi_F(\psi_2^F)^2+\pi_U(\psi_2^U)^2=0.9a^2+0.1\cdot81a^2=9a^2=1$ gives $|\psi_2^F|=1/3,\ |\psi_2^U|=3$ $\checkmark$.

**§4.11.2 — unidirectional excitation.** $c_2^{\text{uni-F}}=-1/3$, $c_2^{\text{uni-U}}=+3$, ratio $=9$ $\checkmark$.

**§4.11.3 — bidirectional identity.** $c_2^{\text{bidir}}(\alpha)=\alpha\psi_2^F+(1-\alpha)\psi_2^U=(\alpha-\pi_F)(\psi_2^F-\psi_2^U)$ verified algebraically (uses $1-\alpha-\pi_U=-(\alpha-\pi_F)$) $\checkmark$. $c_2^{\text{bidir}}(0.5)=(-0.4)(-10/3)=4/3\approx1.33$ $\checkmark$; vanishes iff $\alpha=\pi_F$ $\checkmark$.

**§4.11.4 — tradeoff ratios.** $|c_2^{\text{uni-U}}|/|c_2^{\text{bidir}}|=3/(4/3)=2.25$ $\checkmark$; $|c_2^{\text{bidir}}|/|c_2^{\text{uni-F}}|=(4/3)/(1/3)=4$ $\checkmark$.

**§4.14.2 — destructive case** ($\pi_F=\pi_U=0.5$, $\psi_2^F=-1,\psi_2^U=+1$). Norm $0.5+0.5=1$ $\checkmark$; $c_2^{(A)}=-0.6$, $c_2^{(B)}=+0.6$, equal-weight pool $=0$ $\checkmark$. Partial case $(-0.6,+0.4)\to-0.1$, a sixfold reduction $\checkmark$.

**§4.14.3.4–5 — break-even & worked example.** Three-run $(-0.6,+0.6,-0.6)/3=-0.2$, $\eta_2=0.2/0.6=1/3<1/\sqrt3\approx0.577$ $\to$ net-negative $\checkmark$. Partial $\eta_2=0.1/0.6=1/6<1/\sqrt2\approx0.707$ $\to$ net-negative $\checkmark$. Thresholds $1/\sqrt N$: $N=4\to0.5$, $9\to1/3$, $16\to1/4$ $\checkmark$.

### 2.2 Derivations checked (sound)

- **Similarity transform** $(\mathcal P_\tau\rho)/\pi=\mathcal K_\tau(\rho/\pi)$ (§2.3): correct; follows directly from detailed balance $p_\tau(y\mid x)\pi(x)=p_\tau(x\mid y)\pi(y)$.
- **Density-ratio recursion** $r(x,t+\tau)=(\mathcal K_\tau r(\cdot,t))(x)$ and expansion $r=\sum_{i\ge2}c_i\lambda_i^{t/\tau}\psi_i$ (§2.6): correct; the $-1$ shift removes $\psi_1\equiv1$ and yields $\langle r,1\rangle_\pi=0$.
- **Implied timescales** $\lambda_i(\tau)=e^{-\tau/\tau_i}$ from the semigroup/Chapman–Kolmogorov functional equation (§2.3): correct.
- **KL as Lyapunov function** (§2.6a): $D_{\mathrm{KL}}(p\Vert\pi)\approx\tfrac12\|r\|_{L^2(\pi)}^2=\tfrac12\sum_{i\ge2}c_i^2\lambda_i^{2t/\tau}$, and $\ln|\lambda_2|=\lim_k\frac1{2k}\ln D_{\mathrm{KL}}(t_k)$: second-order expansion verified numerically (ratio $D_{\mathrm{KL}}/\tfrac12\|r\|^2\to0.9997$) $\checkmark$. The asymptotic-slope identity is correct.
- **KDE noise / $\sqrt\pi$-weighting** (§2.6b): $\mathrm{Var}[\eta]\sim f/(n_{\mathrm{eff}}h^d)$ is standard kernel-density form; $\mathrm{SD}[\hat{\mathbf R}]\approx\sigma_{\mathrm{noise}}/\hat\pi$ and the $\sqrt{\hat\pi}$-weighting reducing it to $\sigma_{\mathrm{noise}}/\sqrt{\hat\pi}$ are algebraically correct, and the geometric reading as the change of measure to $L^2(\pi)$ is sound.
- **Pooling noise** $\sigma^{\mathrm{pool}}\sim\sigma^{(\text{single})}/\sqrt{N_{\mathrm{runs}}}$ and $\mathrm{SNR}\propto|c_i^{\mathrm{pool}}|\sqrt{N_{\mathrm{runs}}}$ (§4.14.3): correct for independent runs of comparable size. The break-even $\sqrt{N_{\mathrm{runs}}}\,\eta_i>1$ follows directly.
- **Prony bound** $N\ge2(n-1)$ to identify $n-1$ exponential modes (§4.12.1): correct in spirit — $2M$ samples to fit $M$ damped exponentials; the rank-$(n-1)$ factorization argument is the right justification.
- **Resampling unbiasedness** (§4.2): the splitting identity (weight conserved, $m$ co-located deltas of weight $w_i/m$) and the merging expectation identity are both correct; merging adds variance, not bias.

---

## 3. Statistical-validity assessment

The statistical architecture is, on the whole, well-posed and the claims are defensible. Observations:

1. **Single-loop joint bootstrap (§1, abstract / companion §21.3).** Replacing a nested $B_{\text{outer}}\times B_{\text{inner}}$ scheme with one replicate loop that (i) draws a Dirichlet weight on the run axis, (ii) block-resamples iterations within runs, (iii) rebuilds and re-runs the pipeline, is statistically coherent. The Dirichlet (Bayesian/Rubin) bootstrap on the run axis is the appropriate device when $N_{\text{runs}}$ is small, and the citation to Mostofian & Zuckerman for log-space CIs on high-log-variance observables is apt. **Caveat (already folded into §4.14):** treating runs as exchangeable resampling units is in tension with the destructive-interference structure, where some run pairs are "compatible" and others cancel; this is correctly flagged as second-order but worth checking empirically.

2. **Block length.** Choosing $b_{\text{block}}$ to exceed the empirical iteration-correlation time (Künsch / Liu–Singh) is the correct way to preserve autocorrelation under resampling. The document should, at validation time, *report* the measured correlation time and the chosen block length so the choice is auditable (noted as a recommendation, not an error).

3. **Kish ESS (companion §21.3a).** Retaining per-state Kish $n_{\text{eff}}=(\sum w)^2/\sum w^2$ as a baseline diagnostic is appropriate; the document correctly treats it as a *baseline* that the block correction then refines (Kish assumes independence).

4. **"Cannot reweight what was not sampled" (§4.7, §4.10).** The coverage-vs-equilibrium-resemblance distinction (§4.14.6) and the support-coverage requirement are stated correctly and are the right framing; the claim that spectral objects are reconstructible from *local* operator action under *coverage* (not reactive crossings) is sound under the stated reversibility/ergodicity-on-accessible-support assumptions.

5. **Reversibility as a checked assumption.** The document is appropriately careful: it assumes detailed balance, derives the real-spectrum/self-adjoint structure from it, and routes finite-sample violations to a reported detailed-balance residual rather than presuming exactness. This is good practice.

No statistical claim was found to be invalid. The principal statistical *risk* the document itself identifies — destructive pooling interference biasing detectability (not the estimates) and complicating the run-axis bootstrap — is correctly characterized and now carries explicit caveats.

---

## 4. Internal-consistency and cross-reference issues

> **Revision status (this pass).** All items in this section have since been addressed in `WEeDS_Background_reorganized.md`: the abstract summation index was corrected to $\sum_{i\ge2}$ (I1); the two-vs-three weight-type count was reconciled (I2); the dangling $H$-theorem reference was repaired (I4); the wording/typos were fixed (I7); and **every in-text cross-reference carrying a specific section number was removed from the prose** (figure `\S` labels retained), which also moots the RevVAMPnet §15/§19 mismatch (I3) and supersedes the remap-table approach of (I6). The brief's introductory roadmap, inadvertently dropped in the first merge, was restored as a lead-in to the Chapter 4 seeding/pooling analysis. The descriptions below are retained as a record of what was found.



These require author attention. None is a computational error; they are notation/labeling/reference inconsistencies, several of which predate the reorganization.

**I1 — Abstract summation index (likely typo).** The abstract writes the spectral expansion as
$p(x,t_k)/\pi(x)-1=\sum_{i\ge 1}c_i\lambda_i^{t_k/\tau_{\mathrm{WE}}}\psi_i(x)$ with the sum starting at $i\ge 1$. Everywhere else in the document the sum starts at $i\ge 2$, because the $-1$ shift removes the stationary mode $\psi_1\equiv1$ (with $\lambda_1=1$). Recommend changing the abstract to $\sum_{i\ge 2}$ for consistency with §2.6, §2.6a, §4.9, and §4.14.

**I2 — "Two" vs "three" walker weight types.** §1 (line ~56) states *"Three walker weight types are kept strictly distinct"* (raw $w_i^{(k)}$, RiteWeight $\hat w_i^\pi$, and per-iteration normalized $\tilde w_i^{(k)}$), while §3.2 opens *"Two distinct walker-weight types appear…"* (treating $\tilde w_i^{(k)}$ as a derived normalization of Type 1). Both are internally defensible, but the count should be reconciled. Suggest: in §3.2, state explicitly that there are two *primary* types plus one *derived* per-iteration normalization, matching §1's count of three. (The section heading was made neutral — "Walker weight types and their roles" — during reorganization.)

**I3 — RevVAMPnet section number: §15 vs §19.** §4.8 (formerly §5.5) refers to *"the RevVAMPnet tertiary route (§15)"*, whereas the figure caption, §2.5, §3.1, §3.2, and the program statement all use **§19** (e.g., §19.2, §19.4). Recommend standardizing to §19 (companion document).

**I4 — Dangling references to §2.7 / §2.8.** §2.6a (line ~432) says the $H$-theorem is *"strengthened to the quantitative de Bruijn form in §2.7 below,"* but there is no §2.7 in this document (the original ToC listed "2.7 Wasserstein gradient-flow structure" and "2.8 Identifiability conditions," neither of which was ever written). Either author §2.7 or soften the reference. This predates the reorganization.

**I5 — Density-ratio symbol $r$ (resolved by the merge, noted for awareness).** The source was internally inconsistent: §2.3/§2.6/§2.6a define $r(x,t)=p/\pi-1$ (shift included), whereas brief §1.3 defined $r=p/\pi$ and subtracted unity separately ($\langle r-1,1\rangle_\pi=0$). Because §1.3 was a duplicate of §2.6a and was dropped, the canonical convention $r\equiv p/\pi-1$ now holds throughout. No action needed unless §1.3's phrasing is reintroduced elsewhere.

**I6 — In-text cross-references to renumbered sections (action required).** Section *headings* were fully renumbered, but in-text cross-references were **not** mechanically rewritten, because the brief's reference tokens are overloaded against the companion pipeline's numbering and cannot be disambiguated automatically without risking corruption. Specifically:

- `§1.8` (x4) should now read **§2.6b** (unambiguous — there is no other §1.8).
- `§7.x` tokens are overloaded: inside the former brief they mean the destructive-interference material (now **§4.14.x**, e.g. former §7.3.4 $\to$ §4.14.3.4), but `§7`/`§7.1`/`§7.2` also legitimately reference the companion *density-ratio* stage. Each occurrence must be read in context.
- `§4.5` (x5) is overloaded: it means either the former Doc-A §4.5 (seeding, now **§4.6**) or the former brief §4.5 (now **§4.11.5**).
- `§6` is overloaded: brief-internal pooling (now **§4.13**) vs companion RiteWeight.
- bare `§2`, `§3`, `§5.x` similarly: brief-internal (now §4.9 / §4.10 / §4.12.x) vs theory/architecture/companion.

Use the map in §1.1 above to update these. An orientation note was added after the Table of Contents stating that references to pipeline stages §5–§30 (tICA, RiteWeight, density-ratio, DMD, kEDMD, PCCA+, GP-FES, VAMPnet, validation, limitations) refer to the **companion Pipeline document**; this is the correct reading for the un-rewritten forward references that genuinely point outside this file (e.g., §9.4, §13, §21.3, §28-L2).

**I7 — Minor wording/typos (pre-existing, not corrected).**

- Abstract: *"we develop and approach"* $\to$ "we develop an approach."
- Brief intro: *"splitting and merging probabilty weights"* $\to$ "probability."
- Figure caption: *"Section numbers refer to chapter each stage is derived"* is ungrammatical; suggest "Section numbers refer to the chapter in which each stage is derived."
- §4.11.6: *"diverse bstates"* — "bstates" (WESTPA basis-state files) is used as jargon; fine internally, but define on first use for external readers.
- §3.2 reads *"appear in the first half of pipeline"* $\to$ "of the pipeline."

---

## 5. Bottom line

The mathematics is sound: every worked example checks out exactly, and the principal derivations (similarity transform, density-relaxation expansion, KL/Lyapunov expansion, KDE-noise and $\sqrt\pi$-weighting, pooling SNR and break-even, Prony bound, resampling unbiasedness) are correct. The statistical design is coherent and its main risk (destructive pooling interference, and its interaction with the run-axis bootstrap) is correctly identified and now carries explicit caveats. The outstanding items are reference/labeling consistency (I1–I7), the most consequential being the abstract's summation index (I1), the two-vs-three weight-type count (I2), the RevVAMPnet §15/§19 mismatch (I3), and the systematic cross-reference update implied by the renumbering (I6). Issues I1–I7 have since been resolved in the current document (see §6): the abstract summation index, the weight-type count, the RevVAMPnet label, and the wording items were corrected, and the in-text section-number cross-references were removed wholesale rather than remapped, which is why the §1.1 old→new map is retained here only as a historical record of the restructuring.


---

## 6. Revision log — alignment with the current document

This section records the changes made after the initial reorganization so that the report matches the document as delivered. None removed technical content; paragraph splitting was verified to preserve non-whitespace content exactly.

**Front matter and formatting.** The Table of Contents was moved to the very front (first content after the title block) at the author's request, and the stray thematic-break rules that produced the horizontal line and excess whitespace beneath the title were removed. The title is set to `\large` and the subtitle to `\normalsize`; the document body is set to `fontsize: 10pt`; and `titling` with `\setlength{\droptitle}{-3.5em}` pulls the title block upward to reduce the top whitespace. The YAML header — which an external formatter had repeatedly flattened (un-nesting `latex_engine`/`extra_dependencies` from under `pdf_document` and at one point duplicating `fontsize`, which silently dropped the `header-includes` block and caused a `\coloneqq` "undefined control sequence" failure) — was restored to correct nesting. Section 1 was retitled "Methodological Statement and Scope."

**Cross-references (resolves I3, supersedes I6).** Every in-text cross-reference carrying a specific section number was removed from the prose; the figure's internal `\S` stage labels were retained. This moots the RevVAMPnet §15/§19 mismatch and replaces the remap-table approach of I6.

**Corrections (I1, I2, I4, I7).** The abstract summation index was corrected to begin at $i\ge 2$; the two-vs-three weight-type wording was reconciled; the dangling $H$-theorem ("de Bruijn form in §2.7") reference was repaired; and the wording/typo items were fixed.

**Expanded weighted-ensemble methodology (§4.1–§4.2).** The WESTPA stage was expanded to cover the actual WE algorithm and theory in more detail: the von Neumann splitting / Russian-roulette origins and the Huber–Kim binning refinement; the framing of WE as unbiased importance sampling in trajectory space; the bias-free-splitting / variance-bounded-merging decomposition and the variance-optimal allocation and recent mathematical-developments literature; and the steady-state (NESS, Hill-relation flux $\mathrm{MFPT}=1/J_{\mathrm{SS}}$) versus equilibrium (WEED) operating modes, with the pipeline's exclusive use of the equilibrium mode motivated by the coverage requirement. New citations were added accordingly (Huber–Kim, Zhang–Jasnow–Zuckerman, Suárez et al., Zwier et al., Russo et al., Bhatt–Zhang–Zuckerman, Aristoff–Zuckerman, Aristoff et al., Ryu et al., Zuckerman–Chong, Copperman–Zuckerman, and the RiteWeight references of Kania et al. and Otten et al.).

**Symbol standardization.** Configuration space, previously written as both $\Omega$ and $\mathcal{X}$, is unified to $\Omega$ throughout (the feature set $\mathcal{F}$ is left untouched); the thermal energy is uniformly $k_B T$; and "Schmid-projected" is hyphenated consistently.

**Readability and math display.** Twenty-nine long paragraphs were split at genuine sentence boundaries (abbreviations, decimals, and inline math protected), and the one extensive inline set-builder — the WE-corpus definition $\mathcal{C}=\{\dots\}$ — was promoted to display math.

**Citations and reference list.** All citations were renumbered contiguously [1]–[25] in order of first appearance, and a numbered reference list was appended. Bibliographic details were drawn from the companion *Foundational Literature* guide; the block-/Bayesian-bootstrap classics (Künsch 1989; Liu & Singh 1992; Rubin 1981) are not in that guide and were supplied from their standard sources.
