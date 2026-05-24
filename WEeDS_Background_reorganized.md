---
title: "\\large WEeDS:: Weighted Ensemble equilibrium-mode Dispersion Seeding"
subtitle: "\\normalsize WE Equilibrium Dynamics (WEED) and multi-basin seeding to drive coverage of conformational space"
fontsize: 10pt
output:
  pdf_document:
    keep_tex: true
    latex_engine: lualatex
    extra_dependencies: ["mathtools", "amsmath", "amssymb", "tikz", "booktabs", "longtable"]
monofont: "Menlo"
header-includes:
 - \usepackage{amsmath}
 - \usepackage{amssymb}
 - \usepackage{mathtools}
 - \usepackage{tikz}
 - \usepackage{booktabs}
 - \usepackage{longtable}
 - \usepackage{pdflscape}
 - \usetikzlibrary{positioning, arrows.meta, shapes.geometric, calc}
 - \definecolor{primaryfill}{RGB}{232,240,247}
 - \definecolor{primaryborder}{RGB}{45,85,130}
 - \definecolor{coreborder}{RGB}{20,55,95}
 - \definecolor{corefill}{RGB}{210,225,240}
 - \definecolor{valfill}{RGB}{245,243,238}
 - \definecolor{valborder}{RGB}{120,110,95}
 - \definecolor{outfill}{RGB}{248,242,232}
 - \definecolor{outborder}{RGB}{150,110,60}
---

*The object of interest is not an individual trajectory — it is the time-evolving probability distribution of the weighted ensemble itself.*

------------------------------------------------------------------------

## Table of Contents

-   [0. Abstract](#0-abstract)
-   [1. Program Statement and Scope](#1-program-statement-and-scope)
-   [2. Theoretical Foundations](#2-theoretical-foundations)
-   2.1 The free energy surface and why it is hard to compute
-   2.2 Molecular dynamics as a stochastic process
-   2.3 The transfer operator and its spectrum
-   2.4 Metastability and the timescale gap
-   2.5 The VAMP variational principle and the Galerkin viewpoint
-   2.5.1 Estimator equivalence and complementarity (DR-DMD, kEDMD, VAMPnet)
-   2.6 The density-relaxation expansion
-   2.6.1 The amplitude--temporal--spatial decomposition
-   2.6a The density ratio as a natural observable of the Koopman operator
-   2.6b Noise structure of the density-ratio estimator
-   [3. Pipeline Architecture](#3-pipeline-architecture)
-   3.1 Overview and data flow
-   3.2 Walker weight types and their roles
-   3.3 The two RiteWeight configurations
-   [4. Stage 1 --- Weighted Ensemble Sampling (WESTPA)](#4-stage-1--weighted-ensemble-sampling-westpa)
-   4.1 The weighted-ensemble strategy
-   4.2 Resampling rules and probability conservation
-   4.3 The WE corpus as an importance-weighted sample
-   4.4 The non-equilibrium history of individual walkers
-   4.5 Equilibrium WE and the WEED driver
-   4.6 Progress-coordinate choice and bidirectional seeding
-   4.7 Progress-coordinate independence and its structural limits
-   4.8 What WESTPA produces
-   4.9 Mode excitation, precisely
-   4.10 Coverage of supp(pi) versus barrier crossings
-   4.11 The coverage--excitation tradeoff and bidirectional seeding
-   4.12 Identifiability conditions and the window of observability
-   4.13 Pooling multiple WE runs
-   4.14 Destructive interference in pooling
-   4.15 Operational recommendations

------------------------------------------------------------------------

> **Scope and cross-references.** This is the self-contained *Background* document for the WEeDS pipeline; it covers the theoretical foundations, the pipeline architecture, and the weighted-ensemble sampling stage, including the bidirectional-seeding and pooling analysis. References in the text to downstream pipeline stages denote stages developed in full in the **companion Pipeline document**.

------------------------------------------------------------------------

## 0. Abstract

Unlike traditional methods that require treatment of the unconverged nature enhanced sampling methods due to the sampling problem and rare event timescales, we develop an approach that leverages the non-equilibrium structure of weighted-ensemble (WE) output as the *signal source* for transfer-operator spectral decomposition, rather than treating it as bias to be corrected away. Existing WE analysis pipelines — history-augmented Markov state models, standard reweighting-to-equilibrium, neural operator learning on reweighted trajectory pairs — all treat the iteration-indexed sequence of transient distributions $\{p(x, t_k)\}$ as an obstacle to equilibrium analysis. The density-ratio construction inverts this: the transient distributions are exactly the data the spectral expansion $p(x,t_k)/\pi(x) - 1 = \sum_{i\ge 2} c_i\,\lambda_i^{t_k/\tau_{\mathrm{WE}}}\,\psi_i(x)$ consumes. WE's non-equilibrium character is what makes the decomposition identifiable; longer brute-force equilibrium MD cannot produce this signal because it does not sample the transient relaxation modes by which the operator spectrum is exposed.

A single WE run on a generous, deliberately over-complete feature set $\mathcal{F}$ produces transient distributions whose density-ratio decomposition yields transfer-operator eigenfunctions $\{\psi_i\}$ and eigenvalues $\{\lambda_i\}$ directly — no grid discretization, no microstate clustering, no neural-network training, no manual progress-coordinate choice beyond the initial feature-set specification. The recovered eigenfunctions are smooth analytical functions on $\mathcal{F}$, evaluable at arbitrary new walker configurations via the kernel representation. A Gaussian-process posterior variance field in eigenfunction coordinates, produced by the dynamical-kernel GP over RiteWeight-reweighted walker observations, simultaneously tracks eigenfunction resolution, density-ratio noise, and free-energy uncertainty in a single scalar field. The pipeline developed here also supplies the structural commitments — the kernel-representation eigenfunctions (evaluable at arbitrary new walker configurations) and the dynamical-kernel GP posterior variance (a single scalar uncertainty oracle) — that make *iterative coordinate-discovery* tractable: the recovered eigenfunctions can serve as progress coordinates of subsequent WE runs, and the GP posterior variance supplies the bin-density criterion in the refined coordinates, closing the sampling–analysis loop around the transfer-operator spectrum. The method offers a substantial step toward data driven progress coordinate discovery and refinement, targeting the 'dark art' of collective variable identification.

The approach treats the WE approach as a method for driving adequate coverage of conformational space on both sides of the free energy barrier in contrast to the traditional path sampling approach which drives sampling along a single reaction path seeded from a single basin toward a target state. The present method relies on WE Equilibrium Dynamics (WEED) to generate a corpus of weighted walkers at each iteration that are unbiased estimates of the probability density as the system relaxes toward the stationary equilibrium distribution. Dynamics remain unbiased and decoupled from binning along chosen progress coordinates, as in traditional WE approaches. The central methodological reframing shifts the goal of WE runs as a rare event path sampling strategy to a massively parallel method for dispersing simulations across conformational space without the explicit requirement of barrier crossing observation to extract mechanistic, kinetic, and thermodynamic information from WE results.

------------------------------------------------------------------------

## 1. Methodological Statement and Scope

We briefly describe the larger computational pipeline for estimating time-resolved and equilibrium free energy surfaces (FES) from weighted ensemble (WE) molecular-dynamics simulations via recovery of the slow eigenfunctions of the system's transfer operator. The pipeline is organized around the observation that the FES is most cleanly defined not on an externally imposed collective coordinate but on the leading transfer-operator eigenfunctions themselves: the slow eigenfunctions are, by construction, the coordinates on which metastability is most sharply expressed, and the barrier regions that confound externally imposed coordinates are replaced by interpolation regions between eigenfunction basin assignments. WE sampling (WESTPA 2.0, Russo et al.\textsuperscript{[8]}) supplies the raw corpus; every downstream stage is an inverse problem targeting one of three objects — the stationary distribution $\pi$, the transfer-operator spectrum $\{(\lambda_i,\psi_i)\}$, or the FES on the eigenfunction embedding — for which the identifiability conditions and estimator consistency are derived. The present document defines the WE method and addresses the ontological stance with regard to total probability density as the object of interest.

The **primary pipeline** is

$$
\text{WESTPA} \;\rightarrow\; \text{IW-tICA} \;\rightarrow\; \text{thermodynamic RiteWeight} \;\rightarrow\; \text{kernel-evaluated density-ratio matrix } \mathbf{R} \;
$$ $$
\rightarrow\; \pi\text{-weighted Exact DMD} \;\rightarrow\; \text{PCCA+} \;\rightarrow\; \text{kinetic RiteWeight} \;\rightarrow\; \text{GP FES with dynamical kernel},
$$

The pipeline is organized around **three spectral estimators of the same transfer operator**, each consuming a different data structure derived from the WE corpus. The primary estimator is kernel-evaluated DR-DMD, which consumes iteration-indexed density snapshots. The secondary estimator is kernel Extended Dynamic Mode Decomposition (kEDMD) applied to RiteWeight-reweighted trajectory pairs $\{(x_i^{(k,\mathrm{start})}, x_i^{(k,\mathrm{end})}, \hat{w}_i^{\pi})\}$ aggregated across all walkers and iterations, with the importance weight $\hat{w}_i^{\pi}$ applied to the pair endpoint in the Wu–Nüske–Paul\textsuperscript{[22]} convention. kEDMD operates in the same tICA coordinate space and uses the same isotropic Matérn-$5/2$ kernel with the same bandwidth $h$ as the primary route, so the two routes address the same Koopman operator in the same Sobolev RKHS $\mathbb{H}_h = H^\sigma(\Omega)$ at order $\sigma = 5/2 + d_{\mathrm{tIC}}/2$ (Bold 2025\textsuperscript{[117]} / Hertel 2026\textsuperscript{[118]} unconditional Koopman-invariance certificate) and must agree under reversibility and sampling convergence to within bootstrap CI. Agreement on eigenvalues to within bootstrap confidence intervals, together with $L^2(\pi)$ eigenfunction inner products $|\langle \hat{\psi}_i^{\mathrm{DR}}, \hat{\psi}_i^{\mathrm{kEDMD}} \rangle_{\pi}| \ge 0.9$, is the primary spectral consistency check. Disagreement is diagnostic: eigenvalue disagreement points to convergence failure in one or both routes, and eigenfunction disagreement with matching eigenvalues points to degenerate or near-degenerate modes that require Hankel-DMD or longer simulation. The tertiary estimator is a RevVAMPnet operating on full-feature trajectory pairs without tICA preprocessing; it provides an independent nonlinear reference spectrum not confined to the tIC subspace. Three-way agreement (DR-DMD $\approx$ kEDMD $\approx$ VAMPnet) is the strongest spectral validation criterion; disagreement between VAMPnet and the two tIC-space routes is diagnostic of tICA linearity limitations or incomplete tIC coverage of the slow manifold.

The pipeline is grounded in two complementary theoretical anchors. The unifying transfer-operator framework of Klus et al.\textsuperscript{[29]} places tICA, DMD, kEDMD, and VAMP as special cases of one generalized eigenvalue problem targeting the Koopman/transfer operator spectrum, which justifies the three-way spectral cross-validation architecture and motivates the shared RKHS in which all spectral estimators operate. The off-equilibrium Koopman estimation framework of Wu, Nüske, Paul, Klus, Koltai & Noé\textsuperscript{[22]} provides the specific algorithmic templates: Algorithm 1 proves that the non-symmetrized covariance estimator used in importance-weighted tICA is unbiased and consistent on off-equilibrium data, removing any equilibrium-WE-convergence prerequisite for the tICA stage, and Algorithm 3 is the explicit-basis template that the kEDMD secondary route extends to the RKHS, with RiteWeight stationary weights $\hat{w}_i^{\pi}$ substituting for the internally computed Koopman-reweighting weights of Wu et al.'s Algorithm 2. Together these two anchors fix the methodological content of every spectral stage in this document and identify the pipeline as a **non-equilibrium WE extension of an established off-equilibrium spectral-estimation framework**, rather than a collection of ad hoc estimators. 

Together these make a **self-consistent coordinate-discovery loop** algorithmically tractable: the eigenfunctions from WE run $n$ become the progress coordinates for run $n+1$, and the GP variance field from run $n$ determines the bin scheme in the new coordinates. Successive runs refine eigenfunction map, stationary distribution, and transient-ensemble relaxation spectrum *together*, terminating at mutual consistency rather than at convergence of any single observable. 

Uncertainty quantification in the pipeline is executed as a **single-loop joint bootstrap on the pooled corpus**. Per-walker observation variance for the GP-FES training data is derived by propagating the pointwise kernel-evaluated density-ratio mean-squared error through the $-\ln$ transformation (the delta method). Distributional uncertainty on every other reported observable — eigenvalues, implied timescales, eigenfunction overlaps, MFPTs, rates, committors, barrier heights — is estimated by a single replicate loop in which each replicate (i) draws a Dirichlet weight vector $\boldsymbol{\alpha}^{(b)} \sim \mathrm{Dir}(1,\ldots,1)$ on the $N_{\mathrm{runs}}$ run axis (Rubin\textsuperscript{[127]}; Mostofian & Zuckerman\textsuperscript{[11]}), (ii) block-resamples the iteration sequence within each run with block length $b_{\mathrm{block}}$ chosen to exceed the empirical iteration-correlation time (Künsch\textsuperscript{[78]}; Liu–Singh\textsuperscript{[79]}), (iii) rebuilds the pool with bootstrap run + iteration weights, and (iv) re-executes the full pipeline (RiteWeight refit, density-ratio reassembly, $\pi$-weighted DMD, kEDMD, PCCA+, weighted MSM, GP-FES) on the rebuilt pool. Headline percentile CIs are taken in log space for high-log-variance observables — following Mostofian & Zuckerman\textsuperscript{[11]} whose analysis established that the standard percentile bootstrap produces anti-conservative or unphysical lower CIs on small-sample high-log-variance observables and that the Dirichlet-weighted Bayesian bootstrap on the run axis is the correct alternative — and in arithmetic space otherwise. Variance is decomposed post-hoc into inter-run and intra-run contributions by importance-weighted reweighting of the same $B$ replicates toward the marginal on each axis. Per-state Kish effective sample size is retained as a baseline independent-sample diagnostic and as the analytical motivation for the correlation-aware block correction. The previous nested two-level scheme (separate $B_{\mathrm{outer}} \times B_{\mathrm{inner}}$ pipeline executions) is superseded; its statistical content is preserved in a single replicate loop on the pooled corpus.

The remainder of the document develops the WE equilibrium mode approach and discusses the use of a novel bidirectional seeding scheme to improve state space coverage while maintaining sufficient excitation of slow modes for eigenfunction extraction via density ratio dynamic mode decomposition.

------------------------------------------------------------------------

## 2. Theoretical Foundations

### 2.1 The free energy surface and why it is hard to compute

The central object of interest is the **free energy surface** (FES), also called the potential of mean force, along a set of collective coordinates $q(x)$: $$
F(q) = -k_B T \ln p(q) + \text{const},
$$ where $p(q)$ is the marginal probability density of the collective coordinate under the Boltzmann distribution and $k_B T$ is the thermal energy. For a system with Hamiltonian $H(x)$ over a configuration space $\mathcal{X}$, the Boltzmann distribution is $$
\pi(x) = Z^{-1} \exp\!\left(-\frac{H(x)}{k_B T}\right), \qquad Z = \int_{\mathcal{X}} \exp\!\left(-\frac{H(x)}{k_B T}\right) dx,
$$ and the marginal of $q$ is obtained by integrating out all degrees of freedom not captured by $q$: $$
p(q) = \int_{\mathcal{X}} \delta(q(x') - q)\,\pi(x')\,dx'.
$$ Substituting, $$
F(q) = -k_B T \ln \int_{\mathcal{X}} \delta(q(x') - q)\,\exp\!\left(-\frac{H(x')}{k_B T}\right) dx' + \text{const},
$$ which carries the intuitive interpretation that the free energy at $q$ is proportional to the negative logarithm of the total statistical weight of all configurations that map to $q$. Free-energy minima are densely populated; free-energy maxima are exponentially rare, and their statistical weight is swamped by that of the minima in any finite sample drawn from $\pi$.

The difficulty in computing $F(q)$ from unbiased molecular dynamics is therefore immediate. The Boltzmann distribution is sharply concentrated near free-energy minima, and for a barrier of height $\Delta F^{\ddagger}$ the probability of a spontaneous barrier-crossing event in brute-force molecular dynamics at temperature $T$ is proportional to $\exp(-\Delta F^{\ddagger}/k_B T)$. For the biomolecular systems this pipeline targets — RNA loops, peptide folding, small-molecule conformational rearrangements — barriers of $10$–$30\,k_B T$ are typical, which means spontaneous crossings occur $e^{10}$ to $e^{30}$ times less frequently than intra-basin fluctuations. The mean first-passage time over such a barrier is measured in microseconds to milliseconds; individual MD timesteps are femtoseconds. Brute-force sampling is therefore computationally intractable for all but the smallest barriers, and every rare-event method that has been proposed in the last four decades — umbrella sampling, metadynamics, replica exchange, milestoning, forward-flux sampling, weighted ensemble — can be read as a different strategy for amplifying the rare configurations that carry the signal about $F(q)$ at and above the barrier.

Weighted ensemble (WE) resampling, the sampling instrument for this pipeline, exploits the rare-event-amplification problem without biasing the underlying dynamics: walkers propagate under exact, unbiased MD, and the rare-event amplification is achieved purely through a resampling rule that equalizes walker populations across progress-coordinate bins while preserving the total probability carried by each bin. The resulting walker corpus is an importance-weighted sample of the transient distribution $p(x, t_k)$ at each WE iteration $k$, not a sample of the Boltzmann distribution $\pi$; recovering $\pi$ from this importance-weighted transient sample is the job of the downstream stages of the pipeline, beginning with tICA featurization and thermodynamic RiteWeight and continuing through the density-ratio spectral analysis that ultimately defines the eigenfunction coordinate system on which the FES is constructed.

### 2.2 Molecular dynamics as a stochastic process

The molecular systems targeted by this pipeline are modeled as continuous-time Markov processes on a configuration space $\Omega \subseteq \mathbb{R}^{3N}$, where $N$ is the number of atoms retained after any reduction to a coarse-grained representation (solute heavy atoms, coarse-grained beads, or a selected subset of internal coordinates). The canonical dynamical model is the overdamped Langevin (Itô diffusion) equation $$
dx_t = -\nabla U(x_t)\,dt + \sqrt{2 \beta^{-1}}\,dW_t,
$$ where $U: \Omega \to \mathbb{R}$ is the potential of mean force at inverse temperature $\beta = 1/(k_B T)$, $W_t$ is a standard $\dim(\Omega)$-dimensional Wiener process, and the friction tensor is implicitly absorbed into the rescaling of time. This is the overdamped limit appropriate for biomolecular simulations in explicit solvent after integrating out the momenta; it captures the essential features of solvent-averaged biomolecular dynamics on timescales long compared to momentum decorrelation and is the dynamical assumption under which the transfer-operator spectral theory is formulated. The infinitesimal generator of this process is $$
\mathcal{L} = -\nabla U \cdot \nabla + \beta^{-1} \Delta,
$$ acting on smooth observables, and the forward Kolmogorov (Fokker–Planck) equation governing the evolution of probability densities $\rho_t(x)$ is $$
\frac{\partial \rho_t}{\partial t} = \mathcal{L}^{*} \rho_t, \qquad \mathcal{L}^{*} \rho = \nabla \cdot (\rho\,\nabla U) + \beta^{-1} \Delta \rho,
$$ where $\mathcal{L}^{*}$ is the $L^2$-adjoint of $\mathcal{L}$. The stationary distribution $\pi$ defined above satisfies $\mathcal{L}^{*}\pi = 0$; under mild regularity conditions on $U$ (confinement at infinity and sufficient smoothness) it is the unique invariant measure of the Fokker–Planck semigroup.

Underdamped Langevin and Newtonian-with-thermostat dynamics, which are the integrators actually used in the atomistic MD simulations driving the WE engine, reduce to the same effective transfer operator in the appropriate limit and on the coarse-grained time axis relevant to slow-process identification, so the spectral framework developed here applies to them as well. The relevant timescale separation is between the friction-induced momentum decorrelation time (picoseconds) and the slow conformational modes we seek (nanoseconds to microseconds and longer). Above the former and below the latter, the overdamped Langevin description is faithful; the Koopman and Perron–Frobenius operators are then the correct linear operators governing the evolution of observables and densities.

Three distributions are distinguished throughout the document and are not interchangeable. The **stationary distribution** $\pi(x)$ is the equilibrium target. The **transient distribution** $p(x, t_k)$ is the distribution of walker positions produced by WE at iteration $k$; under the Zhang–Jasnow–Zuckerman\textsuperscript{[2]} unbiasedness result, $p(x, t_k)$ is the true distribution at physical time $t_k = k\,\tau_{\mathrm{WE}}$ starting from the WE initial distribution $p(x, 0)$, and the raw WE iteration weights $w_i^{(k)}$ are unbiased estimators of this distribution in the sense $$
\mathbb{E}\!\left[\sum_i w_i^{(k)} f(x_i^{(k)})\right] = \int_{\Omega} f(x)\,p(x, t_k)\,dx
$$ for any bounded test function $f$. The **aggregate WE sampling distribution** $q(x)$, defined as the $\tau_{\mathrm{WE}}$-weighted empirical distribution of all walker positions pooled across all iterations, is the distribution that RiteWeight reweights against $\pi$ to produce the stationary weights $\hat{w}_i^{\pi} \approx \pi(x_i)/q(x_i)$. The distinction matters because different pipeline stages consume different distributions: tICA covariance estimators are consistent estimators of $q$-weighted covariances (which, via the Wu–Nüske–Paul Algorithm 1 unbiasedness argument, recover the correct non-symmetrized Koopman projection regardless of the mismatch between $q$ and $\pi$); the density-ratio numerator at iteration $k$ uses $p(x, t_k)$; the density-ratio denominator and the kEDMD endpoint reweighting use $\pi$.

Two structural properties of the dynamics are assumed throughout the spectral theory and must be verified empirically wherever possible.

1.  **Reversibility** is equivalent to the detailed-balance condition on the transition density $p_\tau(y \mid x)$, $$
    p_\tau(y \mid x)\,\pi(x) = p_\tau(x \mid y)\,\pi(y) \qquad \text{for all } x, y \in \Omega,\; \tau \ge 0,
    $$ which is guaranteed by construction for the overdamped Langevin process with gradient drift above and holds to excellent approximation for thermostatted atomistic MD. Reversibility implies that the transfer operator is self-adjoint in $L^2(\pi)$ and has a real spectrum, which is the regime in which DMD and kEDMD eigenvalues can be interpreted as real decay rates and eigenfunctions as orthonormal slow-mode coordinates.

2.  **Ergodicity**, i.e. the absence of non-trivial invariant sets that partition $\Omega$ into dynamically disconnected regions, guarantees that every slow mode can in principle be excited by some initial condition and observed in the WE corpus; its failure on practical simulation timescales produces the metastability that motivates WE sampling in the first place. Within each connected component of the free-energy landscape reachable under WE resampling, the dynamics restricted to that component are reversible and ergodic on the pipeline's lag time $\tau_{\mathrm{WE}}$, and the spectral theory applies component-wise.

Throughout this document, $x \in \Omega$ denotes a full-dimensional configuration, $\chi(x) \in \mathbb{R}^{F}$ denotes a feature vector, $q(x) \in \mathbb{R}^{d}$ with $d = 2$–$4$ denotes the tICA projection to slow-mode coordinates (notation chosen to match the kernel-argument convention), and $\Psi(x) = (\psi_2(x), \ldots, \psi_{r+1}(x)) \in \mathbb{R}^{r}$ denotes the eigenfunction-coordinate embedding recovered from DR-DMD. The physical lag time is $\tau = \tau_{\mathrm{WE}}$ throughout, and WE iteration $k$ corresponds to physical time $t_k = k\,\tau_{\mathrm{WE}}$. Where confusion between the sampling distribution and the tICA projection is possible, the latter is written $q(\cdot)$ with an explicit argument and the former is left unargumented.

### 2.3 The transfer operator and its spectrum

The dynamical content of the stochastic process is captured, for fixed lag time $\tau > 0$, by two adjoint linear operators on function spaces over $\Omega$. The **Koopman (transfer) operator** $\mathcal{K}_\tau$ acts on observables $f \in L^2(\pi)$ by $$
(\mathcal{K}_\tau f)(x) = \mathbb{E}[f(x_\tau) \mid x_0 = x] = \int_{\Omega} f(y)\,p(y \mid x; \tau)\,dy,
$$ where $p(y \mid x; \tau)$ is the transition density from $x$ to $y$ over lag $\tau$. The **Perron–Frobenius operator** $\mathcal{P}_\tau$ acts on densities $\rho \in L^1(\Omega)$ by $$
(\mathcal{P}_\tau \rho)(y) = \int_{\Omega} p(y \mid x; \tau)\,\rho(x)\,dx,
$$ propagating a probability density forward in time. The Koopman operator propagates observables *pulled back* along trajectories; the Perron–Frobenius operator propagates probability densities *pushed forward* along trajectories. Both act on the same transition kernel but on dual spaces, and the pair $(\mathcal{K}_\tau, \mathcal{P}_\tau)$ contains the same dynamical information.

For reversible dynamics with stationary distribution $\pi$, the two operators are related by a similarity transform with weight $\pi$: if $$r(x, t) = p(x, t)/\pi(x) - 1$$ denotes the density ratio (the central object), then a direct calculation using detailed balance yields $$
(\mathcal{P}_\tau \rho)(y)/\pi(y) = (\mathcal{K}_\tau (\rho/\pi))(y).
$$ The Perron–Frobenius action on a density $\rho$ is equivalent, after division by $\pi$, to the Koopman action on the same density divided by $\pi$. This identity is the theoretical bridge between the density-ratio viewpoint (in which the data matrix $\mathbf{R}$ contains evaluations of $r(x, t)$) and the Koopman viewpoint (in which DMD, kEDMD, and VAMP are spectral estimators). The two viewpoints are mathematically equivalent under reversibility; the DR-DMD primary route uses the density-ratio viewpoint because WE produces density snapshots rather than trajectory pairs, while the kEDMD and VAMPnet routes use the Koopman viewpoint because they consume trajectory pairs directly. Both extract the same spectrum from the same underlying operator.

Detailed balance, $$
p_\tau(x' \mid x)\,\pi(x) = p_\tau(x \mid x')\,\pi(x'),
$$ implies that $\mathcal{K}_\tau$ is **self-adjoint** with respect to the $\pi$-weighted inner product on $L^2(\pi)$, $$
\langle f, g \rangle_{\pi} \coloneqq \int_{\Omega} f(x)\,g(x)\,\pi(x)\,dx = \mathbb{E}_{\pi}[f(x)\,g(x)],
$$ that is, $\langle \mathcal{K}_\tau f, g \rangle_{\pi} = \langle f, \mathcal{K}_\tau g \rangle_{\pi}$ for all $f, g \in L^2(\pi)$. By the spectral theorem for self-adjoint operators, $\mathcal{K}_\tau$ admits a spectral decomposition $$
(\mathcal{K}_\tau f)(x) = \sum_{i=1}^{\infty} \lambda_i(\tau)\,\langle \psi_i, f \rangle_{\pi}\,\psi_i(x), \qquad \langle \psi_i, \psi_j \rangle_{\pi} = \delta_{ij},
$$ where the $\{\psi_i\}$ form a complete orthonormal basis in $L^2(\pi)$ and $\lambda_1 = 1 \ge |\lambda_2| \ge |\lambda_3| \ge \cdots$ are real eigenvalues with $\lambda_1 = 1$ (since $\mathcal{K}_\tau \mathbf{1} = \mathbf{1}$).

Each eigenfunction $\psi_i$ together with its eigenvalue $\lambda_i$ describes a distinct slow relaxation mode. The leading eigenfunction $\psi_1 \equiv 1$ with eigenvalue $\lambda_1 = 1$ is the trivial stationary mode, corresponding to the equilibrium distribution itself, which is fixed by $\mathcal{K}_\tau$. The sub-leading mode $\psi_2$ with eigenvalue $|\lambda_2| < 1$ is the slowest non-trivial mode, corresponding to the slowest conformational transition; for a system with two long-lived conformations (folded versus unfolded), $\psi_2(x) \approx +c$ in one basin and $\psi_2(x) \approx -c$ in the other, and its sign partitions configuration space into the two basins. The next mode $\psi_3$ with $|\lambda_3| \le |\lambda_2|$ corresponds to the second-slowest transition (an alternative folding pathway, a subdomain rearrangement, or a tertiary contact rearrangement, depending on the system).

The eigenvalue $\lambda_i(\tau)$ decays exponentially with lag time $\tau$ according to the Chapman–Kolmogorov property of the semigroup, $$
\mathcal{K}_{\tau_1 + \tau_2} = \mathcal{K}_{\tau_1} \circ \mathcal{K}_{\tau_2} \implies \lambda_i(\tau_1 + \tau_2) = \lambda_i(\tau_1)\,\lambda_i(\tau_2).
$$ The unique solution to this functional equation is $\lambda_i(\tau) = e^{-\tau/\tau_i}$, where $$
\tau_i = -\frac{\tau}{\ln |\lambda_i(\tau)|}
$$ is the **implied (relaxation) timescale** associated with mode $i$. Crucially, $\tau_i$ is independent of the lag time $\tau$ used to compute it; this is the Markov property, and the relaxation rates are intrinsic to the dynamics rather than to the observation window. Empirically verifying lag-time independence of $\tau_i$ across a range of $\tau$ is the central convergence diagnostic for any spectral estimator of $\mathcal{K}_\tau$ — including DR-DMD, kEDMD, and VAMPnet — and is the content of the implied-timescales test discussed.

Substituting the spectral decomposition into the transfer operator yields the eigen-expansion of the transition density, $$
p_\tau(x' \mid x) = \pi(x')\,\sum_{i=1}^{\infty} \lambda_i(\tau)\,\psi_i(x)\,\psi_i(x'),
$$ which expresses the transition density as a sum of separable terms, each decaying exponentially at rate $1/\tau_i$. For large $\tau$, only the stationary term $\pi(x')$ survives and the system has lost memory of its initial condition — it has equilibrated. The same spectral decomposition, applied to a time-evolved density $p(x, t) = (\mathcal{P}_t p_0)(x)$ via the similarity transform above, yields the density-relaxation expansion — the operational identity that the DR-DMD primary route inverts.

**Operating prerequisites for the spectral framework.** The derivation above rests on three operating prerequisites that the pipeline assumes hold for every system it analyses, and which the validation panel audits empirically rather than presuming.

*(i) Reversibility.* The self-adjointness of $\mathcal{K}_\tau$ on $L^2(\pi)$ used to derive the orthonormal spectral decomposition (337) presumes detailed balance. For molecular dynamics under standard stochastic thermostats (Langevin, Andersen, Bussi-Donadio-Parrinello), reversibility is the operating assumption; on such systems the weighted-MSM validator constructs a single reversible MLE $\hat T_{\mathrm{rev}}$, and finite-sample violations of detailed balance in the RiteWeight-weighted count statistic are reported as the detailed-balance residual $\overline{\mathrm{DB}}$ for diagnostic action rather than absorbed by maintaining a parallel non-reversible matrix. Driven systems with explicit external forces or constant-gradient boundary conditions are out of scope for the present methods document — their transfer operators are not self-adjoint under any choice of $L^2$ inner product, and their analysis routes to the dedicated Extensions document where the dual-$T$ apparatus (parallel construction of $\hat T_{\mathrm{rev}}$ and $\hat T_{\mathrm{nonrev}}$) is specified.

*(ii) Markovianity at lag* $\tau_{\mathrm{WE}}$. The Chapman–Kolmogorov property used to derive the implied-timescale relation (348) requires the dynamics at the working lag $\tau_{\mathrm{WE}}$ to be approximately Markovian — equivalently, the gap between the slowest mode retained by the spectral truncation and the fastest mode discarded must be sufficient that memory effects below $\tau_{\mathrm{WE}}$ have decayed. The implied-timescales lag-scan and the Chapman–Kolmogorov test are the diagnostics for this assumption; failure of either flags the lag as too short and the Hankel-DMD extension or a longer $\tau_{\mathrm{WE}}$ as the remediation.

*(iii) Ergodicity within the accessible support of* $\pi$. The spectral expansion (337) and the time-lagged covariance identities are stated for $\mathcal{K}_\tau$ acting on $L^2(\pi)$ with $\pi$ the *true* equilibrium distribution. Finite WE coverage may sample only an accessible subset of $\operatorname{supp}(\pi)$ within the simulation horizon — landscape regions separated from the corpus by barriers higher than the WE budget can cross will be absent from the empirical $\hat\pi$, and the recovered spectrum is the spectrum of $\mathcal{K}_\tau$ restricted to the accessible support, not of the full equilibrium operator. The WESTPA convergence prerequisite and the initial-condition completeness limitation record the operational consequences. None of the three prerequisites is treated as automatically satisfied; together they constitute the standing-assumption layer beneath the rest of the methods document.

The generalized-eigenvalue-problem formulation of Klus et al.\textsuperscript{[29]} unifies the spectral estimators used throughout the pipeline: tICA, DMD, kEDMD, and VAMP all solve a single generalized eigenvalue problem $$
\hat{C}_{0\tau}\,v = \lambda\,\hat{C}_{00}\,v
$$ for differently chosen basis sets and inner products, with each estimator targeting a finite-rank approximation to $\mathcal{K}_\tau$ in the span of its basis. This unification is the mathematical foundation for the three-way spectral cross-validation architecture and is given explicit treatment.

### 2.4 Metastability and the timescale gap

A system is **metastable** if there is a large gap in the eigenvalue spectrum of $\mathcal{K}_\tau$ between $\lambda_n$ and $\lambda_{n+1}$, $$
1 = \lambda_1 > |\lambda_2| \approx \cdots \approx |\lambda_n| \;\gg\; |\lambda_{n+1}| \ge \cdots,
$$ which corresponds to a **timescale gap** $\tau_2 \approx \cdots \approx \tau_n \gg \tau_{n+1}$. Physically, the system has $n$ long-lived conformational basins (metastable states) that interconvert on timescales of order $\tau_2$, and within each basin the system equilibrates much faster, on timescale of order $\tau_{n+1}$. This is exactly the situation for biomolecular folding: the molecule has a few distinct states (folded, unfolded, and small numbers of intermediates) with long lifetimes, and rapid intra-basin fluctuations on much shorter timescales. The spectral gap at position $n$ is the quantitative marker of a clean metastable decomposition into $n$ states; the leading $n - 1$ non-trivial eigenfunctions $\psi_2, \ldots, \psi_n$ are approximately piecewise constant on the basins and encode the basin assignments via their sign and magnitude.

The size of the spectral gap $g = |\lambda_n| - |\lambda_{n+1}|$ is a harder constraint for DR-DMD than it is for VAMPnet. VAMPnet optimizes the VAMP-2 score and automatically ranks modes by their contribution to this variational criterion; a small gap produces a graceful degradation in which the recovered modes are linear combinations of the true modes but the leading $n$-dimensional subspace is still correctly identified. DR-DMD, by contrast, relies on the SVD rank selection within the Schmid projected algorithm to identify the slow subspace; when the gap is small, the SVD singular values corresponding to $\psi_n$ and $\psi_{n+1}$ cluster numerically, the DMD eigenvalues $\hat{\lambda}_n$ and $\hat{\lambda}_{n+1}$ appear as a near-degenerate pair (or, in the non-symmetrized estimator, as a complex-conjugate pair with small imaginary part), and the recovered eigenfunctions mix between the two modes. The operational consequence is that DR-DMD requires either a clean spectral gap ($g \gtrsim 0.1$ in typical biomolecular applications) or the **Hankel-DMD** extension to recover near-degenerate modes reliably. The three-way cross-validation architecture provides a mechanism for detecting this failure: disagreement between DR-DMD eigenvalues and kEDMD or VAMPnet eigenvalues in the near-degenerate regime is diagnostic and triggers the Hankel-DMD fallback.

For the GP-FES stage, the relevant spectral gap is $|\lambda_{r+1}| - |\lambda_{r+2}|$, where $r$ is the number of eigenfunctions retained as coordinates for the FES. A clean gap at this position defines the dimensionality of the slow manifold and is the empirical criterion for choosing $r$; implied-timescale plots identifying the ITS plateau \textsuperscript{[53]} are the standard diagnostic. In the absence of a clean gap, the FES is constructed on the largest $r$ for which the first $r$ eigenvalues agree across the three spectral routes, which is both a dimensionality criterion and a consistency check in one. The metastability assumption is thus not merely qualitative scene-setting: it determines the dimensionality of the eigenfunction coordinate embedding, the rank truncation in the DMD SVD, the number of macrostates in the PCCA+ assignment, and the size of the GP regression problem.

### 2.5 The VAMP variational principle and the Galerkin viewpoint

A key theoretical result underlies every spectral estimator in the pipeline: the slow eigenvalues of the Koopman operator $\mathcal{K}_\tau$ can be approximated by a variational principle, and every finite-rank estimator of the spectrum — tICA, DMD, kEDMD, VAMPnet — can be cast as a **Galerkin projection** of $\mathcal{K}_\tau$ onto a finite-dimensional function basis (Klus et al.\textsuperscript{[29]}; Wu & Noé\textsuperscript{[23]}; Koltai et al.\textsuperscript{[124]} for the propagation-error optimal-low-rank framing of the same Galerkin projection on non-equilibrium data, which classifies the WE pipeline within the regime-(iii) reversible-with-non-stationary-data taxonomy of that reference). The variational approach for Markov processes (VAMP) recasts Koopman-eigenfunction estimation as the optimization of a scalar score over a finite-dimensional function basis; the Galerkin viewpoint recasts it as a generalized eigenvalue problem on basis-expansion coefficients. Both views are equivalent under the assumptions of reversibility and sufficient basis richness, and the equivalence is the mathematical foundation of the three-way spectral cross-validation architecture.

**The generalized eigenvalue problem.** Given any set of $n$ trial functions $\chi_1, \ldots, \chi_n: \Omega \to \mathbb{R}$ assembled into a vector $\boldsymbol{\chi}(x) = (\chi_1(x), \ldots, \chi_n(x))^\top$, the best $n$-dimensional linear approximation to the leading eigenvalues of $\mathcal{K}_\tau$ is given by the generalized eigenvalue problem $$
\hat{C}_{0\tau}\,v = \lambda\,\hat{C}_{00}\,v,
$$ where the instantaneous and time-lagged covariance matrices are $$
\hat{C}_{00} = \mathbb{E}_{\pi}[\boldsymbol{\chi}(x_t)\,\boldsymbol{\chi}(x_t)^\top], \qquad \hat{C}_{0\tau} = \mathbb{E}_{\pi}[\boldsymbol{\chi}(x_t)\,\boldsymbol{\chi}(x_{t+\tau})^\top].
$$ The eigenvalues $\hat{\lambda}_1 \ge \hat{\lambda}_2 \ge \cdots \ge \hat{\lambda}_n$ of this generalized eigenvalue problem are the best $n$-dimensional approximations to the leading eigenvalues of $\mathcal{K}_\tau$; the corresponding eigenvectors $v_i$ supply the basis-expansion coefficients for the estimated eigenfunctions, $\hat{\psi}_i(x) = v_i^\top \boldsymbol{\chi}(x)$. As the trial-function space becomes richer — a deeper or wider neural network, a denser explicit dictionary, a finer kernel — the approximation converges to the true eigenvalues and eigenfunctions. This is the content of the variational principle: every finite-rank estimator under-approximates the true spectrum, and the estimator whose basis most closely spans the slow-mode subspace delivers the tightest approximation.

**The VAMP-2 score.** For the same basis $\{\chi_i\}_{i=1}^{n}$ and its Koopman-transformed counterpart $\{g_i\}_{i=1}^{n} = \{\mathcal{K}_\tau \chi_i\}$, the VAMP-2 score is $$
\mathcal{V}_2(\{\chi_i\}) = \sum_{i=1}^{n} \hat{\sigma}_i^{\,2},
$$ where $\hat{\sigma}_i$ are the singular values of the reweighted cross-correlation matrix $\hat{C}_{00}^{-1/2}\hat{C}_{0\tau}\hat{C}_{\tau\tau}^{-1/2}$ computed from trajectory pairs. The variational principle states that $\mathcal{V}_2$ is maximized, over all $n$-dimensional function bases in $L^2(\pi)$, when the basis spans the subspace of the top-$n$ Koopman eigenfunctions; the maximum value equals $\sum_{i=1}^{n} \lambda_i^{\,2}$, the sum of squared Koopman eigenvalues. This is the theoretical guarantee that underlies VAMPnet: neural networks trained to maximize $\mathcal{V}_2$ converge to bases that span the slow subspace of the Koopman operator, with the recovered eigenfunctions obtainable by singular-value decomposition of the post-training correlation matrix.

**The Galerkin-projection unification.** The density-ratio decomposition, tICA, EDMD, kEDMD, and VAMPnet are all instances of the same underlying operation: Galerkin projection of the transfer operator onto a finite basis. They differ in the choice of basis, the inner product under which the projection is taken, and the numerical route by which the resulting generalized eigenvalue problem is solved. The choices determine the tradeoffs, not the mathematical framework.

-   **tICA.** The basis consists of the input features $\chi(x)$ themselves (a linear basis in the feature space). The inner product is the $q$-weighted empirical inner product, which, via the Wu–Nüske–Paul Algorithm 1 unbiasedness result, converges to the $\pi$-weighted inner product without requiring equilibrium convergence of the WE corpus. tICA is the fastest, least expressive estimator; it is exact for linear dynamics and defines the geometric reduction that the other estimators build on.

-   **DR-DMD.** The basis consists of kernel evaluations centered at the walker reference set $\{x_j\}_{j=1}^{M}$ with $\pi$-weighting applied to the SVD (Schmid projected variant). The inner product is $L^2(\pi)$ with $\pi$ estimated by the kernel-evaluated sum weighted by $\hat{w}_i^{\pi}$. The Galerkin generalized eigenvalue problem becomes the density-ratio SVD/DMD problem on $\mathbf{R}$, which recovers $(\lambda_i, \psi_i)$ directly from density snapshots without requiring trajectory-pair data.

-   **kEDMD.** The basis is the Sobolev RKHS $\mathbb{H}_h = H^\sigma(\Omega)$ at order $\sigma = 5/2 + d_{\mathrm{tIC}}/2$, generated by the same isotropic Matérn-$5/2$ kernel $K_h$ used by DR-DMD, evaluated on the trajectory-pair endpoints. The inner product is the RKHS inner product, and the generalized eigenvalue problem is solved in the kernel matrix representation via the Wu–Nüske–Paul Algorithm 3 template with RiteWeight importance weights substituted for the Koopman-reweighting weights of Wu et al.'s Algorithm 2. kEDMD has infinite-dimensional basis expressiveness without neural-network training and operates in the same RKHS as DR-DMD, which is the mathematical basis for the primary/secondary cross-validation.

-   **VAMPnet.** The basis consists of the outputs of a neural-network featurizer $\boldsymbol{\chi}_{\theta}(x)$, optimized to maximize $\mathcal{V}_2$. The basis is nonlinear, learned from data, and can capture structure invisible to any fixed physical coordinate set or fixed kernel basis. VAMPnet is the most expressive estimator and the one whose structural assumptions are weakest; it is also the most computationally expensive and the one whose output is the hardest to interpret.

-   **EDMD.** The basis is an explicit dictionary of fixed functions (polynomials, radial basis functions, Hermite functions, etc.). EDMD is a classical approach that sits between tICA (linear in features) and kEDMD (infinite-dimensional RKHS). In this pipeline it appears only as a pedagogical waypoint; the primary explicit-basis route is the kernel variant.

All four estimators solve the same variational problem — finding the finite-dimensional subspace that best approximates the slow modes of $\mathcal{K}_\tau$ under the VAMP-2 score — in different function spaces. The density-ratio route has the advantage of simplicity and interpretability (no training, no explicit dictionary, minimal hyperparameters beyond the kernel bandwidth $h$) but the disadvantage of being confined to the tIC subspace used to build the walker reference set. VAMPnet has the opposite tradeoffs: maximal expressiveness at the cost of training expense and interpretability. kEDMD sits between the two: infinite-dimensional RKHS expressiveness with no training, at the cost of scaling with corpus size and sensitivity to the bandwidth choice.

DR-DMD does not itself optimize a variational objective — it solves a linear algebra problem (a generalized eigenvalue problem on the density-ratio matrix with $\pi$-weighting) rather than a non-convex optimization — and therefore does not inherit the graceful-degradation property of variational methods under poor basis choice. Where VAMPnet automatically identifies the best $n$-dimensional subspace within whatever function class its neural architecture spans, DR-DMD recovers whatever the data and the kernel-evaluated representation allow, with the quality of the recovery controlled by the kernel bandwidth, the SVD rank threshold, and the coverage of configuration space by the walker reference set. The VAMP-2 score's role is to quantify how close the DR-DMD output is to the variational optimum after the fact, not to drive the estimation.

The propagation-error optimal-low-rank framing of Koltai, Wu, Noé & Schütte\textsuperscript{[124]} (their Algorithm 1 and the time-lagged canonical correlation algorithm) is the explicit estimator-side companion to the unifying generalized eigenvalue problem of Klus, Nüske, Koltai, Wu, Kevrekidis, Schütte & Noé\textsuperscript{[29]}: under propagation error in $L^2(\pi)$, the rank-$k$-optimal Koopman approximant coincides with the projection on the dominant $k$ singular modes of the whitened operator (their Theorem A1, an Eckart–Young–Mirsky theorem for compact operators). This grounds the $\pi$-weighted Schmid-projected DMD primary route in a first-principles optimality argument rather than as a heuristic rank-truncation rule, and locates the DR-DMD, kEDMD, and tICA estimators as Galerkin projections of the same propagation-error-optimal estimator under different basis choices.

#### 2.5.1 Estimator equivalence and complementarity (DR-DMD, kEDMD, VAMPnet)

Both views recover the same operator. The Klus framework places tICA, DMD, kEDMD, and VAMP as Galerkin projections of $\mathcal{K}_\tau$ onto different bases; the three-way cross-validation architecture (DR-DMD, kEDMD, RevVAMPnet) checks that the three estimators agree, which is the strongest validation short of analytic ground truth.

The views differ in their data requirements and failure modes:

-   **DR-DMD** reads spectrum from temporal decay of iteration-indexed snapshots. Requires snapshot sequence spanning enough iterations for slow modes to decay observably, and initial distribution to project nontrivially onto modes (the excitation requirement).

-   **kEDMD** reads spectrum from time-lagged pair correlations. Requires coverage of $\mathrm{supp}(\pi)$ by pooled pairs and does not have an explicit excitation requirement at the same algorithmic level.

The two routes can fail in different ways on the same corpus. A corpus with poor mode excitation in the initial condition can still support kEDMD recovery if the pairs span enough of $\mathrm{supp}(\pi)$. The three-way agreement check catches both failure types.

### 2.6 The density-relaxation expansion

The density-relaxation expansion is the theoretical core of the DR-DMD primary route. Its content is that the time evolution of any transient density under the reversible dynamics admits an exact decomposition into a stationary component plus a sum of separable spatial-temporal modes, each mode the product of a transfer-operator eigenfunction and an exponentially decaying amplitude. The structure of this decomposition — separability of space and time, exponential decay in time, spatial form fixed by the eigenfunctions of $\mathcal{K}_\tau$ — is exactly the structure that dynamic mode decomposition is designed to recover from snapshot data. The significance of the expansion for this pipeline is not merely motivational: it establishes that WE's iteration-indexed density output contains, in principle, the full slow-mode spectral information of the underlying dynamics, and that this information can be extracted by any estimator that can invert a separable snapshot sequence. DR-DMD is one such estimator; kEDMD on trajectory pairs is another, targeting the same spectrum through different data.

Let $p(x, t)$ denote a transient density evolving under the Perron–Frobenius semigroup, so that $p(\cdot, t + \tau) = \mathcal{P}_\tau\,p(\cdot, t)$ for all $t$. Under the reversibility assumption, $\pi$ is invariant ($\mathcal{P}_\tau\,\pi = \pi$), and the ratio $r(x, t) = p(x, t)/\pi(x) - 1$ evolves under the Koopman operator by the identity $(\mathcal{P}_\tau \rho)(x)/\pi(x) = (\mathcal{K}_\tau(\rho/\pi))(x)$ derived. Substituting $\rho = p(\cdot, t)$ and using $\mathcal{K}_\tau\,\mathbf{1} = \mathbf{1}$, $$
r(x, t + \tau) = \frac{p(x, t + \tau)}{\pi(x)} - 1 = \frac{(\mathcal{P}_\tau\,p(\cdot, t))(x)}{\pi(x)} - 1
$$ $$
(\mathcal{K}_\tau (p(\cdot, t)/\pi))(x) - 1 = (\mathcal{K}_\tau (r(\cdot, t) + 1))(x) - 1 = (\mathcal{K}_\tau\,r(\cdot, t))(x),
$$ so the density-ratio function satisfies the exact recursion $r(x, t + \tau) = (\mathcal{K}_\tau\,r(\cdot, t))(x)$. The shift by $-1$ that defines $r$ removes the trivial stationary mode $\psi_1 \equiv 1$ from the evolution, which is essential: without the shift, the evolution is dominated by the stationary mode and the slow non-trivial modes appear as small corrections, whereas with the shift the evolution is purely in the non-trivial subspace and the spectral structure of the slow modes is exposed directly.

Expanding $r(\cdot, 0)$ in the eigenbasis of $\mathcal{K}_\tau$ as $r(\cdot, 0) = \sum_{i \ge 2} c_i\,\psi_i$ with modal coefficients $$
c_i = \langle r(\cdot, 0), \psi_i \rangle_{\pi} = \int_{\Omega}\!\left(\frac{p(x, 0)}{\pi(x)} - 1\right)\psi_i(x)\,\pi(x)\,dx = \mathbb{E}_{p(\cdot, 0)}[\psi_i] - \mathbb{E}_{\pi}[\psi_i],
$$ (the sum starts at $i = 2$ because $\psi_1 \equiv 1$ and $\langle r(\cdot, 0), 1 \rangle_{\pi} = \int(p - \pi) = 0$ by normalization of both distributions), and applying the recursion $t/\tau$ times gives $$
r(x, t) = (\mathcal{K}_\tau^{\,t/\tau}\,r(\cdot, 0))(x) = \sum_{i \ge 2} c_i\,\lambda_i^{\,t/\tau}\,\psi_i(x),
$$ which is the **density-relaxation expansion** in its operational form. Equivalently, the density itself admits the representation $$
p(x, t) = \pi(x)\!\left[1 + \sum_{i \ge 2} c_i\,\lambda_i^{\,t/\tau}\,\psi_i(x)\right],
$$ in which the temporal evolution is explicit: each mode decays as $\lambda_i^{\,t/\tau} = \exp(-t/\tau_i)$ with characteristic timescale $\tau_i = -\tau/\ln|\lambda_i|$, and the slowest non-trivial mode ($i = 2$) dominates the long-time approach to equilibrium.

The operational consequence for DR-DMD is that if the density ratio $r(x, t_k)$ is sampled at a set of reference points $\{x_j\}_{j=1}^{M}$ and iterations $t_k = k\tau$ for $k = 1, \ldots, N$, the resulting matrix $\mathbf{R} \in \mathbb{R}^{M \times N}$ with entries $R_{jk} = r(x_j, t_k)$ decomposes as $$
R_{jk} = \sum_{i \ge 2} c_i\,\lambda_i^{\,k}\,\psi_i(x_j) \;+\; \text{higher-order modes that decay faster,}
$$ where $\lambda_i^{\,k} = \exp(-k\tau/\tau_i)$ is the per-iteration decay factor. This is exactly the rank-reduced separable form that Exact DMD \textsuperscript{[26]} is designed to recover from snapshot matrices: columns of $\mathbf{R}$ are samples of the evolving density ratio at successive times; rows are samples at fixed reference points; the DMD algorithm extracts the eigenvalues $\lambda_i$ from the temporal decay and the eigenfunction values $\psi_i(x_j)$ from the spatial structure. The rank of $\mathbf{R}$ (in the noise-free limit and neglecting fast modes already relaxed by iteration $k = 1$) equals the number of slow modes present with non-zero modal coefficient $c_i$; the effective rank in the presence of sampling noise is set by the spectral gap and the signal-to-noise ratio.

Three **identifiability conditions** determine which modes can be recovered from a given WE experiment, and all three are non-negotiable inputs to the algorithm design.

(1) First, the mode must be **excited** in the initial distribution, i.e. $c_i = \mathbb{E}_{p(\cdot, 0)}[\psi_i] - \mathbb{E}_{\pi}[\psi_i] \ne 0$. A mode with $c_i = 0$ is invisible to DR-DMD regardless of how long the simulation runs, because the transient density is identical to what it would be in the absence of that mode. In practice, WE initial distributions localized on a single metastable basin excite the sub-leading mode $\psi_2$ (which distinguishes basins) strongly and higher modes weakly; this determines the order in which modes are recovered as simulation length increases. A deliberate initial-condition design in which walkers are seeded from multiple basins and multiple high-energy transition regions increases excitation of higher modes and is recommended where practicable; the bidirectional seeding strategy is designed to address this condition directly.

(2) Second, the mode must be **resolved** in the data — its decay timescale $\tau_i$ must be long enough compared to the sampling lag $\tau_{\mathrm{WE}}$ that $\lambda_i^{\,k}$ varies across the observed iteration range $k = 1, \ldots, N$ (otherwise the mode appears either as a constant — already equilibrated — or as noise — already decayed) and short enough compared to $N\tau_{\mathrm{WE}}$ that the total simulation captures a meaningful fraction of its decay (otherwise $\lambda_i^{\,k}$ is numerically indistinguishable from $1$ across all observed iterations). The window of observability is thus $$
    N\tau_{\mathrm{WE}} \gg \tau_i \gg \tau_{\mathrm{WE}},
    $$ which translates into a concrete constraint on the number of WE iterations required to recover a mode of given timescale.

(3) Third, the mode must be **separable from noise and from neighboring modes** — the spectral gap $|\lambda_i| - |\lambda_{i+1}|$ and the signal-to-noise ratio in the sampled ratio matrix must both be large enough that the DMD rank selection can distinguish the mode from its neighbors and from the noise floor.

The density-relaxation expansion makes the non-equilibrium character of WE operational rather than incidental. A WE simulation converged to equilibrium would have $p(x, t) = \pi(x)$ for all iterations, $r(x, t) = 0$ identically, and $\mathbf{R} = 0$; no spectral information could be extracted because there is nothing to decompose. It is precisely the **unconverged**, transient character of the WE walker ensemble — the fact that the ensemble is still relaxing toward $\pi$ and therefore the density ratio $r(\cdot, t_k)$ carries a non-zero signal in the sub-leading eigenbasis — that makes DR-DMD possible.

This is the opposite dependency from equilibrium-based methods such as reversible MSMs, which require convergence to $\pi$ within each metastable basin to estimate transition probabilities; DR-DMD exploits the distance from equilibrium as signal. The raw WE iteration weights $w_i^{(k)}$, the RiteWeight stationary weights $\hat{w}_i^{\pi}$, and the kernel-evaluated density-ratio construction are the finite-data realization of this theoretical picture: $w_i^{(k)}$ estimates the transient distribution $p(\cdot, t_k)$ in the numerator, $\hat{w}_i^{\pi}$ estimates the stationary distribution $\pi$ in the denominator, and the kernel sum evaluates $r(x_j, t_k)$ at the walker reference set.

Two further theoretical notes complete the foundation. First, the expansion above assumes reversibility; under irreversible dynamics the spectrum of $\mathcal{K}_\tau$ may be complex and the eigenfunctions are not orthogonal in $L^2(\pi)$, which complicates the spectral-estimator convergence theory. *The pipeline assumes reversibility and uses the emergence of predominantly real eigenvalues from the DMD output as a consistency check* — eigenvalues with significant imaginary parts signal either convergence failure or a genuine departure from the reversibility assumption (for example, through non-equilibrium driving that has not equilibrated). Second, the expansion is pointwise in $x$ and assumes $p(x, 0)/\pi(x) \in L^2(\pi)$; initial distributions with support strictly smaller than $\pi$ (for example, point masses or deltas on low-volume sets) violate this assumption formally but are handled correctly in practice because the kernel evaluation smooths singular initial distributions, and the density-ratio shift $-1$ ensures the data matrix has zero row mean in the large-$N$ limit.

#### 2.6.1 The amplitude--temporal--spatial decomposition

The spectral expansion

$$
\mathbf{R}_{jk} \;=\; \sum_{i \geq 2} c_i\, \lambda_i^k\, \psi_i(x_j) + \text{noise}
$$

decomposes the signal into three structurally distinct pieces, each contributing differently to what DMD recovers:

-   **Amplitude** $c_i$: scalar set by the initial condition (excitation coefficient). Determines whether mode $i$ is present in the data at all.

-   **Temporal structure** $\lambda_i^k$: exponential decay set by the operator. Determines the eigenvalue that DMD reads.

-   **Spatial structure** $\psi_i(x_j)$: eigenfunction values set by the operator. Determines the eigenfunction that DMD reads.

DMD's algorithmic separation of these three pieces is the content of the Schmid-Tu construction: the SVD of the snapshot matrix separates rank (signal vs. noise), the linear regression on the temporal axis extracts $\lambda_i$, and the spatial structure $\psi_i$ is recovered from the right singular vectors of the truncated SVD. **The eigenvalue** $\lambda_i$ and the eigenfunction shape $\psi_i$ are properties of the operator that do not depend on $c_i$; the amplitude $c_i$ sets only whether the mode is detectable, not what its eigenvalue or eigenfunction look like when detected.

This decomposition is the structural principle that resolves several of the apparent paradoxes in pooling: pooling can attenuate $c_i^{\mathrm{pool}}$ without affecting $\lambda_i$ or $\psi_i$, so destructive interference reduces detectability but does not bias the recovered values; pooling reduces noise floor without affecting $\lambda_i$ or $\psi_i$, so pooling improves estimation precision; the only failure mode of pooling is amplitude attenuation falling below the noise floor. The phrase "DMD reads the changes between iterations" captures this correctly: DMD reads the temporal ratio $\lambda_i = \mathbf{R}_{j, k+1}/\mathbf{R}_{jk}$ in each spatial mode, and that ratio is invariant under multiplication of the signal by any constant amplitude factor, attenuated or not.

### 2.6a The Density Ratio as a Natural Observable of the Koopman Operator

The density ratio $r(x, t) = p(x, t)/\pi(x) - 1$ is not merely a convenient normalisation — it is the natural observable of the Koopman operator under the reversibility and ergodicity assumptions. Three facts support this identification, completing the case developed for why the density-relaxation expansion is the correct operational object for a spectral estimator.

**Linearity under** $\mathcal{K}_\tau$. The recursion derived,

$$
r(x, t + \tau) \;=\; (\mathcal{K}_\tau\, r(\cdot, t))(x),
$$

shows that $r$ evolves under the Koopman operator directly, with no residual shift or inhomogeneous term. This is not true of $p(x, t)$ itself, which evolves under the Perron–Frobenius operator with the stationary distribution $\pi$ as a fixed point; it is also not true of $\log p$ or $p - \pi$, which satisfy nonlinear or inhomogeneous evolution equations. The ratio $p/\pi - 1$ is the unique pointwise nonlinear transformation of $p$ (up to affine rescaling) that linearises the evolution under $\mathcal{K}_\tau$, which is the operational content of the similarity-transform identity.

$L^2(\pi)$ geometry. The $\pi$-weighted inner product $\langle f, g\rangle_\pi = \int f(x) g(x) \pi(x)\,dx$ is the inner product under which $\mathcal{K}_\tau$ is self-adjoint, and the eigenfunctions $\{\psi_i\}$ form an orthonormal basis in this inner product. The density ratio satisfies

$$
\langle r(\cdot, t), 1\rangle_\pi \;=\; \int \left(\frac{p(x, t)}{\pi(x)} - 1\right) \pi(x)\,dx \;=\; \int p(x, t)\,dx - \int \pi(x)\,dx \;=\; 1 - 1 \;=\; 0,
$$

i.e. it is orthogonal to the stationary mode in $L^2(\pi)$ by construction. This means the density ratio lives entirely in the non-trivial subspace $\mathrm{span}(\psi_2, \psi_3, \ldots)$, and its expansion $r(x, t) = \sum_{i \ge 2} c_i\,\lambda_i^{t/\tau}\,\psi_i(x)$ is the spectral expansion of a vector in that subspace. The subtraction of $1$ thus aligns the density ratio with the subspace the DMD estimator is designed to recover.

**Kullback–Leibler divergence as a Lyapunov function.** The distance from equilibrium, quantified by the Kullback–Leibler divergence

$$
D_{\mathrm{KL}}(p(\cdot, t)\,\Vert\,\pi) \;=\; \int_\Omega p(x, t)\,\log\frac{p(x, t)}{\pi(x)}\,dx,
$$

is a monotonically decreasing function of $t$ under reversible dynamics, with $D_{\mathrm{KL}} \to 0$ as $t \to \infty$ — the standard $H$-theorem for the Fokker–Planck semigroup, strengthened to the quantitative de Bruijn form. Expanding to second order in the small-amplitude regime where $|r(x, t)| \ll 1$ uniformly,

$$
D_{\mathrm{KL}}(p(\cdot, t)\,\Vert\,\pi) \;\approx\; \tfrac{1}{2}\int r(x, t)^2\,\pi(x)\,dx \;=\; \tfrac{1}{2}\,\lVert r(\cdot, t)\rVert_{L^2(\pi)}^2 \;=\; \tfrac{1}{2}\sum_{i \ge 2} c_i^2\,\lambda_i^{2t/\tau},
$$

so the KL divergence is, to leading order, the squared $L^2(\pi)$ norm of the density ratio, and its decay is dominated by the slowest mode: $D_{\mathrm{KL}}(t) \sim c_2^2\, \lambda_2^{2t/\tau}/2$ as $t \to \infty$. This identification — the density ratio as the natural small-amplitude variable around equilibrium, with $\lVert r\rVert_{L^2(\pi)}^2$ as the quadratic approximation of KL — grounds the burn-in criterion and noise-floor analysis: the empirical $D_{\mathrm{KL}}(p(\cdot, t_k)\,\Vert\,\hat\pi)$ curve computed from kernel-evaluated densities is used to identify the iteration at which the linear-response regime is entered, and $\lambda_2$ extracted by DR-DMD should match the asymptotic slope

$$
\ln|\lambda_2| \;=\; \lim_{k \to \infty}\, \frac{1}{2k}\,\ln D_{\mathrm{KL}}(t_k)
$$

as an independent consistency check. These three facts — linearity of $r$ under $\mathcal{K}_\tau$, $\pi$-orthogonality of $r$ to the stationary mode, and the quadratic identification of the KL Lyapunov functional with the $L^2(\pi)$ norm of $r$ — identify the density ratio as the object that a spectral estimator should consume under the reversibility assumption. The kernel-evaluated density-ratio matrix is the finite-data realisation of this object, constructed via pointwise kernel sums at walker reference positions rather than via grid-based KDE; the $\pi$-weighting inside the Schmid projected DMD is the finite-sample implementation of the $L^2(\pi)$ inner product; and the three-way cross-validation uses precisely the $L^2(\pi)$ inner product to compute eigenfunction overlaps between the DR-DMD, kEDMD, and VAMPnet routes.

### 2.6b Noise structure of the density-ratio estimator

The signal-to-noise analysis underlying every subsequent claim about mode recoverability rests on a specific property of the density-ratio construction that is worth stating explicitly. The estimator

$$
\hat{\mathbf{R}}_{jk} \;=\; \frac{\hat{p}(x_j, t_k)}{\hat{\pi}(x_j)} - 1
$$

is a ratio of two kernel density estimates, each with its own sampling fluctuation. Let

$$
\hat{p}(x_j, t_k) \;=\; p(x_j, t_k) + \eta_p(x_j, t_k), \qquad \hat{\pi}(x_j) \;=\; \pi(x_j) + \eta_\pi(x_j)
$$

where $\eta_p$ and $\eta_\pi$ are the kernel density estimation errors with variance of standard kernel-density form,

$$
\mathrm{Var}[\eta_p(x_j, t_k)] \;\sim\; \frac{p(x_j, t_k)}{n_{\mathrm{eff}} h^d}, \qquad \mathrm{Var}[\eta_\pi(x_j)] \;\sim\; \frac{\pi(x_j)}{n_{\mathrm{eff}}^\pi h^d}
$$

where $n_{\mathrm{eff}}$ is the effective sample size of the kernel sum at $x_j$ at iteration $k$, $n_{\mathrm{eff}}^\pi$ the effective sample size of the pooled $\pi$-estimator, and $h$ the kernel bandwidth.

Propagating to the ratio via first-order Taylor expansion,

$$
\hat{\mathbf{R}}_{jk} - \mathbf{R}_{jk} \;\approx\; \frac{\eta_p(x_j, t_k)}{\pi(x_j)} - \frac{p(x_j, t_k)}{\pi(x_j)^2}\, \eta_\pi(x_j).
$$

In the regime where $p(x_j, t_k) \sim \pi(x_j)$ (everywhere except early iterations in regions far from the initial support), both terms scale as $1/\pi(x_j)$. Defining a single effective noise level $\sigma_{\mathrm{noise}}(x_j)$ as the standard deviation of $\hat{p}(x_j, t_k)$ at fixed $x_j$,

$$
\mathrm{SD}[\hat{\mathbf{R}}_{jk}] \;\approx\; \frac{\sigma_{\mathrm{noise}}(x_j)}{\hat{\pi}(x_j)}.
$$

**The noise in the density-ratio estimator scales inversely with the local stationary density.** This is the structural reason that high-free-energy regions (small $\pi$, near barriers and in tails) are noisy in $\hat{\mathbf{R}}$ even when the underlying $\hat{p}$ is well-sampled in absolute terms: dividing by a small denominator amplifies the relative noise. It is also the structural justification for the $\pi$-weighting in the Schmid projected DMD,

$$
\tilde{\mathbf{R}}_{jk} \;=\; \sqrt{\hat{\pi}(x_j)}\, \hat{\mathbf{R}}_{jk},
$$

which exactly cancels the leading $1/\hat{\pi}$ scaling of the noise standard deviation, leaving a residual noise of approximately uniform amplitude across the support:

$$
\mathrm{SD}[\tilde{\mathbf{R}}_{jk}] \;\approx\; \frac{\sigma_{\mathrm{noise}}(x_j)}{\sqrt{\hat{\pi}(x_j)}}.
$$

The $\pi$-weighting trades a mild bias near boundaries (where the $\sqrt{\hat{\pi}}$ factor still suppresses signal in addition to noise) for a much better-conditioned numerical problem in the SVD that DMD applies. The geometric interpretation is that $\sqrt{\hat{\pi}}$-weighting is the change-of-measure to the $L^2(\pi)$ inner product on $\mathbf{R}$ itself, matching the inner product in which $\mathcal{K}_\tau$ is self-adjoint under reversibility.

For the signal-to-noise analysis of subsequent sections, the salient fact is that the noise in $\hat{\mathbf{R}}_{jk}$ at reference point $x_j$ is

$$
\sigma_{\mathrm{noise}}^{\mathbf{R}}(x_j) \;\approx\; \frac{\sigma_{\mathrm{noise}}(x_j)}{\hat{\pi}(x_j)}
$$

with $\sigma_{\mathrm{noise}}(x_j)$ the kernel-density-estimator standard deviation of $\hat{p}(x_j, t_k)$. This is the denominator of the SNR criterion.

------------------------------------------------------------------------

## 3. Pipeline Architecture

```{=latex}
\begin{landscape}
\begin{figure}[p]
\centering
\resizebox{!}{0.85\textheight}{%
\begin{tikzpicture}[
  node distance=0.85cm and 1.0cm,
  >={Latex[length=2.2mm,width=1.8mm]},
  every node/.style={font=\footnotesize},
  block/.style={
    rectangle, draw=primaryborder, line width=0.5pt, rounded corners=2pt,
    minimum width=3.6cm, minimum height=0.95cm, align=center,
    fill=primaryfill, font=\footnotesize
  },
  core/.style={
    rectangle, draw=coreborder, line width=1.2pt, rounded corners=2pt,
    minimum width=4.0cm, minimum height=1.25cm, align=center,
    fill=corefill, font=\footnotesize
  },
  valblock/.style={
    rectangle, draw=valborder, line width=0.4pt, rounded corners=2pt,
    minimum width=3.6cm, minimum height=0.95cm, align=center,
    fill=valfill, font=\footnotesize, dash pattern=on 3pt off 1.5pt
  },
  output/.style={
    rectangle, draw=outborder, line width=0.4pt, rounded corners=2pt,
    minimum width=3.4cm, minimum height=0.9cm, align=center,
    fill=outfill, font=\footnotesize\itshape
  },
  rolelabel/.style={font=\scriptsize\itshape, inner sep=1pt},
  tierlabel/.style={font=\scriptsize\itshape, text=gray!65, align=right},
  primaryarrow/.style={draw, ->, line width=0.7pt, color=primaryborder},
  corearrow/.style={draw, ->, line width=1.0pt, color=coreborder},
  valarrow/.style={draw, ->, line width=0.4pt, color=valborder, dashed},
  outarrow/.style={draw, ->, line width=0.6pt, color=outborder}
]

% --- Tier 1: shared upstream pipeline ---
\node[block] (we) {WESTPA corpus \;(\S4)\\$\{w_i^{(k)},\, x_i^{(k)}\}$};
\node[block, right=of we] (tica) {IW-tICA \;(\S5)\\tIC coordinates $q(x)$};
\node[block, right=of tica] (rw) {Thermo.\ RiteWeight \;(\S6)\\stationary weights $\hat{w}_i^{\pi}$};

\draw[primaryarrow] (we) -- (tica);
\draw[primaryarrow] (tica) -- (rw);

% --- Tier 2: three parallel spectral estimators (independent, no cross-feeds) ---
\node[core, below=2.4cm of tica] (dmd)
  {density-ratio matrix \;(\S7)\\
   $R_{jk}=\hat{p}(x_j,t_k)/\hat{\pi}(x_j)-1$\\[1pt]
   \textbf{$\pi$-weighted Exact DMD} \;(\S9)\\
   $\Rightarrow \{\hat{\lambda}_i^{\mathrm{DR}},\,\hat{\psi}_i^{\mathrm{DR}}\}$};

\node[valblock, left=1.5cm of dmd] (kedmd)
  {kEDMD \;(\S10)\\
   $\{\hat{\lambda}_i^{\mathrm{kEDMD}},\,\hat{\psi}_i^{\mathrm{kEDMD}}\}$};

\node[valblock, right=1.5cm of dmd] (vamp)
  {RevVAMPnet \;(\S19)\\
   $\{\hat{\lambda}_i^{\mathrm{VAMP}},\,\hat{\psi}_i^{\mathrm{VAMP}}\}$};

\node[rolelabel, text=valborder, below=2pt of kedmd.south] {secondary};
\node[rolelabel, text=coreborder, below=2pt of dmd.south] {\textbf{primary}};
\node[rolelabel, text=valborder, below=2pt of vamp.south] {tertiary};

% Inputs to DR-DMD (primary path)
\draw[corearrow] (tica.south) -- (dmd.north);
\draw[corearrow] (rw.south) -- ++(0,-0.45) -| ([xshift=10mm]dmd.north);

% Inputs to kEDMD (tICA coords + RW endpoint weighting from corpus)
\draw[valarrow] (tica.south west) .. controls +(down:0.8cm) and +(up:1.0cm) .. (kedmd.north);
\draw[valarrow] (we.south) -- ++(0,-0.45) -| ([xshift=-6mm]kedmd.north);

% Inputs to RevVAMPnet (raw feature space directly from WE corpus)
\draw[valarrow] (we.south east) .. controls +(1.0cm,-1.6cm) and +(up:1.6cm) .. (vamp.north);

% --- Tier 3: spectral agreement (validates parallel estimators) ---
\node[block, below=1.6cm of dmd, minimum width=4.8cm] (agree)
  {Spectral agreement check \;(\S\S10,19)\\
   compare $\{\hat{\lambda}_i,\hat{\psi}_i\}$ across estimators\\
   eigenvalue spectra, mode overlap};

\draw[corearrow] (dmd.south) -- (agree.north);
\draw[valarrow] (kedmd.south) |- ([yshift=2mm]agree.west);
\draw[valarrow] (vamp.south) |- ([yshift=2mm]agree.east);

% --- Tier 4: three parallel eigenfunction consumers ---
\node[block, below=1.4cm of agree] (gp)
  {GP FES \;(\S13)\\
   dynamical kernel $k_{\mathrm{dyn}}$\\
   thermodynamic observables};

\node[output, left=1.0cm of gp, minimum width=3.6cm, minimum height=1.25cm] (kde)
  {Visualization KDE \;(\S14)\\
   $\hat{p}_{\mathrm{viz}}(\Psi,t_k)$\\
   diagnostic in $\Psi$-coords};

\node[block, right=1.0cm of gp] (pcca)
  {PCCA+ \;(\S12) + Kinetic RW \;(\S11)\\
   soft states, MFPTs, committors\\
   kinetic observables};

\draw[corearrow] (agree.south) -- (gp.north);
\draw[outarrow]  (agree.south) -| (kde.north);
\draw[corearrow] (agree.south) -| (pcca.north);

% Weighted-MSM validator: hangs off PCCA+ specifically (kinetic check)
\node[valblock, right=0.7cm of pcca, minimum width=2.8cm] (msm)
  {Weighted-MSM\\validator \;(\S12)\\CK test, det.\ balance};
\draw[valarrow] (pcca.east) -- (msm.west);

% --- Tier 5: combined ensemble observables ---
\node[output, below=1.4cm of gp, minimum width=4.6cm] (obs)
  {Ensemble observables \;(\S17)\\
   thermodynamic + kinetic};

\draw[outarrow] (gp.south) -- (obs.north);
\draw[outarrow] (pcca.south) |- ([yshift=2mm]obs.east);

% --- Tier annotations on left margin ---
\coordinate (leftmargin) at ([xshift=-0.4cm]kedmd.west);
\node[tierlabel, anchor=east] at (leftmargin |- we.center)
  {\begin{tabular}{r}shared\\inputs\end{tabular}};
\node[tierlabel, anchor=east] at (leftmargin |- dmd.center)
  {\begin{tabular}{r}parallel\\spectral\\estimators\end{tabular}};
\node[tierlabel, anchor=east] at (leftmargin |- agree.center)
  {\begin{tabular}{r}spectral\\agreement\end{tabular}};
\node[tierlabel, anchor=east] at (leftmargin |- gp.center)
  {\begin{tabular}{r}eigenfunction\\consumers\end{tabular}};
\node[tierlabel, anchor=east] at (leftmargin |- obs.center)
  {observables};

\end{tikzpicture}%
}
\caption{\textbf{Pipeline architecture.} Shared upstream stages — the WESTPA corpus is projected onto IW-tICA coordinates, and a thermodynamic RiteWeight pass recovers stationary weights $\hat{w}_i^\pi$. Three \emph{parallel} spectral estimators operate independently on the same upstream corpus. The \textbf{primary} route constructs the density-ratio matrix $R_{jk}=\hat{p}(x_j,t_k)/\hat{\pi}(x_j)-1$ and extracts transfer-operator eigenfunctions via $\pi$-weighted Exact DMD. The \textbf{secondary} route (kEDMD) operates in the same tICA coordinate space with RiteWeight endpoint weighting; the \textbf{tertiary} route (RevVAMPnet) operates directly on the full feature space. The three estimator outputs are compared at the spectral agreement check. Three sibling consumers of the validated eigenfunctions $\{\hat{\psi}_i\}$: The Gaussian-process free-energy surface with dynamical kernel $k_{\mathrm{dyn}}$ produces thermodynamic observables (free energies, equilibrium populations with calibrated uncertainty); PCCA+ with kinetic RiteWeight produces kinetic observables (soft macrostates, MFPTs, committors), with the weighted-MSM validator running CK and detailed-balance tests on the resulting kinetic construction. The visualization KDE produces $\hat{p}_{\mathrm{viz}}(\Psi,t_k)$ as a diagnostic of density evolution in eigenfunction coordinates. Section numbers refer to the chapter in which each stage is derived.}
\label{fig:pipeline}
\end{figure}
\end{landscape}
```

### 3.1 Overview and data flow

It is useful to fix the overall shape of the pipeline in a single figure. The pipeline is organized as one primary sequential path, one secondary parallel spectral cross-validation route, and one tertiary parallel validation route, all operating on the same WE simulation corpus. The primary path is the density-ratio DMD route that this document is about; the secondary and tertiary routes exist to provide independent estimates of the same spectral object against which the primary route is validated. The entire pipeline is post-hoc in the sense that no stage downstream of WESTPA modifies the trajectory data or influences the dynamics — every subsequent estimator is a deterministic or stochastic function of the walker ensemble emitted by the WESTPA run.

The primary path begins at a WESTPA simulation that produces an iteration-indexed walker ensemble with weights $\{(x_i^{(k,\text{start})}, x_i^{(k,\text{end})}, w_i^{(k)})\}_{i,k}$. The walker configurations are passed through a featurization $\chi: \Omega \to \mathbb{R}^D$ of moderate dimension $D$ and then through importance-weighted tICA using the raw WE iteration weights $w_i^{(k)}$ in the covariance estimators; the output is the tICA projection $q: \Omega \to \mathbb{R}^d$ with $d = 2$–$4$. Thermodynamic RiteWeight is then run in tIC space on the pooled walker ensemble to produce stationary weights $\hat{w}_i^\pi$ on the existing walkers. The kernel-evaluated density-ratio matrix $\mathbf{R} \in \mathbb{R}^{M \times N}$ is constructed on a walker reference set $\{x_j\}_{j=1}^M$ in tIC space, using a single isotropic Matérn-$5/2$ kernel of bandwidth $h$ for both numerator and denominator. A $\pi$-weighted Exact DMD with Schmid projection is applied to $\mathbf{R}$ to produce eigenvalues $\{\hat\lambda_i^{\mathrm{DR}}\}$ and continuous-valued eigenfunctions $\{\hat\psi_i^{\mathrm{DR}}\}$ represented as kernel sums over the reference set. State assignment via PCCA+ proceeds in eigenfunction coordinates and feeds two parallel branches: kinetic RiteWeight, which produces the pipeline's reportable kinetic observables (rate matrix, MFPTs, committors, reactive fluxes) directly from the RiteWeight-reweighted soft-membership flux; and the weighted-MSM validator, which constructs the reversible-MLE transition matrix $\hat T_{\mathrm{rev}}$ whose stationary distribution, implied timescales, and Chapman–Kolmogorov / detailed-balance residuals feed the validation suite as numerical comparators against the kinetic-RiteWeight values. The GP-FES estimator is fit in the DR-DMD eigenfunction coordinate system $\Psi = (\hat\psi_2^{\mathrm{DR}}, \hat\psi_3^{\mathrm{DR}}, \ldots)$ on the per-walker training observations $\hat F_i = -k_BT \ln \hat\pi(x_i)$ supplied by RiteWeight, using a dynamical kernel whose spectral decay is tied to the recovered eigenvalues; a visualization-only KDE movie in the same $\Psi$-coordinates is produced for qualitative interpretation.

The secondary route operates in parallel on the same WE corpus. The trajectory-pair representation $\{(x_i^{(k,\text{start})}, x_i^{(k,\text{end})})\}_{i,k}$ is aggregated across all walkers and all iterations, reweighted via the RiteWeight stationary weight $\hat{w}_i^\pi$ applied to the endpoint in the Wu–Nüske–Paul\textsuperscript{[22]} convention, and fed to kernel Extended DMD using the **same** tICA projection $q(\cdot)$ and the **same** isotropic Matérn-$5/2$ kernel $K_h$ as the primary route. The tertiary route trains a RevVAMPnet on the same trajectory pairs, but using the full feature vector $\chi(x)$ without tICA preprocessing, with the RiteWeight weights entering an importance-weighted VAMP-2 loss. All three routes target the same Koopman spectrum of the same lag-$\tau_{\mathrm{WE}}$ transfer operator, and their outputs are compared in the validation stage.

The tIC coordinate space serves three architectural roles simultaneously. First, it is the input space to thermodynamic RiteWeight's random-projection cell partition: the full feature space $\chi(\Omega) \subset \mathbb{R}^D$ is too high-dimensional for RiteWeight's piecewise-constant density estimate to be tractable, but the $d$-dimensional tIC subspace is small enough to admit fine partitions while preserving the slow-mode structure that RiteWeight's random-projection discretization is designed to exploit. Second, it is the coordinate space in which the Matérn-$5/2$ kernel bandwidth $h$ is chosen and in which all kernel evaluations are performed for both the density-ratio construction and the kEDMD Gram matrix — this shared kernel choice is the mathematical content of the statement that DR-DMD and kEDMD operate in the same reproducing kernel Hilbert space, and it is what allows direct comparison of the eigenfunction outputs of the two routes. Third, the tICA projection provides a linear-Koopman reference eigenvalue spectrum, in the sense that the tICA eigenvalues themselves are linear Koopman estimates in the Klus framework and therefore serve as a lower bound on the eigenvalue magnitude achievable by the nonlinear DR-DMD and kEDMD routes; any slow-mode eigenvalue the nonlinear routes return that is smaller in magnitude than the corresponding tICA eigenvalue flags an inconsistency that must be traced to either insufficient rank, over-regularization, or a misspecified bandwidth.

The kernel-evaluated density-ratio construction is explicit and is the computational core of the primary route. At each WE iteration $k$, $$
\hat{p}(x_j, t_k) = \sum_i w_i^{(k)}\, K_h\bigl(q(x_j) - q(x_i^{(k,\text{end})})\bigr), \qquad \hat{\pi}(x_j) = \sum_i \hat{w}_i^\pi\, K_h\bigl(q(x_j) - q(x_i)\bigr),
$$ with $\{x_j\}_{j=1}^M$ the walker reference set — typically the union of walker endpoints across all iterations, possibly subsampled. The density-ratio matrix entries are $R_{jk} = \hat{p}(x_j, t_k) / \hat{\pi}(x_j) - 1$, with the $-1$ shift removing the equilibrium offset so that the matrix decays to zero as the ensemble relaxes (cf.). No spatial grid is constructed at any stage: the data matrix is walker-reference $\times$ WE-iteration, not grid $\times$ iteration, so the method scales with the number of walkers rather than exponentially with the intrinsic dimension $d$ of the tIC space. The kernel is evaluated pointwise at walker reference positions, and the DMD eigenfunctions recovered from $\mathbf{R}$ are represented as kernel sums over the reference points — a Nadaraya-Watson-style interpolation that yields continuous-valued eigenfunctions evaluable at any point in tIC space without grid lookups.

The RKHS unification is the mathematical content of this shared-kernel choice and is what distinguishes the present pipeline from earlier DMD-on-WE efforts that used incompatible kernels at different stages. The isotropic Matérn-$5/2$ kernel $K_h$ on tIC space defines the Sobolev RKHS $\mathbb{H}_h = H^\sigma(\Omega)$ at order $\sigma = 5/2 + d_{\mathrm{tIC}}/2$ (Bold 2025\textsuperscript{[117]} / Hertel 2026\textsuperscript{[118]} unconditional Koopman-invariance certificate). The kernel-evaluated density-ratio construction produces point evaluations of functions in $\mathcal{H}_h$; the kEDMD Gram matrix is the Gram matrix of feature maps in the same $\mathcal{H}_h$; the GP-FES posterior mean in eigenfunction coordinates is a function in an RKHS induced by the dynamical kernel $k_{\mathrm{dyn}}$, which is itself a specific weighted combination of the eigenfunctions in $\mathcal{H}_h$. All spectral estimators and all downstream regressors therefore operate on the same function space, and cross-validation between routes compares objects in a common mathematical framework rather than objects defined by incompatible representations.

### 3.2 Walker weight types and their roles

Two primary walker-weight types --- together with one derived per-iteration normalization, which makes up the three weights enumerated in the Program Statement --- appear in the first half of the pipeline, each with a specific theoretical grounding and a specific set of pipeline stages at which it is consumed. The weights are not interchangeable, and using the wrong weight type at any stage corrupts the corresponding estimator; accordingly, the two types are tracked under distinct symbols throughout the document and are never mixed within a single estimator.

**Type 1 — Raw WE iteration weights** $w_i^{(k)}$. These are the per-walker weights produced by the WESTPA splitting-and-merging algorithm at iteration $k$, normalized so that $\sum_i w_i^{(k)} = 1$ for each iteration by WE construction. Under the Zhang–Jasnow–Zuckerman\textsuperscript{[2]} unbiasedness theorem, the raw weights provide an unbiased estimator of the transient distribution at iteration $k$: for any bounded test function $f$, $$
\mathbb{E}\left[\sum_i w_i^{(k)}\, f(x_i^{(k)})\right] = \int_{\Omega} f(x)\, p(x, t_k)\, dx.
$$ Raw WE weights enter four pipeline estimators. First, they enter the importance-weighted tICA covariance matrices pooled across all iterations $k$ — this is the mechanism by which tICA accounts for the non-equilibrium character of the WE corpus without requiring WE to have converged to $\pi$. Second, they enter the density-ratio numerator at each iteration $k$, where $\hat{p}(x_j, t_k) = \sum_i w_i^{(k)} K_h(q(x_j) - q(x_i^{(k,\text{end})}))$ is the kernel-smoothed estimator of $p(\cdot, t_k)$ evaluated at reference point $x_j$. Third, their per-iteration normalizations $\tilde{w}_i^{(k)} = w_i^{(k)} / \sum_j w_j^{(k)}$ (trivially equal to $w_i^{(k)}$ under WE normalization, but logically distinguished) drive the visualization-only KDE in recovered eigenfunction coordinates. Fourth — and only in specific contexts — they can enter the walker reference set construction via weighted subsampling if a reference-set reduction is performed. Raw WE weights are not valid estimators of $\pi$; using them where a stationary estimator is required is a category error that produces biased results.

**Type 2 — RiteWeight stationary weights** $\hat{w}_i^\pi$. These are produced by the thermodynamic RiteWeight procedure run in the tIC coordinate space on the aggregate walker ensemble. Theoretically, $\hat{w}_i^\pi \approx \pi(x_i) / q_{\mathrm{agg}}(x_i)$, where $q_{\mathrm{agg}}(x)$ is the aggregate WE sampling distribution across the whole corpus pooled across all iterations (not to be confused with the tICA projection, which shares notation but is always written with an argument $q(\cdot)$). After RiteWeight correction, the weighted walker ensemble $\{(x_i, \hat{w}_i^\pi)\}$ is an unbiased-in-expectation estimator of $\pi$-weighted integrals: $$
\mathbb{E}[\sum_i \hat{w}_i^\pi f(x_i)] \approx \int f(x) \pi(x)\,dx
$$. RiteWeight weights enter seven pipeline stages: First, the density-ratio denominator, $$
\hat\pi(x_j) = \sum_i \hat{w}_i^\pi K_h(q(x_j) - q(x_i))
$$. Second, the kEDMD trajectory-pair endpoint reweighting, with the Wu–Nüske–Paul\textsuperscript{[22]} convention that $\hat{w}_i^\pi$ is applied to the endpoint of the $i$-th trajectory pair. Third, the $\pi$-weighted SVD projection inside the Schmid projected Exact DMD, which weights each row of the density-ratio matrix by $\sqrt{\hat\pi(x_j)}$ before SVD. Fourth, the $L^2(\pi)$ orthonormality checks on recovered eigenfunctions, $\langle \hat\psi_i, \hat\psi_j \rangle_\pi \approx \delta_{ij}$. Fifth, the weighted MSM count matrix. Sixth, the GP-FES training observations, where $\hat{w}_i^\pi$ weights the contribution of each walker to the negative log-density estimator that serves as GP training target. Seventh, the importance-weighted VAMP-2 loss in the tertiary RevVAMPnet route. RiteWeight weights are not valid estimators of the transient distribution $p(x, t_k)$; using them in the density-ratio numerator or in per-iteration visualization KDE is a category error.

### 3.3 The two RiteWeight configurations

RiteWeight is run in two distinct configurations at different pipeline stages, with different inputs, different outputs, and different physical interpretations. The two configurations share the same underlying algorithmic machinery — a random-projection cell partition of the walker cloud in a low-dimensional coordinate space, followed by a linear system on these cells that produces corrected walker weights — but they differ in whether a state decomposition is provided as input and in what the output weights estimate. It is important to state at the outset that neither configuration consumes, requires, or reweights walker genealogy or parent/child lineage; RiteWeight is a stationary-reweighting estimator that operates on the walker cloud as a set of configurations with assigned iteration-weights, and it is not history-augmented in any algorithmic sense.

**Configuration A — Thermodynamic RiteWeight.** Run once, after tICA has produced the tIC coordinate space, on the full walker ensemble pooled across all iterations. No state definition is required. The input is the pooled walker set $\{x_i\}$ in tIC coordinates together with the raw WE iteration weights $w_i^{(k)}$ and the random-projection cell partition of tIC space. The algorithm constructs a piecewise-constant density-ratio estimator on these cells by solving a linear system whose unknowns are per-cell reweighting factors, and applies the resulting factors to the walker weights to yield per-walker stationary weights $\hat{w}_i^\pi$. The output weights satisfy the normalization $\sum_i \hat{w}_i^\pi = 1$ and the equilibrium expectation identity $$
\mathbb{E}[\sum_i \hat{w}_i^\pi f(x_i)] \approx \int f(x) \pi(x)\,dx
$$ for bounded $f$. The primary consumers are the density-ratio denominator $$
\hat\pi(x_j) = \sum_i \hat{w}_i^\pi K_h(q(x_j) - q(x_i))
$$, the kEDMD endpoint reweighting, and the $\pi$-weighted inner products inside the Schmid-projected DMD. Thermodynamic RiteWeight does not produce kinetic observables such as mean first-passage times or state-to-state rates; it produces only a stationary reweighting of the walker cloud.

**Configuration B — Kinetic RiteWeight.** Run after state definitions are available, either in tICA coordinates (before DR-DMD, as a first-pass estimate) or in DMD eigenfunction coordinates (after DR-DMD, as the primary estimate); both are constructed and compared as a consistency check. The input is a set of $n$ soft state assignments $\{\chi_s(x)\}_{s=1}^n$, typically from PCCA+ applied to the DMD eigenfunctions, together with the raw WE trajectory data. The output is state-to-state mean first-passage times (MFPTs), state-wise transition rates, and committor functions — kinetic observables not available from the thermodynamic configuration. Kinetic RiteWeight is run last in the pipeline and consumes the stationary DR-DMD eigenfunctions together with the raw WE trajectory data to produce committor and rate estimates. In the kinetic configuration, the RiteWeight linear system's unknowns are the state-to-state transition rates and the committor values on each cell, rather than the per-walker stationary reweighting factors of the thermodynamic configuration.

------------------------------------------------------------------------

## 4. Stage 1 — Weighted Ensemble Sampling (WESTPA)

### 4.1 The weighted ensemble strategy

Weighted ensemble (WE) sampling is a rigorous statistical-mechanics resampling algorithm that enhances the sampling of rare events and of regions of configuration space that are not favorably populated under unbiased dynamics, **without modifying the underlying equations of motion**. It achieves this by maintaining a population of trajectories ("walkers") that run in parallel, each carrying an explicit probability weight, and periodically resampling the walker ensemble according to a progress coordinate (PC) to maintain coverage across the PC range while preserving the correct statistical weights of configurations. No bias potential is ever introduced: the walkers propagate under exact, unbiased molecular dynamics; the resampling step only adjusts the partitioning of probability mass across walker-delta measures and never alters the force field or the equations of motion.

The foundational property proved by Zhang, Jasnow and Zuckerman\textsuperscript{[2]} is that the raw WE walker weights $w_i^{(k)}$ at iteration $k$ are unbiased estimators of the transient distribution $p(x, t_k)$ of the underlying unbiased dynamics at physical time $t_k = k\tau_{\mathrm{WE}}$, where $\tau_{\mathrm{WE}}$ is the WE iteration lag time. This unbiasedness holds for each iteration **individually**, regardless of whether the WE run has reached a steady state, and is the theoretical foundation on which all downstream stages of the pipeline rest. The point deserves emphasis because it is structurally different from the assumptions underlying conventional time-average MD: WE does not require ergodic convergence of any single walker, and it does not require the walker cloud as a whole to have converged to $\pi$; it requires only that the dynamics within each WE propagation segment are exact unbiased MD, and that the resampling rule conserve total weight in expectation.

The WE protocol proceeds iteratively. At iteration $k = 0$, the PC range of interest is partitioned into bins $\{B_\alpha\}_{\alpha=1}^{N_{\mathrm{bins}}}$, and an initial population of walkers is seeded across the bins with assigned weights $w_i^{(0)}$ summing to unity (or to the user-specified total probability when combining ensembles). Each walker is a physical configuration $x_i^{(k)} \in \mathcal{X}$ carrying a weight $w_i^{(k)} \in [0, 1]$. The interpretation is that the walker represents a Dirac delta at $x_i^{(k)}$ in the transient distribution with measure $w_i^{(k)}$, so that $$
p(x, t_k) \;\approx\; \sum_i w_i^{(k)}\, \delta(x - x_i^{(k)})
$$ in the weak sense, and expectations are computed by $\mathbb{E}_{p(\cdot,t_k)}[f] \approx \sum_i w_i^{(k)} f(x_i^{(k)})$ for any bounded test function $f$.

Each WE iteration consists of three stages executed in sequence: **propagation** of each walker independently under the unbiased dynamics for a fixed lag time $\tau_{\mathrm{WE}}$ (typically 10–100 ps of MD integration, depending on the system and the slowest timescales of interest); **binning** of the resulting endpoint configurations according to their PC values into the user-specified set of bins; and **resampling** within each bin via the Huber-Kim splitting and merging protocol to drive the walker count in each occupied bin toward a target value (typically 4–8 walkers per bin). The resampling conserves total probability exactly and preserves the unbiasedness of walker weights for representing $p(x, t_k)$: no walker configurations are modified, no biasing forces are applied, and no reweighting is imposed by hand. Only the partitioning of probability across walker-delta measures is adjusted, so that underpopulated bins receive additional walkers (via splitting of existing walkers in those bins) and overpopulated bins are thinned (via merging of low-weight walkers). Splitting and merging decisions are local to each bin and are independent of the dynamics in every other bin.

The output of an $N$-iteration WE run is the set of walker trajectories $\{x_i^{(k)}, w_i^{(k)}\}$ for $k = 0, 1, \ldots, N-1$, together with the per-walker parent-child lineage structure recording which walker at iteration $k$ produced which walker(s) at iteration $k+1$. The lineage structure permits reconstruction of trajectory pairs $(x_i^{(k, \mathrm{start})}, x_i^{(k, \mathrm{end})})$ with the physical lag time $\tau_{\mathrm{WE}}$ for each walker at each iteration, and these pairs — aggregated across all walkers and all iterations and reweighted by the RiteWeight stationary weights — supply the input to the kEDMD secondary spectral route and the RevVAMPnet tertiary route. The transient distribution snapshots $p(x, t_k)$ indexed by iteration $k$, taken in their entirety across all iterations, supply the input to the density-ratio primary spectral route (–).

The WE iteration count $N$ is therefore a dual-purpose parameter. In the conventional WE analysis, $N$ sets the physical horizon of the simulation; in the present pipeline, it additionally sets **the number of columns of the density-ratio matrix** $\mathbf{R} \in \mathbb{R}^{M \times N}$, which is the data matrix consumed by the DMD stage. Each WE iteration produces exactly one measurement of the evolving density $p(x, t_k)$, and the iteration count $N$ directly determines the number of snapshots available to the DMD algorithm. Sizing $N$ correctly — specifically, ensuring $N$ exceeds the longest implied timescale of interest by several multiples and provides enough columns for the SVD rank truncation to resolve $n-1$ slow modes — is therefore a pipeline design decision that couples the WE runtime to the spectral resolution of the downstream DMD. Typical values of $N$ in this pipeline are 200–2000 iterations, with the lower end adequate for systems with a single slow timescale and two basins, and the upper end required for multi-basin systems with multiple slow timescales and fine spectral-gap resolution.

WESTPA 2.0\textsuperscript{[8]} is the reference open-source implementation of the WE algorithm used in this pipeline. The software exposes the propagation, binning, and resampling stages as modular components and supports a wide range of PC definitions, adaptive binning schemes, and parallel execution backends. The pipeline makes no algorithmic modifications to WESTPA itself; all novelty lies in the downstream analysis stages that consume the WESTPA output.

### 4.2 Resampling rules and probability conservation

The Huber-Kim splitting and merging protocol\textsuperscript{[1]} is the standard resampling rule used by WESTPA and by the pipeline. Given a target walker count $n_{\mathrm{target}}$ per occupied bin and the current set of walkers $\{(x_i, w_i)\}$ occupying a given bin after propagation, the rule proceeds as follows.

**Splitting.** If the current walker count in the bin $n_{\mathrm{current}} < n_{\mathrm{target}}$, the highest-weight walker in the bin is selected and split into $m$ copies, where $m$ is chosen to bring $n_{\mathrm{current}}$ closer to $n_{\mathrm{target}}$. Each of the $m$ copies retains the configuration $x_i$ and is assigned weight $w_i / m$. The split is exact: total weight in the bin is conserved because $$
m \cdot \frac{w_i}{m} + \sum_{j \neq i} w_j \;=\; w_i + \sum_{j \neq i} w_j \;=\; \sum_j w_j^{\mathrm{before}},
$$ and each split copy evolves independently from that iteration forward. If $n_{\mathrm{current}} + (m-1)$ still falls short of $n_{\mathrm{target}}$, the next-highest-weight walker is split next, and the procedure repeats until the target count is reached. Splitting never modifies a configuration and introduces no bias; it merely refines the walker-delta representation of $p(x, t_k)$ in the bin from $w_i \delta(x - x_i)$ into $m$ co-located deltas each of weight $w_i / m$, which represent the same transient distribution.

**Merging.** If $n_{\mathrm{current}} > n_{\mathrm{target}}$, the two lowest-weight walkers in the bin are selected for merging. Denote their weights and configurations by $(w_a, x_a)$ and $(w_b, x_b)$. One of the two configurations is chosen at random with probability proportional to its weight: $$
P(\text{keep } x_a) \;=\; \frac{w_a}{w_a + w_b}, \qquad P(\text{keep } x_b) \;=\; \frac{w_b}{w_a + w_b}.
$$ The retained configuration is assigned the summed weight $w_a + w_b$; the other configuration is discarded. The procedure repeats until $n_{\mathrm{current}} = n_{\mathrm{target}}$. Merging conserves total bin weight exactly. It also preserves the expectation of any bounded observable $f$: $$
\mathbb{E}\bigl[w_{\mathrm{merged}}\, f(x_{\mathrm{merged}})\bigr] \;=\; \frac{w_a}{w_a + w_b}(w_a + w_b)\, f(x_a) + \frac{w_b}{w_a + w_b}(w_a + w_b)\, f(x_b) \;=\; w_a f(x_a) + w_b f(x_b),
$$ which is the unmerged contribution. Merging introduces variance (one configuration is lost at each merge), but no bias; the unbiasedness property of WE walker weights is preserved exactly across splitting and merging.

The formal statement \textsuperscript{[2]} is that for any bounded observable $f$ and any iteration $k$, $$
\mathbb{E}\!\left[\sum_i w_i^{(k)}\, f(x_i^{(k)})\right] \;=\; \int_\Omega f(x)\, p(x, t_k)\, dx,
$$ where the outer expectation is over the randomness of the WE resampling and the underlying dynamics. This identity holds iteration by iteration; it does not require the WE run to have reached a stationary distribution, and it is the foundation on which the pipeline's use of raw WE weights $w_i^{(k)}$ in the density-ratio numerator $\hat{p}(x_j, t_k) = \sum_i w_i^{(k)} K_h(q(x_j) - q(x_i^{(k, \mathrm{end})}))$ and in the importance-weighted tICA covariance estimators rests. No equilibrium-convergence prerequisite is invoked for these uses.

### 4.3 The WE corpus as an importance-weighted sample

After $N$ WE iterations the simulation produces a **segment corpus** $$
\mathcal{S} \;=\; \bigl\{(x_i^{(k, \mathrm{start})},\; x_i^{(k, \mathrm{end})},\; w_i^{(k)})\bigr\}_{i, k},
$$ where $x_i^{(k, \mathrm{start})}$ and $x_i^{(k, \mathrm{end})}$ are the initial and final configurations of the $i$-th walker's propagation at iteration $k$, and $w_i^{(k)}$ is its per-iteration weight. The interpretation of the weights follows directly from the Zhang-Jasnow-Zuckerman theorem: for any measurable set $A \subset \Omega$ and for any bounded observable $f$, the per-iteration weighted sum $$
\sum_{i: x_i^{(k, \mathrm{start})} \in A} w_i^{(k)} \;\approx\; \mathrm{Pr}\bigl[X_{t_k} \in A\bigr], \qquad \sum_i w_i^{(k)}\, f(x_i^{(k, \mathrm{start})}) \;\approx\; \mathbb{E}_{p(\cdot, t_k)}[f],
$$ provides an unbiased estimator of the transient probability mass of $A$ and of the transient expectation of $f$ at physical time $t_k$. The identity applies equally to instantaneous and time-lagged observables: because the dynamics within each segment are exact unbiased MD, the pair $(x_i^{(k, \mathrm{start})}, x_i^{(k, \mathrm{end})})$ is a valid sample of the joint process at lag $\tau_{\mathrm{WE}}$ weighted by $w_i^{(k)}$, and $$
\sum_i w_i^{(k)}\, f(x_i^{(k, \mathrm{start})})\, g(x_i^{(k, \mathrm{end})}) \;\approx\; \mathbb{E}\bigl[f(X_{t_k})\, g(X_{t_k + \tau_{\mathrm{WE}}})\bigr].
$$ This property — that $w_i^{(k)}$ correctly weights both instantaneous and time-lagged expectations of the transient ensemble — is what enables every downstream importance-weighted estimator in the pipeline.

The full WE corpus pooled across iterations is *not*, however, a sample from the equilibrium distribution $\pi$. It is a sample from an aggregate distribution $q_{\mathrm{agg}}(x)$ that reflects the history of WE splitting and merging and is generally non-equilibrium. Writing the pooled walker-delta density, $$
q_{\mathrm{agg}}(x) \;=\; \frac{1}{N} \sum_{k=0}^{N-1} p(x, t_k),
$$ the pooled density is the time-averaged transient distribution over the WE run. For a WE simulation that has reached steady state, $p(x, t_k) \to \pi_{\mathrm{WE\text{-}SS}}(x)$ and $q_{\mathrm{agg}}(x) \approx \pi_{\mathrm{WE\text{-}SS}}$, but the WE steady state $\pi_{\mathrm{WE\text{-}SS}}$ is biased by the PC binning and by the splitting/merging dynamics, and generally differs from the unbiased Boltzmann $\pi$. For a WE run that has not reached steady state, $q_{\mathrm{agg}}(x)$ is a genuine time average of transient distributions and does not approximate any stationary distribution at all.

The pipeline treats this pooled distribution correctly in two complementary ways. First, any quantity that is *native* to the transient ensemble — individual $p(x, t_k)$ snapshots, time-averaged features of the walker trajectories, non-equilibrium correlations — is computed directly from the raw WE weights $w_i^{(k)}$, exploiting the Zhang-Jasnow-Zuckerman unbiasedness per iteration. The density-ratio numerator, the importance-weighted tICA covariance estimators, and the visualization KDE per iteration all fall in this category. Second, any quantity that requires an equilibrium reference — stationary expectations, equilibrium free energies, the density-ratio denominator, kEDMD endpoint reweighting — is computed from the RiteWeight-corrected stationary weights $\hat{w}_i^\pi$, which estimate $\pi(x_i) / q_{\mathrm{agg}}(x_i)$ and thereby convert the pooled walker set into an importance-weighted sample from $\pi$. Specifically, the identity $$
\sum_i \hat{w}_i^\pi\, f(x_i) \;\approx\; \mathbb{E}_\pi[f]
$$ holds to RiteWeight accuracy for any bounded $f$, so $\hat{w}_i^\pi$ supplies the equilibrium oracle wherever one is required.

This two-weight discipline — raw $w_i^{(k)}$ for transient-native quantities, RiteWeight $\hat{w}_i^\pi$ for equilibrium-reference quantities — is a structural feature of the pipeline and is enforced throughout the methods document. Mixing the two within a single estimator leads to bias (using raw WE weights where an equilibrium reference is required undercounts rare basins whose walkers have not been up-weighted to their correct $\pi$-share; using RiteWeight weights where a transient snapshot is required erases the time-resolution that the iteration index provides and collapses the density-ratio matrix to zero). The weight-type mapping is the authoritative reference for which weight enters which stage, and (symbol table) lists the correct weight type for every equation in the document.

### 4.4 The non-equilibrium history of individual walkers

The Zhang-Jasnow-Zuckerman per-iteration unbiasedness is a statement about the walker *ensemble* at iteration $k$, not about the history of any individual walker. It does not imply that any given walker trajectory is a sample from equilibrium dynamics; it implies only that the weighted sum over all walkers at iteration $k$ correctly reproduces the transient distribution at time $t_k$. The history of an individual walker is, in general, strongly non-equilibrium: a walker that was born (via splitting) at WE iteration $k' < k$ began from a configuration selected by the resampling rule — typically a position near a bin boundary or near its parent's location — which is not drawn from the equilibrium conditional distribution within that bin, and the walker's subsequent trajectory carries the history of its initial condition through however many intervening resampling events have occurred. Two walkers that are both located in the same metastable basin at iteration $k$ may have entirely different ancestries, and their conditional-on-basin distributions of configuration and of configuration derivatives are not the equilibrium conditional distributions.

The key observation is that every estimator in this pipeline operates on quantities that the Zhang-Jasnow-Zuckerman theorem directly certifies — either (i) the raw per-iteration weighted expectations $\sum_i w_i^{(k)} f(x_i^{(k)})$, which are unbiased estimators of $\mathbb{E}_{p(\cdot, t_k)}[f]$ irrespective of individual walker histories, or (ii) the RiteWeight-corrected weighted expectations $\sum_i \hat{w}_i^\pi f(x_i)$, which are unbiased estimators of $\mathbb{E}_\pi[f]$ after the RiteWeight reweighting has corrected the pooled corpus back to $\pi$. No estimator in the pipeline requires the equilibrium conditional distribution within a bin, and no estimator requires a walker's history to be an equilibrium trajectory. The non-equilibrium history of individual walkers is irrelevant to the estimators as constructed, because the estimators are designed to consume only the weighted-sum-over-walkers quantities that the theorem certifies.

Concretely, the density-ratio numerator $$
\hat{p}(x_j, t_k) = \sum_i w_i^{(k)} K_h(q(x_j) - q(x_i^{(k, \mathrm{end})}))
$$ is a kernel-smoothed estimator of $p(\cdot, t_k)$ evaluated at reference point $x_j$; its unbiasedness follows from the Zhang-Jasnow-Zuckerman theorem applied to the test function $f(y) = K_h(q(x_j) - q(y))$, which is bounded and does not depend on walker history. The density-ratio denominator $$
\hat{\pi}(x_j) = \sum_i \hat{w}_i^\pi K_h(q(x_j) - q(x_i))
$$ is a kernel-smoothed estimator of $\pi$ at $x_j$; its unbiasedness follows from the RiteWeight identity $$
\sum_i \hat{w}_i^\pi f(x_i) \approx \mathbb{E}_\pi[f]
$$ applied to the same test function. The importance-weighted tICA covariance estimators are per-iteration Zhang-Jasnow-Zuckerman estimators of $\mathbb{E}_{p(\cdot, t_k)}[\chi \chi^\top]$ and $\mathbb{E}[\chi(X_{t_k}) \chi(X_{t_k + \tau})^\top]$ pooled across iterations. The kEDMD endpoint reweighting uses $\hat{w}_i^\pi$ on the endpoint of each trajectory pair and is a valid estimator of the reweighted Koopman covariance (Wu–Nüske–Paul\textsuperscript{[22]} Algorithm 3). In every case, the estimator consumes walker weights and walker configurations, and the Zhang-Jasnow-Zuckerman or RiteWeight unbiasedness is invoked to certify the estimator, with no appeal to individual walker histories.

The haMSM was designed to correct biases that arise when an estimator implicitly assumes equilibrium conditional distributions within bins. The pipeline avoids such estimators entirely, replacing them with either kernel-smoothed density estimators (which avoid bins altogether and use weighted sums that the Zhang-Jasnow-Zuckerman theorem directly certifies) or RiteWeight-reweighted trajectory-pair estimators (which carry the equilibrium-conditional correction as the RiteWeight weights themselves). The net effect is that the non-equilibrium history of individual walkers does not enter the correctness argument for any pipeline estimator, and no history-augmented state space is required. Reweighting the pooled corpus with RiteWeight produces stationary weights that, combined with the raw iteration weights for transient-native estimators, are sufficient to extract the transfer-operator spectrum and to construct the FES.

### 4.5 Equilibrium WE and the WEED driver

The pipeline is configured to run WESTPA in its **equilibrium (WEED) mode**: no recycling boundary condition is imposed, so the WE stationary distribution coincides with the unbiased Boltzmann $\pi$ rather than with a steady-state non-equilibrium flux. The WEED driver (weighted-ensemble equilibrium dynamics\textsuperscript{[3]}) is the variant of WE designed for this setting. Walkers evolve under unbiased dynamics, splitting and merging proceed according to the Huber-Kim protocol described, and no configurations are recycled back to a source state upon reaching a sink. The asymptotic distribution of the walker ensemble is therefore $\pi$, and the iteration-indexed walker set samples transient distributions $p(\cdot, t_k)$ that relax toward $\pi$ rather than toward a non-equilibrium steady state.

This configuration is the correct pairing for the density-ratio primary route because $R_{jk} = \hat{p}(x_j, t_k) / \hat{\pi}(x_j) - 1$ is defined relative to $\pi$, and the RiteWeight-corrected denominator $\hat{\pi}$ is an estimator of the true equilibrium density. A NESS-WE recycling rule would instead concentrate walkers along the reactive flux tube and leave the non-reactive portions of $\operatorname{supp}(\pi)$ underpopulated, producing a RiteWeight estimator of a different stationary distribution and biasing the density-ratio denominator. The kEDMD secondary route has the same requirement: its eigenvalue guarantees ($|\lambda_i| \leq 1$, real spectrum under reversibility) require the reweighted samples to represent $\pi$, not a NESS.

The WEED driver is compatible with bidirectional seeding: multiple initial basins are treated as independent walker populations evolved in parallel under the same equilibrium WE rules and then pooled. The pooled raw-weight distribution $q_{\mathrm{agg}}(x) = N^{-1} \sum_k p(x, t_k)$ converges to $\pi$ only in the long-time limit, so the two-stage workflow — tICA with raw WE weights first, then RiteWeight in tIC space — remains the mechanism by which equilibrium estimates are recovered at finite simulation time. In this sense, the WEED driver provides the correct asymptotic target but does not obviate the need for explicit reweighting; RiteWeight provides the finite-time correction that converts the pooled walker cloud into an importance-weighted sample from $\pi$ before any of the equilibrium-reference estimators are applied.

The practical consequence of the WEED+bidirectional-seeding+RiteWeight protocol is that the full pipeline is well-defined at any WE iteration count $N \geq N_{\min}$, where $N_{\min}$ is the minimum count required to resolve the longest implied timescale of interest (cf. on sizing $N$) and to achieve adequate RiteWeight effective sample size in each basin. The pipeline does not require the WE run to have converged in any strong sense; it requires only that the corpus cover $\operatorname{supp}(\pi)$ and that $N$ be large enough to resolve the slow modes in the iteration index. This is a substantial relaxation of the conventional MD requirement of ergodic convergence and is the operational reason why the pipeline is feasible for systems with timescales far beyond the reach of brute-force MD.

### 4.6 Progress-coordinate choice and bidirectional seeding

The progress coordinate choice is left to the analyst and reflects prior knowledge of the system. Any scalar or low-dimensional function of configuration $\phi: \mathcal{X} \to \mathbb{R}^{d_\phi}$ can serve as a PC; common choices include interatomic distances, RMSD from a reference structure, chemical or physical order parameters, radius of gyration, number of native contacts, or low-dimensional projections onto previously computed collective variables. For RNA systems specifically, practical choices include end-to-end distance $r_{ee}$, radius of gyration $R_g$, heavy-atom RMSD from a reference folded structure, number of native base pairs, or a low-dimensional projection from a pre-trained autoencoder. For protein folding, backbone RMSD, fraction of native contacts, or a one-dimensional projection onto a previously identified slow collective variable are typical. For ligand binding and unbinding, ligand-receptor distance, pocket exit coordinate, or a binding-pose RMSD are standard. The pipeline makes no assumption that the PC coincides with the slow manifold or with any Koopman eigenfunction; indeed, PC independence of the downstream estimators is engineered precisely to avoid conditioning the spectral output on the PC choice. The PC determines **sampling coverage** but not **analysis quality** — any unbiased progress coordinate will eventually produce a corpus from which the correct dynamics can be extracted, given sufficient iterations.

The only requirement the analysis places on the PC is that it drive the WE walkers into coverage of the configuration regions of interest, including all metastable basins of interest and the transition regions between them. A PC that is strongly correlated with the slowest mode of interest will achieve coverage quickly; a PC that is nearly orthogonal to the slowest mode will waste walkers on rapid intra-basin fluctuations and may fail to achieve basin-to-basin coverage within a feasible walker budget. Adaptive progress coordinates — in which the PC is updated between WE runs based on preliminary analysis of early iterations — can dramatically accelerate coverage of barrier regions. Although the PC choice does not affect the asymptotic spectral recovery, it does affect the **convergence rate** of both the WE coverage and the downstream RiteWeight and DMD stages. A PC that is weakly correlated with the slowest conformational transitions will require many more WE iterations to achieve coverage of $\operatorname{supp}(\pi)$ than a PC that is strongly correlated. The recovered DR-DMD eigenfunctions from an early analysis round can themselves serve as improved PCs for subsequent WE rounds, closing the loop between sampling and analysis.

**Bidirectional seeding** is recommended for systems with two or more metastable basins and is, in fact, the default protocol for this pipeline. Walkers are initialized in each basin separately, with equal initial weight per basin (e.g. $1/N$ per walker for an $N$-walker initial ensemble), and the WE protocol is run outward from each basin. The resulting two (or more) walker ensembles are combined via the WESTPA ensemble-combination machinery, and the pooled corpus covers both basins and the intervening transition region from the outset. This is particularly important for the present pipeline because the density-ratio DMD and kEDMD routes require coverage of the transition region to resolve the slow modes that interpolate between basins; a one-sided WE run that covers only a single basin will produce a density-ratio matrix with negligible variation in the transition region and will fail to resolve the slow modes regardless of how many WE iterations are run.

Bidirectional seeding does not, of course, represent the equilibrium distribution: the relative populations of the two basins in the initial ensemble are artificially set by the seeding ratio (typically 50/50 or 1/$n$ per basin for $n$ basins), not by the Boltzmann-correct populations, which are typically unknown in advance. The non-equilibrium initialization is corrected by the RiteWeight stage, which reweights the pooled walker corpus back to $\pi$ regardless of the initial basin populations. The importance-weighting machinery throughout the pipeline — raw WE weights for transient observables, RiteWeight weights for equilibrium observables — correctly absorbs this initialization convention, and the final FES output is invariant to the choice of seeding ratios, provided only that every basin of interest is seeded with a non-vanishing initial population so that coverage is achieved. The 50/50 default is also safe for the downstream modal-excitation condition: because the excitation coefficient $c_2(0)$ is proportional to the deviation of the seeding ratio from the equilibrium basin ratio, the 50/50 allocation gives $c_2(0) \ne 0$ precisely when the equilibrium is asymmetric — the generic case. An asymmetric allocation is activated as a remediation only when a preliminary modal-excitation audit flags a near-symmetric equilibrium.

For systems with more than two metastable basins (e.g. native, major misfolded, and unfolded states for a protein; several binding poses for a ligand), the strategy generalizes: seed walkers from all known or suspected basins, with initial weights reflecting equal allocation rather than equilibrium populations. The WE resampling rule, the RiteWeight reweighting, and the importance-weighting discipline in the spectral estimators handle the rest. In practice, identifying which basins to seed is itself a substantive prior step in the analysis, and missing a basin in the seed set will result in that basin being absent from the final FES — this is the coverage limitation common to all reweighting methods, and the architectural response is to expand the seed set and re-run WE when a missing basin is suspected.

### 4.7 Progress-coordinate independence and its structural limits

WE sampling uses a progress coordinate (PC) to drive splitting and merging; the choice of PC affects sampling efficiency and which regions of configuration space are explored. The pipeline's downstream estimators are designed to be as independent as possible of the PC choice, so that the quality of the recovered spectrum and FES is not dominated by the analyst's choice of a good PC. The degree of PC independence varies by pipeline stage, and the structural limits are made explicit here so that the reader can calibrate expectations and interpret diagnostic failures correctly.

Thermodynamic RiteWeight is exactly PC-independent in construction: its random-projection cell partition depends only on the aggregate walker ensemble and the tIC coordinate space, not on the PC that drove the WE splitting/merging. If two WE runs on the same system with different PCs produce walker ensembles that cover the same regions of configuration space, RiteWeight applied to the pooled tIC coordinates will return weights consistent with the same underlying $\pi$. The only residual PC dependence is through coverage: a PC that fails to drive walkers into a metastable basin results in absent walkers in that basin, which no reweighting procedure can recover. This "cannot reweight what was not sampled" constraint is a fundamental property shared by all reweighting methods, and it is not a failure of RiteWeight specifically — it is the mathematical content of the statement that probability cannot be manufactured from zero samples.

The density-ratio DMD primary route and the kEDMD secondary route are fully PC-independent in the same respect as RiteWeight: both operate in tIC space, and the tIC subspace is derived from the WE data using raw WE iteration weights in the importance-weighted covariance estimators. If the WE PC leads to severely incomplete coverage of the slow manifold — for example, if the PC is chosen perpendicular to the true slowest mode — the tICA projection will not resolve the slowest mode, and both DR-DMD and kEDMD will fail to recover it, because both operate in a subspace that does not contain the mode. This is a structural architectural property of the pipeline rather than a caveat: coverage of the tIC subspace is required for both the primary and the secondary spectral routes. In practice, tIC subspace coverage is verified by checking that the leading tICA eigenvalues stabilize as the WE simulation lengthens, that the dominant tICs display multi-modal distributions aligned with known or suspected metastable basins, and that the projection of walker density onto the leading tICs reveals a clear spectral gap; makes these checks concrete.

The tertiary RevVAMPnet route operates in the full feature space without tICA preprocessing and therefore does not share this tIC-subspace limitation. Disagreement between the two tIC-space routes (DR-DMD and kEDMD) and the full-feature VAMPnet route is diagnostic of tICA linearity or tIC coverage limitations: if VAMPnet finds a slow mode that is invisible to DR-DMD and kEDMD, the slow manifold has nonlinear structure that the tICA basis cannot capture with the current feature set. The architectural response in this case is to expand the feature set, revisit the tICA stage with the expanded features, and re-run the primary and secondary routes; the tertiary route thus serves as the detector for when this diagnostic intervention is required. The converse situation — DR-DMD and kEDMD agree on a slow mode that VAMPnet misses — is typically a sign that VAMPnet has under-trained or that the VAMPnet network architecture is too small to resolve the mode; the architectural response is to increase network capacity or training epochs, not to discard the tIC-space result. The companion Pipeline document's validation framework codifies these three-way agreement patterns into concrete decision rules.

The overall picture is that the pipeline is PC-independent at the thermodynamic layer (stationary weights, thermodynamic quantities integrated against $\pi$) and PC-dependent-through-coverage at the spectral and kinetic layers, with the coverage dependence made visible by the convergence diagnostics. A reader interpreting a pipeline output should read the coverage and convergence diagnostics first — they tell the reader whether the spectral claims of the pipeline are backed by the data or whether the downstream estimates are inheriting a sampling gap that no amount of sophisticated reweighting can close.

#### 4.8 What WESTPA produces

Stage 1 output is the **WE corpus** $\mathcal{C} = \{(x_i^{(k, \mathrm{start})}, x_i^{(k, \mathrm{end})}, w_i^{(k)}): i = 1, \ldots, n_k;\; k = 1, \ldots, N\}$, together with auxiliary data (walker trajectories, bin assignments, parent pointers, the WESTPA HDF5 archive) that are used for diagnostics but do not enter the spectral pipeline proper. Three distinct objects are extracted downstream:

1.  **Trajectory-pair corpus.** The paired start–end configurations $\{(x_i^{(k, \mathrm{start})}, x_i^{(k, \mathrm{end})}, \hat w_i^\pi)\}$ — with stationary weights applied after Stage 2B — feed the kEDMD cross-validation route and the RevVAMPnet tertiary route. The lag of each pair is exactly $\tau_{\mathrm{WE}}$, which becomes the physical lag of the DMD and kEDMD eigenvalue problems.
2.  **Walker-endpoint corpus with iteration labels.** The configurations $\{(x_i^{(k, \mathrm{end})}, w_i^{(k)})\}_{i, k}$ with iteration index $k$ retained are the primary input to the density-ratio kernel assembly of Stage 3; the iteration labels index the columns of $\mathbf{R}$, and the weights $w_i^{(k)}$ weight the kernel sum in the numerator (4.1a).
3.  **Pooled walker-configuration corpus.** The pooled set $\{x_i\}_i$ — the union of all walker positions across all iterations and across start/end, without iteration labels — feeds both the importance-weighted tICA estimator (Stage 2A, with raw weights pooled into the weighted covariance) and RiteWeight (Stage 2B, which reweights this pooled set to $\pi$). The tICA output is the projection $q$; RiteWeight's output is $\hat w_i^\pi$.

All three extracted objects derive from the same underlying WE corpus and inherit the unbiasedness guarantee (5.2). This is the last point in the pipeline at which "the data" is a single, unambiguous object. Every subsequent stage is an algebraic transformation of one of the three extracted corpora, which is why the pipeline's correctness can be audited stage by stage rather than end to end.

Bidirectional seeding and pooling of independent WE runs into a single corpus require careful treatment of the statistical and methodological assumptions behind splitting and merging probability weights, of what mode excitation means under the distribution view of simulations, and of the need for thorough coverage of conformational space. The subsections that follow establish five points for WEeDS protocol design: (i) the asymmetric advance of bidirectionally seeded populations is not a structural problem; (ii) pooling multiple runs with diverse initial conditions is the correct strategy and is what the pipeline is designed to consume; (iii) the framing of "mode excitation" requires care, because it does not mean what the rate-method intuition suggests; (iv) there is a genuine and underappreciated tradeoff between coverage and excitation that bidirectional seeding navigates; and (v) there is a genuine and underappreciated risk that destructive interference between runs can attenuate specific modes in the pooled DR-DMD signal --- a risk that kEDMD on the same corpus is immune to, and that the three-way cross-validation architecture catches.

### 4.9 Mode excitation, precisely

#### 4.9.1 The Koopman eigenvalue equation and its conditional-expectation form

The Koopman operator acts on observables by

$$
\bigl(\mathcal{K}_\tau f\bigr)(x) \;=\; \mathbb{E}\bigl[f(X_{t+\tau}) \,\big|\, X_t = x\bigr].
$$

The eigenvalue problem

$$
\mathcal{K}_\tau\, \psi_i \;=\; \lambda_i\, \psi_i
$$

states that propagating $\psi_i$ forward by lag $\tau$ multiplies it by $\lambda_i$. The eigenfunctions are global properties of the operator on $\mathrm{supp}(\pi)$; the eigenvalues are scalars that the same propagation step yields at every $x$.

#### 4.9.2 What "excitation" means in the snapshot view

In the snapshot view, mode $i$ contributes to $\mathbf{R}_{jk}$ only when its coefficient $c_i \neq 0$. The coefficient

$$
c_i \;=\; \mathbb{E}_{p(\cdot,0)}[\psi_i] - \mathbb{E}_\pi[\psi_i]
$$

is determined by the initial distribution alone; it does not depend on subsequent dynamics. A mode with $c_i = 0$ is invisible to DR-DMD regardless of how long the simulation runs, because the transient density is identical to what it would be in the absence of that mode.

The excitation requirement is therefore a statement about the *initial condition*: the initial walker distribution must differ from $\pi$ in a way that has nonzero projection onto each mode the analyst intends to recover. This is the algorithmic content of Condition 2 of the methods doc.

#### 4.9.3 What "excitation" does not mean

The excitation requirement does not mean that:

-   The system must be "driven" or "perturbed" during the simulation.
-   Flux across barriers is required during the run.
-   Mode amplitudes grow with time, or the simulation "stirs up" slow modes.

Once the initial condition has $c_i \neq 0$, mode $i$ contributes to the density-ratio matrix at every subsequent iteration with amplitude

$$
c_i\, \lambda_i^k, \qquad |\lambda_i| < 1
$$

which decays as $k$ grows. The decay is the signal. The simulation's role is to allow this decay to be sampled across the iteration range, not to excite anything further. Mode excitation is an initial-condition property; mode observation is a sampling property.

#### 4.9.4 Why long equilibrium MD trajectories do not have this requirement

A long equilibrium MD trajectory does not have an explicit excitation requirement at the algorithmic level. The pair-based estimators built on such a trajectory read the spectrum from time-lagged correlations of the form

$$
C_\tau(f, g) \;=\; \mathbb{E}_\pi\bigl[f(X_t)\, g(X_{t+\tau})\bigr].
$$

Under reversibility and ergodic sampling, this satisfies the spectral decomposition

$$
C_\tau(f, g) \;=\; \sum_i \lambda_i\, \langle f, \psi_i \rangle_\pi\, \langle g, \psi_i \rangle_\pi.
$$

The information about $\lambda_i$ and $\psi_i$ is in correlations of fluctuations, not in any deviation of the empirical distribution from $\pi$.

The reason this difference exists is structural: the long-trajectory estimator works in the pair-based view, where the excitation requirement is automatically satisfied by sufficient coverage of $\mathrm{supp}(\pi)$ in the pair distribution. DR-DMD works in the snapshot view, where the same physical information is accessed through temporal decay rather than pair correlations, and the snapshot-view access route surfaces the initial-condition dependence explicitly.

Critically, this is not "DR-DMD requires something more than long-trajectory MSM analysis." Both are extracting the same spectrum of the same operator. The snapshot view exposes the initial-condition dependence as an explicit algorithmic requirement; the pair view absorbs it into the sampling-coverage condition. In a properly ergodic long trajectory, the initial condition has been forgotten by the time the analysis window is selected, and the trajectory's empirical distribution is approximately $\pi$. In a finite-iteration WE run with seeded initial conditions, the initial condition is the dominant feature of early snapshots and remains visible until the slowest mode has decayed appreciably. The information content is the same; the algorithmic surface is different.

#### 4.9.5 Why this matters for biomolecular systems

For systems where the slowest timescale exceeds the wall-clock simulation budget, the long-trajectory analysis is unavailable. A microsecond MD trajectory of a folding-scale RNA does not ergodically sample $\pi$ --- it spends most of its time in one basin, and the pair-based estimators on this trajectory have not seen $\pi$ on the slow-mode-relevant scales. The trajectory's empirical distribution is close to the conditional distribution within one basin, not to $\pi$, and tICA on it will recover within-basin slow modes (which are fast on the relevant timescale) rather than the inter-basin slow mode.

WE is the workaround for this regime. The WE corpus is deliberately non-ergodic: walkers are seeded in chosen regions, the resampling rule preserves probability mass, and the iteration-indexed snapshots track the relaxation of the deliberately non-equilibrium initial condition toward $\pi$. The snapshot view of DR-DMD is well-matched to this data structure, and the excitation requirement is a feature, not a bug: the snapshot view's algorithmic surface explicitly handles the non-equilibrium initial condition that WE is built to exploit.

#### 4.9.6 Two dual sampling strategies for the same operator

The framing that has emerged is that long-trajectory MD and WE-snapshot DR-DMD are *dual sampling strategies* for hitting the same operator $\mathcal{K}_\tau$:

-   **Strategy A (equilibrium-fluctuation reading)**: Run a long trajectory at or near equilibrium. The system perpetually fluctuates away from $\pi$ and relaxes back; those fluctuations are small in amplitude but constantly present. Read the spectrum from time-lagged correlations of fluctuations. Pay in trajectory length; gain in algorithmic simplicity (no reweighting, no excitation tracking). Works when the slowest timescale is short enough to be sampled by ergodic exploration.

-   **Strategy B (large-perturbation snapshot reading)**: Seed walkers in a deliberately far-from-$\pi$ initial distribution. The system relaxes toward $\pi$ over the simulation window; the iteration-indexed density snapshots track this relaxation. Read the spectrum from temporal decay of the density-ratio matrix. Pay in non-equilibrium bookkeeping (RiteWeight, importance-weighting, excitation tracking); gain in feasibility for high-barrier systems where Strategy A is intractable. Works whenever the initial perturbation projects nontrivially onto the modes of interest.

The signal magnitudes are dual:

-   Strategy A: signal amplitude per pair $\sim$ equilibrium fluctuation amplitude $\sim \sqrt{\mathrm{Var}_\pi[\psi_i]} \sim 1$ in $L^2(\pi)$ units; total signal $\sim N_{\mathrm{pair}}^{\mathrm{MD}} \cdot 1$ scaled by $\sqrt{N_{\mathrm{pair}}}$ in the noise.

-   Strategy B: signal amplitude per snapshot $\sim c_i$, the initial-condition perturbation magnitude, which can be much larger than $1$ when the initial distribution is far from $\pi$; total signal $\sim N_{\mathrm{snap}} \cdot c_i$ at each snapshot, decaying as $\lambda_i^k$ across snapshots.

For folding-scale problems, $c_i$ can be order unity or larger and $\lambda_i^k$ decays slowly over many iterations, so the per-snapshot signal is strong. Strategy A's per-pair signal is also order unity in $L^2(\pi)$, but Strategy A's *coverage* of $\mathrm{supp}(\pi)$ fails on the timescale of folding events --- the trajectory never visits the rare basin --- so the slow-mode information is absent regardless of pair counts.

This dual framing is structurally important because it locates DR-DMD as the natural estimator for a sampling regime (large initial perturbation, short observation window with non-equilibrium relaxation) that long-trajectory analysis was never designed for. The excitation requirement is not an additional complication; it is the algorithmic surface of the strategy that makes folding-scale spectral analysis feasible at all.

### 4.10 Coverage of $\operatorname{supp}(\pi)$ versus observed barrier crossings

The WE configuration chosen for this pipeline — equilibrium WE via the WEED driver with bidirectional seeding — is the *correct* simulation protocol for a spectral-recovery analysis precisely because the analysis depends on coverage of $\operatorname{supp}(\pi)$ rather than on observation of reactive trajectories. This distinction is structural to the method and is stated here as a formal part of the simulation-protocol specification rather than as a motivational remark; misunderstanding it has been a recurring source of confusion about WE-based spectral methods, and the coverage-versus-crossing distinction is the single most important operational point that makes to the reader.

**Claim.** Let $\mathcal{C} = \{(x_i^{(k,\mathrm{start})}, x_i^{(k,\mathrm{end})}, w_i^{(k)})\}_{i, k}$ be the WE corpus. The density-ratio entries $$
\mathbf{R}_{jk} \;=\; \frac{\hat{p}(x_j, t_k)}{\hat{\pi}(x_j)} - 1 \;=\; \frac{\sum_i w_i^{(k)}\, K_h\bigl(q(x_j) - q(x_i^{(k,\mathrm{end})})\bigr)}{\sum_i \hat{w}_i^\pi\, K_h\bigl(q(x_j) - q(x_i)\bigr)} - 1
$$ are functionals of the walker configurations $\{x_i^{(k,\mathrm{end})}\}$, the walker positions $\{x_i\}$ entering $\hat{\pi}$, and the weights $w_i^{(k)}$ and $\hat{w}_i^\pi$. They are *not* functionals of whether any individual walker trajectory connects two basins in a single $\tau_{\mathrm{WE}}$ interval. Analogously, the kEDMD Gram-matrix entries $G_{ij} = K_h(q(x_i) - q(x_j))$ and the kEDMD reweighted covariance matrices $\hat{C}_0, \hat{C}_\tau$ are sums over trajectory pairs $(x_i^{(k,\mathrm{start})}, x_i^{(k,\mathrm{end})})$ evaluated at lag $\tau_{\mathrm{WE}}$; they are *not* sums over reactive events.

**Consequence.** A walker that samples *local* action of the transfer operator over $\tau_{\mathrm{WE}}$ in basin $A$, without ever reaching basin $B$, contributes three distinct quantitative pieces of information to the spectral estimators: (i) a term $w_i^{(k)} K_h(q(x_j) - q(x_i^{(k,\mathrm{end})}))$ to every column $k$ of $\hat{p}$ at every reference point $x_j$ in the kernel neighborhood of the walker's endpoint, which updates every entry $\mathbf{R}_{jk}$ for which $x_j$ is within a few bandwidths of the walker's endpoint; (ii) a term $\hat{w}_i^\pi K_h(q(x_j) - q(x_i))$ to $\hat{\pi}(x_j)$ at every reference point in the kernel neighborhood, which enters the denominator of *every* column of $\mathbf{R}$; (iii) a trajectory-pair contribution to the kEDMD Gram matrix and reweighted covariances, discretizing the transfer operator's local action near $x_i^{(k,\mathrm{start})}$. None of these contributions requires the walker to have crossed the barrier.

**Analytical content.** The spectral decomposition writes the transfer operator's eigenvalues and eigenfunctions as objects defined pointwise on $\operatorname{supp}(\pi)$, weighted by the stationary density. Equivalently, the Klus et al.\textsuperscript{[29]} generalized eigenvalue formulation recasts the spectral problem as a statement about covariance operators $C_0$ and $C_\tau$ built from samples drawn from $\pi$, where the samples need only **cover** the support — they need not trace reactive trajectories between basins. The spectral objects (eigenvalues $\lambda_i$, eigenfunctions $\psi_i$, implied timescales $\tau_i$) are *global* properties of the transfer operator that are reconstructible from sufficient *local* samples of the operator's action, where "sufficient" means coverage of $\operatorname{supp}(\pi)$ adequate to resolve $\hat{\pi}$ and $\hat{p}(\cdot, t_k)$ at the reference set, not observation of any specific class of reactive trajectories. This is the central mathematical content of the spectral-recovery formulation and the structural reason why WE's coverage-oriented sampling is the right input to this pipeline.

**Simulation-protocol consequence.** Bidirectional seeding is prescribed not as a heuristic for increasing crossing frequency but as the simulation-protocol mechanism that **guarantees coverage of** $\operatorname{supp}(\pi)$ on both sides of every barrier in the initial basin set. The walker budget required is then the budget required to achieve adequate reference-set density and RiteWeight effective sample size in each basin, which is substantially smaller than the budget required to observe one or more transition events at NESS-WE rates. The WEED driver without recycling is the correct pairing because its stationary distribution is $\pi$ itself, so the walker corpus asymptotically covers $\operatorname{supp}(\pi)$ with the correct local density; a NESS-WE recycling rule would instead concentrate walkers near the reactive flux tube and underpopulate non-reactive regions of $\operatorname{supp}(\pi)$, violating the coverage requirement that the density-ratio construction rests on.

**Failure diagnostic.** When the spectral output is unreliable, the correct diagnostic is coverage failure, not crossing-count failure. The primary coverage diagnostics are the RiteWeight per-basin effective sample size, the reference-set density $\hat{\pi}(x_j)$ at the walker reference points, and the mode-excitation criterion. A corpus that passes these diagnostics is sufficient for spectral recovery even if it contains zero reactive trajectories; a corpus that fails them cannot be rescued by observing a small number of reactive trajectories. This reversal of the conventional MD-based diagnostic hierarchy is a direct consequence of the coverage-oriented analytical structure.

**Relation to reactive-trajectory methods.** Transition-path sampling, milestoning, and string methods are organized around harvesting reactive events and estimating rates from the reactive flux. These methods have different information requirements and a different relationship to the force-field-generated dynamics than spectral-recovery methods have, and they should be understood as a separate class of methods rather than as direct competitors. The present pipeline does not replace them; it consumes a different kind of data (coverage of $\operatorname{supp}(\pi)$) and produces a different kind of output (a continuous-coordinate transfer-operator spectrum with a GP posterior on the FES in eigenfunction coordinates). The choice of equilibrium WE with bidirectional seeding, rather than NESS WE with reactive-flux recycling, is the simulation-protocol operationalization of this distinction.

**Validation and consistency diagnostics.** Stage (Stage 1 — WESTPA) is audited by the following diagnostics, with full specifications in the noted sections. *In-section diagnostics* (computed alongside the stage output): Hamiltonian-identity precondition; coverage of $\operatorname{supp}(\pi)$. *Core validation program* (computed for every pipeline invocation): WE convergence prerequisite; multi-run pipeline comparison (cross-run bandwidth $h^*$ stability and per-run coverage of the pooled basin set). *Additional statistical tools* (campaign-specific or deployed on demand): convergence with simulation length on the pooled corpus. *Limitations constraining this stage*: L2 initial-condition completeness; L7 RiteWeight regularization dependence (when WE coverage failure produces RiteWeight failure). The pipeline does not certify this stage's output as publication-ready until the in-section diagnostics and the core program both pass; the tools are deployed when one of those checks flags a problem or when publication-grade extension is required.

#### 4.10.1 Empirical observation: asymmetric basin advance

When the folded-seeded population of a bidirectional WE run remains spatially quiet while the unfolded-seeded population advances into the barrier face, this is not a failure of the pipeline. Local sampling of $\mathcal{K}_\tau$ inside the folded basin is being collected at every iteration by every folded-seeded walker; the density-ratio entries at folded-basin reference points evolve according to the same temporal structure as those at unfolded-basin reference points, because the eigenvalues $\lambda_i$ are global properties of the operator.

The asymmetric advance affects only one specific aspect: spatial resolution of eigenfunctions in the barrier region. For reference points $x_j$ in the transition region, both the numerator $\hat{p}(x_j, t_k)$ and the denominator $\hat{\pi}(x_j)$ require walker contributions within kernel range. If walkers approach the transition region from one side only, the kernel-smoothed density at barrier reference points is supported by one-sided contributions, giving higher sampling variance at those points than two-sided coverage would. This shows up as wider uncertainty bands on eigenfunction values near the transition --- not as loss of identifiability and not as bias in eigenvalue estimates.

For PCCA+ state assignment, GP-FES construction in eigenfunction coordinates, and most downstream observables, basin-interior eigenfunction values are what matters. Barrier-region eigenfunction structure is interesting scientifically but not essential for the pipeline to function.

### 4.11 The coverage--excitation tradeoff and bidirectional seeding

Bidirectional seeding solves the coverage problem at the cost of partially reducing the excitation amplitude for the dominant slow mode. This tradeoff is rarely articulated in the WE literature and deserves explicit statement.

#### 4.11.1 The eigenfunction normalization for asymmetric equilibria

Consider a UUCG-like system with $\pi_F \gg \pi_U$. For concreteness, take

$$
\pi_F = 0.9, \qquad \pi_U = 0.1.
$$

The eigenfunction $\psi_2$ that distinguishes the two basins satisfies the orthogonality

$$
\mathbb{E}_\pi[\psi_2] \;=\; \pi_F\, \psi_2^F + \pi_U\, \psi_2^U \;=\; 0,
$$

which gives

$$
\psi_2^U \;=\; -\frac{\pi_F}{\pi_U}\, \psi_2^F \;=\; -9\, \psi_2^F.
$$

With $L^2(\pi)$ normalization $\langle \psi_2, \psi_2 \rangle_\pi = 1$:

$$
\pi_F\, (\psi_2^F)^2 + \pi_U\, (\psi_2^U)^2 \;=\; 1
$$

$$
0.9\, (\psi_2^F)^2 + 0.1 \cdot 81\, (\psi_2^F)^2 \;=\; 9\, (\psi_2^F)^2 \;=\; 1.
$$

So

$$
|\psi_2^F| \;=\; 1/3, \qquad |\psi_2^U| \;=\; 3.
$$

The eigenfunction has small magnitude in the high-probability basin and large magnitude in the rare basin. This is a generic feature when basin populations are asymmetric: rare regions carry small probability but large eigenfunction values.

Choose the sign convention $\psi_2^F = -1/3$, $\psi_2^U = +3$.

#### 4.11.2 Unidirectional excitation coefficients

For unidirectional seeding entirely in the folded basin, the initial distribution puts all walkers at $\psi_2^F = -1/3$:

$$
c_2^{\text{uni-F}} \;=\; \psi_2^F - 0 \;=\; -1/3.
$$

For unidirectional seeding entirely in the unfolded basin:

$$
c_2^{\text{uni-U}} \;=\; \psi_2^U - 0 \;=\; 3.
$$

The unidirectional unfolded seeding excites $\psi_2$ much more strongly than the unidirectional folded seeding, by

$$
\frac{|c_2^{\text{uni-U}}|}{|c_2^{\text{uni-F}}|} \;=\; \frac{3}{1/3} \;=\; 9,
$$

because the eigenfunction's magnitude is concentrated in the rare basin.

But unidirectional folded seeding has only the folded basin covered; the unfolded basin acquires walkers only through rare folding events, which at UUCG barriers occur on timescales far exceeding any feasible budget. The denominator $\hat{\pi}(x_j)$ at unfolded-basin reference points is zero or kernel-tail-only, and the spectral analysis fails categorically on the unfolded side. Symmetrically, unidirectional unfolded seeding fails on the folded side.

#### 4.11.3 Bidirectional excitation

Bidirectional seeding splits walkers $\alpha$ folded / $(1-\alpha)$ unfolded. The excitation coefficient is

$$
c_2^{\text{bidir}}(\alpha) \;=\; \alpha\, \psi_2^F + (1-\alpha)\, \psi_2^U - 0.
$$

A general identity: subtracting and adding $\mathbb{E}_\pi[\psi_2] = \pi_F \psi_2^F + \pi_U \psi_2^U$ inside the expression,

$$
c_2^{\text{bidir}}(\alpha) \;=\; (\alpha - \pi_F)\, \psi_2^F + ((1-\alpha) - \pi_U)\, \psi_2^U \;=\; (\alpha - \pi_F)\, (\psi_2^F - \psi_2^U).
$$

This shows that $c_2 = 0$ if and only if $\alpha = \pi_F$ --- the measure-zero coincidence where the seeding ratio matches the equilibrium basin ratio. The 50/50 default ($\alpha = 0.5$) gives

$$
c_2^{\text{bidir}}(0.5) \;=\; (0.5 - 0.9)(-1/3 - 3) \;=\; (-0.4)(-10/3) \;=\; 4/3 \;\approx\; 1.33.
$$

#### 4.11.4 The tradeoff quantified

Comparing magnitudes for $\pi_F = 0.9$:

| Protocol            | $c_2$   | Coverage of $\mathrm{supp}(\pi)$ |
|:--------------------|:--------|:---------------------------------|
| Unidirectional U    | $+3.00$ | Unfolded only (folded fails)     |
| Bidirectional 50/50 | $+1.33$ | Both basins                      |
| Unidirectional F    | $-0.33$ | Folded only (unfolded fails)     |

Bidirectional seeding *loses excitation amplitude relative to unidirectional unfolded seeding* by a factor of

$$
\frac{|c_2^{\text{uni-U}}|}{|c_2^{\text{bidir}}|} \;=\; \frac{3}{4/3} \;=\; 2.25
$$

and *gains amplitude relative to unidirectional folded seeding* by a factor of 4, while *guaranteeing both-basin coverage from iteration 1* in either case.

The structural meaning: bidirectional seeding trades some excitation amplitude for the coverage guarantee. The trade is asymmetric in a useful way:

-   **Coverage failure is identifiability failure** (a mode is unrecoverable on the uncovered region).
-   **Excitation reduction is estimation-quality degradation** (the mode is still recoverable, with wider confidence intervals).

Trading identifiability for CI width is strongly favorable when identifiability failure is otherwise certain, which is the regime for high-barrier biomolecular systems on accessible simulation budgets.

#### 4.11.5 The asymmetric-equilibrium effect, generalized

The above exposes a structural feature worth stating: for asymmetric equilibria, the eigenfunction's magnitude is concentrated in the rare basin, and unidirectional seeding *of the rare basin* gives the strongest excitation of the dominant slow mode. This follows from $L^2(\pi)$ normalization: rare regions carry small $\pi$ but large $|\psi|$, and an initial distribution localized there has a large $\mathbb{E}[\psi]$ that the equilibrium $\mathbb{E}_\pi[\psi] = 0$ does not match.

In practical protocol design, this means a mix of bidirectional and unidirectional unfolded-side runs can strengthen $\psi_2$ excitation in the pooled corpus relative to pure bidirectional, without sacrificing coverage. See for the pooled-corpus analysis.

#### 4.11.6 Higher modes

For higher modes $\psi_3, \psi_4, \ldots$, the analogous excitation coefficients depend on the spatial detail of the initial distribution within each basin, not just on basin-occupancy ratios. Higher modes typically distinguish substructures within basins or alternative pathways between basins; their eigenfunctions have more complex sign patterns than $\psi_2$'s simple folded/unfolded split.

The implication: diverse bstates (WESTPA basis-state structures used to seed walkers) within each basin are required to excite higher modes, where "diverse" means spread across the basin in directions aligned with the modes' eigenfunction variations. A bstate set clustered tightly within an unfolded basin (e.g., all unfolded seeds at one local region) projects onto $\psi_2$ well, because both basins are represented, but projects onto $\psi_3$ poorly if $\psi_3$ distinguishes substructures within the unfolded basin and all seeds sit at the same value of $\psi_3$. This is the operational reason for prioritizing within-basin bstate diversity, separate from the bidirectional/unidirectional decision.

> **Caveat (higher-mode cancellation, open).** Whether higher modes ($\psi_3, \psi_4, \ldots$) have predictable cancellation patterns under specific within-basin seed designs is an empirical question that the dipeptide validation can address. A controlled experiment varying within-basin spread and tracking $c_3^{\mathrm{pool}}$ and $c_4^{\mathrm{pool}}$ across pooled corpora would establish operational principles for higher-mode-targeted run design.

> **Caveat (asymmetric-equilibrium generalization, tentative).** The observation that unidirectional rare-basin seeding maximizes $\psi_2$ excitation may generalize beyond UUCG: for ligand binding/unbinding with strongly asymmetric equilibria, preferentially seeding the dissociated state, and for protein folding with strongly stable native states, preferentially seeding the unfolded ensemble, would both be expected to strengthen the dominant slow mode. This is suggestive rather than established, and indicates that the 50/50 default is not always optimal.

### 4.12 Identifiability conditions and the window of observability

Beyond excitation, three further conditions determine which modes are recoverable from a finite WE experiment.

Recovering $n - 1$ eigenmodes from the density-ratio matrix $\mathbf{R}$ requires several conditions to hold simultaneously. Where a mode fails to satisfy one of them, the pipeline is expected to return a specific, identifiable failure mode rather than silently mis-report. This enumeration is exhaustive for the DR-DMD primary route; analogous conditions for kEDMD (aggregation across walkers and iterations) and VAMPnet (neural-network expressivity and training convergence) are stated for those routes in the companion Pipeline document, which develops these conditions further in their pipeline-specific contexts.

#### 4.12.1 Number of time points (Prony bound); mode excitation

**Condition 1: number of time points.** At minimum, $N \ge 2(n-1)$ WE iterations are required to separate $n - 1$ modes — this is a Prony-type argument: the density-ratio matrix must admit a rank-$(n-1)$ factorisation in its non-trivial subspace, and fewer than $2(n-1)$ columns cannot uniquely identify that many distinct exponential rates. Substantial oversampling, $N \gg n$, is needed for numerical stability, particularly when eigenvalues are closely spaced. For biomolecular systems with the typical target $n = 2$–$5$ metastable states, $N \ge 50$–$100$ iterations is a practical lower bound; the validation pipeline reports the effective rank from the SVD spectrum as an empirical check.

At minimum,

$$
N \;\geq\; 2(n-1)
$$

WE iterations are required to separate $n-1$ modes. This is a Prony-type argument: the density-ratio matrix must admit a rank-$(n-1)$ factorization in its non-trivial subspace, and fewer than $2(n-1)$ columns cannot uniquely identify that many distinct exponential rates. Substantial oversampling, $N \gg n$, is needed for numerical stability, particularly when eigenvalues are closely spaced.

For typical biomolecular targets $n = 2$–$5$ metastable states, $N \geq 50$–$100$ iterations is a practical lower bound.

**Condition 2: mode excitation.** The coefficient

$$
c_i \;=\; \langle p(\cdot, 0)/\pi - 1,\, \psi_i\rangle_\pi \;=\; \mathbb{E}_{p(\cdot, 0)}[\psi_i] - \mathbb{E}_\pi[\psi_i]
$$

must be nonzero for each eigenfunction to be recoverable. A mode with $c_i = 0$ is invisible to DR-DMD regardless of how long the simulation runs, because the transient density is identical to what it would be in the absence of that mode. In practice, localised initial distributions (walkers started from a single conformation or a narrow structural ensemble) excite the sub-leading mode $\psi_2$ (which distinguishes basins) strongly and higher modes weakly; this determines the order in which modes are recovered as simulation length increases. A seeded initial distribution

$$
p(\cdot, 0) = \sum_\alpha W_\alpha\, \delta(\cdot - y_\alpha)
$$ has $c_i = \sum_\alpha W_\alpha\, \psi_i(y_\alpha)$

which vanishes only in the coincidental case $\sum_\alpha W_\alpha\, \psi_i(y_\alpha) = 0$ — a measure-zero condition on the allocation weights $\{W_\alpha\}$ given fixed seed positions $\{y_\alpha\}$. For two-basin seeding with weights $(\alpha, 1-\alpha)$, this reduces to $\alpha = \pi_F$ in the delta-seed limit — i.e., the seeding ratio exactly matching the equilibrium basin ratio — which for asymmetric equilibria ($\pi_F \ne 1/2$) is not satisfied by the default 50/50 allocation. Symmetric 50/50 seeding therefore nulls $\psi_2$ only on systems whose equilibrium basin populations are themselves 50/50. Using multiple initial conditions (e.g. bidirectional seeding from folded and unfolded basins, or seeding from multiple high-energy transition regions) generically excites all modes and is strongly recommended where practicable. This is the theoretical ground for the bidirectional-seeding strategy.

#### 4.12.2 Temporal coverage and the window of observability; temporal spacing

**Condition 3: temporal coverage and the window of observability.** The time span of the WE simulation must extend from $t \sim 0$ (where fast modes are still active and detectable) to $t \gg \tau_2$ (long enough for the slowest mode to have decayed appreciably). For each mode $i$ to be resolved, its decay factor $\lambda_i^k$ must vary meaningfully across the observed iteration range $k = 1, \ldots, N$: if the simulation is too short for the slowest mode to decay significantly, its eigenvalue cannot be estimated — the temporal profile is approximately constant over the observation window, and the corresponding spatial mode is indistinguishable from the stationary distribution. Conversely, if the simulation is much longer than $\tau_i$ ($\tau_i \ll \tau_{\mathrm{WE}}$), the mode has already decayed to within noise across most of the corpus and contributes no signal beyond an early transient.

For each mode $i$ to be resolved, its decay factor $\lambda_i^k$ must vary meaningfully across the observed iteration range $k = 1, \ldots, N$. The **window of observability** for mode $i$ is:

$$
N\tau_{\mathrm{WE}} \;\gg\; \tau_i \;\gg\; \tau_{\mathrm{WE}},
$$

which translates into a concrete constraint on the number of WE iterations required to recover a mode of given timescale. This is the analog of the implied-timescale plateau test of Noé & Fischer\textsuperscript{[53]}: modes whose timescales are much longer than the total simulation time are not identifiable; modes whose timescales are comparable to a single iteration are not resolved.

**Condition 4: temporal spacing.** For a Prony-type inverse problem, exponentially spaced time points $t_k = t_0\, \alpha^k$ are generally optimal, as they provide resolution across multiple decades of decay rates. The WE data structure, however, produces equally spaced iterations $t_k = k\,\tau_{\mathrm{WE}}$ — an artefact of the WE resampling cadence that is beyond the pipeline's control. Equally spaced points are acceptable but less efficient: they oversample the fast modes and undersample the slow modes. In practice this is the regime in which DR-DMD operates, and the effective consequence is that the long-timescale end of the spectrum has the tightest per-decay-time sampling, which is beneficial for the slowest modes but requires the total simulation length to be long enough for the slow modes to decay. Hankel-DMD partially compensates by constructing stacked delay embeddings that effectively interpolate between the observed iterations.

#### 4.12.3 Spatial resolution and kernel bandwidth as a coupled condition

**Condition 5: spatial resolution.** The reference set $\{x_j\}_{j=1}^M$ must resolve the eigenfunction structure — $M$ must be large enough that the eigenfunctions are not aliased by the sparsity of the reference set. In conventional grid-based density-ratio estimation this would impose $M \sim 50$–$100$ per dimension (so $M \sim 2500$–$10\,000$ for a two-dimensional slow manifold), which quickly becomes intractable in higher dimensions. The kernel-evaluated density-ratio construction sidesteps this constraint: the reference set is drawn from the walker positions themselves (or a subsample thereof), so $M$ scales with the number of walkers (typically $10^3$–$10^5$) rather than with a dimensional-grid product, and the density-ratio evaluation is pointwise at walker positions rather than on an imposed grid. The effective resolution is set by the kernel bandwidth $h$ in tICA space, and the pipeline's only free hyperparameter at this stage is $h$, which is chosen by the standard bias–variance trade-off.

The kernel bandwidth $h$ is the parameter that simultaneously controls both the numerator $\hat{p}(x_j, t_k)$ and the denominator $\hat{\pi}(x_j)$ of the density-ratio matrix, and its choice creates a coupled tradeoff that is worth surfacing as a structural feature, not just a hyperparameter to be optimized.

The two competing pressures on $h$:

-   $h$ too small: The kernel sums $\hat{p}$ and $\hat{\pi}$ become spiky and noisy. In regions where walker density is low (high free energy), the kernel sums fluctuate strongly between zero (no walker within bandwidth) and singleton-walker contributions, producing high-variance estimates of both numerator and denominator. The density ratio inherits the variance of both, and the $1/\hat{\pi}$ noise scaling becomes severe in low-density regions.

-   $h$ too large: The kernel sums over-smooth the spatial structure of the eigenfunctions. The eigenfunction $\psi_i$ has interesting structure on the scale of basin sizes and barrier widths; if $h$ exceeds those scales, the kernel evaluation washes out the eigenfunction's variation at the relevant scale, and the recovered $\hat{\psi}_i$ is a smeared version of the true $\psi_i$ that may fail to satisfy the eigenfunction sign-structure or orthogonality checks against kEDMD.

The Sobolev-RKHS Koopman-invariance certificate of Bold 2025 / Hertel 2026 provides a theoretical bandwidth choice that places the kernel in the Sobolev space $H^\sigma(\Omega)$ at the order

$$
\sigma \;=\; 5/2 + d_{\mathrm{tIC}}/2
$$

for the Matérn-$5/2$ kernel in $d_{\mathrm{tIC}}$-dimensional tICA coordinates. This is a theoretically grounded choice but it is one choice out of a family, and the empirical question of whether this specific $h$ delivers acceptable SNR on a given finite corpus is what the alanine dipeptide validation answers.

The structural coupling of $h$ across numerator and denominator means that $h$ tuning is not a degree of freedom that can be exercised separately for the two; a single bandwidth must serve both. The practical consequence: if alanine dipeptide validation finds that the Sobolev-optimal $h$ is too small for adequate denominator support in barrier regions (high noise floor) or too large for adequate spatial resolution (smeared eigenfunctions), the response is to vary $h$ as a calibration parameter and choose the value that maximizes recovery of known PyEMMA modes, accepting that this empirical choice may not coincide exactly with the theoretical optimum.

#### 4.12.4 Noise amplification in the tails; spectral gap

**Condition 6: noise amplification in the tails.** The division by $\pi(x)$ in the density ratio amplifies statistical noise in regions where $\pi(x) \approx 0$ (high-free-energy regions). This is mathematically inevitable — the density ratio $p/\pi$ can take arbitrarily large values at barriers — and it introduces noise into the eigenfunction estimates in precisely the regions where the eigenfunctions have the most interesting structure (they change sign near saddle points). Regularisation strategies include thresholding ($\hat\pi(x) \to \max(\hat\pi(x), \pi_{\min})$), $\pi$-weighting of the DMD (promoting $\mathbf{R}$ to $\tilde{\mathbf{R}}_{jk} = \sqrt{\hat\pi(x_j)}\,R_{jk}$, which is the $\sqrt{\pi}$-weighted inner product corresponding to the $L^2(\pi)$ geometry), and restricting the reference set to regions where $\hat\pi(x) > \pi_{\mathrm{threshold}}$ and accepting that the eigenfunction values in the tails are extrapolated via the kernel representation. The pipeline defaults to the $\pi$-weighted SVD within the Schmid projected DMD with a mild $\hat\pi$-floor; the noise-floor derivation quantifies the resulting bias–variance trade-off via a delta-method expansion.

The noise floor in $\mathbf{R}$ near barriers, where $\pi(x)$ is small, follows the structural property of the density-ratio estimator developed in: the noise standard deviation scales as $\sigma_{\mathrm{noise}}(x_j)/\hat{\pi}(x_j)$, which diverges as $\hat{\pi}(x_j) \to 0$. The $\pi$-weighting

$$
\tilde{\mathbf{R}}_{jk} \;=\; \sqrt{\hat{\pi}(x_j)}\, \mathbf{R}_{jk}
$$

corresponds to the $L^2(\pi)$ inner product geometry and suppresses the divergence at the cost of mild bias near boundaries.

**Condition 7: spectral gap.** As discussed, DR-DMD requires a clean gap $|\lambda_n| - |\lambda_{n+1}| \gtrsim 0.1$ to recover $n - 1$ slow modes cleanly. In the absence of a clean gap, near-degenerate modes are mixed in the DMD output; the consistency check is that VAMPnet and kEDMD recover the same near-degenerate pair (or that Hankel-DMD separates them). This is the theoretical condition that underpins the dimensionality-selection procedure and the three-way cross-validation.

The spectral gap for mode $i$ is

$$
g_i \;=\; |\lambda_i| - |\lambda_{i+1}|.
$$

DR-DMD relies on SVD rank selection within the Schmid projected algorithm. When $g_i$ is small, the SVD singular values corresponding to $\psi_i$ and $\psi_{i+1}$ cluster numerically, the DMD eigenvalues appear as a near-degenerate pair, and the recovered eigenfunctions mix between modes. The operational consequence: DR-DMD requires either $g_i \gtrsim 0.1$ for typical biomolecular applications, or the Hankel-DMD extension to recover near-degenerate modes reliably. The three-way cross-validation catches this failure when DR-DMD disagrees with kEDMD or VAMPnet.

**Summary.** The seven conditions above are jointly sufficient to guarantee recovery of $n - 1$ slow modes from a WE corpus via DR-DMD. In practice, Conditions 1, 3, and 7 are the binding ones: the simulation length and the spectral gap are the two dominant constraints, with mode excitation (Condition 2) addressed by bidirectional seeding and temporal spacing (Condition 4) accepted as a suboptimality of the WE data structure rather than a limiting factor. Conditions 5 and 6 — spatial resolution and noise amplification — are mitigated by the kernel-evaluated construction and the $\pi$-weighting. Where one of these conditions fails, the validation framework returns a specific diagnostic that maps back to the failing condition: short simulation $\to$ ITS non-plateau (Condition 3); missing modes $\to$ bidirectional-seeding recheck (Condition 2); near-degenerate pair $\to$ Hankel-DMD fallback (Condition 7). The density-ratio-specific and DMD-specific refinements of these conditions are developed.

> **Caveat (bandwidth--pooling interaction).** The kernel-bandwidth tradeoff of this subsection interacts with the pooling break-even criterion (\S4.14.3.4). A larger bandwidth $h$ reduces per-cell noise (more walkers within range) at the cost of smearing eigenfunction structure; a smaller $h$ preserves structure at higher per-cell noise. Because the noise reduction from pooling $N_{\mathrm{runs}}$ runs can substitute for the noise reduction from a larger bandwidth, a pooled corpus may tolerate a smaller $h$ than any single-run analysis would. Whether this is a useful degree of freedom is an empirical question for the dipeptide validation.

### 4.13 Pooling multiple WE runs: pooled-corpus architecture and reweighting

The pipeline consumes a **pooled corpus** across $N_{\mathrm{runs}}$ independent WE simulations. The pool is the formal input to RiteWeight, kernel-evaluated density-ratio assembly, $\pi$-weighted Exact DMD, kEDMD, kinetic RiteWeight, PCCA+ and weighted-MSM construction, and the GP free energy surface. Per-run analyses are *diagnostic projections* of the pooled headline, not parallel meta-analytic estimates that get combined into an alternative CI. The pooled-corpus architecture is what allows the spectral and statistical machinery to operate at the largest possible effective sample size, and it is the structural feature on which the joint-loop bootstrap rests.

**Pooled walker weight specification.** Let $W^{(\alpha)} \ge 0$ be the analyst-assigned weight for run $\alpha = 1, \ldots, N_{\mathrm{runs}}$ with $\sum_\alpha W^{(\alpha)} = 1$. The pooled walker weight at iteration $k$ of run $\alpha$ is $$
\tilde{w}_i^{(\alpha,k)} \;=\; W^{(\alpha)}\, w_i^{(\alpha,k)},
$$ where $w_i^{(\alpha,k)}$ is the raw WE weight for walker $i$ at iteration $k$ of run $\alpha$. The pooled weights replace the per-run $w_i^{(k)}$ throughout every downstream stage; the kernel-evaluated density-ratio numerator, the importance-weighted tICA covariance estimators, and the kEDMD reweighted covariances all see the pool as a single reweighted sample.

**Per-iteration sum and unbiasedness.** Because each run satisfies $\sum_i w_i^{(\alpha,k)} = 1$ (WE weight conservation) and $\sum_\alpha W^{(\alpha)} = 1$, the pooled per-iteration weight sum is $$
\sum_\alpha \sum_i \tilde{w}_i^{(\alpha,k)} \;=\; \sum_\alpha W^{(\alpha)} \;=\; 1
$$ at every iteration $k$. By the Zhang–Jasnow–Zuckerman unbiasedness theorem applied to each run independently and the linearity of expectation, $$
\mathbb{E}\!\left[\sum_\alpha \sum_i \tilde{w}_i^{(\alpha,k)}\, f\!\bigl(x_i^{(\alpha,k,\mathrm{end})}\bigr)\right] \;=\; \sum_\alpha W^{(\alpha)} \int f(x)\, p^{(\alpha)}(x, t_k)\, dx
$$ for any bounded measurable $f$, where $p^{(\alpha)}(x, t_k)$ is the transient distribution from run $\alpha$'s initial condition. Under Schemes where all runs share the same initial-condition distribution (Scheme 2 with iid-equivalent seeding across runs), $p^{(\alpha)} = p$ is independent of $\alpha$ and the pooled estimator $\hat{p}_{\mathrm{pool}}(x, t_k)$ is unbiased for $p(x, t_k)$. Under Scheme 1 with mixed seeding types, the pool is unbiased for the *mixture transient distribution* $\sum_\alpha W^{(\alpha)}\, p^{(\alpha)}(x, t_k)$; the spectral expansion is linear and the recovered eigenfunctions and eigenvalues remain those of the same transfer operator, with mixture-weighted excitation coefficients $c_i^{\mathrm{mix}} = \sum_\alpha W^{(\alpha)}\, c_i^{(\alpha)}$.

**Run-and-iteration label bookkeeping.** Every walker in the pool carries a triple $(\alpha, k, i)$ identifying its source run $\alpha$, iteration $k$, and walker index $i$ within iteration $k$ of run $\alpha$. The labels are required at three points downstream: (i) the density-ratio matrix assembly groups walkers by iteration $k$ to form column $k$ of $\mathbf{R}$, summing across all runs $\alpha$ at fixed $k$; (ii) the joint-loop bootstrap preserves the run-label $\alpha$ when applying the Dirichlet draws and preserves the iteration label $k$ when block-resampling within run $\alpha$; (iii) the multi-run pipeline comparison reruns the full pipeline on each run $\alpha$ in isolation, requiring the run label to project the pool back to its constituent runs. The label triples are stored alongside the WE corpus on disk and propagate through every stage.

**Iteration-count handling: truncation to** $N_{\min}$. When the pooled runs have differing iteration counts $\{N_\alpha\}_{\alpha=1}^{N_{\mathrm{runs}}}$, the density-ratio matrix is built from columns $k = 1, \ldots, N_{\min}$ where $N_{\min} = \min_\alpha N_\alpha$. Iterations $k > N_{\min}$ in any longer run are *not* included in the density-ratio matrix and are *not* used by the spectral pipeline (DR-DMD, kEDMD, RevVAMPnet) or by the time-resolved visualization. This truncation enforces uniform per-column effective sample size across the matrix, preserves the standard DMD theory's uniform-quality assumption (Tu 2014, Schmid 2010, Klus 2018), keeps the singular-value spectrum and Gavish–Donoho rank threshold interpretable without per-column heterogeneity correction, and yields a bootstrap CI of uniform statistical weight at every retained eigenvalue. The surplus iterations of longer runs are *not* discarded from the corpus: they remain available to RiteWeight, which is iteration-blind and pools across all walkers from all iterations of all runs to estimate $\hat{\pi}$. Surplus iterations *cannot* contribute to any iteration-resolved object — the time-resolved visualization, the per-iteration density-ratio columns, the per-iteration mode-coefficient series(c) — because the methodological discipline is that no figure or movie shows data that was not included in the eigenfunction extraction. The recommended analyst protocol is to plan WE campaigns with comparable target iteration counts so the truncation discards little. When the surplus is large (any run with $N_\alpha > 1.5\,N_{\min}$), the truncation is flagged in the campaign report as a sampling-efficiency limitation and the analyst is advised to extend the shorter runs before re-running the analysis.

**Three pooling schemes.** Three seeding configurations are documented; the analyst chooses based on the tractability of single-basin coverage rather than on statistical preference among the schemes.

-   **Scheme 1 — Single-basin pooled.** $N_F$ folded-start runs and $N_U$ unfolded-start runs pooled at $(W_F, W_U)$ with $W_F + W_U = 1$. Each constituent run individually fails to cover the opposite basin; the pool covers both. The per-run weight $W^{(\alpha)}$ is set to $W_F / N_F$ for folded-start runs and $W_U / N_U$ for unfolded-start runs, so that the *total* contribution of folded-start runs to the pool is $\sum_{\alpha \in F} W^{(\alpha)} = W_F$ and the *total* contribution of unfolded-start runs is $\sum_{\alpha \in U} W^{(\alpha)} = W_U$, with $W_F + W_U = 1$ recovering the per-iteration sum-to-one property. When tractable, Scheme 1 is preferred because it admits an asymmetry-sweep validator: thermodynamic observables (FES, populations) are invariant under sweeps of $(W_F, W_U)$ while kinetic observables re-weight predictably.
-   **Scheme 2 — Multi-run bidirectional pooled (recommended default).** $N_{\mathrm{runs}} \ge 5$ independent runs, each itself bidirectionally seeded. Pooled at uniform weights $W^{(\alpha)} = 1/N_{\mathrm{runs}}$; reduces to $\tilde{w}_i^{(\alpha,k)} = w_i^{(\alpha,k)} / N_{\mathrm{runs}}$. Default choice for systems where single-basin coverage is intractable.
-   **Scheme 3 — Mixed pool.** Both run types present, weighted at analyst-chosen $\{W^{(\alpha)}\}$ summing to unity. Catch-all configuration for diagnostic and incremental campaigns.

**Hamiltonian-identity precondition.** Pooling validity requires identical force field, identical solvent and ion molecule counts, and identical box geometry across all $N_{\mathrm{runs}}$ runs in the pool. Independent solvent and ion *placements* within those fixed counts are the source of run independence and are required (and sufficient) for the runs to constitute independent samples of the same equilibrium distribution. Violations of the Hamiltonian-identity precondition (e.g., different solvent box sizes, different ion concentrations, different protonation states between runs) invalidate the pool: the reweighted sample no longer represents a single $\pi$, and every downstream estimator inherits a systematic bias that the validation framework is not designed to detect. The Hamiltonian-identity precondition is verified at the simulation-protocol definition stage and is distinct from the per-run runtime convergence prerequisite.

**Hamiltonian-identity audit checklist for archival pooling.** The Hamiltonian-identity precondition above is not generically satisfied by archival WE corpora aggregated post-hoc; the analyst must verify it explicitly before any pooled estimator is computed. The operational checklist comprises five items:

(i) **parameter-file hash equality** across all $\alpha$ — bytewise checksum match of the force-field parameter files (parm7, prmtop, XML, itp, or equivalent), including any small-molecule parameter overrides, ligand parameter sets, and modified-residue libraries;

(ii) **solvent-and-ion molecule-count equality** — exact integer match of water-model molecule counts (TIP3P, TIP4P-EW, OPC, SPC/E), Na$^{+}$ / Cl$^{-}$ / K$^{+}$ / Mg$^{2+}$ counts, and any buffer-species counts;

(iii) **box-geometry equality** — identical periodic-cell vectors at $t = 0$, with NPT-equilibrated runs deemed compatible only if their post-equilibration mean cell volumes agree to within $0.5\%$;

(iv) **thermostat / barostat protocol equality** — identical integrator family, timestep, thermostat coupling constant, barostat coupling constant, and temperature / pressure setpoints;

(v) **progress-coordinate definition equality** — identical bin definitions, identical recycling boundaries, and identical WESTPA driver settings (target walker count per bin, splitting / merging thresholds).

Runs failing any of (i)–(v) are excluded from the pool; the pipeline does not silently pool runs from incompatible Hamiltonians or operationally divergent simulation protocols. Failed checks and the per-run weights $W^{(\alpha)}$ that would have been assigned to excluded runs are reported in the SI audit table that accompanies every campaign-level publication. The checklist is run once at the simulation-protocol definition stage for prospectively designed campaigns and once at the corpus-aggregation stage for archival pools.

**Single-point failure mode flagged here.** If RiteWeight fails on the pool (weight degeneracy, ESS collapse, regression entering a numerically pathological regime), the headline is biased and the per-run pipeline reruns *will all show the same bias* because they share the analyst-pipeline. The weight-distribution + ESS diagnostics, the per-state weight degeneracy check, and the L7 pre-pipeline acceptance gate are the only protections; all three apply to the pool. The multi-run consistency check is *not* a check on RiteWeight failure.

#### 4.13.1 The shared iteration index and the pooled density-ratio matrix

When $N_{\mathrm{runs}}$ independent WE runs each producing $N_{\mathrm{iter}}$ iterations of length $\tau_{\mathrm{WE}}$ are pooled, the iteration index $k$ is shared across runs. The pooled transient distribution at iteration $k$ is

$$
\hat{p}_{\mathrm{pool}}(x_j, t_k) \;=\; \sum_{r=1}^{N_{\mathrm{runs}}} \alpha_r \sum_{i \in \text{run } r} w_i^{(k,r)}\, K_h\!\bigl(q(x_j) - q(x_i^{(k,r,\mathrm{end})})\bigr),
$$

where $\alpha_r$ is the run weight (typically $1/N_{\mathrm{runs}}$, or set by the Dirichlet bootstrap). Each run contributes its iteration-$k$ snapshot at the same physical time $t_k = k\tau_{\mathrm{WE}}$.

The pooled stationary distribution $\hat{\pi}_{\mathrm{pool}}(x_j)$ comes from RiteWeight applied to the pooled aggregate corpus. The density-ratio matrix is

$$
\mathbf{R}_{jk}^{\mathrm{pool}} \;=\; \frac{\hat{p}_{\mathrm{pool}}(x_j, t_k)}{\hat{\pi}_{\mathrm{pool}}(x_j)} - 1.
$$

This is one matrix, $M \times N_{\mathrm{iter}}$. DR-DMD operates on this single pooled matrix; temporal decay is read along the shared iteration axis.

The structural constraint: pooled runs must share the same $\tau_{\mathrm{WE}}$ and (after truncation if necessary) the same iteration count. Different progress coordinates across runs are allowed because the analysis space is the tICA-derived feature space, not PC space. Different bstates and different seeding ratios are allowed; their effect on excitation is what develops.

#### 4.13.2 Functional equivalence to a long ergodic trajectory

A pooled bidirectional WE corpus contains approximately

$$
N_{\mathrm{runs}} \cdot N_{\mathrm{iter}} \cdot N_{\mathrm{walkers}}
$$

trajectory pairs of length $\tau_{\mathrm{WE}}$, concentrated in regions of $\mathrm{supp}(\pi)$ that PCs drove walkers into, weighted by RiteWeight to look like draws from $\pi$. A long equilibrium MD trajectory of total length $T$ at lag $\tau_{\mathrm{WE}}$ contains

$$
N_{\mathrm{pair}}^{\mathrm{MD}} \;=\; T/\tau_{\mathrm{WE}}
$$

pairs distributed according to $\pi$ automatically, concentrated where $\pi$ puts them.

If totals are matched, the pooled WE corpus is more informative per pair about slow modes, because WE pairs are concentrated in transition regions where eigenfunctions have interesting structure, while long-trajectory pairs are concentrated in basins where eigenfunctions are flat. This is the precise sense in which pooled WE samples "effective microseconds" of dynamics in a fraction of the wall-clock time --- not literal continuous time, but microseconds-equivalent in spectral-information content.

The differences are statistical efficiency and bookkeeping, not information content:

-   WE pays for RiteWeight reweighting; the long trajectory does not.
-   WE provides natural snapshot structure for DR-DMD; the long trajectory does not.
-   WE concentrates statistical effort in transition regions; the long trajectory undersamples them.
-   Multiple independent WE runs are statistically independent; chunks of a single long trajectory are not.

All favor pooled WE for biomolecular folding-scale problems; the only reason to prefer a long trajectory is when the system's slowest timescale is short enough that ergodic sampling is feasible.

### 4.14 Destructive interference: when pooling attenuates modes

The pooling discussion has so far emphasized benefits. There is, however, a genuine and underappreciated risk that pooling can attenuate specific modes in the pooled DR-DMD signal if runs' excitation coefficients for those modes have opposite signs.

#### 4.14.1 The linearity of pooled excitation

Each run $r$ contributes its initial distribution $p^{(r)}(x, 0)$ to the pooled initial distribution

$$
p_{\mathrm{pool}}(x, 0) \;=\; \sum_r \alpha_r\, p^{(r)}(x, 0).
$$

The pooled excitation coefficient for mode $i$ is linear in the run-level excitation coefficients:

$$
c_i^{\mathrm{pool}} \;=\; \mathbb{E}_{p_{\mathrm{pool}}(\cdot,0)}[\psi_i] - \mathbb{E}_\pi[\psi_i] \;=\; \sum_r \alpha_r\, c_i^{(r)}.
$$

This is exact; it follows from linearity of expectation and the construction of $p_{\mathrm{pool}}$ as a mixture.

Two consequences:

(i) **Pooling cannot manufacture excitation not present in any individual run.** If every $c_i^{(r)} = 0$, then $c_i^{\mathrm{pool}} = 0$. The pooled corpus inherits the modes its constituent runs excited.

(ii) **Pooling can produce destructive cancellation** when individual runs' excitation coefficients have opposite signs, because the weighted average can be smaller in magnitude than the largest individual $|c_i^{(r)}|$. In the worst case, $c_i^{\mathrm{pool}} = 0$ even when every $|c_i^{(r)}|$ is substantial.

#### 4.14.2 A concrete destructive case

Take a symmetric equilibrium $\pi_F = \pi_U = 0.5$ for clean arithmetic. The eigenfunction normalization gives $\psi_2^F = -1$, $\psi_2^U = +1$.

Run A is seeded 80% folded / 20% unfolded:

$$
c_2^{(A)} \;=\; 0.8(-1) + 0.2(+1) \;=\; -0.6.
$$

Run B is seeded 20% folded / 80% unfolded:

$$
c_2^{(B)} \;=\; 0.2(-1) + 0.8(+1) \;=\; +0.6.
$$

Equal-weighted pooling:

$$
c_2^{\mathrm{pool}} \;=\; 0.5(-0.6) + 0.5(+0.6) \;=\; 0.
$$

The pooled corpus has zero $\psi_2$ signal in DR-DMD even though both individual runs have substantial signal.

Cancellation can also be partial. Two runs with coefficients $c_2^{(A)} = -0.6$ and $c_2^{(B)} = +0.4$ pool to

$$
c_2^{\mathrm{pool}} \;=\; 0.5(-0.6) + 0.5(+0.4) \;=\; -0.1,
$$

a sixfold reduction in signal that may push the mode below the noise floor for spectral recovery.

#### 4.14.3 The signal-to-noise structure of pooling

The structural property — that the noise in $\hat{\mathbf{R}}_{jk}$ at reference point $x_j$ scales as $\sigma_{\mathrm{noise}}(x_j)/\hat{\pi}(x_j)$ — combines with the linearity of pooled excitation to give the precise SNR criterion for whether mode $i$ is recoverable from a pooled corpus.

##### 4.14.3.1 The signal in the pooled density-ratio matrix

After pooling, the density-ratio matrix retains the spectral form

$$
\mathbf{R}_{jk}^{\mathrm{pool}} \;=\; \sum_{i \geq 2} c_i^{\mathrm{pool}}\, \lambda_i^k\, \psi_i(x_j) + \text{noise},
$$

with $c_i^{\mathrm{pool}} = \sum_r \alpha_r c_i^{(r)}$. The contribution of mode $i$ to the signal at reference point $x_j$ and iteration $k$ has magnitude

$$
S_i(x_j, k) \;=\; \bigl|c_i^{\mathrm{pool}}\bigr| \cdot |\lambda_i|^k \cdot \bigl|\psi_i(x_j)\bigr|.
$$

Each factor has a structurally distinct origin: $|c_i^{\mathrm{pool}}|$ depends on pooling, $|\lambda_i|^k$ on the iteration index, $|\psi_i(x_j)|$ on the reference point.

##### 4.14.3.2 The noise in the pooled density-ratio matrix

By, the per-cell noise standard deviation of the pooled density-ratio estimator is

$$
\sigma_{\mathrm{noise}}^{\mathbf{R},\mathrm{pool}}(x_j) \;\approx\; \frac{\sigma_{\mathrm{noise}}^{\mathrm{pool}}(x_j)}{\hat{\pi}_{\mathrm{pool}}(x_j)}
$$

where $\sigma_{\mathrm{noise}}^{\mathrm{pool}}(x_j)$ is the kernel-density-estimator standard deviation of the pooled $\hat{p}_{\mathrm{pool}}$ at the reference point. For $N_{\mathrm{runs}}$ independent runs with comparable per-run sample sizes contributing to the kernel sum at $x_j$, the effective sample size scales as $N_{\mathrm{runs}}$, and the per-cell noise scales as

$$
\sigma_{\mathrm{noise}}^{\mathrm{pool}}(x_j) \;\sim\; \frac{\sigma_{\mathrm{noise}}^{(\mathrm{single})}(x_j)}{\sqrt{N_{\mathrm{runs}}}}.
$$

The denominator $\hat{\pi}_{\mathrm{pool}}$ is similarly improved by pooling, but the leading $1/\hat{\pi}$ scaling is invariant: more sampling improves $\hat{\pi}$ precision but does not change its mean value.

##### 4.14.3.3 The SNR criterion

The SNR for detecting mode $i$ at reference point $x_j$ and iteration $k$ is the ratio of signal to noise:

$$
\mathrm{SNR}_i(x_j, k) \;=\; \frac{S_i(x_j, k)}{\sigma_{\mathrm{noise}}^{\mathbf{R},\mathrm{pool}}(x_j)} \;\approx\; \frac{\bigl|c_i^{\mathrm{pool}}\bigr|\, |\lambda_i|^k\, \bigl|\psi_i(x_j)\bigr|\, \hat{\pi}_{\mathrm{pool}}(x_j)}{\sigma_{\mathrm{noise}}^{\mathrm{pool}}(x_j)}.
$$

Substituting the $N_{\mathrm{runs}}$ scaling of the noise:

$$
\mathrm{SNR}_i(x_j, k) \;\sim\; \frac{\bigl|c_i^{\mathrm{pool}}\bigr| \cdot \sqrt{N_{\mathrm{runs}}}}{\sigma_{\mathrm{noise}}^{(\mathrm{single})}(x_j)/[\hat{\pi}(x_j) |\psi_i(x_j)| |\lambda_i|^k]}.
$$

The right-hand denominator is independent of pooling; the dependence on pooling is concentrated in the numerator factor $|c_i^{\mathrm{pool}}| \cdot \sqrt{N_{\mathrm{runs}}}$.

##### 4.14.3.4 The pooling break-even criterion

Define the cancellation factor

$$
\eta_i \;=\; \frac{\bigl|c_i^{\mathrm{pool}}\bigr|}{\max_r \bigl|c_i^{(r)}\bigr|} \;\in\; [0, 1].
$$

$\eta_i = 1$ means no cancellation (the pooled signal equals the strongest individual run's signal, up to weighting), $\eta_i = 0$ means total cancellation. Pooling $N_{\mathrm{runs}}$ runs is net-beneficial for mode $i$ when the noise reduction exceeds the signal attenuation:

$$
\sqrt{N_{\mathrm{runs}}} \cdot \eta_i \;>\; 1.
$$

Equivalently, pooling is net-beneficial when

$$
\eta_i \;>\; \frac{1}{\sqrt{N_{\mathrm{runs}}}}.
$$

This is the break-even criterion. For $N_{\mathrm{runs}} = 4$, pooling helps if $\eta_i > 1/2$; for $N_{\mathrm{runs}} = 9$, pooling helps if $\eta_i > 1/3$; for $N_{\mathrm{runs}} = 16$, pooling helps if $\eta_i > 1/4$. Pooling more runs raises the bar on how badly the cancellation can hurt before pooling becomes net-negative, but it does not eliminate the risk: a mode with $\eta_i$ close to zero is below threshold regardless of $N_{\mathrm{runs}}$.

##### 4.14.3.5 Worked example

For the destructive case: Run A has $c_2^{(A)} = -0.6$, Run B has $c_2^{(B)} = +0.6$, equal-weighted pooling gives $c_2^{\mathrm{pool}} = 0$, so $\eta_2 = 0$. Pooling A and B is catastrophic for $\psi_2$: the signal is gone regardless of how much noise reduction $\sqrt{N_{\mathrm{runs}}} = \sqrt{2}$ provides. Adding more bidirectional runs to the pool helps only if their excitation coefficients consistently have one sign — adding a third run with $c_2^{(C)} = -0.6$ gives $c_2^{\mathrm{pool}} = (-0.6 + 0.6 - 0.6)/3 = -0.2$ and $\eta_2 = 0.2/0.6 = 1/3$, which is below the $1/\sqrt{3} \approx 0.577$ break-even for three runs, so pooling all three is still net-negative for $\psi_2$ relative to using Run C alone.

For the partial-cancellation case: $c_2^{(A)} = -0.6$, $c_2^{(B)} = +0.4$ gives $c_2^{\mathrm{pool}} = -0.1$ and $\eta_2 = 0.1/0.6 = 1/6$. The break-even threshold for $N_{\mathrm{runs}} = 2$ is $1/\sqrt{2} \approx 0.707$. Since $1/6 \ll 1/\sqrt{2}$, pooling is net-negative for $\psi_2$ — Run A alone gives better $\psi_2$ recovery than Runs A and B pooled.

##### 4.14.3.6 The asymmetry that makes the criterion practically tractable

The break-even criterion $\eta_i > 1/\sqrt{N_{\mathrm{runs}}}$ can be checked *before* spending the GPU time on pooled analysis, because the excitation coefficients $c_i^{(r)}$ are functionals of the initial walker distributions and the eigenfunctions $\psi_i$. The eigenfunctions are estimated from individual-run analyses or from a preliminary pooled analysis; the initial walker distributions are known from the simulation protocol. The cancellation factor $\eta_i$ is therefore a *pre-pooling diagnostic* that can be computed from each candidate run set without re-executing the spectral analysis.

The operational protocol this enables: estimate $\hat{\psi}_i^{(r)}$ from each individual run (or from a preliminary pool), compute $\hat{c}_i^{(r)} = \sum_j w_j^{(0,r)} \hat{\psi}_i^{(r)}(x_j^{(0,r)})$ for each run, compute the pooled $\hat{c}_i^{\mathrm{pool}}$ for the candidate run set, and check $\hat{\eta}_i$ against $1/\sqrt{N_{\mathrm{runs}}}$. If the diagnostic passes, run the pooled analysis; if it fails, identify which runs are contributing destructive cancellation and either re-weight them, exclude them from the pool, or commission additional runs with seeding chosen to restore $\eta_i$ above threshold.

#### 4.14.4 The asymmetry between DR-DMD and kEDMD

A crucial structural feature: **kEDMD on the pooled corpus does not suffer from destructive interference.** kEDMD aggregates trajectory pairs across all runs and reweights to $\pi$; each pair contributes the same to covariance matrices regardless of which run it came from. There is no sign that depends on run-level seeding. Pooling for kEDMD is pure noise reduction.

This asymmetry is why the three-way cross-validation architecture is essential. When DR-DMD attenuates a specific mode via destructive pooling interference, kEDMD on the same pooled corpus does not. Disagreement between DR-DMD and kEDMD eigenvalues at a specific mode is diagnostic of destructive pooling interference, and the operational response is to inspect run-level $c_i^{(r)}$ coefficients and either re-pool with different weights or accept kEDMD as the primary estimate for the affected mode.

#### 4.14.5 What attenuation does and does not do

Destructive-interference attenuation reduces amplitude $c_i^{\mathrm{pool}}$ but does not affect temporal structure $\lambda_i^k$ or spatial structure $\psi_i(x_j)$. The pooled density-ratio matrix retains the form

$$
\mathbf{R}_{jk}^{\mathrm{pool}} \;=\; \sum_{i \geq 2} c_i^{\mathrm{pool}}\, \lambda_i^k\, \psi_i(x_j) + \text{noise},
$$

with each $c_i^{\mathrm{pool}}$ possibly reduced. Eigenvalues are still encoded in temporal decay; eigenfunctions in spatial structure. What changes is the SNR per mode: a mode whose amplitude has been cancelled below the noise floor is no longer recoverable in practice, even though the temporal and spatial information is nominally still present.

#### 4.14.6 Coverage versus equilibrium-resemblance

A final precision on a point that can be intuitively confusing. Pooling broader coverage at each iteration index might seem to imply that the pooled snapshot $\hat{p}_{\mathrm{pool}}(\cdot, t_k)$ "looks closer to $\pi$" than any individual run's snapshot, and that this should reduce the density-ratio signal. This intuition is partly right and partly wrong.

Coverage of $\mathrm{supp}(\pi)$ and equilibrium-resemblance are different conditions:

-   **Coverage** means every region where $\pi > 0$ has walkers in it (a support condition).
-   **Equilibrium-resemblance** means the proportions of walker mass across regions match the proportions of $\pi$ (a distributional condition).

Pooling runs with diverse initial conditions but *similar* seeding ratios produces broad coverage and *the same* proportions as any individual run. The pooled snapshot is no closer to $\pi$ than any individual snapshot. No signal loss.

Pooling runs with *different* seeding ratios that average to $\pi$ produces broad coverage and proportions close to $\pi$ --- which is exactly the destructive-interference scenario. There is signal loss, but the loss is the cancellation of $c_i^{\mathrm{pool}}$ rather than an "averaging to equilibrium" effect at the snapshot level. The same cancellation would happen even if the pooled snapshot looked nothing like $\pi$ at any specific iteration.

The relevant quantity for spectral recovery is always the projection

$$
c_i^{\mathrm{pool}} \;=\; \bigl\langle p_{\mathrm{pool}}(\cdot, 0)/\pi - 1,\; \psi_i \bigr\rangle_\pi.
$$

Whether the snapshot at $t = 0$ "looks like $\pi$" in some overall visual sense is irrelevant. Whether it projects onto each $\psi_i$ with substantial $c_i^{\mathrm{pool}}$ is what matters.

> **Caveat (heterogeneous coverage).** The break-even criterion $\eta_i > 1/\sqrt{N_{\mathrm{runs}}}$ is derived under the idealization that the per-cell kernel-density-estimator noise scales as $1/\sqrt{N_{\mathrm{runs}}}$ across pooled runs with comparable per-run effective sample sizes at the reference points of interest. When runs have heterogeneous coverage --- e.g. a unidirectional unfolded run contributes much more to $\hat{p}$ at unfolded reference points than a unidirectional folded run does --- the noise reduction is spatially nonuniform and the criterion becomes a per-reference-point quantity rather than a global one. The global criterion should then be read as a screening tool that may be conservative in some regions and permissive in others.

> **Caveat (interaction with the run-axis bootstrap).** The uncertainty protocol's Dirichlet bootstrap on the run axis treats runs as exchangeable resampling units. The destructive-interference structure implies that some pairs of runs are *compatible* (consistent-sign $c_i^{(r)}$) and others are not, so the bootstrap distribution may misrepresent the relevant sampling variability when the run set has heterogeneous excitation patterns: it can under-report uncertainty when resamples weight destructively-cancelling runs equally and over-report it when resamples preferentially exclude one of two cancelling runs. The effect is expected to be second-order but should be checked on the alanine-dipeptide validation system.

### 4.15 Operational recommendations

#### 4.15.1 Run-set design for the UUCG campaign

The framework above implies specific operational choices:

**Continue the existing bidirectional WE configuration to completion.** The 50-iteration data already supports the conceptual picture; running to the planned iteration count completes the corpus for the dominant slow mode.

**Set up a second run with substantially more diverse bstates.** The clustered seed problem (folded seeds at one tight position, unfolded seeds at another tight position) limits higher-mode excitation. Diversification within each basin across the two runs is the lever for unlocking higher-mode recovery in the pooled corpus. Strategies:

-   Pull bstates from different time windows of source trajectories rather than random sampling from one trajectory segment.
-   Use a diversity criterion (pairwise RMSD threshold) on candidate frames.
-   Supplement with frames from higher-temperature simulations if available.

**Consider adding one unidirectional unfolded-side run with diverse bstates.** This strengthens $\psi_2$ excitation in the pooled corpus through the asymmetric-equilibrium mechanism:

$$
c_2^{\mathrm{pool}} \;=\; \alpha_{\mathrm{bidir}}\, c_2^{\text{bidir}} + \alpha_{\mathrm{uni-U}}\, c_2^{\text{uni-U}}
$$

with both terms of consistent sign (positive in the $\pi_F = 0.9$ case), so the sum strengthens. The unidirectional run's walkers also improve within-basin coverage of unfolded substructure.

**Do not add unidirectional folded-side runs** unless there is a specific reason. They partially cancel $\psi_2$ in the pooled corpus per.

**Vary the progress coordinate across runs.** Different PCs drive walkers into different regions and reduce coordinate-driven sampling bias.

**Seed transition-region intermediates** (e.g., partially-formed-stem structures) in at least one run if such structures can be constructed. These directly excite the eigenfunctions in their nodal regions and improve spatial resolution of $\psi_2$ near the barrier, addressing the asymmetric-advance observation.

**For misfolded conformations**, distinguish:

(i) Near-native substates that broaden the folded seed cloud and help higher modes;

(ii) Topologically distinct misfolds that constitute scientific scope expansion and require Boltzmann-weight justification;

(iii) High-energy near-transition structures that improve barrier-region resolution.

Category (iii) has the cleanest pipeline rationale and is the safest addition.

#### 4.15.2 Alanine dipeptide validation

**Use unidirectional or single-basin seeding as the primary validation protocol.** Dipeptide barriers are low enough that coverage is not the binding constraint, and matching the PyEMMA ground truth's typical generation protocol gives the cleanest comparison.

**Optionally run one bidirectional dipeptide WE as a controlled comparison artifact.** This provides a controlled within-system comparison of single-seed vs. bidirectional DR-DMD recovery against PyEMMA ground truth. The honest answer (no bias, because RiteWeight corrects basin populations exactly) is cheaper to demonstrate than to argue.

#### 4.15.3 Pooled-corpus analysis protocol

**After each new WE run completes**, re-run the pooled DR-DMD analysis with the previous corpus plus the new run, and check whether eigenvalues for $\psi_2$ and $\psi_3$ strengthen or weaken relative to the previous pool. Strengthening confirms compatibility; weakening flags destructive interference and triggers inspection of run-level $c_i^{(r)}$.

**Cross-validate DR-DMD against kEDMD on each pooled corpus.** Disagreement at a specific mode flags destructive pooling interference for that mode. Use kEDMD as the fallback estimate for any mode where DR-DMD fails the agreement check.

**Track pooled excitation coefficients** $c_i^{\mathrm{pool}}$ as a corpus-level diagnostic. A summary table reporting $c_2^{(r)}, c_3^{(r)}, c_4^{(r)}$ for each run and $c_2^{\mathrm{pool}}, c_3^{\mathrm{pool}}, c_4^{\mathrm{pool}}$ for the pool is the natural representation. The pooling break-even criterion of — that pooling is net-beneficial for mode $i$ when

$$
\eta_i \;=\; \frac{\bigl|c_i^{\mathrm{pool}}\bigr|}{\max_r \bigl|c_i^{(r)}\bigr|} \;>\; \frac{1}{\sqrt{N_{\mathrm{runs}}}}
$$

— is the operational test. Modes that fail the criterion are better recovered from a subset of the pool (typically the strongest-excitation single run) or by re-weighting the pool to suppress destructive contributions; modes that pass the criterion benefit from pooling and should be reported from the pooled analysis.
