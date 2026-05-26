#!/usr/bin/env python3
"""
02_analyze_bstate_pcoords.py
============================================================================
Assemble the full progress-coordinate table for every WESTPA basis state in
03_bstates_final/, characterise the diversity of the seed set, and recommend
concrete bstate subsets for the bidirectional-coverage and perturbation run
schemes described in README.md and WEeDS_Background_reorganized.md.

DESIGN PRINCIPLE (per user instruction):
    Seed-set RECOMMENDATIONS are derived ONLY from the measured progress
    coordinates, never from the folded/intermediate/unfolded directory labels.
    Basins ("folded side" / "unfolded side") and the transition band are
    assigned by k-means + a continuous folded<->unfolded margin in standardized
    pcoord space. The original labels are reported alongside for reference and
    for a label-vs-pcoord cross-check, but they do not drive any selection.

Run from the bstates repo root.

INPUT  (produced by 01_compute_raw_distances.sh):
    04_analyze/raw/struct_X/pcoord_candidates.dat   rmsd_global/stem/loop, chi_G9, Rg, d_e2e
    04_analyze/raw/struct_X/loop_hbonds_raw.dat     6 loop H-bond distances
    04_analyze/raw/struct_X/stem_hbonds_raw.dat     14 stem WC H-bond distances
    04_analyze/raw/struct_X/mindist.dat             G1-C14 native-contact min distance
    03_bstates_final/struct_mapping.log             struct_X -> bstate_id (reference only)

OUTPUT (written to 04_analyze/results/):
    bstate_pcoords.csv        one row per bstate, all pcoords + data-driven basin
    bstate_pcoords.dat        same table, whitespace-aligned text
    recommendations.csv       recommended struct sets per scheme
    analysis.log              human-readable log: data + diversity + recommendations
    fig_*.png                 diversity / pcoord figures

Q-value conventions replicate compute_Q.1.py exactly (4 loop H-bonds, 14 stem
H-bonds, 3.5 A WC cutoff, 3.7 A base-phosphate cutoff, sigmoid + hard switches).
============================================================================
"""
import os
import sys
import contextlib
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
HERE      = os.path.dirname(os.path.abspath(__file__))      # .../bstates (script lives at repo root)
RAW_ROOT  = os.path.join(HERE, "04_analyze", "raw")         # per-bstate cpptraj output (incl. NMR ref)
RESULTS   = os.path.join(HERE, "04_analyze", "results")     # tables, figures, log
MAP_LOG   = os.path.join(HERE, "03_bstates_final", "struct_mapping.log")

STEM_CUTOFF = 3.5   # A, WC H-bond cutoff
LOOP_CUTOFF = 3.5   # A, loop H-bond cutoff
BPHO_CUTOFF = 3.7   # A, base-phosphate cutoff

# Features defining the folded<->unfolded basin axis (the slow-mode-like
# coordinate). Loop metrics are deliberately excluded here: the UUCG loop stays
# largely formed while the stem melts, so it is a poor basin discriminator.
BASIN_FEATS = ["rmsd_global", "rmsd_stem", "mindist", "Q_stem_sig",
               "n_stem_bp", "rog", "d_e2e"]
# Features used for within-group diversity (farthest-point sampling). Loop
# metrics ARE included here: within-basin loop/substructure spread feeds the
# higher modes (S4.11.6).
DIV_FEATS   = ["rmsd_global", "mindist", "rmsd_stem", "rmsd_loop",
               "Q_stem_sig", "Q_loop_sig", "rog", "d_e2e"]

BASIN_COLORS = {"folded_side": "#1b7837", "transition": "#762a83", "unfolded_side": "#b35806"}
BASIN_ORDER  = ["folded_side", "transition", "unfolded_side"]
CAT_ORDER    = ["folded", "intermediate", "unfolded"]
CAT_COLORS   = {"folded": "#1b7837", "intermediate": "#762a83", "unfolded": "#b35806"}


# ---------------------------------------------------------------------------
def sigmoid_switch(d, d0, n=6, m=12):
    """
    Smooth sigmoid switching function for native contacts.
    Returns ~1 when d < d0, ~0 when d >> d0.
    Standard form from Best, Hummer, Eaton (2013):

        Q_ij = (1 - (d/d0)^n) / (1 - (d/d0)^m)

    A uniform 3.5 A (3.7 A for base-phosphate) cutoff for d0 is simple and
    standard, and matches compute_Q.1.py.
    """
    r = d / d0
    # Avoid numerical issues at r = 1
    r = np.where(np.abs(r - 1.0) < 1e-6, 1.0 + 1e-6, r)
    return (1.0 - r**n) / (1.0 - r**m)


def read_row(path):
    """Read the single data row from a one-frame cpptraj .dat file."""
    arr = np.loadtxt(path, comments="#")
    return np.atleast_1d(arr)


def load_mapping(path):
    """struct_index (int) -> (bstate_id, original_label) ; reference only."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            sidx = int(parts[0].replace("struct_", ""))
            bid  = parts[1]
            lab  = bid.rsplit("_", 1)[0]
            mapping[sidx] = (bid, lab)
    return mapping


def compute_pcoords(sidx):
    """Return a dict of all pcoords for one struct index."""
    sdir = os.path.join(RAW_ROOT, f"struct_{sidx}")

    cand = read_row(os.path.join(sdir, "pcoord_candidates.dat"))
    # cols: frame, rmsd_global, rmsd_stem, rmsd_loop, chi_G9, Rg, d_e2e
    rmsd_global, rmsd_stem, rmsd_loop, chi_G9, rog, d_e2e = cand[1:7]

    mind = read_row(os.path.join(sdir, "mindist.dat"))
    mindist = mind[1]

    loop = read_row(os.path.join(sdir, "loop_hbonds_raw.dat"))
    # cols: frame, d1..d6 ; use d1,d3,d4,d6 (4 contacts) as in compute_Q.1.py
    d_loop = loop[[1, 3, 4, 6]]
    loop_cutoffs = np.array([LOOP_CUTOFF, LOOP_CUTOFF, LOOP_CUTOFF, BPHO_CUTOFF])

    stem = read_row(os.path.join(sdir, "stem_hbonds_raw.dat"))
    d_stem = stem[1:15]   # 14 distances

    # hard contacts
    loop_hard = (d_loop <= loop_cutoffs).astype(float)
    stem_hard = (d_stem <= STEM_CUTOFF).astype(float)
    Q_loop_hard = loop_hard.sum() / 4.0
    Q_stem_hard = stem_hard.sum() / 14.0
    loop_contact_count = loop_hard.sum()

    # sigmoid contacts (smoother, matches compute_Q.1.py)
    loop_sig = np.array([sigmoid_switch(d_loop[i], loop_cutoffs[i]) for i in range(4)])
    stem_sig = sigmoid_switch(d_stem, STEM_CUTOFF)
    Q_loop_sig = loop_sig.sum() / 4.0
    Q_stem_sig = stem_sig.sum() / 14.0

    # number of stem base pairs (GC: >=2/3 hbonds; AU: 2/2)
    bp = np.zeros(5)
    bp[0] = stem_hard[0:3].sum()  >= 2   # G1-C14
    bp[1] = stem_hard[3:6].sum()  >= 2   # G2-C13
    bp[2] = stem_hard[6:9].sum()  >= 2   # C3-G12
    bp[3] = stem_hard[9:11].sum() >= 2   # A4-U11
    bp[4] = stem_hard[11:14].sum() >= 2  # C5-G10
    n_stem_bp = bp.sum()

    # G9 syn check (UUCG signature: chi in the syn window)
    g9_syn = float(-25 <= chi_G9 <= 115)

    return dict(rmsd_global=rmsd_global, rmsd_stem=rmsd_stem, rmsd_loop=rmsd_loop,
                mindist=mindist, Q_stem_sig=Q_stem_sig, Q_loop_sig=Q_loop_sig,
                Q_stem_hard=Q_stem_hard, Q_loop_hard=Q_loop_hard,
                n_stem_bp=n_stem_bp, loop_contact_count=loop_contact_count,
                chi_G9=chi_G9, g9_syn=g9_syn, rog=rog, d_e2e=d_e2e)


# ---------------------------------------------------------------------------
def fps(Z, k):
    """
    Farthest-point sampling on already-standardized rows Z (no internal
    standardization, so the caller controls the scale). Returns LOCAL indices
    (into Z) of k maximally-spread points.

    Greedy: start from the point nearest the centroid, then repeatedly add the
    point with the largest minimum distance to the already-chosen set. This is
    the diversity criterion called for in S4.11.6 / S4.15.1 (spread within a
    basin to excite higher modes).
    """
    n = Z.shape[0]
    k = min(k, n)
    if n == 0:
        return []
    start = int(np.argmin(np.linalg.norm(Z - Z.mean(0), axis=1)))
    chosen = [start]
    mind = np.linalg.norm(Z - Z[start], axis=1)
    while len(chosen) < k:
        nxt = int(np.argmax(mind))
        if nxt in chosen:
            break
        chosen.append(nxt)
        mind = np.minimum(mind, np.linalg.norm(Z - Z[nxt], axis=1))
    return chosen


def mean_pairwise_dist(Z):
    """Mean pairwise Euclidean distance of (already-standardized) rows Z.

    A scalar spread metric: higher = more diverse seed cloud.
    """
    if len(Z) < 2:
        return 0.0
    d = [np.linalg.norm(Z[i] - Z[j]) for i in range(len(Z)) for j in range(i + 1, len(Z))]
    return float(np.mean(d))


def kmeans2(X, n_init=25, iters=200, seed=0):
    """Minimal 2-cluster k-means (Lloyd) with multiple restarts. Returns labels, centroids."""
    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_init):
        c = X[rng.choice(len(X), 2, replace=False)].copy()
        lab = np.zeros(len(X), dtype=int)
        for _ in range(iters):
            d = np.linalg.norm(X[:, None, :] - c[None, :, :], axis=2)
            new = d.argmin(1)
            cc = np.array([X[new == g].mean(0) if (new == g).any() else c[g] for g in range(2)])
            if np.array_equal(new, lab) and np.allclose(cc, c):
                c = cc
                break
            lab, c = new, cc
        inertia = sum(((X[lab == g] - c[g]) ** 2).sum() for g in range(2))
        if best is None or inertia < best[0]:
            best = (inertia, lab.copy(), c.copy())
    return best[1], best[2]


# ===========================================================================
def main():
    os.makedirs(RESULTS, exist_ok=True)
    mapping = load_mapping(MAP_LOG)
    sidxs = sorted(mapping)

    pc = [compute_pcoords(s) for s in sidxs]
    bid   = [mapping[s][0] for s in sidxs]
    label = [mapping[s][1] for s in sidxs]   # reference only

    keys = ["rmsd_global", "rmsd_stem", "rmsd_loop", "mindist",
            "Q_stem_sig", "Q_loop_sig", "Q_stem_hard", "Q_loop_hard",
            "n_stem_bp", "loop_contact_count", "chi_G9", "g9_syn", "rog", "d_e2e"]
    M = {k: np.array([p[k] for p in pc]) for k in keys}
    n = len(sidxs)

    # ---- data-driven basin assignment (label-free, from pcoords only) -------
    # Two basins come from k-means(2) on the standardized basin features; the
    # continuous folded<->unfolded reaction coordinate is the centroid margin.
    Xb = np.column_stack([M[k] for k in BASIN_FEATS])
    Zb = (Xb - Xb.mean(0)) / (Xb.std(0) + 1e-12)
    klab, cent = kmeans2(Zb)
    # folded cluster = the one with lower mean rmsd_global
    rmsd_col = BASIN_FEATS.index("rmsd_global")
    folded_cluster = int(np.argmin([cent[g][rmsd_col] for g in range(2)]))
    cF, cU = cent[folded_cluster], cent[1 - folded_cluster]
    dF = np.linalg.norm(Zb - cF, axis=1)
    dU = np.linalg.norm(Zb - cU, axis=1)
    margin = (dU - dF) / (dU + dF + 1e-12)        # +1 fully folded, -1 fully unfolded
    fold_frac = 0.5 * (margin + 1.0)              # 1 folded ... 0 unfolded (reported)
    in_folded_cluster = (klab == folded_cluster)

    # fold_frac is sharply BIMODAL (a folded cloud and an unfolded cloud with an
    # empty gap between), so the k-means dividing surface itself is nearly
    # unpopulated -- the true barrier is not a thin margin band but the set of
    # partially-formed-stem structures held at an intermediate terminal MinDist.
    # Define the transition/barrier set purely from pcoords: a partly-formed stem
    # (1 <= n_stem_bp <= 4) AND an intermediate G1-C14 MinDist -- neither the
    # closed native pair (~2 A) nor a fully separated terminus (>~25 A). The
    # thresholds are read off the bimodal MinDist distribution (folded core
    # tops out ~2.1 A; the deep-unfolded cloud starts ~24 A).
    MD_LO, MD_HI = 3.0, 25.0
    is_transition = ((M["n_stem_bp"] >= 1) & (M["n_stem_bp"] <= 4)
                     & (M["mindist"] > MD_LO) & (M["mindist"] < MD_HI))

    # transition takes priority; the rest fall to their k-means basin.
    basin = np.where(is_transition, "transition",
             np.where(in_folded_cluster, "folded_side", "unfolded_side"))

    # ---- standardized diversity space (global) -----------------------------
    Xd = np.column_stack([M[k] for k in DIV_FEATS])
    Zd = (Xd - Xd.mean(0)) / (Xd.std(0) + 1e-12)

    # ---- write table -------------------------------------------------------
    cols = (["struct_index", "bstate_id", "orig_label", "basin", "fold_frac", "margin"]
            + keys)
    rows = []
    for i, s in enumerate(sidxs):
        rows.append([s, bid[i], label[i], basin[i], fold_frac[i], margin[i]]
                    + [M[k][i] for k in keys])

    csv_path = os.path.join(RESULTS, "bstate_pcoords.csv")
    with open(csv_path, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            out = []
            for v in r:
                if isinstance(v, str):
                    out.append(v)
                elif isinstance(v, (int, np.integer)):
                    out.append(f"{v}")
                else:
                    out.append(f"{v:.4f}")
            fh.write(",".join(out) + "\n")

    dat_path = os.path.join(RESULTS, "bstate_pcoords.dat")
    with open(dat_path, "w") as fh:
        fh.write(f"# {'idx':>4s} {'bstate_id':>16s} {'orig_label':>13s} {'basin':>13s} "
                 f"{'foldFr':>7s} {'margin':>7s} {'rmsdGlob':>9s} {'rmsdStem':>9s} "
                 f"{'rmsdLoop':>9s} {'mindist':>9s} {'Qstem_s':>8s} {'Qloop_s':>8s} "
                 f"{'nBP':>4s} {'chiG9':>9s} {'Rg':>8s} {'d_e2e':>8s}\n")
        for r in rows:
            d = dict(zip(cols, r))
            fh.write(f"  {d['struct_index']:>4d} {d['bstate_id']:>16s} {d['orig_label']:>13s} "
                     f"{d['basin']:>13s} {d['fold_frac']:>7.3f} {d['margin']:>7.3f} "
                     f"{d['rmsd_global']:>9.4f} {d['rmsd_stem']:>9.4f} {d['rmsd_loop']:>9.4f} "
                     f"{d['mindist']:>9.4f} {d['Q_stem_sig']:>8.4f} {d['Q_loop_sig']:>8.4f} "
                     f"{d['n_stem_bp']:>4.0f} {d['chi_G9']:>9.3f} {d['rog']:>8.4f} "
                     f"{d['d_e2e']:>8.4f}\n")

    # index helpers
    pos_of = {s: i for i, s in enumerate(sidxs)}
    idx2bid = {s: bid[pos_of[s]] for s in sidxs}
    idx2lab = {s: label[pos_of[s]] for s in sidxs}

    def members(mask):
        return [sidxs[i] for i in range(n) if mask[i]]

    def select(struct_pool, k):
        """FPS pick k structs from a pool (list of struct indices)."""
        pos = [pos_of[s] for s in struct_pool]
        loc = fps(Zd[pos], k)
        return [struct_pool[j] for j in loc]

    def spread(struct_pool):
        pos = [pos_of[s] for s in struct_pool]
        return mean_pairwise_dist(Zd[pos])

    folded_pool      = members(basin == "folded_side")
    unfolded_pool    = members(basin == "unfolded_side")
    transition_pool  = members(basin == "transition")
    unfolded_leaning = members(margin < 0.0)        # for the perturbation run

    # =======================================================================
    # SCHEME SELECTIONS  (purely pcoord-driven)
    # =======================================================================
    N_BASIN = 8
    A_folded   = select(folded_pool,   N_BASIN)
    A_unfolded = select(unfolded_pool, N_BASIN)
    B_trans    = select(transition_pool, min(6, len(transition_pool)))
    C_unfolded = select(unfolded_leaning, min(10, len(unfolded_leaning)))
    D_inter    = select(transition_pool, min(12, len(transition_pool)))

    # =======================================================================
    # FIGURES
    # =======================================================================
    def scatter_basin(ax, x, y, annotate=False, highlight=None):
        for b in BASIN_ORDER:
            m = basin == b
            ax.scatter(x[m], y[m], s=40, alpha=0.8, label=b,
                       color=BASIN_COLORS[b], edgecolor="k", linewidth=0.3)
        if highlight:
            hp = [pos_of[s] for s in highlight]
            ax.scatter(x[hp], y[hp], s=180, facecolors="none",
                       edgecolors="red", linewidths=1.6, label="selected")
        if annotate:
            for i in range(len(x)):
                ax.annotate(str(sidxs[i]), (x[i], y[i]), fontsize=5,
                            alpha=0.6, xytext=(2, 2), textcoords="offset points")

    # Fig 1: WE 2D pcoord (RMSD vs MinDist), colored by DATA-DRIVEN basin
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    scatter_basin(ax, M["rmsd_global"], M["mindist"], annotate=True)
    ax.set_xlabel("RMSD to NMR folded (Å)")
    ax.set_ylabel("G1–C14 native-contact min distance (Å)")
    ax.set_title("Basis states in the WE 2D progress coordinate\n"
                 "colored by data-driven basin (k-means + margin, labels NOT used)")
    ax.legend(title="pcoord basin")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig_pcoord_RMSD_MinDist.png"), dpi=160)
    plt.close(fig)

    # Fig 2: Q_stem vs Q_loop
    fig, ax = plt.subplots(figsize=(7, 5.5))
    scatter_basin(ax, M["Q_stem_sig"], M["Q_loop_sig"])
    ax.set_xlabel("Q_stem (fraction of 14 WC H-bonds)")
    ax.set_ylabel("Q_loop (fraction of 4 UUCG H-bonds)")
    ax.set_title("Native-contact space — the loop stays formed while the stem melts")
    ax.legend(title="pcoord basin")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig_Qstem_Qloop.png"), dpi=160)
    plt.close(fig)

    # Fig 3: per-coordinate distributions, by ORIGINAL label (diversity overview)
    panel = ["rmsd_global", "mindist", "rmsd_stem", "rmsd_loop",
             "Q_stem_sig", "Q_loop_sig", "n_stem_bp", "rog", "d_e2e"]
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    for ax, name in zip(axes.ravel(), panel):
        allv = M[name]
        bins = np.linspace(allv.min(), allv.max(), 15)
        for c in CAT_ORDER:
            m = np.array([ll == c for ll in label])
            ax.hist(allv[m], bins=bins, alpha=0.6, label=c, color=CAT_COLORS[c])
        ax.set_title(name); ax.set_ylabel("count")
    axes.ravel()[0].legend(fontsize=8)
    fig.suptitle("Distribution of each progress coordinate across the 70 basis states "
                 "(original labels)", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig_pcoord_distributions.png"), dpi=160)
    plt.close(fig)

    # Fig 4: selected seeds for Scheme B marked on the WE pcoord
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    scatter_basin(ax, M["rmsd_global"], M["mindist"],
                  highlight=A_folded + A_unfolded + B_trans)
    ax.set_xlabel("RMSD to NMR folded (Å)")
    ax.set_ylabel("G1–C14 native-contact min distance (Å)")
    ax.set_title("Scheme B recommended seeds (red rings)\n"
                 "diverse folded + unfolded + transition, chosen by farthest-point sampling")
    ax.legend(title="pcoord basin")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig_recommended_seeds.png"), dpi=160)
    plt.close(fig)

    # Fig 5: Rg vs end-to-end, sized by n_stem_bp
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for b in BASIN_ORDER:
        m = basin == b
        ax.scatter(M["rog"][m], M["d_e2e"][m], s=30 + 40 * M["n_stem_bp"][m],
                   alpha=0.75, label=b, color=BASIN_COLORS[b], edgecolor="k", linewidth=0.3)
    ax.set_xlabel("Radius of gyration (Å)")
    ax.set_ylabel("End-to-end distance G1(O6)–C14(O3') (Å)")
    ax.set_title("Global compaction (marker size ∝ number of stem base pairs)")
    ax.legend(title="pcoord basin")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "fig_Rg_e2e.png"), dpi=160)
    plt.close(fig)

    # =======================================================================
    # recommendations CSV
    # =======================================================================
    rec_csv = os.path.join(RESULTS, "recommendations.csv")
    with open(rec_csv, "w") as fh:
        fh.write("scheme,role,struct_index,bstate_id,orig_label,basin,fold_frac\n")
        def w(scheme, role, structs):
            for s in structs:
                fh.write(f"{scheme},{role},{s},{idx2bid[s]},{idx2lab[s]},"
                         f"{basin[pos_of[s]]},{fold_frac[pos_of[s]]:.3f}\n")
        w("A_bidirectional", "folded_seed",   A_folded)
        w("A_bidirectional", "unfolded_seed", A_unfolded)
        w("B_bidir+transition", "folded_seed",     A_folded)
        w("B_bidir+transition", "unfolded_seed",   A_unfolded)
        w("B_bidir+transition", "transition_seed", B_trans)
        w("C_unidirectional_unfolded", "unfolded_seed", C_unfolded)
        w("D_intermediate_only", "transition_seed", D_inter)

    # =======================================================================
    # LOG
    # =======================================================================
    log_path = os.path.join(RESULTS, "analysis.log")
    with open(log_path, "w") as fh, contextlib.redirect_stdout(fh):
        print("=" * 78)
        print("BASIS-STATE PROGRESS-COORDINATE ANALYSIS  —  UUCG hairpin (2KOC), 14 nt")
        print("=" * 78)
        print(f"Basis states analysed : {n}")
        print(f"Reference structure   : 2KOCFolded_NMR (folded native)")
        print(f"WE 2D progress coord  : (RMSD-to-NMR, G1-C14 native-contact MinDist)")
        print("Pcoord defs follow 1_we_pcoord_distances.in / get_pcoord.cpptraj;")
        print("Q-values follow compute_Q.1.py (sigmoid + 3.5/3.7 A hard cutoffs).")
        print()
        print("SELECTION POLICY: recommendations use ONLY measured pcoords. Basins are")
        print("assigned by k-means(2) + a continuous folded<->unfolded margin in")
        print("standardized pcoord space. Directory labels are shown for reference only.")
        print()

        # data-driven basin vs original label cross-check
        print("-" * 78)
        print("DATA-DRIVEN BASIN ASSIGNMENT  (label-free; from pcoords only)")
        print("-" * 78)
        print("  basins  = k-means(2) on standardized basin features")
        print(f"  transition = partly-formed stem (1<=nBP<=4) AND {MD_LO:.0f} < MinDist < {MD_HI:.0f} A")
        print(f"  folded_side   : {len(folded_pool):3d}   unfolded_side : {len(unfolded_pool):3d}"
              f"   transition : {len(transition_pool):3d}")
        print()
        print("  Cross-tab  (rows = original label, cols = pcoord basin):")
        print(f"    {'':<14s}{'folded_side':>13s}{'transition':>13s}{'unfolded_side':>15s}")
        for c in CAT_ORDER:
            r = {b: sum(1 for i in range(n) if label[i] == c and basin[i] == b)
                 for b in BASIN_ORDER}
            print(f"    {c:<14s}{r['folded_side']:>13d}{r['transition']:>13d}{r['unfolded_side']:>15d}")
        print()
        print("  => Many 'intermediate'-labelled seeds are pcoord-folded or pcoord-unfolded;")
        print("     they broaden the folded / unfolded seed clouds exactly as intended,")
        print("     and only the genuinely in-between ones form the transition band.")
        print()

        print("-" * 78)
        print("SUMMARY STATISTICS  (mean +/- std  [min, max])  by data-driven basin")
        print("-" * 78)
        for name in ["rmsd_global", "rmsd_stem", "mindist", "Q_stem_sig",
                     "Q_loop_sig", "n_stem_bp", "rog", "d_e2e"]:
            v = M[name]
            print(f"  {name:<14s} ALL          {v.mean():8.3f} +/- {v.std():7.3f}  "
                  f"[{v.min():7.3f}, {v.max():7.3f}]")
            for b in BASIN_ORDER:
                vb = v[basin == b]
                if len(vb):
                    print(f"  {'':<14s} {b:<13s} {vb.mean():8.3f} +/- {vb.std():7.3f}  "
                          f"[{vb.min():7.3f}, {vb.max():7.3f}]")
            print()

        print("-" * 78)
        print("SEED-CLOUD DIVERSITY  (mean pairwise dist in standardized pcoord space)")
        print("-" * 78)
        print("  Higher = more diverse. S4.11.6: within-basin spread is what excites the")
        print("  higher modes psi_3, psi_4; clustered seeds project poorly onto them.")
        for b in BASIN_ORDER:
            pool = members(basin == b)
            print(f"    {b:<14s}: {spread(pool):.3f}   (n={len(pool)})")
        print()

        print("=" * 78)
        print("RECOMMENDATIONS  (all subsets chosen by farthest-point sampling on pcoords)")
        print("=" * 78)
        print("""
Grounding (README.md + WEeDS_Background_reorganized.md S4.6, S4.10, S4.11, S4.14, S4.15.1):

 * Bidirectional 50/50 seeding is the DEFAULT. RiteWeight corrects the
   non-equilibrium basin populations downstream, so the FES is invariant to the
   seeding ratio provided every basin is seeded (S4.6).
 * The binding constraint is COVERAGE of supp(pi), not observed crossings
   (S4.10). Coverage failure = identifiability failure (a mode is unrecoverable
   on the uncovered region); excitation loss = wider CIs only (S4.11.4).
 * Within-basin DIVERSITY is the lever on higher modes (S4.11.6) -> FPS spread.
 * Unidirectional UNFOLDED (rare-basin) seeding maximally excites psi_2 because
   |psi_2| is concentrated in the rare basin (S4.11.1-2). A folded-side
   unidirectional run partially CANCELS psi_2 in the pool (S4.14) -> avoid.
""")

        def show(role, structs):
            print(f"  {role}  (n={len(structs)})")
            print(f"    structs : {structs}")
            print(f"    bstates : {[idx2bid[s] for s in structs]}")
            print(f"    labels  : {[idx2lab[s] for s in structs]}")
            print(f"    fold_fr : {[round(float(fold_frac[pos_of[s]]),2) for s in structs]}")
            if len(structs) >= 2:
                print(f"    spread  : {spread(structs):.3f}")
            print()

        print("-" * 78)
        print("SCHEME A — Bidirectional 50/50 coverage run (default)")
        print("-" * 78)
        print(f"  {N_BASIN} folded-side + {N_BASIN} unfolded-side seeds, diverse within each basin.")
        print()
        show("folded-side seeds",   A_folded)
        show("unfolded-side seeds", A_unfolded)

        print("-" * 78)
        print("SCHEME B — Bidirectional + transition-region seeds")
        print("-" * 78)
        print("  Scheme A plus partially-formed-stem seeds spanning the barrier; sharpens")
        print("  psi_2's node and the FES along the pathway (S4.10.1, S4.15.1 cat. iii).")
        print()
        show("folded-side seeds",   A_folded)
        show("unfolded-side seeds", A_unfolded)
        show("transition seeds",    B_trans)

        print("-" * 78)
        print("SCHEME C — Unidirectional UNFOLDED perturbation run (Strategy B, S4.9.6)")
        print("-" * 78)
        print("  Diverse far-from-pi unfolded seeds. Strengthens psi_2 in the pool via the")
        print("  asymmetric-equilibrium mechanism; pool WITH A/B and check break-even.")
        print()
        show("unfolded perturbation seeds", C_unfolded)

        # ---- the intermediate-only-run analysis ----
        print("=" * 78)
        print("SCHEME D — Intermediate / barrier-ONLY run: what it does mathematically")
        print("=" * 78)
        print(f"""
Candidate seeds (transition band, FPS):
    structs : {D_inter}
    labels  : {[idx2lab[s] for s in D_inter]}
    fold_fr : {[round(float(fold_frac[pos_of[s]]),2) for s in D_inter]}
    spread  : {spread(D_inter):.3f}

For any slow mode i>=2 the run's excitation coefficient is (S4.9.2)
        c_i^(r) = E_{{p^(r)(.,0)}}[psi_i]            (since E_pi[psi_i] = 0).
psi_2 is the folded/unfolded contrast mode: it has a NODE (psi_2 ~ 0) in the
transition region and large magnitude in the basins (|psi_2^U|=3, |psi_2^F|=1/3
for pi_F=0.9; S4.11.1).

(1) psi_2, standalone.  Seeds sitting on the node give E[psi_2] ~ 0, so
    c_2^(int) ~ 0.  A barrier-only run is the WEAKEST possible exciter of the
    dominant slow mode -- it never lights up psi_2's temporal decay, so lambda_2
    is essentially unreadable from this run alone. (If the transition pool spans
    near-folded to near-unfolded, its -psi_2 and +psi_2 seeds average toward 0:
    same conclusion, by internal cancellation.)

(2) Higher modes psi_3, psi_4, standalone.  These distinguish substructure /
    pathways and have their large-magnitude regions and nodes INSIDE the
    transition / within-basin structure (S4.11.6). A barrier-seeded distribution
    can carry substantial c_3, c_4. The intermediate-only run is therefore a
    psi_2-BLIND, higher-mode-TARGETED probe -- the run you would design to chase
    psi_3/psi_4 that endpoint seeding excites poorly.

(3) Coverage, standalone.  Barrier walkers populate pi-hat and p-hat at
    transition reference points, cutting the 1/pi-hat noise blow-up there
    (S4.14.3.2) and sharpening the SPATIAL resolution of psi_2's node and the FES
    along the pathway (S4.10.1). This is eigenfunction-SHAPE quality, distinct
    from the eigenvalue/amplitude signal in (1).

POOLED WITH AN A/B RUN (linearity S4.14.1; break-even S4.14.3.4):
    c_2^pool = a*c_2^(int) + (1-a)*c_2^(AB) ~ (1-a)*c_2^(AB).
    Equal weight a=0.5, with c_2^(AB) ~ 1.33 (the S4.11.3 50/50 value) and
    c_2^(int) ~ 0:   c_2^pool ~ 0.665,  eta_2 = |c_2^pool|/max_r|c_2^(r)| = 0.5.
    Break-even for N=2 is 1/sqrt(2) ~ 0.707.  Since 0.5 < 0.707, equal-weight
    pooling is NET-NEGATIVE for psi_2: you recover psi_2 BETTER from the A/B run
    alone. The harm is DILUTION (adding zero-psi_2-signal mass), not sign
    cancellation. (Numbers are the doc's pi_F=0.9 illustration; the qualitative
    result -- c_2^int ~ 0, dilution, sub-break-even -- holds generally.)
    Remedies: (a) DOWN-WEIGHT the intermediate run (small a) so it barely dilutes
    c_2 while still giving barrier coverage + higher-mode signal; (b) keep it as a
    SEPARATE corpus analysed for psi_3/psi_4; (c) only merge into the psi_2 pool
    once N_runs is large enough that sqrt(N)*eta_2 > 1.

POOLED WITH A PERTURBATION (uni-U) RUN INSTEAD:
    uni-U has c_2 ~ +3 (SAME sign as A/B) and the intermediate run ~0, so there
    is no psi_2 sign cancellation between them. BUT {{uni-U + intermediate}} ALONE
    leaves the FOLDED basin uncovered -> identifiability failure on the folded
    side (S4.11.2: coverage failure is unrecoverable, not merely noisy). So an
    intermediate+perturbation pool is sound ONLY when folded coverage is supplied
    by a third (A/B) run. The productive 3-run pool is:
        {{ A/B bidirectional  (coverage of BOTH basins, c_2 ~ +1.33),
          uni-U perturbation (c_2 ~ +3, same sign -> STRENGTHENS pooled c_2),
          intermediate       (c_2 ~ 0, down-weighted: barrier coverage + psi_3/psi_4) }}.
    With same-sign c_2 from the first two, eta_2 stays high and rises with N_runs;
    the intermediate run rides along for higher-mode + barrier benefit without
    dragging psi_2 below break-even, provided it is not over-weighted.
    For psi_3/psi_4 the intermediate run's sign is NOT guaranteed to match the
    others, so check eta_3, eta_4 run-by-run after pooling (S4.15.3) before
    trusting pooled higher-mode estimates.

BOTTOM LINE: an intermediate-only run is scientifically interesting but is a
HIGHER-MODE / BARRIER-RESOLUTION tool, NOT a psi_2 tool. Pool it down-weighted
(or analyse it separately); never let it dilute the basin-seeded c_2 at equal
weight; and always pair it with runs that cover BOTH basins so coverage -- the
binding constraint (S4.10) -- is never lost.
""")
        print("-" * 78)
        print("DO NOT: add a unidirectional folded-side run (cancels psi_2, S4.14).")
        print("NEXT  : vary the PC across runs (S4.6); after each new run re-pool and")
        print("        re-check eta_2, eta_3, eta_4 against 1/sqrt(N_runs) (S4.15.3).")
        print("-" * 78)

    # ---- console summary ----
    print(f"Analysed {n} basis states (label-free selection from pcoords).")
    print(f"  basins: folded_side={len(folded_pool)}  transition={len(transition_pool)}"
          f"  unfolded_side={len(unfolded_pool)}")
    for p in ("bstate_pcoords.csv", "bstate_pcoords.dat", "recommendations.csv", "analysis.log"):
        print(f"  -> results/{p}")
    print("  figs: fig_pcoord_RMSD_MinDist.png, fig_Qstem_Qloop.png,")
    print("        fig_pcoord_distributions.png, fig_recommended_seeds.png, fig_Rg_e2e.png")


if __name__ == "__main__":
    sys.exit(main())
