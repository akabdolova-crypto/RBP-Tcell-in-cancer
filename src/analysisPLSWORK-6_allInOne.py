# -*- coding: utf-8 -*-
"""
analysisPLSWORK-6_allInOne.py
--------------------------------
This script augments your existing pipeline with:
  - Cleaned barplots (fixed scales, spacing, p-value subtitle) from the CSVs your pipeline writes
  - A–E follow-ups you requested:
      A) Pre-treatment R vs NR tests within CD8 stages (Mann-Whitney + BH-FDR)
      B) Paired Pre->Post tests within stages (Wilcoxon, fallback to MWU; BH-FDR)
      C) Composition-informed prediction of response (logistic regression w/ CV AUROC)
      D) DPT monotonicity QC per cohort (Spearman rho vs stage order) when available
      E) Extra visualizations: dot-heatmap (RBPs x CD8 stages), transition tables and heatmaps, optional ridgelines

It preserves your working logic. By default, it first runs your core script
('analysisPLSWORK-5_spacedPvals_fixed.py' if present), then executes the posthoc A–E suite.

Usage:
  python3 analysisPLSWORK-6_allInOne.py --mtx <path> --meta <path> --outdir results
  python3 analysisPLSWORK-6_allInOne.py --outdir results --skip-core  # only posthoc A–E
"""

import os, sys, argparse, subprocess, math, warnings
import numpy as np
import pandas as pd

# plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# stats
try:
    from scipy.stats import mannwhitneyu, wilcoxon, gaussian_kde, spearmanr
except Exception:
    mannwhitneyu = None
    wilcoxon = None
    gaussian_kde = None
    spearmanr = None

# sklearn
try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_score
    from sklearn.metrics import roc_auc_score
except Exception:
    StandardScaler = None
    LogisticRegression = None
    StratifiedKFold = None
    cross_val_score = None
    roc_auc_score = None

# ---------------------------
# Utilities
# ---------------------------
def log(msg: str):
    print(msg, flush=True)

def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)
    return p

def safe_slug(x):
    return "".join([c if c.isalnum() or c in "._-" else "_" for c in str(x)])

def bh_fdr(pvals):
    p = np.array(pvals, dtype=float)
    mask = ~np.isnan(p)
    if not np.any(mask):
        return p
    p_eff = p[mask]
    n = p_eff.size
    order = np.argsort(p_eff)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, n+1)
    fdr = p_eff * n / ranks
    # monotone
    fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
    fdr_final = np.empty_like(p_eff)
    fdr_final[order] = np.clip(fdr_sorted, 0, 1)
    out = np.full_like(p, np.nan, dtype=float)
    out[mask] = fdr_final
    return out

def mwu(a, b):
    if mannwhitneyu is None:
        return np.nan
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size < 2 or b.size < 2:
        return np.nan
    try:
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        return np.nan

def wilcoxon_paired(pre_vals, post_vals, pre_ids, post_ids):
    # Pair by intersection of patient IDs; fallback to MWU if not enough pairs.
    pre_ids = list(map(str, pre_ids)); post_ids = list(map(str, post_ids))
    common = sorted(set(pre_ids).intersection(set(post_ids)))
    if wilcoxon is None:
        return np.nan
    if len(common) < 2:
        return np.nan
    s1 = []
    s2 = []
    for pid in common:
        s1.append(pre_vals[pre_ids.index(pid)])
        s2.append(post_vals[post_ids.index(pid)])
    s1 = np.asarray(s1, float); s2 = np.asarray(s2, float)
    try:
        return float(wilcoxon(s1, s2).pvalue)
    except Exception:
        return np.nan

# ---------------------------
# Posthoc barplot helpers (clean spacing, fixed scales)
# ---------------------------
def plot_bar_prepost(ax, xcats, pre_mean, pre_sd, post_mean, post_sd,
                     title, ylabel, ylim=None, subtitle=None):
    w = 0.38
    x = np.arange(len(xcats))
    ax.bar(x - w/2, pre_mean, width=w, yerr=pre_sd, capsize=3, edgecolor="black", label="Pre")
    ax.bar(x + w/2, post_mean, width=w, yerr=post_sd, capsize=3, edgecolor="black", label="Post")
    ax.set_xticks(x); ax.set_xticklabels(xcats, rotation=30, ha="right")
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, pad=16, fontsize=12)
    ax.legend()
    ax.set_axisbelow(True); ax.grid(axis="y", linestyle=":", alpha=0.35)
    ax.tick_params(axis="x", labelsize=10); ax.tick_params(axis="y", labelsize=10)
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(0, ymax*1.15 if ymax>0 else 1.0)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    plt.tight_layout(rect=[0.02, 0.06, 0.98, 0.90])
    if subtitle:
        plt.gcf().text(0.5, 0.935, subtitle, ha="center", va="top", fontsize=10)

# ---------------------------
# Posthoc (cleaned figures) & A–E analyses
# ---------------------------
def posthoc_clean_barplots(results_dir: str, use_cd8_only=True, unify_y_per_gene=True):
    outdir = ensure_dir(os.path.join(results_dir, "barplots_improved"))
    # 1) Global patient means (Pre/Post x Response)
    fname = "rbp_expression_patientMeans_CD8only.csv" if use_cd8_only else "rbp_expression_patientMeans.csv"
    fpath = os.path.join(results_dir, fname)
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        ignore = {"therapy","timepoint","response","patient_id","n_cells"}
        genes = [c for c in df.columns if c not in ignore and pd.api.types.is_numeric_dtype(df[c])]
        for ther, DT in df.groupby("therapy"):
            for gene in genes:
                fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
                ylim_common = None
                if unify_y_per_gene:
                    vals = DT[gene].dropna().values
                    if vals.size:
                        ylim_common = (0, float(np.nanmax(vals))*1.20)
                for j, resp in enumerate(["Responder","Non-responder"]):
                    ax = axes[j]
                    D = DT[DT["response"]==resp]
                    pre = D[D["timepoint"]=="Pre"][gene].dropna().values
                    post= D[D["timepoint"]=="Post"][gene].dropna().values
                    pre_ids = D[D["timepoint"]=="Pre"]["patient_id"].astype(str).tolist()
                    post_ids= D[D["timepoint"]=="Post"]["patient_id"].astype(str).tolist()
                    # pre vs post (paired if possible)
                    p_prepost = wilcoxon_paired(pre.tolist(), post.tolist(), pre_ids, post_ids)
                    if np.isnan(p_prepost):
                        p_prepost = mwu(pre, post)
                    subtitle = f"Pre vs Post {resp}: " + ("n/a" if np.isnan(p_prepost) else f"p={p_prepost:.3g}")
                    cats = ["All patients"]
                    pre_m, pre_s = np.array([np.nanmean(pre) if pre.size else np.nan]), np.array([np.nanstd(pre, ddof=1) if pre.size>=2 else np.nan])
                    post_m,post_s= np.array([np.nanmean(post) if post.size else np.nan]), np.array([np.nanstd(post, ddof=1) if post.size>=2 else np.nan])
                    title = f"{ther} — {resp} — {gene}"
                    ylabel = f"{gene} (patient mean)"
                    plot_bar_prepost(ax, cats, pre_m, pre_s, post_m, post_s, title, ylabel, ylim=ylim_common, subtitle=subtitle)
                plt.suptitle(f"Pre/Post by Response — {gene}", y=0.999, fontsize=12)
                safe_t = safe_slug(ther); safe_g = safe_slug(gene)
                plt.savefig(os.path.join(outdir, f"bar_global_{safe_g}_{safe_t}_PrePost_byResponse_IMPROVED.pdf"), bbox_inches="tight")
                plt.close()

    # 2) Stage fractions (group stats)
    f_fracs = os.path.join(results_dir, "cd8_stageFractions_groupStats.csv")
    if os.path.exists(f_fracs):
        G = pd.read_csv(f_fracs)
        stages = [s for s in G["CD8_stage"].dropna().unique() if s != "Non-CD8"]
        for ther, GT in G.groupby("therapy"):
            for resp, GR in GT.groupby("response"):
                pre = GR[GR["timepoint"]=="Pre"].set_index("CD8_stage")
                post= GR[GR["timepoint"]=="Post"].set_index("CD8_stage")
                cats = [s for s in stages if (s in pre.index) or (s in post.index)]
                pre_m = np.array([pre.loc[s,"mean_fraction"] if s in pre.index else np.nan for s in cats])
                post_m= np.array([post.loc[s,"mean_fraction"] if s in post.index else np.nan for s in cats])
                pre_s = np.array([pre.loc[s,"sd_fraction"] if s in pre.index else np.nan for s in cats])
                post_s= np.array([post.loc[s,"sd_fraction"] if s in post.index else np.nan for s in cats])
                fig, ax = plt.subplots(figsize=(max(7.5, 0.55*len(cats)+2), 4.6))
                title = f"{ther} — {resp}"
                plot_bar_prepost(ax, cats, pre_m, pre_s, post_m, post_s, title, "CD8 stage fraction (mean ± SD)", ylim=(0,1), subtitle=None)
                safe_t = safe_slug(ther); safe_r = safe_slug(resp)
                plt.savefig(os.path.join(outdir, f"bar_stageFractions_{safe_t}_{safe_r}_PrePost_IMPROVED.pdf"), bbox_inches="tight")
                plt.close()

    # 3) RBP by CD8 stage (group stats)
    f_rbpstage = os.path.join(results_dir, "rbp_byCD8stage_groupStats.csv")
    if os.path.exists(f_rbpstage):
        RS = pd.read_csv(f_rbpstage)
        gene_cols = [c for c in RS.columns if c.startswith("mean_")]
        for ther, T in RS.groupby("therapy"):
            for resp, R2 in T.groupby("response"):
                for g in gene_cols:
                    gene = g.replace("mean_","")
                    pre = R2[R2["timepoint"]=="Pre"].set_index("CD8_stage")
                    post= R2[R2["timepoint"]=="Post"].set_index("CD8_stage")
                    stages = sorted(set(pre.index).union(set(post.index)))
                    stages = [s for s in stages if s!="Non-CD8"]
                    pre_m = np.array([pre.loc[s, g] if s in pre.index else np.nan for s in stages])
                    post_m= np.array([post.loc[s, g] if s in post.index else np.nan for s in stages])
                    pre_s = np.array([pre.loc[s, f"sd_{gene}"] if s in pre.index and f"sd_{gene}" in pre.columns else np.nan for s in stages])
                    post_s= np.array([post.loc[s, f"sd_{gene}"] if s in post.index and f"sd_{gene}" in post.columns else np.nan for s in stages])
                    fig, ax = plt.subplots(figsize=(max(7.5, 0.55*len(stages)+2), 4.6))
                    title = f"{ther} — {resp} — {gene}"
                    plot_bar_prepost(ax, stages, pre_m, pre_s, post_m, post_s, title, f"{gene} (patient mean within stage)", ylim=None, subtitle=None)
                    safe_t = safe_slug(ther); safe_r = safe_slug(resp); safe_g = safe_slug(gene)
                    plt.savefig(os.path.join(outdir, f"bar_rbpByStage_{safe_g}_{safe_t}_{safe_r}_PrePost_IMPROVED.pdf"), bbox_inches="tight")
                    plt.close()

# A) pre-treatment R vs NR within CD8 stages (patient-level)
def analysis_A_state_tests(results_dir: str, stages_focus=("T-EFF/CTL","Tex-int","Term-Tex")):
    outdir = ensure_dir(os.path.join(results_dir, "stats_A_state"))
    f = os.path.join(results_dir, "rbp_byCD8stage_patientMeans.csv")
    if not os.path.exists(f):
        log(f"[A] Missing {f}; skipping A-tests.")
        return
    DF = pd.read_csv(f)
    # Keep Pre
    DF = DF[DF["timepoint"]=="Pre"].copy()
    if "response" not in DF or "patient_id" not in DF:
        log("[A] Missing response/patient_id columns.")
        return
    gene_cols = [c for c in DF.columns if c not in ["therapy","timepoint","response","patient_id","CD8_stage","n_cells"] and pd.api.types.is_numeric_dtype(DF[c])]
    rows = []
    for ther, DT in DF.groupby("therapy"):
        for stage in stages_focus:
            DS = DT[DT["CD8_stage"]==stage]
            if DS.empty:
                continue
            for g in gene_cols:
                r = DS[DS["response"]=="Responder"][g].dropna().values
                n = DS[DS["response"]=="Non-responder"][g].dropna().values
                p = mwu(r, n)
                rows.append({"therapy": ther, "stage": stage, "gene": g, "n_R": len(r), "n_NR": len(n), "p_MWU": p})
    if not rows:
        log("[A] No rows; skip.")
        return
    OUT = pd.DataFrame(rows)
    OUT["FDR"] = bh_fdr(OUT["p_MWU"].values)
    OUT.sort_values(["therapy","stage","FDR","gene"], inplace=True)
    OUT.to_csv(os.path.join(outdir, "A_state_R_vs_NR_Pre.csv"), index=False)
    log(f"[A] Saved A_state_R_vs_NR_Pre.csv with {len(OUT)} rows.")

# B) longitudinal Pre->Post within stages, per response (patient-level)
def analysis_B_paired_tests(results_dir: str, stages_focus=("T-EFF/CTL","Tex-int","Term-Tex")):
    outdir = ensure_dir(os.path.join(results_dir, "stats_B_paired"))
    f = os.path.join(results_dir, "rbp_byCD8stage_patientMeans.csv")
    if not os.path.exists(f):
        log(f"[B] Missing {f}; skipping B-tests.")
        return
    DF = pd.read_csv(f)
    gene_cols = [c for c in DF.columns if c not in ["therapy","timepoint","response","patient_id","CD8_stage","n_cells"] and pd.api.types.is_numeric_dtype(DF[c])]
    rows = []
    for ther, DT in DF.groupby("therapy"):
        for resp, DR in DT.groupby("response"):
            for stage in stages_focus:
                S = DR[DR["CD8_stage"]==stage]
                pre = S[S["timepoint"]=="Pre"]
                post= S[S["timepoint"]=="Post"]
                for g in gene_cols:
                    pre_ids = pre["patient_id"].astype(str).tolist()
                    post_ids= post["patient_id"].astype(str).tolist()
                    pre_vals = pre[g].dropna().values
                    post_vals= post[g].dropna().values
                    p = wilcoxon_paired(pre_vals.tolist(), post_vals.tolist(), pre_ids, post_ids)
                    if np.isnan(p):
                        p = mwu(pre_vals, post_vals)
                    rows.append({"therapy": ther, "response": resp, "stage": stage,
                                 "gene": g, "n_Pre": pre_vals.size, "n_Post": post_vals.size,
                                 "p_Wilcoxon_or_MWU": p})
    if not rows:
        log("[B] No rows; skip.")
        return
    OUT = pd.DataFrame(rows)
    OUT["FDR"] = bh_fdr(OUT["p_Wilcoxon_or_MWU"].values)
    OUT.sort_values(["therapy","response","stage","FDR","gene"], inplace=True)
    OUT.to_csv(os.path.join(outdir, "B_paired_PrePost_withinStage.csv"), index=False)
    log(f"[B] Saved B_paired_PrePost_withinStage.csv with {len(OUT)} rows.")

# C) composition-informed prediction (Pre only): CTL/EFF + Term-Tex + RBP module + key RBPs
def analysis_C_prediction(results_dir: str, use_cd8_only=True):
    outdir = ensure_dir(os.path.join(results_dir, "stats_C_prediction"))
    f_frac = os.path.join(results_dir, "cd8_stageFractions_patient_all.csv")
    f_rbp  = os.path.join(results_dir, "rbp_expression_patientMeans_CD8only.csv" if use_cd8_only else "rbp_expression_patientMeans.csv")
    if not os.path.exists(f_frac) or not os.path.exists(f_rbp):
        log("[C] Missing patient-level fractions or patient means; skipping.")
        return
    FR = pd.read_csv(f_frac)  # per patient, many stage columns + therapy, timepoint, response
    RBP = pd.read_csv(f_rbp)

    # Focus Pre only; build features at patient level by therapy
    stages_needed = {"T-EFF/CTL","Term-Tex"}
    stage_cols = [c for c in FR.columns if c not in ["patient_id","response","therapy","timepoint"]]
    for s in stages_needed:
        if s not in stage_cols:
            FR[s] = 0.0  # absent stage -> 0 fraction

    # Keep Pre and only rows present in both files
    FRp = FR[FR["timepoint"]=="Pre"].copy()
    RBPp= RBP[RBP["timepoint"]=="Pre"].copy()
    join_keys = ["therapy","response","patient_id"]
    X = FRp[join_keys + list(stages_needed)].merge(RBPp[join_keys + [c for c in RBPp.columns if c not in ["timepoint","n_cells"]]],
                                                   on=join_keys, how="inner")
    if X.empty:
        log("[C] No joinable Pre rows; skipping.")
        return

    # y labels
    y_map = {"Responder":1, "Non-responder":0}
    if not set(X["response"].unique()).intersection(y_map.keys()):
        log("[C] No valid response labels.")
        return

    # For each therapy, build & evaluate model
    rows_cv = []
    coef_rows = []
    for ther, D in X.groupby("therapy"):
        y = D["response"].map(y_map).values
        # predictors: CTL/EFF, Term-Tex fractions; rbp_module-ish if present; plus individual RBPs
        ignore = set(["therapy","response","patient_id"])
        feat = [c for c in D.columns if c not in ignore and pd.api.types.is_numeric_dtype(D[c])]
        # some scripts might use "rbp_module_pm" as CD8-only mean module
        # keep all numeric predictors but guard dimensionality vs sample size
        Xmat = D[feat].values.astype(float)
        n, p = Xmat.shape
        if n < 8 or StandardScaler is None or LogisticRegression is None:
            rows_cv.append({"therapy": ther, "n_samples": n, "note": "too few or sklearn missing"})
            continue
        scaler = StandardScaler()
        Xs = scaler.fit_transform(Xmat)
        # mild L2 to prevent overfitting when n small
        model = LogisticRegression(penalty="l2", C=1.0, solver="liblinear", max_iter=200)
        # CV AUROC
        k = 5 if n >= 10 else (3 if n >= 6 else 2)
        try:
            cv = StratifiedKFold(n_splits=k, shuffle=True, random_state=1)
            aucs = cross_val_score(model, Xs, y, cv=cv, scoring="roc_auc")
            rows_cv.append({"therapy": ther, "n_samples": n, "n_features": p,
                            "cv_splits": k, "AUROC_mean": float(np.mean(aucs)),
                            "AUROC_sd": float(np.std(aucs))})
            # Fit on full set to export coefficients
            model.fit(Xs, y)
            coef = model.coef_[0]
            for fname, c in zip(feat, coef):
                coef_rows.append({"therapy":ther, "feature":fname, "coef":float(c)})
        except Exception as e:
            rows_cv.append({"therapy": ther, "n_samples": n, "n_features": p,
                            "cv_splits": k, "error": str(e)})

    if rows_cv:
        pd.DataFrame(rows_cv).to_csv(os.path.join(outdir, "C_prediction_CV_AUROC.csv"), index=False)
    if coef_rows:
        pd.DataFrame(coef_rows).to_csv(os.path.join(outdir, "C_prediction_full_fit_coefs.csv"), index=False)
    log("[C] Saved prediction outputs.")

# D) DPT monotonicity QC (per cohort): Spearman(DPT, stage_order)
def analysis_D_dpt_monotonicity(results_dir: str):
    """
    This step requires per-cohort CSVs with cell-level 'CD8_stage' and 'dpt_pseudotime'.
    If your core pipeline stores them as 'results/<COHORT>/cells_with_dpt.csv', we'll use them.
    Otherwise we gracefully skip and still produce a placeholder summary.
    """
    outdir = ensure_dir(os.path.join(results_dir, "stats_D_dpt"))
    # Search for cohort subfolders
    cohorts = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    rows = []
    stage_order = ["Stem-like","T-CM","T-EFF/CTL","Tex-int","Prog-Tex","Term-Tex"]
    rank = {s:i for i,s in enumerate(stage_order)}
    for c in sorted(cohorts):
        f = os.path.join(results_dir, c, "cells_with_dpt.csv")
        if not os.path.exists(f):
            continue
        try:
            DF = pd.read_csv(f)
            DF = DF.dropna(subset=["CD8_stage","dpt_pseudotime"])
            # CD8-only, non-innate
            DF = DF[DF["CD8_stage"]!="Non-CD8"].copy()
            if DF.empty:
                continue
            DF["stage_rank"] = DF["CD8_stage"].map(rank).astype(float)
            # Only rows with known rank
            DF = DF[~DF["stage_rank"].isna()].copy()
            if spearmanr is None or DF.empty:
                continue
            rho, p = spearmanr(DF["stage_rank"].values, DF["dpt_pseudotime"].values)
            rows.append({"cohort": c, "n_cells": DF.shape[0], "spearman_rho": rho, "p_value": p})
            # quick scatter-like plot
            plt.figure(figsize=(4.0,3.2))
            plt.scatter(DF["stage_rank"]+np.random.uniform(-0.07,0.07, size=DF.shape[0]), DF["dpt_pseudotime"], s=3, alpha=0.35)
            plt.xlabel("CD8 stage (ordinal)"); plt.ylabel("DPT")
            plt.title(f"{c} — Spearman rho={rho:.2f}, p={p:.2g}")
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"DPT_monotonicity_{safe_slug(c)}.pdf"))
            plt.close()
        except Exception as e:
            rows.append({"cohort": c, "error": str(e)})
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(outdir, "DPT_monotonicity_summary.csv"), index=False)
        log("[D] Saved DPT monotonicity summary.")
    else:
        log("[D] No cell-level DPT CSVs found; produce them in the core to enable this QC.")

# E) extra visualizations (dot-heatmap; stage transitions; optional ridgelines)
def analysis_E_visuals(results_dir: str):
    outdir = ensure_dir(os.path.join(results_dir, "visuals_extras"))
    # Dot heatmap: RBPs x CD8 stages (means). Use groupStats as an omnibus summary.
    f = os.path.join(results_dir, "rbp_byCD8stage_groupStats.csv")
    if os.path.exists(f):
        RS = pd.read_csv(f)
        # Average across therapy/timepoint/response to one mean per stage
        gene_cols = [c for c in RS.columns if c.startswith("mean_")]
        if gene_cols:
            M = RS.groupby("CD8_stage")[gene_cols].mean()
            stages = [s for s in M.index.tolist() if s!="Non-CD8"]
            M = M.loc[stages]
            genes = [g.replace("mean_","") for g in gene_cols]
            V = M.values  # rows=stages, cols=mean_gene
            # Dot heatmap
            fig, ax = plt.subplots(figsize=(max(8, 0.5*len(genes)+2), max(3.5, 0.45*len(stages)+1)))
            vmax = np.nanpercentile(V, 98)
            vmin = np.nanpercentile(V, 2)
            # normalize size by mean magnitude
            vnorm = (V - vmin) / max(1e-9, (vmax - vmin))
            for i, s in enumerate(stages):
                for j, g in enumerate(genes):
                    size = 50 + 250 * float(np.clip(vnorm[i,j], 0, 1))
                    ax.scatter(j, i, s=size, edgecolor="black", linewidth=0.3)
            ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=60, ha="right")
            ax.set_yticks(range(len(stages))); ax.set_yticklabels(stages)
            ax.set_title("RBPs × CD8 stages (mean expression across cohorts)", pad=12)
            ax.invert_yaxis()
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, "dotheatmap_RBPs_by_CD8stage.pdf"), bbox_inches="tight"); plt.close()

    # Stage transitions: table and heatmaps using per-patient fractions Pre/Post
    fT = os.path.join(results_dir, "cd8_stageFractions_patient_all.csv")
    if os.path.exists(fT):
        DF = pd.read_csv(fT)
        # choose dominant stage per patient at each timepoint
        stage_cols = [c for c in DF.columns if c not in ["patient_id","response","therapy","timepoint"]]
        rows = []
        for (ther, resp, pid), sub in DF.groupby(["therapy","response","patient_id"]):
            pre = sub[sub["timepoint"]=="Pre"]
            post= sub[sub["timepoint"]=="Post"]
            if pre.empty or post.empty:
                continue
            s_pre = pre[stage_cols].T[pre[stage_cols].T.idxmax()].index[0] if pre[stage_cols].size else None
            s_post= post[stage_cols].T[post[stage_cols].T.idxmax()].index[0] if post[stage_cols].size else None
            if s_pre and s_post:
                rows.append({"therapy":ther,"response":resp,"patient_id":pid,"pre_stage":s_pre,"post_stage":s_post})
        if rows:
            T = pd.DataFrame(rows)
            ensure_dir(os.path.join(outdir, "transitions"))
            T.to_csv(os.path.join(outdir, "transitions", "per_patient_pre_post_stage.csv"), index=False)
            # heatmap counts per therapy/response
            for (ther, resp), grp in T.groupby(["therapy","response"]):
                tab = pd.crosstab(grp["pre_stage"], grp["post_stage"]).astype(int)
                stages = sorted(set(tab.index).union(tab.columns))
                tab = tab.reindex(index=stages, columns=stages, fill_value=0)
                plt.figure(figsize=(max(5, 0.4*len(stages)+2), max(4, 0.4*len(stages)+2)))
                im = plt.imshow(tab.values, cmap="Blues")
                plt.colorbar(im, fraction=0.046, pad=0.04)
                plt.xticks(range(len(stages)), stages, rotation=60, ha="right")
                plt.yticks(range(len(stages)), stages)
                plt.title(f"Stage transitions Pre→Post — {ther} — {resp}", pad=12)
                for i in range(len(stages)):
                    for j in range(len(stages)):
                        v = tab.values[i,j]
                        if v>0:
                            plt.text(j, i, str(v), ha="center", va="center", fontsize=8)
                plt.tight_layout()
                plt.savefig(os.path.join(outdir, "transitions", f"transition_heatmap_{safe_slug(ther)}_{safe_slug(resp)}.pdf"), bbox_inches="tight")
                plt.close()

    # Optional: ridgelines of exhaustion score (requires patient-level exhaustion CSV)
    fEx = os.path.join(results_dir, "exhaustion_patientMeans.csv")
    if os.path.exists(fEx) and gaussian_kde is not None:
        EX = pd.read_csv(fEx)
        for ther, D in EX.groupby("therapy"):
            plt.figure(figsize=(8, 5.5))
            y0 = 0
            order = [("Pre","Responder"), ("Pre","Non-responder"), ("Post","Responder"), ("Post","Non-responder")]
            for tp, resp in order:
                vals = D[(D["timepoint"]==tp) & (D["response"]==resp)]["exhaustion_score"].dropna().values
                if vals.size >= 5:
                    kde = gaussian_kde(vals)
                    xs = np.linspace(np.percentile(vals, 1), np.percentile(vals, 99), 200)
                    ys = kde(xs); ys = ys / ys.max() * 0.8
                    plt.fill_between(xs, y0, y0+ys, alpha=0.35)
                    plt.text(xs[0], y0+0.85, f"{tp} {resp} (n={vals.size})", fontsize=9, va="bottom")
                    y0 += 1.0
                elif vals.size > 0:
                    plt.text(np.median(vals), y0+0.4, f"{tp} {resp} (n={vals.size})", fontsize=9, ha="center")
                    y0 += 1.0
            plt.yticks([]); plt.xlabel("Exhaustion score (patient mean)")
            plt.title(f"Ridgelines: exhaustion score — {ther}")
            plt.tight_layout()
            plt.savefig(os.path.join(outdir, f"ridgeline_exhaustion_{safe_slug(ther)}.pdf"), bbox_inches="tight")
            plt.close()

# ---------------------------
# Orchestration
# ---------------------------
def run_core_if_present(args):
    # Try your latest core script; fall back to previous name variants if needed.
    candidates = [
        "analysisPLSWORK-5_spacedPvals_fixed.py",
        "analysisPLSWORK-5_spacedPvalues.py",
        "analysisPLSWORK-5_spacedPvals.py",
    ]
    core = None
    for c in candidates:
        if os.path.exists(c):
            core = c; break
    if core is None:
        log("[core] No core script found; skipping core run and proceeding to posthoc A–E.")
        return
    # assemble command
    cmd = ["python3", core, "--outdir", args.outdir]
    if args.mtx:  cmd += ["--mtx", args.mtx]
    if args.meta: cmd += ["--meta", args.meta]
    if args.rbps and len(args.rbps)>0:
        cmd += ["--rbps"] + args.rbps
    log("[core] Running: " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        log(f"[core] Core script exited with error: {e}. Continuing with posthoc A–E with whatever outputs exist.")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mtx", default=None, help="matrix file")
    ap.add_argument("--meta", default=None, help="meta file")
    ap.add_argument("--outdir", default="results", help="output directory")
    ap.add_argument("--skip-core", action="store_true", help="skip running core pipeline; only run posthoc A–E")
    ap.add_argument("--rbps", nargs="*", default=None, help="optional: list of RBPs")
    ap.add_argument("--use-cd8-only", action="store_true",
                    help="prefer CD8-only patient means for global barplots/prediction if present (default True)", default=True)
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # 0) run core if not skipped
    if not args.skip_core:
        run_core_if_present(args)

    # 1) posthoc cleaned barplots (scales & spacing)
    log("[posthoc] Generating cleaned barplots…")
    posthoc_clean_barplots(args.outdir, use_cd8_only=args.use_cd8_only, unify_y_per_gene=True)

    # 2) A–E analyses
    log("[A] State-resolved tests (Pre, R vs NR, within CD8 stages)…")
    analysis_A_state_tests(args.outdir)

    log("[B] Paired/longitudinal (Pre->Post within stage, per response)…")
    analysis_B_paired_tests(args.outdir)

    log("[C] Composition-informed prediction (Pre only)…")
    analysis_C_prediction(args.outdir, use_cd8_only=args.use_cd8_only)

    log("[D] DPT monotonicity QC…")
    analysis_D_dpt_monotonicity(args.outdir)

    log("[E] Extra visualizations (dot-heatmap; transitions; ridgelines if available)…")
    analysis_E_visuals(args.outdir)

    log("Done. See results/ subfolders for outputs.")

if __name__ == "__main__":
    main()
