

# Run with: python3 analysisBLYAT.py
# Sade-Feldman melanoma scRNA-seq – RBP-centric analysis for CTL→Tex
# Adds: (i) UMAP pages for ALL 11 RBPs, (ii) statistics across CTL CD8+,
#       pre-Tex, Tex-int, Term-Tex states using exhaustion score,
#       (iii) per-patient tests + heatmaps.

import os, re, time, csv, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy.sparse as sp
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
def _sanitize_uns_for_h5ad(adata: ad.AnnData):
    """Convert any tuples/sets in .uns to JSON/HDF5-friendly types (lists)."""
    def _conv(x):
        if isinstance(x, tuple):
            return [_conv(v) for v in x]
        if isinstance(x, set):
            return sorted([_conv(v) for v in x])
        if isinstance(x, dict):
            return {k: _conv(v) for k, v in x.items()}
        if isinstance(x, list):
            return [_conv(v) for v in x]
        return x
    # mutate in place
    for k in list(adata.uns.keys()):
        adata.uns[k] = _conv(adata.uns[k])



warnings.filterwarnings("ignore")
sc.settings.verbosity = 2

# ---------------------------
# Paths / Config
# ---------------------------
MTX  = "data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"
META = "data/GSE120575_patient_ID_single_cells.txt"
OUTDIR = "results"

HVG = 2000
N_CELLS_CAP = None     # e.g. 8000 if RAM limited per cohort
SUBSET_CTL = False     # global prefilter to CD8 T (heuristic). Stats use their own CD8 gate.

THERAPIES  = ["anti-CTLA4", "anti-PD1", "anti-CTLA4+PD1"]
TIMEPOINTS = ["Pre", "Post"]

EXHAUSTION_GENES = ["PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","TOX","EOMES","CXCL13","ENTPD1","LAYN"]
CORE_RBPS        = ["ZFP36","ZFP36L1","ZFP36L2","ZC3H12A","RC3H1","RC3H2","PCBP1","ELAVL1","HNRNPLL","TIA1","TIA1B"]
DEFAULT_CORE_RBPS = CORE_RBPS.copy()
EXPECTED_PATIENTS = 32

# ---------------------------
# tiny logger
# ---------------------------
def ensure_outdir():
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "run.log"), "a") as f:
        f.write(f"\n=== New run @ {time.ctime()} ===\n")
def log(msg: str):
    print(msg, flush=True)
    with open(os.path.join(OUTDIR, "run.log"), "a") as f:
        f.write(msg + "\n")

# ---------------------------
# Normalization helpers
# ---------------------------
def canon_patient_id(pid: str) -> str:
    if not isinstance(pid, str): return "UNK"
    s = pid.strip().upper()
    m = re.search(r"P(\d+)$", s) or re.search(r"(?:PRE|POST)[\-_ ]?P(\d+)", s)
    return f"P{int(m.group(1))}" if m else ("UNK" if s in ("", "NA", "NONE") else s)

def normalize_therapy(t: str) -> str:
    if not isinstance(t, str): return "UNK"
    s = t.strip().lower().replace("–","-").replace("—","-").replace(" ", "")
    s = s.replace("pd-1","pd1").replace("ctla-4","ctla4").replace("ipilimumab","ctla4")
    s = s.replace("nivolumab","pd1").replace("pembrolizumab","pd1")
    has_pd1, has_ctla4 = ("pd1" in s), ("ctla4" in s)
    if has_pd1 and has_ctla4: return "anti-CTLA4+PD1"
    if has_pd1: return "anti-PD1"
    if has_ctla4: return "anti-CTLA4"
    return "UNK"

def normalize_timepoint(tp: str) -> str:
    s = str(tp).strip().lower()
    if s.startswith("pre") or "baseline" in s: return "Pre"
    if s.startswith("post") or "on-treatment" in s or "ontreatment" in s or "during" in s: return "Post"
    if re.search(r"(?:^|[^a-z])(pre)[^a-z]*p\d+", s): return "Pre"
    if re.search(r"(?:^|[^a-z])(post)[^a-z]*p\d+", s): return "Post"
    return "UNK"

def normalize_response_label(s: str) -> str:
    s = str(s).strip().lower()
    if s in {"responder","r","yes","response","cr","pr","partial","complete"}: return "Responder"
    if s in {"non-responder","nonresponder","nr","no","sd","pd","progression"}: return "Non-responder"
    return "UNK"

# ---------------------------
# Read matrix header (cells only)
# ---------------------------
def read_header_cells(mtx_path: str):
    with open(mtx_path, "r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    gene_col = header[0]
    cells = header[1:]
    cell_to_idx = {c: i+1 for i, c in enumerate(cells)}  # 0 is gene column
    with open(os.path.join(OUTDIR, "header_first50.txt"), "w") as fh:
        for c in cells[:50]: fh.write(c + "\n")
    return gene_col, cells, cell_to_idx

# ---------------------------
# META parser (GEO SAMPLES table)
# ---------------------------
HEADER_KEYS = {
    "sample_id_like": ["title"],
    "patient_like":   ["characteristics: patinet id", "characteristics: patient id", "characteristics: patient"],
    "response_like":  ["characteristics: response", "response"],
    "therapy_like":   ["characteristics: therapy", "therapy", "treatment", "drug", "checkpoint"],
}
SID_REGEX = re.compile(r"^[A-Za-z0-9]+_[A-Za-z0-9]+_[A-Za-z0-9]+(?:_[A-Za-z0-9]+)?$")

def split_meta_line(line: str) -> list[str]:
    if "\t" in line: return line.rstrip("\n").split("\t")
    try:
        row = next(csv.reader([line]));  row = [c.strip() for c in row]
        if len(row) > 1: return row
    except Exception: pass
    return [p for p in re.split(r"\s{2,}", line.strip()) if p]

def build_header_map(header_cells: list[str]) -> dict:
    col_map = {}
    low = [h.strip().lower() for h in header_cells]
    for want, keys in HEADER_KEYS.items():
        idx = None
        for k in keys:
            for j, col in enumerate(low):
                if col.startswith(k): idx = j; break
            if idx is not None: break
        col_map[want] = idx
    return col_map

def parse_samples_table(meta_path: str) -> list[dict]:
    rows, in_table, header, col_map = [], False, None, None
    with open(meta_path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            s = raw.strip().strip("\ufeff")
            if not s or s.startswith(("#","//","\"#")): continue
            if not in_table and s.lower().startswith("sample name"):
                header = split_meta_line(raw);  in_table = True
                col_map = build_header_map(header);  continue
            if in_table and s.startswith("Sample "):
                parts = split_meta_line(raw)
                if header and len(parts) < len(header): parts += [""]*(len(header)-len(parts))
                sid_val = ""
                if col_map and col_map.get("sample_id_like") is not None:
                    sid_val = parts[col_map["sample_id_like"]]
                if not sid_val:
                    for v in parts:
                        if SID_REGEX.match(v.strip()): sid_val = v.strip(); break
                if not sid_val:
                    toks = [p for p in parts if p];  sid_val = toks[1] if len(toks) >= 2 else ""
                patient_field = parts[col_map["patient_like"]] if col_map and col_map.get("patient_like") is not None else ""
                response_field = parts[col_map["response_like"]] if col_map and col_map.get("response_like") is not None else ""
                therapy_field  = parts[col_map["therapy_like"]]  if col_map and col_map.get("therapy_like")  is not None else ""
                rows.append(dict(
                    sid=(sid_val or "").strip().replace(" ", "_"),
                    patient_id_raw=patient_field.strip(),
                    timepoint_raw=patient_field.strip(),
                    response_raw=(response_field or "").strip(),
                    therapy_raw=(therapy_field or "").strip(),
                ))
    if not rows: raise RuntimeError("Failed to parse GEO SAMPLES table; check META file.")
    return rows

def canonicalize_sample_id(s: str) -> str:
    s2 = s.replace("-", "_").strip()
    parts = s2.split("_")
    if parts: parts[-1] = re.sub(r"([0-9]+)[A-Za-z]*$", r"\1", parts[-1])
    return "_".join(parts)

def token_key(s: str, n: int) -> str:
    toks = s.split("_")
    return "_".join(toks[:n]) if len(toks) >= n else s

def build_sample_maps(sample_rows: list[dict]) -> tuple[dict, dict, dict]:
    map3, map4, full = {}, {}, {}
    for r in sample_rows:
        sid_full = canonicalize_sample_id(r["sid"])
        if not sid_full: continue
        info = dict(
            patient_id = canon_patient_id(r["patient_id_raw"]),
            timepoint  = normalize_timepoint(r["timepoint_raw"]),
            response   = normalize_response_label(r["response_raw"]),
            therapy    = normalize_therapy(r["therapy_raw"]),
            sid_full   = sid_full,
        )
        k3 = token_key(sid_full, 3);  k4 = token_key(sid_full, 4)
        map3.setdefault(k3, info);  map4.setdefault(k4, info);  full.setdefault(sid_full, info)
    return map3, map4, full

def map_cell_to_sample(cid: str, map3: dict, map4: dict, full: dict) -> dict | None:
    cid_c = canonicalize_sample_id(cid)
    if cid_c in full: return full[cid_c]
    k4 = token_key(cid_c, 4)
    if k4 in map4: return map4[k4]
    k3 = token_key(cid_c, 3)
    if k3 in map3: return map3[k3]
    return None

def load_meta_per_cell_from_samples(meta_path: str, header_cells: list[str]) -> pd.DataFrame:
    rows = parse_samples_table(meta_path)
    map3, map4, full = build_sample_maps(rows)
    cov3 = sum(1 for c in header_cells if token_key(canonicalize_sample_id(c), 3) in map3)
    cov4 = sum(1 for c in header_cells if token_key(canonicalize_sample_id(c), 4) in map4)
    log(f"[META] SAMPLES-mode; token coverage: 3-tok={cov3}/{len(header_cells)}, 4-tok={cov4}/{len(header_cells)}")
    out, unmatched = [], []
    for cid in header_cells:
        info = map_cell_to_sample(cid, map3, map4, full)
        if info is None:
            unmatched.append(cid)
            out.append((cid, "UNK","UNK","UNK","UNK",""))
        else:
            out.append((cid, info["patient_id"], info["timepoint"], info["response"], info["therapy"], info["sid_full"]))
    if unmatched:
        with open(os.path.join(OUTDIR, "unmatched_cells.tsv"), "w") as fh:
            for x in unmatched: fh.write(x + "\n")
    meta = pd.DataFrame(out, columns=["cell_id","patient_id","timepoint","response","therapy","sample_id"]).set_index("cell_id")
    return meta

# ---------------------------
# Scanpy helpers
# ---------------------------
def qc_preprocess_keep_all(adata: ad.AnnData, hvg=HVG):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg, flavor="seurat")
    try:
        sc.tl.pca(adata, mask_var="highly_variable", svd_solver="arpack")
    except TypeError:
        sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True)
    n_pcs = int(min(30, adata.obsm["X_pca"].shape[1]))
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata)
    return adata

def ensure_neighbors(adata: ad.AnnData):
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        if "X_pca" not in adata.obsm:
            try:
                sc.tl.pca(adata, mask_var="highly_variable", svd_solver="arpack")
            except TypeError:
                sc.tl.pca(adata, svd_solver="arpack", use_highly_variable=True)
        n_pcs = int(min(30, adata.obsm["X_pca"].shape[1]))
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)

def ensure_exhaustion(adata: ad.AnnData):
    sig = [g for g in EXHAUSTION_GENES if g in adata.var_names]
    if len(sig) >= 3:
        sc.tl.score_genes(adata, sig, score_name="exhaustion_score", use_raw=False)
    else:
        adata.obs["exhaustion_score"] = np.nan

def ensure_rbp_module(adata: ad.AnnData, rbp_list):
    present = [g for g in rbp_list if g in adata.var_names]
    if len(present) >= 3:
        sc.tl.score_genes(adata, present, score_name="rbp_module", use_raw=False)
    else:
        adata.obs["rbp_module"] = np.nan

def subset_cd8_t(adata: ad.AnnData) -> ad.AnnData:
    pos = [g for g in ["CD3D","CD3E","CD8A","CD8B"] if g in adata.var_names]
    if len(pos) < 2:
        log("   [info] CD8 gating skipped (markers missing).")
        return adata
    X = adata[:, pos].X
    X = X.A if hasattr(X, "A") else np.array(X)
    score = X.mean(axis=1)
    adata.obs["cd8_score"] = score
    thr = np.percentile(score, 60)  # keep top 40%
    sub = adata[score >= thr].copy()
    log(f"   [CD8 gate] kept {sub.n_obs}/{adata.n_obs} cells.")
    return sub

def add_pseudotime(adata: ad.AnnData):
    ensure_neighbors(adata)
    sc.tl.diffmap(adata)
    root_idx = 0
    cand = [g for g in ["TCF7","CCR7","IL7R","LEF1"] if g in adata.var_names]
    if cand:
        vals = adata[:, cand].X
        vals = vals.A if hasattr(vals, "A") else np.array(vals)
        scores = np.nan_to_num(vals, nan=0.0).sum(axis=1)
        root_idx = int(np.argmax(scores))
    adata.uns["iroot"] = int(root_idx)
    try:
        sc.tl.dpt(adata, n_dcs=10)
    except Exception:
        try:
            sc.tl.dpt(adata)
        except Exception:
            log("   [warn] DPT failed; continuing without pseudotime.")
    log(f"   DPT present: {'dpt_pseudotime' in adata.obs.columns}")
    return adata

# ---------------------------
# Cohort data loader (column-index pull)
# ---------------------------
def load_cohort_adata_by_index(mtx_path, meta, all_cells_set, cell_to_idx,
                               therapy: str, timepoint: str, n_cells_cap=None):
    mask = (meta["therapy"].astype(str) == therapy) & \
           (meta["timepoint"].astype(str) == timepoint) & \
           (meta["patient_id"].astype(str) != "UNK")
    cells_c = [c for c in meta.index[mask] if c in all_cells_set]
    if not cells_c: return None
    if n_cells_cap is not None and len(cells_c) > n_cells_cap: cells_c = cells_c[:n_cells_cap]
    usecols_idx = [0] + sorted({cell_to_idx[c] for c in cells_c})
    expr_c = pd.read_csv(mtx_path, sep="\t", header=0, index_col=0,
                         usecols=usecols_idx, low_memory=False)
    expr_num = expr_c.apply(pd.to_numeric, errors="coerce")
    dropped_mask = expr_num.isna().all(axis=1)
    if dropped_mask.any():
        tag = f"{therapy.replace('/','-')}_{timepoint}"
        with open(os.path.join(OUTDIR, f"dropped_rows_{tag}.txt"), "w", encoding="utf-8") as fh:
            for idx in expr_num.index[dropped_mask]: fh.write(str(idx) + "\n")
    expr_num = expr_num[~dropped_mask].dropna(axis=1, how="all")
    dfX = expr_num.T.fillna(0.0).astype("float32")  # cells × genes
    adata = ad.AnnData(X=dfX.values)
    adata.obs_names = dfX.index.astype(str)
    adata.var_names = dfX.columns.astype(str)
    adata.obs = adata.obs.join(meta.loc[adata.obs_names], how="left")
    return adata

# ---------------------------
# RBP ranking (for ordering, optional)
# ---------------------------
def rank_rbps(adata: ad.AnnData, rbp_list=None):
    if rbp_list is None: rbp_list = CORE_RBPS
    rbp_list = [g for g in rbp_list if g in adata.var_names]
    if not rbp_list: return pd.DataFrame([])
    ensure_exhaustion(adata)

    # simple co-expression network
    centrality = {g:0.0 for g in rbp_list}
    if len(rbp_list) >= 3:
        Xm = adata[:, rbp_list].X
        Xm = Xm.A if hasattr(Xm,"A") else np.array(Xm)
        R = np.corrcoef(Xm, rowvar=False)
        G = nx.Graph()
        for i, gi in enumerate(rbp_list):
            for j in range(i+1, len(rbp_list)):
                r = R[i, j]
                if np.isfinite(r) and abs(r) >= 0.2:
                    G.add_edge(gi, rbp_list[j], weight=float(abs(r)))
        if len(G): centrality = nx.pagerank(G, alpha=0.85, max_iter=500)

    dpt = adata.obs["dpt_pseudotime"].values if "dpt_pseudotime" in adata.obs.columns else None
    rows = []
    for g in rbp_list:
        x = adata[:, g].X
        x = x.A.squeeze() if hasattr(x, "A") else np.array(x).squeeze()
        r_tex, _ = spearmanr(x, adata.obs["exhaustion_score"].values, nan_policy="omit")
        r_dpt, _ = (spearmanr(x, dpt, nan_policy="omit") if dpt is not None else (np.nan, None))
        rows.append({"rbp": g, "spearman_tex": r_tex, "rho_pseudotime": r_dpt, "centrality": centrality.get(g, 0.0)})
    df = pd.DataFrame(rows)
    if df.empty: return df
    df["score"] = (df["centrality"].fillna(0)*0.2
                   + df["spearman_tex"].abs().fillna(0)*0.5
                   + df["rho_pseudotime"].abs().fillna(0)*0.3)
    return df.sort_values("score", ascending=False).reset_index(drop=True)

# ---------------------------
# Exhaustion-state labeling for CD8+ cells
# ---------------------------
def label_tex_states(adata: ad.AnnData) -> ad.AnnData:
    """Define CD8+ and 4 Tex states using exhaustion-score quartiles among CD8+."""
    # CD8 score (recompute in case SUBSET_CTL=False)
    pos = [g for g in ["CD8A","CD8B"] if g in adata.var_names]
    if len(pos) == 0:
        adata.obs["cd8_score"] = np.nan
        adata.obs["tex_state"] = "Other"
        return adata
    X = adata[:, pos].X
    X = X.A if hasattr(X, "A") else np.array(X)
    adata.obs["cd8_score"] = X.mean(axis=1)
    cd8_thr = np.nanpercentile(adata.obs["cd8_score"].values, 60)  # keep top 40%
    is_cd8 = adata.obs["cd8_score"].values >= cd8_thr
    ex = adata.obs["exhaustion_score"].values
    # handle missing scores
    valid = np.isfinite(ex) & is_cd8
    tex_state = np.array(["Other"]*adata.n_obs, dtype=object)
    if valid.sum() >= 50:
        q1, q2, q3 = np.nanpercentile(ex[valid], [25,50,75])
        tex_state[(valid) & (ex <= q1)] = "CTL_CD8"
        tex_state[(valid) & (ex > q1) & (ex <= q2)] = "pre-Tex"
        tex_state[(valid) & (ex > q2) & (ex <= q3)] = "Tex-int"
        tex_state[(valid) & (ex > q3)] = "Term-Tex"
    adata.obs["tex_state"] = tex_state
    return adata

# ---------------------------
# Per-patient tests across states
# ---------------------------
def patient_state_table(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.var_names]
    rows = []
    for (pid, state), obs_g in adata.obs.groupby(["patient_id","tex_state"]):
        if state not in ["CTL_CD8","pre-Tex","Tex-int","Term-Tex"]: continue
        sub = adata[obs_g.index, :]
        means = {g: float(np.mean(sub[:, g].X)) for g in genes}
        rows.append({"patient_id": pid, "tex_state": state, **means, "n_cells": int(sub.n_obs)})
    return pd.DataFrame(rows)

def state_stats_per_gene(df_pat: pd.DataFrame, gene: str) -> dict:
    # gather vectors per state
    groups = {s: df_pat.loc[df_pat["tex_state"]==s, gene].dropna().values
              for s in ["CTL_CD8","pre-Tex","Tex-int","Term-Tex"]}
    # counts
    counts = {f"n_pat_{s}": int(np.sum(~np.isnan(groups[s]))) for s in groups}
    # Kruskal across available states with >=3 patients
    avail = [v for v in groups.values() if len(v) >= 3]
    kw_p = kruskal(*avail).pvalue if len(avail) >= 2 else np.nan
    # Adjacent pairwise tests
    pairs = [("CTL_CD8","pre-Tex"), ("pre-Tex","Tex-int"), ("Tex-int","Term-Tex")]
    out_p = {}
    for a,b in pairs:
        A, B = groups[a], groups[b]
        if len(A) >= 3 and len(B) >= 3:
            U, p = mannwhitneyu(A, B, alternative="two-sided")
            out_p[f"p_{a}_vs_{b}"] = p
        else:
            out_p[f"p_{a}_vs_{b}"] = np.nan
    # BH-FDR over the three pairwise p-values
    ps = [out_p[k] for k in out_p if np.isfinite(out_p[k])]
    padj_map = {}
    if len(ps) >= 1:
        padj = multipletests(ps, method="fdr_bh")[1]
        j = 0
        for k in out_p:
            if np.isfinite(out_p[k]):
                padj_map[f"q_{k[2:]}"] = float(padj[j]); j += 1
            else:
                padj_map[f"q_{k[2:]}"] = np.nan
    return {"gene": gene, "kw_p": kw_p, **counts, **out_p, **padj_map}


# ---------------------------
# NEW: Patient-averaged barplots (mean +/- SD) and statistical tests
# ---------------------------
def _mk_group_label(tp, resp):
    return f"{tp} | {resp}"

def _group_vectors(df: pd.DataFrame, gene: str, therapy: str):
    """Return value arrays and counts for Pre/Responder, Pre/Non-responder, Post/Responder, Post/Non-responder within a therapy."""
    out = {}
    for tp in ["Pre","Post"]:
        for resp in ["Responder","Non-responder"]:
            mask = (df["therapy"]==therapy) & (df["timepoint"]==tp) & (df["response"]==resp)
            vals = df.loc[mask, gene].dropna().values
            out[(tp, resp)] = vals
    return out

def _paired_overlap(a_df: pd.DataFrame, b_df: pd.DataFrame, gene: str):
    """Return paired arrays aligned by patient_id (intersection)."""
    # requires columns: patient_id, gene
    A = a_df[["patient_id", gene]].dropna()
    B = b_df[["patient_id", gene]].dropna()
    import numpy as _np
    common = _np.intersect1d(A["patient_id"].values, B["patient_id"].values)
    if common.size == 0:
        return _np.array([]), _np.array([])
    A2 = A.set_index("patient_id").loc[common, gene].values
    B2 = B.set_index("patient_id").loc[common, gene].values
    return A2, B2

def _stats_prepost_by_response(df: pd.DataFrame, therapy: str, gene: str) -> dict:
    """Compute group means/SDs/N and p-values for:
       * R vs NR within Pre and within Post (Mann-Whitney)
       * Pre vs Post within R and within NR (paired Wilcoxon if matched patients exist, else Mann-Whitney)."""
    from scipy.stats import mannwhitneyu, wilcoxon
    res = {"therapy": therapy, "gene": gene}
    # group summary
    G = _group_vectors(df, gene, therapy)
    for (tp, resp), vals in G.items():
        lab = _mk_group_label(tp, resp)
        res[f"mean_{lab}"] = float(np.mean(vals)) if vals.size else np.nan
        res[f"sd_{lab}"]   = float(np.std(vals, ddof=1)) if vals.size >= 2 else np.nan
        res[f"n_{lab}"]    = int(vals.size)
    # R vs NR at Pre and Post
    for tp in ["Pre","Post"]:
        A = G[(tp,"Responder")]; B = G[(tp,"Non-responder")]
        if A.size >= 3 and B.size >= 3:
            U, p = mannwhitneyu(A, B, alternative="two-sided")
            res[f"p_R_vs_NR_{tp}"] = float(p)
        else:
            res[f"p_R_vs_NR_{tp}"] = np.nan
    # Pre vs Post within R and within NR
    for resp in ["Responder","Non-responder"]:
        a = df[(df["therapy"]==therapy) & (df["timepoint"]=="Pre")  & (df["response"]==resp)][["patient_id", gene]]
        b = df[(df["therapy"]==therapy) & (df["timepoint"]=="Post") & (df["response"]==resp)][["patient_id", gene]]
        A, B = _paired_overlap(a, b, gene)
        if A.size >= 3 and B.size >= 3:
            try:
                stat, p = wilcoxon(A, B, zero_method="wilcox", alternative="two-sided")
                res[f"p_Pre_vs_Post_{resp}_paired"] = float(p)
                res[f"n_paired_{resp}"] = int(A.size)
            except Exception:
                # fallback to unpaired
                A2 = a[gene].dropna().values; B2 = b[gene].dropna().values
                if A2.size >= 3 and B2.size >= 3:
                    U, p = mannwhitneyu(A2, B2, alternative="two-sided")
                    res[f"p_Pre_vs_Post_{resp}_unpaired"] = float(p)
                    res[f"n_unpaired_{resp}"] = int(min(A2.size, B2.size))
                else:
                    res[f"p_Pre_vs_Post_{resp}_unpaired"] = np.nan
                    res[f"n_unpaired_{resp}"] = int(min(A2.size if A2.size else 0, B2.size if B2.size else 0))
        else:
            # not enough matched patients; do unpaired if possible
            A2 = a[gene].dropna().values; B2 = b[gene].dropna().values
            if A2.size >= 3 and B2.size >= 3:
                U, p = mannwhitneyu(A2, B2, alternative="two-sided")
                res[f"p_Pre_vs_Post_{resp}_unpaired"] = float(p)
                res[f"n_unpaired_{resp}"] = int(min(A2.size, B2.size))
            else:
                res[f"p_Pre_vs_Post_{resp}_unpaired"] = np.nan
                res[f"n_unpaired_{resp}"] = int(min(A2.size if A2.size else 0, B2.size if B2.size else 0))
    return res

def _plot_bar_prepost_by_response(df: pd.DataFrame, therapy: str, gene: str, outdir: str):
    """Barplot: mean +/- SD across patients for 4 groups: Pre|Responder, Pre|Non-responder, Post|Responder, Post|Non-responder."""
    import matplotlib.pyplot as plt
    groups = [("Pre","Responder"), ("Pre","Non-responder"), ("Post","Responder"), ("Post","Non-responder") ]
    stats = _stats_prepost_by_response(df, therapy, gene)
    means, sds, ns, labels = [], [], [], []
    for (tp, resp) in groups:
        lab = _mk_group_label(tp, resp)
        means.append(stats.get(f"mean_{lab}", np.nan))
        sds.append(stats.get(f"sd_{lab}", np.nan))
        ns.append(stats.get(f"n_{lab}", 0))
        labels.append(lab.replace(" | ", "\n"))
    x = np.arange(len(groups))
    width = 0.65
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.bar(x, means, width=width, yerr=sds, capsize=3, edgecolor="black")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(f"{gene} (patient mean, log1p norm.)")
    ax.set_title(f"{gene} — {therapy}: mean +/- SD per patient")
    # Add n below ticks
    for i, n in enumerate(ns):
        ax.text(i, 0.02, f"n={n}", ha="center", va="bottom", transform=ax.get_xaxis_transform(), fontsize=8)
    # Summarize p-values in a subtitle line
    p_parts = []
    for tp in ["Pre","Post"]:
        p = stats.get(f"p_R_vs_NR_{tp}", np.nan)
        p_parts.append("R vs NR {}: ".format(tp) + ("p={:.3g}".format(p) if np.isfinite(p) else "n/a"))
    for resp in ["Responder","Non-responder"]:
        key_p = f"p_Pre_vs_Post_{resp}_paired" if f"p_Pre_vs_Post_{resp}_paired" in stats else "p_Pre_vs_Post_{}_unpaired".format(resp)
        p = stats.get(key_p, np.nan)
        p_parts.append("Pre vs Post {}: ".format(resp) + ("p={:.3g}".format(p) if np.isfinite(p) else "n/a"))
    ax.text(0.5, 1.02, " | ".join(p_parts), transform=ax.transAxes, ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    import re as _re
    safe_gene = _re.sub(r"[^A-Za-z0-9_.-]+", "_", gene)
    safe_therapy = _re.sub(r"[^A-Za-z0-9_.-]+", "_", therapy)
    out_pdf = os.path.join(outdir, f"barplot_{safe_gene}_{safe_therapy}_PrePost_byResponse.pdf")
    plt.savefig(out_pdf, bbox_inches="tight"); plt.close()
    # Save stats CSV
    out_csv = os.path.join(outdir, f"stats_{safe_gene}_{safe_therapy}_PrePost_byResponse.csv")
    pd.DataFrame([stats]).to_csv(out_csv, index=False)


def make_prepost_response_barplots(df_pat_all: pd.DataFrame, genes: list[str], outdir: str):
    """Entry point: builds mean+/-SD barplots and stats for each therapy and each gene.
    Also writes an aggregated CSV per therapy with BH-FDR per comparison type."""
    if df_pat_all is None or df_pat_all.empty:
        return
    therapies = sorted([t for t in df_pat_all["therapy"].dropna().unique() if t != "UNK"])
    from statsmodels.stats.multitest import multipletests as _mt
    for therapy in therapies:
        stats_rows = []
        for gene in genes:
            if gene not in df_pat_all.columns:
                continue
            # generate figure + per-gene CSV
            _plot_bar_prepost_by_response(df_pat_all, therapy, gene, outdir)
            # collect stats for FDR
            st = _stats_prepost_by_response(df_pat_all, therapy, gene)
            stats_rows.append(st)
        if stats_rows:
            S = pd.DataFrame(stats_rows)
            # unify keys for Pre vs Post tests (paired/unpaired)
            for resp in ["Responder","Non-responder"]:
                paired_key = f"p_Pre_vs_Post_{resp}_paired"
                unpaired_key = f"p_Pre_vs_Post_{resp}_unpaired"
                S[f"p_Pre_vs_Post_{resp}"] = S[paired_key].fillna(S[unpaired_key]) if paired_key in S.columns and unpaired_key in S.columns else S.get(paired_key, S.get(unpaired_key, np.nan))
            # Apply BH FDR per comparison type if column exists
            def _bh(col):
                if col in S.columns:
                    p = S[col].astype(float).values
                    mask = np.isfinite(p)
                    q = np.full_like(p, np.nan, dtype=float)
                    if mask.sum() >= 2:
                        q[mask] = _mt(p[mask], method="fdr_bh")[1]
                    S[col.replace("p_","q_")] = q
            for col in ["p_R_vs_NR_Pre", "p_R_vs_NR_Post", "p_Pre_vs_Post_Responder", "p_Pre_vs_Post_Non-responder"]:
                _bh(col)
            # Save aggregated
            safe_therapy = re.sub(r"[^A-Za-z0-9_.-]+", "_", therapy)
            S.to_csv(os.path.join(outdir, f"stats_summary_{safe_therapy}_PrePost_byResponse.csv"), index=False)

def _get_connectivities_with_self(adata: ad.AnnData):
    ensure_neighbors(adata)
    A = adata.obsp.get("connectivities", None)
    if A is None:
        raise RuntimeError("KNN connectivities missing; cannot smooth markers.")
    # add self-loops
    A = A.tocsr()
    I = sp.eye(A.shape[0], format="csr")
    return (A + I)

def _smooth_gene(adata: ad.AnnData, gene: str, A=None):
    if A is None:
        A = _get_connectivities_with_self(adata)
    n = adata.n_obs
    if gene not in adata.var_names:
        return np.full(n, np.nan, dtype=float)
    x = adata[:, gene].X
    if hasattr(x, "A"): x = x.A
    x = np.asarray(x).reshape(n, -1)
    if x.shape[1] != 1:
        x = x[:, 0].reshape(n, 1)
    x = x.ravel()
    denom = np.maximum(np.array(A.sum(1)).ravel(), 1e-9)
    sm = (A @ x) / denom
    return np.asarray(sm).ravel()

def _quantiles(v: np.ndarray, qlo=0.30, qhi=0.70):
    vv = v[np.isfinite(v)]
    if vv.size == 0:
        return (np.nan, np.nan)
    return (np.quantile(vv, qlo), np.quantile(vv, qhi))

def _make_thresholds(adata: ad.AnnData, genes: list[str], A=None):
    thr = {}
    for g in genes:
        v = _smooth_gene(adata, g, A=A)
        adata.obs[g + "_sm"] = v
        lo, hi = _quantiles(v, 0.30, 0.70)
        thr[g] = (lo, hi)
    return thr

def _HIGH(adata: ad.AnnData, g: str):
    lo, hi = adata.uns.get(g + "_thr", (np.nan, np.nan))
    v = adata.obs.get(g + "_sm", pd.Series(np.nan, index=adata.obs_names)).values
    if not np.isfinite(hi): return np.zeros(adata.n_obs, dtype=bool)
    return v >= hi

def _LOW(adata: ad.AnnData, g: str):
    lo, hi = adata.uns.get(g + "_thr", (np.nan, np.nan))
    v = adata.obs.get(g + "_sm", pd.Series(np.nan, index=adata.obs_names)).values
    if not np.isfinite(lo): return np.zeros(adata.n_obs, dtype=bool)
    return v <= lo

def _MID(adata: ad.AnnData, g: str):
    lo, hi = adata.uns.get(g + "_thr", (np.nan, np.nan))
    v = adata.obs.get(g + "_sm", pd.Series(np.nan, index=adata.obs_names)).values
    if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        return np.zeros(adata.n_obs, dtype=bool)
    return (v > lo) & (v < hi)

def _NOT_HIGH_OPTIONAL(adata: ad.AnnData, g: str):
    # If gene missing/undetermined => don't penalize (return True)
    lo, hi = adata.uns.get(g + "_thr", (np.nan, np.nan))
    if not np.isfinite(hi): return np.ones(adata.n_obs, dtype=bool)
    v = adata.obs.get(g + "_sm", pd.Series(np.nan, index=adata.obs_names)).values
    return ~(v >= hi)



# ---------------------------
# CD8 stage marker dictionary and display orders (added to avoid NameError)
# ---------------------------
# Minimal, dataset-agnostic marker panel used by label_cd8_stages().
# Adjust/extend if your dataset uses alternative symbols.
CD8_STAGE_MARKERS = {
    "CD3D": "CD3D", "CD3E": "CD3E",
    "CD8A": "CD8A", "CD8B": "CD8B",
    "SELL": "SELL", "CCR7": "CCR7",          # naive/central memory axis
    "TCF7": "TCF7", "IL7R": "IL7R",          # stem-like axis
    "CXCR6": "CXCR6", "ITGAE": "ITGAE", "CD69": "CD69",  # tissue residency / TRM
    "GZMB": "GZMB", "PRF1": "PRF1",          # effector
    "TOX": "TOX", "EOMES": "EOMES", "CXCL13": "CXCL13", "ENTPD1": "ENTPD1"  # exhaustion axis
}

# Full display order used in composition plots. Only stages present in the data are shown.
CD8_STAGE_ORDER_ALL = [
    "Non-CD8",
    "Naive/T_SCM", "T-CM", "Stem-like",
    "T-EFF/CTL", "T-EM/TEMRA",
    "Tex-int", "Prog-Tex", "Term-Tex",
    "T-RM", "TRM-CD69only", "Other_CD8",
    "CTL_CD8"
]

# Main axis order for monotonicity checks along the exhaustion trajectory.
CD8_STAGE_ORDER_MAIN = ["Stem-like", "T-CM", "T-EFF/CTL", "Tex-int", "Prog-Tex", "Term-Tex"]

# Optional simplified order (not strictly required by the script but referenced in docs).
CD8_STAGE_ORDER_SIMPLE = ["CTL_CD8", "Tex-int", "Term-Tex"]


# ---------------------------
# Utility: ensure CD8_stage pandas.Categorical has unique, ordered categories
# ---------------------------
def _sanitize_cd8_stage_categories(adata: ad.AnnData):
    try:
        if "CD8_stage" in adata.obs.columns:
            s = adata.obs["CD8_stage"].astype(str)
            cats_pref = [c for c in CD8_STAGE_ORDER_ALL if c in set(s.values)]
            seen = set(); cats_unique = []
            for c in cats_pref:
                if c not in seen:
                    cats_unique.append(c); seen.add(c)
            if not cats_unique:
                cats_unique = list(pd.unique(s.values))
            adata.obs["CD8_stage"] = pd.Categorical(s, categories=cats_unique, ordered=True)
    except Exception as e:
        log(f"   [warn] sanitize CD8_stage categories failed: {e}")

def label_cd8_stages(adata: ad.AnnData, qlo=0.30, qhi=0.70) -> ad.AnnData:
    """Assign adata.obs['CD8_stage'] using smoothed marker gating, with CD3+CD8+ gating and extended states."""
    # Smooth & thresholds for all relevant markers
    A = _get_connectivities_with_self(adata)
    genes = [g for g in CD8_STAGE_MARKERS.values()]
    thr = _make_thresholds(adata, genes, A=A)
    for g, (lo, hi) in thr.items():
        adata.uns[g + "_thr"] = [float(lo), float(hi)]

    n = adata.n_obs

    # --- CD3+CD8+ T-cell mask (smoothed thresholds) ---
    has_cd3 = _HIGH(adata, "CD3D") | _HIGH(adata, "CD3E")
    cd8_pos = _HIGH(adata, "CD8A") | _HIGH(adata, "CD8B")
    is_cd8_t = has_cd3 & cd8_pos
    adata.obs["is_CD8_T"] = is_cd8_t.astype(bool)

    # --- Overlays (orthogonal flags) ---
    overlay_cycling = _HIGH(adata,"MKI67") | _HIGH(adata,"TOP2A") | _HIGH(adata,"PCNA") | _HIGH(adata,"TK1") | _HIGH(adata,"HMGB2")
    adata.obs["CD8_overlay_cycling"] = overlay_cycling.astype(bool)

    # Innate/unconventional overlays
    overlay_MAIT = _HIGH(adata,"KLRB1") & (_HIGH(adata,"TRAV1-2") | _HIGH(adata,"SLC4A10")) & (_HIGH(adata,"IL7R") | _HIGH(adata,"RORC"))
    overlay_iNKT = _HIGH(adata,"ZBTB16") & (_HIGH(adata,"TRAV10") | _HIGH(adata,"TRBV25-1") | _HIGH(adata,"KLRB1"))
    overlay_gdT  = _HIGH(adata,"TRDC") | _HIGH(adata,"TRGC1") | _HIGH(adata,"TRGC2") | _HIGH(adata,"TRDV2") | _HIGH(adata,"TRGV9")
    overlay_NKlike = (_HIGH(adata,"NKG7") & _HIGH(adata,"GNLY")) & ( _HIGH(adata,"KIR2DL1") | _HIGH(adata,"KIR2DL3") | _HIGH(adata,"KIR3DL1") | _HIGH(adata,"FGFBP2") | _HIGH(adata,"XCL1") | _HIGH(adata,"KLRD1") | _HIGH(adata,"KLRK1") )
    overlay_DP = _HIGH(adata,"CD4") & (_HIGH(adata,"CD8A") | _HIGH(adata,"CD8B"))
    overlay_Treglike = _HIGH(adata,"FOXP3") & _HIGH(adata,"IKZF2") & _HIGH(adata,"CTLA4") & _HIGH(adata,"IL2RA")

    adata.obs["CD8_overlay_MAIT"] = overlay_MAIT.astype(bool)
    adata.obs["CD8_overlay_iNKT"] = overlay_iNKT.astype(bool)
    adata.obs["CD8_overlay_gammadelta"] = overlay_gdT.astype(bool)
    adata.obs["CD8_overlay_NKlike"] = overlay_NKlike.astype(bool)
    adata.obs["CD8_overlay_DP"] = overlay_DP.astype(bool)
    adata.obs["CD8_overlay_Treglike"] = overlay_Treglike.astype(bool)

    # --- Stage gates (apply only inside CD8 T compartment) ---
    in_cd8 = is_cd8_t

    # New/extended states
    is_TRM  = in_cd8 & _HIGH(adata,"CD69") & _HIGH(adata,"ITGAE")
    is_TRM_CD69only = in_cd8 & _HIGH(adata,"CD69") & (~_HIGH(adata,"ITGAE"))
    is_TCM  = in_cd8 & _HIGH(adata,"SELL") & _HIGH(adata,"CCR7") & _NOT_HIGH_OPTIONAL(adata,"PDCD1")
    is_NAIVE = in_cd8 & _HIGH(adata,"CCR7") & _HIGH(adata,"SELL") & (_HIGH(adata,"IL7R") | _HIGH(adata,"LEF1") | _HIGH(adata,"LTB")) & _HIGH(adata,"TCF7") & _LOW(adata,"PDCD1") & _LOW(adata,"CD44")
    is_STEM = in_cd8 & _HIGH(adata,"TCF7") & _HIGH(adata,"SELL") & _LOW(adata,"CXCR6") & _NOT_HIGH_OPTIONAL(adata,"PDCD1")

    eff_support = _HIGH(adata,"IFNG") | _HIGH(adata,"TNF") | _HIGH(adata,"PRF1") | _HIGH(adata,"GZMB") | _HIGH(adata,"NKG7")
    is_TEFF = in_cd8 & _LOW(adata,"TCF7") & _LOW(adata,"SELL") & _HIGH(adata,"CXCR6") & eff_support

    # T-EM / TEMRA-like
    is_TEM = in_cd8 & _LOW(adata,"SELL") & _LOW(adata,"CCR7") & (_HIGH(adata,"GZMK") | _HIGH(adata,"PRF1") | _HIGH(adata,"GZMB")) & (_HIGH(adata,"KLRG1") | _HIGH(adata,"ZEB2") | _HIGH(adata,"CX3CR1") | _LOW(adata,"IL7R"))

    # Exhaustion spectrum
    is_TERM = in_cd8 & _HIGH(adata,"PDCD1") & _HIGH(adata,"CD44") & _HIGH(adata,"HAVCR2") & _HIGH(adata,"LAG3") & _LOW(adata,"TCF7")
    is_PROG = in_cd8 & _LOW(adata,"TCF7") & _LOW(adata,"SLAMF6") & _LOW(adata,"CXCR5") & _MID(adata,"PDCD1") & _HIGH(adata,"CD44")
    is_TEXINT = in_cd8 & (_MID(adata,"PDCD1") | _HIGH(adata,"PDCD1")) & (_MID(adata,"TCF7")) & ( _MID(adata,"HAVCR2") | _MID(adata,"LAG3") | (~_HIGH(adata,"HAVCR2")) | (~_HIGH(adata,"LAG3")) ) & (_HIGH(adata,"GZMK") | _MID(adata,"GZMK") | _MID(adata,"TOX"))

    # Initialize labels: Non-CD8 by default
    labels = np.array(["Non-CD8"] * n, dtype=object)
    def set_where(mask, name):
        idx = np.where(mask & (labels == "Non-CD8"))[0]
        labels[idx] = name

    # Precedence (specific → general; avoid Tex-int eating Term/Prog)
    set_where(is_TRM,           "T-RM")
    set_where(is_TRM_CD69only,  "TRM-CD69only")
    set_where(is_TERM,          "Term-Tex")
    set_where(is_PROG,          "Prog-Tex")
    set_where(is_TEXINT,        "Tex-int")
    set_where(is_TEFF,          "T-EFF/CTL")
    set_where(is_TEM,           "T-EM/TEMRA")
    set_where(is_STEM,          "Stem-like")
    set_where(is_TCM,           "T-CM")
    set_where(is_NAIVE,         "Naive/T_SCM")

    # Remaining CD8 T cells that didn't match explicit gates
    rem_cd8 = np.where(in_cd8 & (labels == "Non-CD8"))[0]
    labels[rem_cd8] = "Other_CD8"

    cats = []
    for _c in (["Non-CD8"] + CD8_STAGE_ORDER_ALL):
        if _c not in cats:
            cats.append(_c)
    adata.obs["CD8_stage"] = pd.Categorical(labels, categories=cats, ordered=True)
    return adata

def compute_oriented_dpt(adata: ad.AnnData, exclude_label="T-RM", stages_col="CD8_stage") -> ad.AnnData:
    """Compute DPT oriented from Stem-like/T-CM → Tex-int → Prog-Tex → Term-Tex on CD3+CD8+ cells, excluding side branches."""
    adata.obs["DPT"] = np.nan

    # Ensure labels and CD8 mask
    if stages_col not in adata.obs.columns or "is_CD8_T" not in adata.obs.columns:
        log("   [info] Ensuring CD8 stage labels and CD8 mask.")
        adata = label_cd8_stages(adata)
        _sanitize_cd8_stage_categories(adata)

    labels = adata.obs[stages_col].astype(str).values
    cd8_mask = adata.obs["is_CD8_T"].astype(bool).values

    # Exclude innate-like overlays from the Tex axis if present
    innate_flags = []
    for col in ["CD8_overlay_MAIT","CD8_overlay_iNKT","CD8_overlay_gammadelta","CD8_overlay_NKlike","CD8_overlay_DP","CD8_overlay_Treglike"]:
        if col in adata.obs.columns:
            innate_flags.append(adata.obs[col].astype(bool).values)
    is_innate_like = np.any(np.vstack(innate_flags), axis=0) if len(innate_flags) else np.zeros(adata.n_obs, dtype=bool)

    keep = cd8_mask & (~is_innate_like) & (labels != exclude_label) & (labels != "TRM-CD69only") & (labels != "Other_CD8")

    if keep.sum() < 200:
        log(f"   [warn] Too few CD8 cells for stable DPT after exclusions (n={int(keep.sum())}); falling back to global dpt_pseudotime if available.")
        if "dpt_pseudotime" in adata.obs.columns:
            adata.obs["DPT"] = adata.obs["dpt_pseudotime"].values
        return adata

    sub = adata[keep].copy()
    try:
        ensure_neighbors(sub); sc.tl.diffmap(sub)
    except Exception as e:
        log(f"   [warn] diffmap failed: {e}; falling back to global dpt_pseudotime if available.")
        if "dpt_pseudotime" in adata.obs.columns:
            adata.obs["DPT"] = adata.obs["dpt_pseudotime"].values
        return adata

    # Root among Stem-like/T-CM with lowest PDCD1_sm
    root_pool = np.where(sub.obs[stages_col].isin(["Stem-like","T-CM"]))[0]
    if root_pool.size == 0:
        root_pool = np.arange(sub.n_obs)
    pd1 = sub.obs.get("PDCD1_sm", pd.Series(np.zeros(sub.n_obs), index=sub.obs_names)).values
    sub.uns["iroot"] = int(root_pool[np.argmin(pd1[root_pool])])

    sc.tl.dpt(sub)
    adata.obs.loc[sub.obs_names, "DPT"] = sub.obs["dpt_pseudotime"].values

    # Monotonicity along main axis now includes Tex-int
    stage_means = (sub.obs.groupby(stages_col, observed=True)["dpt_pseudotime"].mean())
    seq = [s for s in CD8_STAGE_ORDER_MAIN if s in stage_means.index]
    vals = stage_means.reindex(seq).values
    monotonic = bool(np.all(np.diff(vals[np.isfinite(vals)]) >= -1e-9))
    log("   DPT means (CD8-only, no innate-like): " + ", ".join([f"{k}={stage_means[k]:.3f}" for k in seq if k in stage_means]))
    log(f"   DPT monotonic along {' → '.join(seq)}: {'YES' if monotonic else 'NO'}")
    if not monotonic:
        log("   [note] DPT not monotonic; consider excluding tiny islands or increasing n_neighbors (e.g., 30–40).")
    return adata
def plot_cd8_stage_suite(adata: ad.AnnData, tag: str, outdir: str):
    """UMAP by CD8_stage, DPT-by-stage plots, stacked composition bars."""
    # UMAP colored by CD8_stage
    if "CD8_stage" in adata.obs.columns and "X_umap" in adata.obsm:
        try:
            sc.pl.umap(adata[adata.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in adata.obs.columns else adata, color=["CD8_stage"], show=False)
            plt.savefig(os.path.join(outdir, f"umap_{tag}_CD8_stage.pdf"), bbox_inches="tight"); plt.close()
        except Exception as e:
            log(f"   [warn] UMAP by CD8_stage failed: {e}")

    # DPT means by stage
    if "DPT" in adata.obs.columns and "CD8_stage" in adata.obs.columns:
        means = (adata.obs.groupby("CD8_stage", observed=True)["DPT"].mean().dropna())
        log("   DPT means by CD8_stage (all): " + ", ".join([f"{k}={v:.3f}" for k,v in means.items()]))
        means.to_csv(os.path.join(outdir, f"DPT_means_by_CD8_stage_{tag}.csv"))
        # Monotonicity check on main axis (ignore T-RM/Other)
        seq = [s for s in CD8_STAGE_ORDER_MAIN if s in means.index]
        vals = means.loc[seq].values if len(seq)>0 else np.array([])
        is_monotonic = bool(np.all(np.diff(vals) >= -1e-9)) if vals.size else False
        with open(os.path.join(outdir, f"DPT_monotonicity_{tag}.txt"), "w") as fh:
            fh.write(f"Order: {seq}\nMeans: {list(map(float, vals))}\nMonotonic_non_decreasing: {is_monotonic}\n")
        if not is_monotonic:
            log("   [note] DPT means are not monotonic across main stages. Consider excluding small islands or increasing n_neighbors.")

        # Violin plot (monotonicity check)
        try:
            import seaborn as sns
            order = [s for s in CD8_STAGE_ORDER_ALL if s in adata.obs["CD8_stage"].cat.categories]
            dfp = adata.obs[["CD8_stage","DPT"]].copy()
            dfp = dfp[np.isfinite(dfp["DPT"].values)]
            if not dfp.empty:
                plt.figure(figsize=(8,4))
                sns.violinplot(data=dfp, x="CD8_stage", y="DPT", order=order, cut=0, inner="box")
                plt.title(f"DPT by CD8_stage — {tag}")
                plt.xticks(rotation=30, ha="right")
                plt.tight_layout()
                plt.savefig(os.path.join(outdir, f"plot_DPT_by_CD8_stage_violin_{tag}.pdf"))
                plt.close()

                # Bar plot (means with 95% CI)
                plt.figure(figsize=(8,4))
                sns.barplot(data=dfp, x="CD8_stage", y="DPT", order=order, estimator=np.mean, ci=95)
                plt.title(f"DPT mean by CD8_stage — {tag}")
                plt.xticks(rotation=30, ha="right")
                plt.tight_layout()
                plt.savefig(os.path.join(outdir, f"plot_DPT_by_CD8_stage_bar_{tag}.pdf"))
                plt.close()
        except Exception as e:
            log(f"   [warn] violin/bar DPT plots failed: {e}")

    # Stacked bars per patient × timepoint × response
    if "patient_id" in adata.obs.columns and "CD8_stage" in adata.obs.columns:
        obs = adata.obs.copy()
        # Keep only stages of interest (optional)
        obs["CD8_stage_simple"] = obs["CD8_stage"].astype(str)
        # Composition by patient
        comp = (obs.groupby(["patient_id","CD8_stage_simple"], observed=True).size()
                   .unstack(fill_value=0))
        if comp.shape[0] >= 1:
            prop = comp.div(comp.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
            # plot stacked bar
            try:
                ax = prop.reindex(columns=[c for c in CD8_STAGE_ORDER_ALL if c in prop.columns]).loc[:, [c for c in CD8_STAGE_ORDER_ALL if c in prop.columns]].plot(kind="bar", stacked=True, figsize=(12,4), legend=True)
                ax.set_ylabel("Proportion")
                ax.set_title(f"CD8_stage composition by patient — {tag}")
                plt.tight_layout()
                plt.savefig(os.path.join(outdir, f"stackedbar_CD8_stage_by_patient_{tag}.pdf"))
                plt.close()
            except Exception as e:
                log(f"   [warn] stacked bar by patient failed: {e}")

        # Composition by response (aggregated within cohort tag)
        if "response" in obs.columns:
            comp2 = (obs.groupby(["response","CD8_stage_simple"], observed=True).size()
                       .unstack(fill_value=0))
            if comp2.shape[0] >= 1:
                prop2 = comp2.div(comp2.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
                try:
                    ax = prop2.reindex(columns=[c for c in CD8_STAGE_ORDER_ALL if c in prop2.columns]).loc[:, [c for c in CD8_STAGE_ORDER_ALL if c in prop2.columns]].plot(kind="bar", stacked=True, figsize=(6,4), legend=True)
                    ax.set_ylabel("Proportion")
                    ax.set_title(f"CD8_stage composition by response — {tag}")
                    plt.tight_layout()
                    plt.savefig(os.path.join(outdir, f"stackedbar_CD8_stage_by_response_{tag}.pdf"))
                    plt.close()
                except Exception as e:
                    log(f"   [warn] stacked bar by response failed: {e}")

# ---------------------------
# Per-patient CD8_stage fractions and RBP-by-stage means + exports
# ---------------------------
def _per_patient_stage_fractions(adata: ad.AnnData) -> pd.DataFrame:
    if "CD8_stage" not in adata.obs.columns:
        return pd.DataFrame()
    obs = adata.obs[["patient_id","response","CD8_stage"]].copy()
    comp = (obs.groupby(["patient_id","response","CD8_stage"], observed=True).size()
              .unstack(fill_value=0))
    frac = comp.div(comp.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0).reset_index()
    ther = adata.obs["therapy"].mode(dropna=False).iat[0] if "therapy" in adata.obs.columns else "UNK"
    tp   = adata.obs["timepoint"].mode(dropna=False).iat[0] if "timepoint" in adata.obs.columns else "UNK"
    frac["therapy"] = ther; frac["timepoint"] = tp
    return frac

def _per_patient_rbp_means_by_stage(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.var_names]
    if not genes or "CD8_stage" not in adata.obs.columns:
        return pd.DataFrame()
    rows = []
    for (pid, resp, stage), idx in adata.obs.groupby(["patient_id","response","CD8_stage"]).groups.items():
        sub = adata[list(idx), :]
        means = {g: float(np.mean(sub[:, g].X)) for g in genes}
        rows.append({"patient_id": pid, "response": resp, "CD8_stage": stage, **means, "n_cells": int(sub.n_obs)})
    df = pd.DataFrame(rows)
    ther = adata.obs["therapy"].mode(dropna=False).iat[0] if "therapy" in adata.obs.columns else "UNK"
    tp   = adata.obs["timepoint"].mode(dropna=False).iat[0] if "timepoint" in adata.obs.columns else "UNK"
    df["therapy"] = ther; df["timepoint"] = tp
    return df

def _export_stage_fraction_group_stats(stage_fracs_all: list[pd.DataFrame], outdir: str):
    if not stage_fracs_all: return
    df = pd.concat(stage_fracs_all, ignore_index=True)
    df.to_csv(os.path.join(outdir, "cd8_stageFractions_patient_all.csv"), index=False)
    stage_cols = [c for c in df.columns if c not in ["patient_id","response","therapy","timepoint"]]
    L = df.melt(id_vars=["therapy","timepoint","response","patient_id"], value_vars=stage_cols,
                var_name="CD8_stage", value_name="fraction")
    g = (L.groupby(["therapy","timepoint","response","CD8_stage"], observed=True)
           .agg(mean_fraction=("fraction","mean"),
                sd_fraction=("fraction", lambda x: float(np.std(x, ddof=1)) if len(x)>=2 else np.nan),
                n_patients=("fraction","count"))
           .reset_index())
    g.to_csv(os.path.join(outdir, "cd8_stageFractions_groupStats.csv"), index=False)

def _export_rbp_by_stage_group_stats(rbp_by_stage_all: list[pd.DataFrame], outdir: str):
    if not rbp_by_stage_all: return
    df = pd.concat(rbp_by_stage_all, ignore_index=True)
    df.to_csv(os.path.join(outdir, "rbp_byCD8stage_patientMeans.csv"), index=False)
    gene_cols = [c for c in df.columns if c not in ["therapy","timepoint","response","patient_id","CD8_stage","n_cells"]]
    rows = []
    for (ther, tp, resp, stage), sub in df.groupby(["therapy","timepoint","response","CD8_stage"], observed=True):
        n = sub["patient_id"].nunique()
        rec = {"therapy": ther, "timepoint": tp, "response": resp, "CD8_stage": stage, "n_patients": int(n)}
        for gname in gene_cols:
            vals = sub[gname].dropna().values
            rec[f"mean_{gname}"] = float(np.mean(vals)) if vals.size else np.nan
            rec[f"sd_{gname}"]   = float(np.std(vals, ddof=1)) if vals.size>=2 else np.nan
        rows.append(rec)
    stats = pd.DataFrame(rows)
    stats.to_csv(os.path.join(outdir, "rbp_byCD8stage_groupStats.csv"), index=False)

def _plot_stage_fraction_grouped_bars(stage_fracs_all: list[pd.DataFrame], outdir: str):
    if not stage_fracs_all: return
    df = pd.concat(stage_fracs_all, ignore_index=True)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    stage_cols = [c for c in df.columns if c not in ["patient_id","response","therapy","timepoint"]]
    for therapy, dft in df.groupby("therapy"):
        for resp, dfr in dft.groupby("response"):
            pre = dfr[dfr["timepoint"]=="Pre"]; post = dfr[dfr["timepoint"]=="Post"]
            mean_pre = pre[stage_cols].mean(axis=0); sd_pre = pre[stage_cols].std(axis=0, ddof=1); n_pre = pre["patient_id"].nunique()
            mean_post= post[stage_cols].mean(axis=0); sd_post= post[stage_cols].std(axis=0, ddof=1); n_post= post["patient_id"].nunique()
            x = np.arange(len(stage_cols)); w = 0.38
            fig, ax = plt.subplots(figsize=(max(7.0, 0.55*len(stage_cols)+2), 4.6))
            ax.bar(x-w/2, mean_pre.values, width=w, yerr=sd_pre.values, capsize=3, edgecolor="black", label="Pre")
            ax.bar(x+w/2, mean_post.values, width=w, yerr=sd_post.values, capsize=3, edgecolor="black", label="Post")
            ax.set_xticks(x); ax.set_xticklabels(stage_cols, rotation=30, ha="right")
            ax.set_ylabel("CD8_stage fraction (mean +/- SD)")
            ax.set_title(f"{therapy} — {resp}")
            ax.legend(); ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
            plt.tight_layout()
            safe_t = re.sub(r"[^A-Za-z0-9_.-]+","_", str(therapy)); safe_r = re.sub(r"[^A-Za-z0-9_.-]+","_", str(resp))
            plt.savefig(os.path.join(outdir, f"barplot_CD8stageFractions_{safe_t}_{safe_r}_PrePost.pdf"), bbox_inches="tight"); plt.close()

def _plot_rbp_by_stage_bars(rbp_by_stage_all: list[pd.DataFrame], genes: list[str], outdir: str):
    if not rbp_by_stage_all: return
    df = pd.concat(rbp_by_stage_all, ignore_index=True)
    import matplotlib.pyplot as plt
    stage_order = [s for s in CD8_STAGE_ORDER_ALL if s in df["CD8_stage"].unique().tolist()]
    for therapy, dft in df.groupby("therapy"):
        for resp, dfr in dft.groupby("response"):
            for gene in genes:
                if gene not in df.columns: continue
                rows = []
                for stage in stage_order:
                    pre = dfr[(dfr["timepoint"]=="Pre") & (dfr["CD8_stage"]==stage)][gene].dropna().values
                    post= dfr[(dfr["timepoint"]=="Post") & (dfr["CD8_stage"]==stage)][gene].dropna().values
                    rows.append((stage, np.mean(pre) if pre.size else np.nan, np.std(pre, ddof=1) if pre.size>=2 else np.nan,
                                       np.mean(post) if post.size else np.nan, np.std(post, ddof=1) if post.size>=2 else np.nan))
                if not rows: continue
                stages = [r[0] for r in rows]
                mean_pre = np.array([r[1] for r in rows]); sd_pre = np.array([r[2] for r in rows])
                mean_post= np.array([r[3] for r in rows]); sd_post= np.array([r[4] for r in rows])
                x = np.arange(len(stages)); w = 0.38
                fig, ax = plt.subplots(figsize=(max(7.0, 0.55*len(stages)+2), 4.6))
                ax.bar(x-w/2, mean_pre, width=w, yerr=sd_pre, capsize=3, edgecolor="black", label="Pre")
                ax.bar(x+w/2, mean_post, width=w, yerr=sd_post, capsize=3, edgecolor="black", label="Post")
                ax.set_xticks(x); ax.set_xticklabels(stages, rotation=30, ha="right")
                ax.set_ylabel(f"{gene} (mean +/- SD, patient-level)")
                ax.set_title(f"{therapy} — {resp} — {gene}")
                ax.legend(); plt.tight_layout()
                safe_t = re.sub(r"[^A-Za-z0-9_.-]+","_", str(therapy)); safe_r = re.sub(r"[^A-Za-z0-9_.-]+","_", str(resp)); safe_g = re.sub(r"[^A-Za-z0-9_.-]+","_", str(gene))
                plt.savefig(os.path.join(outdir, f"barplot_RBP_byCD8stage_{safe_g}_{safe_t}_{safe_r}_PrePost.pdf"), bbox_inches="tight"); plt.close()

def main():
    ensure_outdir()
    log(f"Scanpy version: {sc.__version__}")
    # [1/8]
    log("[1/8] Reading matrix header…")
    gene_col, all_cells, cell_to_idx = read_header_cells(MTX)
    all_cells_set = set(all_cells)
    log(f"Header read: {len(all_cells)} cells in matrix")

    # [2/8] META
    log("[2/8] Loading META and deriving per-cell META…")
    meta = load_meta_per_cell_from_samples(META, all_cells)
    meta["patient_id"] = meta["patient_id"].apply(canon_patient_id)
    meta["therapy"]    = meta["therapy"].apply(normalize_therapy)
    meta["timepoint"]  = meta["timepoint"].apply(normalize_timepoint)
    meta["response"]   = meta["response"].apply(normalize_response_label)
    keep = meta["therapy"].isin(THERAPIES) & meta["timepoint"].isin(TIMEPOINTS) & (meta["patient_id"]!="UNK")
    meta = meta.loc[keep].copy()
    uniq_patients = sorted(meta["patient_id"].unique(),
                           key=lambda x: int(re.search(r"\d+$", x).group()) if re.search(r"\d+$", x) else 9999)
    log(f"Per-cell META rows: {len(meta):,}; unique patients: {len(uniq_patients)}; therapies: {meta['therapy'].nunique()}")
    with open(os.path.join(OUTDIR, "patients_list.txt"), "w") as fh:
        for p in uniq_patients: fh.write(p + "\n")


    # [3/8] cohorts present
    present = (meta.assign(n=1).groupby(["therapy","timepoint"])["n"].count().reset_index())
    present = present[present["n"] > 0]
    cohorts = list(zip(present["therapy"], present["timepoint"]))
    log("[3/8] Cohorts: " + ", ".join([f"{t}/{tp}" for t,tp in cohorts]))

    # collectors
    stage_fracs_all = []  # per-patient CD8_stage fractions
    rbp_by_stage_all = []  # per-patient RBP means by CD8_stage
    all_rankings = []
    # cell-level overview (optional)
    all_means_cells = []
    # patient-level overview for cross-cohort heatmap
    all_means_pat   = []

    for idx,(ther,tp) in enumerate(cohorts, start=1):
        tag = f"{ther.replace('/','-')}_{tp}"
        log(f"[4/8] Cohort {idx}/{len(cohorts)} → {tag}")
        try:
            adata = load_cohort_adata_by_index(MTX, meta, all_cells_set, cell_to_idx,
                                               therapy=ther, timepoint=tp, n_cells_cap=N_CELLS_CAP)
            if adata is None or adata.n_obs == 0:
                log("   (no cells; skipping)");  continue

            if SUBSET_CTL:
                adata = subset_cd8_t(adata)

            log(f"   cells: {adata.n_obs:,} | genes: {adata.n_vars:,}")
            adata = qc_preprocess_keep_all(adata, hvg=HVG)
            ensure_exhaustion(adata)
            ensure_rbp_module(adata, CORE_RBPS)
            adata = add_pseudotime(adata)

            # ---- CD8_stage labeling and oriented DPT + plots
            try:
                adata = label_cd8_stages(adata)
                _sanitize_cd8_stage_categories(adata)
                adata = compute_oriented_dpt(adata, exclude_label="T-RM", stages_col="CD8_stage")
                plot_cd8_stage_suite(adata, tag=tag, outdir=OUTDIR)

                # NEW: per-cohort per-patient CD8_stage fractions & RBP-by-stage patient means
                try:
                    frac_df = _per_patient_stage_fractions(adata)
                    if not frac_df.empty:
                        stage_fracs_all.append(frac_df)
                except Exception as e:
                    log(f"   [warn] collecting stage fractions failed: {e}")
                try:
                    present_rbps = [g for g in CORE_RBPS if g in adata.var_names]
                    rbp_df = _per_patient_rbp_means_by_stage(adata, present_rbps)
                    if not rbp_df.empty:
                        rbp_by_stage_all.append(rbp_df)
                except Exception as e:
                    log(f"   [warn] collecting RBP-by-stage means failed: {e}")

            except Exception as e:
                log(f"   [warn] CD8_stage/DPT suite failed: {e}")

            present_rbps = [g for g in CORE_RBPS if g in adata.var_names]
            log(f"   CORE_RBPS present: {len(present_rbps)} / {len(CORE_RBPS)}")

            # RBP ranking (optional ordering)
            ranking = rank_rbps(adata, rbp_list=CORE_RBPS)
            if not ranking.empty:
                ranking["therapy"]   = ther
                ranking["timepoint"] = tp
                ranking.to_csv(os.path.join(OUTDIR, f"rbp_candidates_{tag}.csv"), index=False)
                all_rankings.append(ranking)

            # Save cohort AnnData
            _sanitize_uns_for_h5ad(adata)
            adata.write_h5ad(os.path.join(OUTDIR, f"melanoma_{tag}.h5ad"), compression="gzip")

            # ---- Context UMAPs
            base_colors = [c for c in ["timepoint","response","exhaustion_score","rbp_module","dpt_pseudotime"]
                           if c in adata.obs.columns]
            if base_colors:
                sc.pl.umap(adata[adata.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in adata.obs.columns else adata, color=base_colors, show=False)
                plt.savefig(os.path.join(OUTDIR, f"umap_{tag}_context.pdf"), bbox_inches="tight"); plt.close()

            # ---- [6/8] RBP visualizations – ALL RBPs, paginated
            log("[6/8] RBP visualizations…")
            # order genes: by ranking when available
            if not ranking.empty:
                ordered = [g for g in ranking["rbp"].tolist() if g in present_rbps]
                ordered += [g for g in present_rbps if g not in ordered]
            else:
                ordered = present_rbps

            # module + exhaustion overview
            sc.pl.umap(adata[adata.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in adata.obs.columns else adata, color=[c for c in ["rbp_module","exhaustion_score"] if c in adata.obs.columns], show=False)
            plt.savefig(os.path.join(OUTDIR, f"umap_{tag}_rbpModule+exhaustion.pdf"), bbox_inches="tight"); plt.close()

            # paginate RBPs in chunks of 6
            def chunked(lst, n=6):
                for i in range(0, len(lst), n): yield lst[i:i+n]
            page = 1
            for chunk in chunked(ordered, n=6):
                sc.pl.umap(adata[adata.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in adata.obs.columns else adata, color=chunk, show=False)
                plt.savefig(os.path.join(OUTDIR, f"umap_{tag}_RBPs_page{page:02d}.pdf"),
                            bbox_inches="tight"); plt.close()
                page += 1

            # split by response
            for grp in ["Responder","Non-responder"]:
                sub = adata[adata.obs["response"] == grp].copy()
                if sub.n_obs >= 50 and len(ordered) > 0:
                    colors_first = (["rbp_module"] if "rbp_module" in sub.obs.columns else []) + ordered[:5]
                    if colors_first:
                        sc.pl.umap(sub[sub.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in sub.obs.columns else sub, color=colors_first, show=False)
                        plt.savefig(os.path.join(OUTDIR, f"umap_{tag}_RBPs_{grp}_page01.pdf"),
                                    bbox_inches="tight"); plt.close()
                    rem = ordered[5:]
                    if rem:
                        p = 2
                        for chunk in chunked(rem, n=6):
                            sc.pl.umap(sub[sub.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in sub.obs.columns else sub, color=chunk, show=False)
                            plt.savefig(os.path.join(OUTDIR, f"umap_{tag}_RBPs_{grp}_page{p:02d}.pdf"),
                                        bbox_inches="tight"); plt.close()
                            p += 1

            # ---- [6b/8] Tex state labeling + UMAP
            adata = label_tex_states(adata)
            if "tex_state" in adata.obs.columns:
                sc.pl.umap(adata[adata.obs["is_CD8_T"].astype(bool)] if "is_CD8_T" in adata.obs.columns else adata, color=["tex_state"], show=False)
                plt.savefig(os.path.join(OUTDIR, f"umap_{tag}_texStates.pdf"), bbox_inches="tight"); plt.close()

            # ---- Quick overviews
            for resp in ["Responder","Non-responder"]:
                sub = adata[adata.obs["response"] == resp].copy()
                if sub.n_obs >= 20 and present_rbps:
                    means = {g: float(np.mean(sub[:, g].X)) for g in present_rbps}
                    means.update(dict(therapy=ther, timepoint=tp, response=resp, n_cells=int(sub.n_obs)))
                    all_means_cells.append(means)

            for pid, obs_g in adata.obs.groupby("patient_id"):
                sub = adata[obs_g.index, :].copy()
                resp = obs_g["response"].mode(dropna=False).iloc[0]
                means = {g: float(np.mean(sub[:, g].X)) for g in present_rbps}
                means.update(dict(therapy=ther, timepoint=tp, response=resp, patient_id=pid, n_cells=int(sub.n_obs)))
                all_means_pat.append(means)

            # ---- [7/8] Per-patient statistics across CTL CD8+ / pre-Tex / Tex-int / Term-Tex
            log("[7/8] RBP statistics across Tex states (per-patient)…")
            df_pat = patient_state_table(adata, present_rbps)
            if not df_pat.empty:
                # write per-cohort patient means per state
                df_pat.to_csv(os.path.join(OUTDIR, f"rbp_TexState_PATIENTMEAN_{tag}.csv"), index=False)
                # stats per gene
                rows = [state_stats_per_gene(df_pat, g) for g in present_rbps]
                de_df = pd.DataFrame(rows)
                if not de_df.empty:
                    de_df.to_csv(os.path.join(OUTDIR, f"rbp_TexState_DE_patients_{tag}.csv"), index=False)
                # heatmap per cohort
                try:
                    import seaborn as sns
                    agg = (df_pat.groupby("tex_state", as_index=False)
                           .agg({**{g:"mean" for g in present_rbps}, "patient_id":"nunique", "n_cells":"sum"}))
                    agg["group"] = agg["tex_state"] + "|n_pat=" + agg["patient_id"].astype(str) + "|n_cells=" + agg["n_cells"].astype(str)
                    M = agg.set_index("group")[present_rbps]
                    Z = (M - M.mean(axis=0)) / (M.std(axis=0, ddof=0) + 1e-9)
                    g = sns.clustermap(Z, cmap="vlag", linewidths=0.25, figsize=(10, 6), col_cluster=False)
                    plt.savefig(os.path.join(OUTDIR, f"heatmap_RBP_by_TexState_PATIENTMEAN_{tag}.pdf"),
                                bbox_inches="tight"); plt.close()
                except Exception as e:
                    log(f"   [warn] heatmap per cohort skipped: {e}")

        except Exception as e:
            log(f"   ✗ Error in cohort {tag}: {e}")
            continue

    # [5/8] Cross-cohort summary (median within-cohort rank)
    log("[5/8] Cross-cohort summary…")
    if all_rankings:
        big = pd.concat(all_rankings, axis=0, ignore_index=True)
        big["within_cohort_rank"] = big.groupby(["therapy","timepoint"])["score"].rank(ascending=False, method="min")
        summary = (big.groupby("rbp")["within_cohort_rank"].median().sort_values().reset_index()
                   .rename(columns={"within_cohort_rank":"median_rank_across_cohorts"}))
        summary.to_csv(os.path.join(OUTDIR, "rbp_candidates_summary_across_cohorts.csv"), index=False)

    # NEW: exports for CD8_stage fractions and RBP-by-stage (patient-level) with group mean+SD
    try:
        _export_stage_fraction_group_stats(stage_fracs_all, OUTDIR)
        _export_rbp_by_stage_group_stats(rbp_by_stage_all, OUTDIR)
        # Optional figures (mean+/-SD bars)
        _plot_stage_fraction_grouped_bars(stage_fracs_all, OUTDIR)
        genes_for_stage_bars = [g for g in CORE_RBPS]
        _plot_rbp_by_stage_bars(rbp_by_stage_all, genes_for_stage_bars, OUTDIR)
    except Exception as e:
        log(f"   [warn] post-cohort stage/RBP aggregation failed: {e}")
        log("   ✓ Saved cross-cohort summary")
    else:
        log("   [info] No rankings collected; skipping summary.")

    # [8/8] Packaging: cross-cohort tables + global heatmap by therapy/timepoint/response
    log("[8/8] Packaging & heatmaps…")
    if all_means_cells:
        df_cells = pd.DataFrame(all_means_cells)
        cols = ["therapy","timepoint","response","n_cells"] + [g for g in CORE_RBPS if g in df_cells.columns]
        df_cells[cols].to_csv(os.path.join(OUTDIR, "rbp_expression_by_therapy_timepoint_response_cells.csv"), index=False)

    if all_means_pat:
        df_pat_all = pd.DataFrame(all_means_pat)
        cols = ["therapy","timepoint","response","patient_id","n_cells"] + [g for g in CORE_RBPS if g in df_pat_all.columns]
        df_pat_all[cols].to_csv(os.path.join(OUTDIR, "rbp_expression_patientMeans.csv"), index=False)
        # --- NEW: RBP module patient-mean (row-wise mean of available RBP genes)
        rbp_cols_present = [g for g in CORE_RBPS if g in df_pat_all.columns]
        if rbp_cols_present:
            df_pat_all = df_pat_all.copy()
            df_pat_all["rbp_module_pm"] = df_pat_all[rbp_cols_present].mean(axis=1)
            # Save enriched table with module
            df_pat_all.to_csv(os.path.join(OUTDIR, "rbp_expression_patientMeans_withModule.csv"), index=False)
            # --- NEW: mean+/-SD barplots + stats per therapy and gene (including module)
            genes_for_bars = rbp_cols_present + ["rbp_module_pm"]
            try:
                make_prepost_response_barplots(df_pat_all, genes_for_bars, OUTDIR)
            except Exception as e:
                log(f"   [warn] barplots pre/post by response failed: {e}")
        try:
            import seaborn as sns
            gdf = (df_pat_all
                   .groupby(["therapy","timepoint","response"], as_index=False)
                   .agg({**{g:"mean" for g in CORE_RBPS if g in df_pat_all.columns}, "patient_id":"nunique","n_cells":"sum"}))
            rbps_present = [g for g in CORE_RBPS if g in gdf.columns]
            if rbps_present:
                gdf["group"] = (gdf["therapy"] + "|" + gdf["timepoint"] + "|" + gdf["response"] +
                                "|n_pat=" + gdf["patient_id"].astype(str) + "|n_cells=" + gdf["n_cells"].astype(str))
                M = gdf.set_index("group")[rbps_present]
                Z = (M - M.mean(axis=0)) / (M.std(axis=0, ddof=0) + 1e-9)
                import seaborn as sns
                g = sns.clustermap(Z, cmap="vlag", linewidths=0.25, figsize=(13, 8), col_cluster=False)
                plt.savefig(os.path.join(OUTDIR, "heatmap_RBP_by_therapy_timepoint_response_PATIENTMEAN.pdf"),
                            bbox_inches="tight"); plt.close()
        except Exception as e:
            log(f"   [warn] global heatmap skipped: {e}")

    log(f"Done. See {OUTDIR}/ for outputs.")

if __name__=="__main__":
    import argparse, json

    parser = argparse.ArgumentParser(description="RBP-centric analysis – Sade-Feldman melanoma")
    parser.add_argument("--mtx", type=str, default=MTX,
                        help=f"Expression matrix path (default: {MTX})")
    parser.add_argument("--meta", type=str, default=META,
                        help=f"META table path (default: {META})")
    parser.add_argument("--outdir", type=str, default=OUTDIR,
                        help=f"Output directory (default: {OUTDIR})")
    parser.add_argument("--subset-ctl", action="store_true", default=SUBSET_CTL,
                        help="Prefilter to CD8 T cells before analysis (heuristic).")
    parser.add_argument("--rbps", type=str, default=None,
                        help="Comma/space/semicolon-separated list or JSON array of gene symbols to use as RBPs (overrides default).")
    parser.add_argument("--rbp-file", type=str, default=None,
                        help="Path to a file containing RBPs (one per line, or first column in CSV/TSV).")
    parser.add_argument("--require-11", action="store_true",
                        help="If set, asserts that exactly 11 RBPs are provided via --rbps/--rbp-file.")
    args = parser.parse_args()

    # Override IO and prefs
    MTX = args.mtx
    META = args.meta
    OUTDIR = args.outdir
    SUBSET_CTL = args.subset_ctl

    # Parse RBPs if provided
    def _clean_sym(x: str) -> str:
        return re.sub(r"[^A-Za-z0-9_.-]+", "", x.strip()).upper()

    selected = []
    if args.rbps:
        s = args.rbps.strip()
        try:
            if s.startswith("["):
                arr = json.loads(s)
                if isinstance(arr, list):
                    selected.extend([_clean_sym(x) for x in arr if str(x).strip()])
            else:
                parts = re.split(r"[\s,;]+", s)
                selected.extend([_clean_sym(x) for x in parts if x.strip()])
        except Exception as e:
            log(f"[warn] Failed to parse --rbps: {e}")
    if args.rbp_file and os.path.exists(args.rbp_file):
        try:
            with open(args.rbp_file, "r", encoding="utf-8") as fh:
                for line in fh:
                    if not line.strip() or line.strip().startswith("#"):
                        continue
                    # If delimited, take the first cell
                    cell = re.split(r"[\s,;\t]+", line.strip())[0]
                    if cell: selected.append(_clean_sym(cell))
        except Exception as e:
            log(f"[warn] Failed to read --rbp-file: {e}")

    # De-duplicate while preserving order
    if selected:
        seen = set()
        sel_unique = []
        for g in selected:
            if g and g not in seen:
                sel_unique.append(g); seen.add(g)
        if args.require_11 and len(sel_unique) != 11:
            raise ValueError(f"--require-11 set, but {len(sel_unique)} RBPs were provided: {sel_unique}")
        CORE_RBPS = sel_unique
        log(f"[args] Using RBPs from args/file (n={len(CORE_RBPS)}): {', '.join(CORE_RBPS)}")
    else:
        CORE_RBPS = DEFAULT_CORE_RBPS.copy()
        log(f"[args] Using default RBPs (n={len(CORE_RBPS)}): {', '.join(CORE_RBPS)}")

    ensure_outdir()
    if not os.path.exists(MTX):  raise FileNotFoundError(f"Expression matrix not found: {MTX}")
    if not os.path.exists(META): raise FileNotFoundError(f"META file not found: {META}")
    main()



    # -*- coding: utf-8 -*-
"""
Patch your working analysis script to:
  - ensure CD8_stage categories are unique & ordered (fixes 'Categorical categories must be unique'),
  - add a CD3+CD8+ mask helper,
  - compute per-patient CD8 stage fractions (CD8-only) + group mean±SD by therapy×timepoint×response,
  - compute per-patient RBP means within CD8 stages + group mean±SD,
  - generate barplots for stage fractions and RBP-by-stage (Pre vs Post, mean±SD),
  - generate Pre/Post × Responder/Non-responder barplots from CD8-only patient means,
  - keep all your existing outputs intact.

Usage:
  python3 apply_cd8_only_patch.py analysisPLSWORK-4_full_corrected.py analysisPLSWORK-5_full_CD8only.py
"""
import sys, re, os
from pathlib import Path

def insert_once(haystack, needle, before=None, after=None):
    if needle.strip() in haystack:
        return haystack, False
    if before and before in haystack:
        i = haystack.find(before)
        return haystack[:i] + needle + haystack[i:], True
    if after and after in haystack:
        i = haystack.find(after) + len(after)
        return haystack[:i] + needle + haystack[i:], True
    return haystack, False

def add_constants(text):
    if "CD8_STAGE_MARKERS" in text:
        return text, False
    block = r"""
# ---------------------------
# CD8 stage marker dictionary and display orders
# ---------------------------
CD8_STAGE_MARKERS = {
    "CD3D":"CD3D","CD3E":"CD3E",
    "CD8A":"CD8A","CD8B":"CD8B",
    "SELL":"SELL","CCR7":"CCR7",
    "TCF7":"TCF7","IL7R":"IL7R",
    "CXCR6":"CXCR6","ITGAE":"ITGAE","CD69":"CD69",
    "GZMB":"GZMB","PRF1":"PRF1",
    "TOX":"TOX","EOMES":"EOMES","CXCL13":"CXCL13","ENTPD1":"ENTPD1"
}
CD8_STAGE_ORDER_ALL = ["Non-CD8","Naive/T_SCM","T-CM","Stem-like","T-EFF/CTL","T-EM/TEMRA","Tex-int","Prog-Tex","Term-Tex","T-RM","TRM-CD69only","Other_CD8","CTL_CD8"]
CD8_STAGE_ORDER_MAIN = ["Stem-like","T-CM","T-EFF/CTL","Tex-int","Prog-Tex","Term-Tex"]
CD8_STAGE_ORDER_SIMPLE = ["CTL_CD8","Tex-int","Term-Tex"]
"""
    text, did = insert_once(text, block, before="def label_cd8_stages(", after="def main():")
    return text, did

def fix_labeler_duplicates(text):
    pat = re.compile(r'adata\.obs\["CD8_stage"\]\s*=\s*pd\.Categorical\(\s*labels\s*,\s*categories\s*=\s*\(\s*\[\s*"Non-CD8"\s*\]\s*\+\s*CD8_STAGE_ORDER_ALL\s*\)\s*,\s*ordered\s*=\s*True\s*\)')
    if not pat.search(text):
        return text, False
    repl = (
        'cats = []\n'
        '    for _c in (["Non-CD8"] + CD8_STAGE_ORDER_ALL):\n'
        '        if _c not in cats:\n'
        '            cats.append(_c)\n'
        '    adata.obs["CD8_stage"] = pd.Categorical(labels, categories=cats, ordered=True)'
    )
    return pat.sub(repl, text), True

def add_sanitizer(text):
    if "_sanitize_cd8_stage_categories" in text:
        return text, False
    block = r"""
# ---------------------------
# Ensure CD8_stage categorical has unique, ordered categories
# ---------------------------
def _sanitize_cd8_stage_categories(adata: ad.AnnData):
    try:
        if "CD8_stage" in adata.obs.columns:
            s = adata.obs["CD8_stage"].astype(str)
            cats = [c for c in CD8_STAGE_ORDER_ALL if c in set(s.values)]
            seen=set(); cats_unique=[]
            for c in cats:
                if c not in seen:
                    cats_unique.append(c); seen.add(c)
            if not cats_unique:
                cats_unique = list(pd.unique(s.values))
            adata.obs["CD8_stage"] = pd.Categorical(s, categories=cats_unique, ordered=True)
    except Exception as e:
        log(f"   [warn] sanitize CD8_stage categories failed: {e}")
"""
    text, did = insert_once(text, block, before="def label_cd8_stages(", after="def main():")
    return text, did

def call_sanitizer_after_label(text):
    lines = text.splitlines()
    out, added = [], 0
    i = 0
    while i < len(lines):
        ln = lines[i]
        out.append(ln)
        if "adata = label_cd8_stages(adata)" in ln:
            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                out.append(lines[j]); j += 1
            has_next = j < len(lines)
            if not has_next or "_sanitize_cd8_stage_categories(adata)" not in lines[j]:
                indent = re.match(r"^(\s*)", ln).group(1) or ""
                out.append(indent + "_sanitize_cd8_stage_categories(adata)")
                added += 1
        i += 1
    return "\n".join(out), (added > 0)

def add_cd8_mask(text):
    if "_mask_cd3_cd8" in text:
        return text, False
    block = r"""
# ---------------------------
# Helper: CD3+CD8+ mask (prefer existing is_CD8_T, else markers and stage!=Non-CD8)
# ---------------------------
def _mask_cd3_cd8(adata: ad.AnnData) -> np.ndarray:
    try:
        if "is_CD8_T" in adata.obs.columns:
            return adata.obs["is_CD8_T"].astype(bool).values
    except Exception:
        pass
    n = adata.n_obs
    stage_ok = np.ones(n, dtype=bool)
    if "CD8_stage" in adata.obs.columns:
        stage_ok = (adata.obs["CD8_stage"].astype(str).values != "Non-CD8")
    cd3 = np.zeros(n, dtype=bool)
    for g in ["CD3D","CD3E"]:
        if g in adata.var_names:
            Xg = adata[:, g].X
            Xg = Xg.A if hasattr(Xg, "A") else Xg
            cd3 |= (np.ravel(Xg) > 0)
    cd8 = np.zeros(n, dtype=bool)
    for g in ["CD8A","CD8B"]:
        if g in adata.var_names:
            Xg = adata[:, g].X
            Xg = Xg.A if hasattr(Xg, "A") else Xg
            cd8 |= (np.ravel(Xg) > 0)
    return (cd3 & cd8 & stage_ok)
"""
    return insert_once(text, block, before="# ---------------------------\n# MAIN", after="def main():")

def ensure_collectors(text):
    if "stage_fracs_all =" in text and "rbp_by_stage_all =" in text:
        return text, False
    return text.replace("    # collectors",
                        "    # collectors\n    stage_fracs_all = []  # CD8 stage fractions\n    rbp_by_stage_all = []  # RBP means per CD8 stage"), True

def add_stage_rbp_helpers(text):
    if "def _per_patient_stage_fractions(" in text and "def _per_patient_rbp_means_by_stage(" in text:
        return text, False
    block = r"""
# ---------------------------
# Per-patient CD8_stage fractions and per-stage RBP means; group mean+/-SD exporters
# ---------------------------

def _per_patient_stage_fractions(adata: ad.AnnData) -> pd.DataFrame:
    if "CD8_stage" not in adata.obs.columns:
        return pd.DataFrame()
    obs = adata.obs[["patient_id","response","CD8_stage"]].copy()
    try:
        mask = _mask_cd3_cd8(adata)
        obs = obs.iloc[np.where(mask)[0]]
    except Exception:
        pass
    obs = obs[obs["CD8_stage"].astype(str) != "Non-CD8"]
    comp = (obs.groupby(["patient_id","response","CD8_stage"], observed=True).size().unstack(fill_value=0))
    frac = comp.div(comp.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0).reset_index()
    ther = adata.obs["therapy"].mode(dropna=False).iat[0] if "therapy" in adata.obs.columns and not adata.obs["therapy"].empty else "UNK"
    tp   = adata.obs["timepoint"].mode(dropna=False).iat[0] if "timepoint" in adata.obs.columns and not adata.obs["timepoint"].empty else "UNK"
    frac["therapy"] = ther; frac["timepoint"] = tp
    return frac

def _per_patient_rbp_means_by_stage(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.var_names]
    if not genes or "CD8_stage" not in adata.obs.columns:
        return pd.DataFrame()
    mask = _mask_cd3_cd8(adata)
    obs_cd8 = adata.obs.loc[np.where(mask)[0]]
    rows = []
    for (pid, resp, stage), idx in obs_cd8.groupby(["patient_id","response","CD8_stage"]).groups.items():
        sub = adata[idx, :]
        means = {g: float(np.mean(sub[:, g].X)) for g in genes}
        rows.append({"patient_id": pid, "response": resp, "CD8_stage": stage, **means, "n_cells": int(sub.n_obs)})
    df = pd.DataFrame(rows)
    ther = adata.obs["therapy"].mode(dropna=False).iat[0] if "therapy" in adata.obs.columns and not adata.obs["therapy"].empty else "UNK"
    tp   = adata.obs["timepoint"].mode(dropna=False).iat[0] if "timepoint" in adata.obs.columns and not adata.obs["timepoint"].empty else "UNK"
    df["therapy"] = ther; df["timepoint"] = tp
    return df

def _export_stage_fraction_group_stats(stage_fracs_all: list[pd.DataFrame], outdir: str):
    if not stage_fracs_all: return
    df = pd.concat(stage_fracs_all, ignore_index=True)
    df.to_csv(os.path.join(outdir, "cd8_stageFractions_patient_all.csv"), index=False)
    stage_cols = [c for c in df.columns if c not in ["patient_id","response","therapy","timepoint"]]
    L = df.melt(id_vars=["therapy","timepoint","response","patient_id"], value_vars=stage_cols,
                var_name="CD8_stage", value_name="fraction")
    g = (L.groupby(["therapy","timepoint","response","CD8_stage"], observed=True)
           .agg(mean_fraction=("fraction","mean"),
                sd_fraction=("fraction", lambda x: float(np.std(x, ddof=1)) if len(x)>=2 else np.nan),
                n_patients=("fraction","count"))
           .reset_index())
    g.to_csv(os.path.join(outdir, "cd8_stageFractions_groupStats.csv"), index=False)

def _export_rbp_by_stage_group_stats(rbp_by_stage_all: list[pd.DataFrame], outdir: str):
    if not rbp_by_stage_all: return
    df = pd.concat(rbp_by_stage_all, ignore_index=True)
    df.to_csv(os.path.join(outdir, "rbp_byCD8stage_patientMeans.csv"), index=False)
    gene_cols = [c for c in df.columns if c not in ["therapy","timepoint","response","patient_id","CD8_stage","n_cells"]]
    rows = []
    for (ther, tp, resp, stage), sub in df.groupby(["therapy","timepoint","response","CD8_stage"], observed=True):
        n = sub["patient_id"].nunique()
        rec = {"therapy": ther, "timepoint": tp, "response": resp, "CD8_stage": stage, "n_patients": int(n)}
        for gname in gene_cols:
            vals = sub[gname].dropna().values
            rec[f"mean_{gname}"] = float(np.mean(vals)) if vals.size else np.nan
            rec[f"sd_{gname}"]   = float(np.std(vals, ddof=1)) if vals.size>=2 else np.nan
        rows.append(rec)
    stats = pd.DataFrame(rows)
    stats.to_csv(os.path.join(outdir, "rbp_byCD8stage_groupStats.csv"), index=False)

def _plot_stage_fraction_grouped_bars(stage_fracs_all: list[pd.DataFrame], outdir: str):
    if not stage_fracs_all: return
    df = pd.concat(stage_fracs_all, ignore_index=True)
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    stage_cols = [c for c in df.columns if c not in ["patient_id","response","therapy","timepoint"]]
    for therapy, dft in df.groupby("therapy"):
        for resp, dfr in dft.groupby("response"):
            pre = dfr[dfr["timepoint"]=="Pre"]; post = dfr[dfr["timepoint"]=="Post"]
            mean_pre = pre[stage_cols].mean(axis=0); sd_pre = pre[stage_cols].std(axis=0, ddof=1); n_pre = pre["patient_id"].nunique()
            mean_post= post[stage_cols].mean(axis=0); sd_post= post[stage_cols].std(axis=0, ddof=1); n_post= post["patient_id"].nunique()
            x = np.arange(len(stage_cols)); w = 0.38
            fig, ax = plt.subplots(figsize=(max(7.0, 0.55*len(stage_cols)+2), 4.6))
            ax.bar(x-w/2, mean_pre.values, width=w, yerr=sd_pre.values, capsize=3, edgecolor="black", label="Pre")
            ax.bar(x+w/2, mean_post.values, width=w, yerr=sd_post.values, capsize=3, edgecolor="black", label="Post")
            ax.set_xticks(x); ax.set_xticklabels(stage_cols, rotation=30, ha="right")
            ax.set_ylabel("CD8_stage fraction (mean +/- SD)")
            ax.set_title(f"{therapy} — {resp}")
            ax.legend()
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5, prune=None))
            plt.tight_layout()
            safe_t = re.sub(r"[^A-Za-z0-9_.-]+","_", str(therapy)); safe_r = re.sub(r"[^A-Za-z0-9_.-]+","_", str(resp))
            plt.savefig(os.path.join(outdir, f"barplot_CD8stageFractions_{safe_t}_{safe_r}_PrePost.pdf"), bbox_inches="tight"); plt.close()

def _plot_rbp_by_stage_bars(rbp_by_stage_all: list[pd.DataFrame], genes: list[str], outdir: str):
    if not rbp_by_stage_all: return
    df = pd.concat(rbp_by_stage_all, ignore_index=True)
    import matplotlib.pyplot as plt
    stage_order = [s for s in CD8_STAGE_ORDER_ALL if s in df["CD8_stage"].unique().tolist() and s!="Non-CD8"]
    for therapy, dft in df.groupby("therapy"):
        for resp, dfr in dft.groupby("response"):
            for gene in genes:
                if gene not in df.columns: continue
                rows = []
                for stage in stage_order:
                    pre = dfr[(dfr["timepoint"]=="Pre") & (dfr["CD8_stage"]==stage)][gene].dropna().values
                    post= dfr[(dfr["timepoint"]=="Post") & (dfr["CD8_stage"]==stage)][gene].dropna().values
                    rows.append((stage, float(pre.mean()) if pre.size else float("nan"),
                                      float(pre.std(ddof=1)) if pre.size>=2 else float("nan"),
                                      float(post.mean()) if post.size else float("nan"),
                                      float(post.std(ddof=1)) if post.size>=2 else float("nan")))
                if not rows: continue
                stages = [r[0] for r in rows]
                import numpy as np
                mean_pre = np.array([r[1] for r in rows]); sd_pre = np.array([r[2] for r in rows])
                mean_post= np.array([r[3] for r in rows]); sd_post= np.array([r[4] for r in rows])
                x = np.arange(len(stages)); w = 0.38
                fig, ax = plt.subplots(figsize=(max(7.0, 0.55*len(stages)+2), 4.6))
                ax.bar(x-w/2, mean_pre, width=w, yerr=sd_pre, capsize=3, edgecolor="black", label="Pre")
                ax.bar(x+w/2, mean_post, width=w, yerr=sd_post, capsize=3, edgecolor="black", label="Post")
                ax.set_xticks(x); ax.set_xticklabels(stages, rotation=30, ha="right")
                ax.set_ylabel(f"{gene} (mean +/- SD, patient-level, CD8 only)")
                ax.set_title(f"{therapy} — {resp} — {gene}")
                ax.legend()
                plt.tight_layout()
                safe_t = re.sub(r"[^A-Za-z0-9_.-]+","_", str(therapy)); safe_r = re.sub(r"[^A-Za-z0-9_.-]+","_", str(resp)); safe_g = re.sub(r"[^A-Za-z0-9_.-]+","_", str(gene))
                plt.savefig(os.path.join(outdir, f"barplot_RBP_byCD8stage_{safe_g}_{safe_t}_{safe_r}_PrePost.pdf"), bbox_inches="tight"); plt.close()
"""
    return insert_once(text, block, before="# ---------------------------\n# MAIN", after="def main():")

def wire_collectors(text):
    anchor = "plot_cd8_stage_suite(adata, tag=tag, outdir=OUTDIR)"
    if anchor not in text or "collecting stage fractions failed" in text:
        return text, False
    inject = r"""
                # Collect per-patient CD8 stage fractions (CD8 only) and per-stage RBP means
                try:
                    frac_df = _per_patient_stage_fractions(adata)
                    if not frac_df.empty:
                        stage_fracs_all.append(frac_df)
                except Exception as e:
                    log(f"   [warn] collecting stage fractions failed: {e}")
                try:
                    rbp_df = _per_patient_rbp_means_by_stage(adata, present_rbps)
                    if not rbp_df.empty:
                        rbp_by_stage_all.append(rbp_df)
                except Exception as e:
                    log(f"   [warn] collecting RBP-by-stage means failed: {e}")
"""
    return text.replace(anchor, anchor + "\n" + inject), True

def add_aggregation_block(text):
    anchor = "[5/8] Cross-cohort summary"
    if anchor not in text or "cd8_stageFractions_groupStats.csv" in text:
        return text, False
    pos = text.find(anchor)
    line_start = text.rfind("\n", 0, pos) + 1
    block = r"""
    # Aggregate and export CD8 stage fractions and RBP-by-stage group stats; draw summary barplots
    try:
        _export_stage_fraction_group_stats(stage_fracs_all, OUTDIR)
        _export_rbp_by_stage_group_stats(rbp_by_stage_all, OUTDIR)
        _plot_stage_fraction_grouped_bars(stage_fracs_all, OUTDIR)
        genes_for_stage_bars = CORE_RBPS.copy()
        if rbp_by_stage_all and "rbp_module" in rbp_by_stage_all[0].columns:
            if "rbp_module" not in genes_for_stage_bars:
                genes_for_stage_bars.append("rbp_module")
        _plot_rbp_by_stage_bars(rbp_by_stage_all, genes_for_stage_bars, OUTDIR)
    except Exception as e:
        log(f"   [warn] post-cohort stage/RBP aggregation failed: {e}")
"""
    return text[:line_start] + block + text[line_start:], True

def add_cd8only_patientmeans_barplots(text):
    anchor = 'df_pat_all[cols].to_csv(os.path.join(OUTDIR, "rbp_expression_patientMeans.csv"), index=False)'
    if anchor not in text or "rbp_expression_patientMeans_CD8only.csv" in text:
        return text, False
    block = r"""
        # CD8-only patient means and barplots (Pre/Post x Response), built from per-stage collector
        try:
            rbp_cols_present = [g for g in CORE_RBPS if g in df_pat_all.columns]
            df_cd8 = None
            if 'rbp_by_stage_all' in globals() and rbp_by_stage_all:
                tmp = pd.concat(rbp_by_stage_all, ignore_index=True)
                df_cd8 = (tmp.groupby(["therapy","timepoint","response","patient_id"], observed=True)[rbp_cols_present]
                             .mean().reset_index())
            if df_cd8 is not None and not df_cd8.empty:
                df_cd8 = df_cd8.copy()
                if rbp_cols_present:
                    df_cd8["rbp_module_pm"] = df_cd8[rbp_cols_present].mean(axis=1)
                df_cd8.to_csv(os.path.join(OUTDIR, "rbp_expression_patientMeans_CD8only.csv"), index=False)
                genes_for_bars = rbp_cols_present + (["rbp_module_pm"] if "rbp_module_pm" in df_cd8.columns else [])
                try:
                    make_prepost_response_barplots(df_cd8, genes_for_bars, OUTDIR)
                except Exception as e:
                    log(f"   [warn] CD8-only barplots failed: {e}")
        except Exception as e:
            log(f"   [warn] CD8-only patient means aggregation failed: {e}")
"""
    return text.replace(anchor, anchor + block), True

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 apply_cd8_only_patch.py <input.py> <output.py>")
        sys.exit(1)
    inp = Path(sys.argv[1]); outp = Path(sys.argv[2])
    text = inp.read_text(encoding="utf-8")

    steps = []
    for fn in [add_constants, fix_labeler_duplicates, add_sanitizer, call_sanitizer_after_label,
               add_cd8_mask, ensure_collectors, add_stage_rbp_helpers, wire_collectors,
               add_aggregation_block, add_cd8only_patientmeans_barplots]:
        text, did = fn(text)
        steps.append((fn.__name__, did))

    outp.write_text(text, encoding="utf-8")
    print("Wrote:", outp)
    print("Steps:", steps)

if __name__ == "__main__":
    main()
