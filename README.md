# Dissecting RNA-binding proteins' implication in exhausted CD8+ T-cells in the cancer context (ONGOING).
- Aim: exploring the RBPs implication in T cell exhaustion using public datasets
- How to run notebooks: run the analysisRBP.py in Terminal (developped and tested on MacOS M1)
- Current version of the script uses the Sade-Feldman (2021) scRNA-seq dataset on 32 melanoma patients

Background: T cell exhaustion is a limitation to durable immune benefit from the immune checkpoint blockade (ICB) in cancer patients. Exhausted CD8+ (Tex) tumour-infiltrating lymphocytes (TILs) progressively lose cytotoxic function and acquire a distinct transcriptional and epigenetic state; PD-1/CTLA-4 blockade only partially reinvigorates these cells and many patients do not respond (Wherry, 2015). Understanding the regulating mechanisms from CD8+ T cells towards Tex is therefore critical for improving the ICB outcome.

Hypothesis: In addition to metabolic and epigenetic control of CD8+ T cell fate, post‑transcriptional regulation by RNA‑binding proteins (RBPs) seems an under‑appreciated layer that can bias CTL→Tex trajectories in some human tumours. Here, we focus on eleven RBPs: ZFP36, ZFP36L1, ZFP36L2, ZC3H12A (Regnase‑1), RC3H1/RC3H2 (Roquin‑1/2), PCBP1, ELAVL1 (HuR), HNRNPLL, TIA1 and TIA1B because they span ARE‑mediated mRNA decay, Regnase–Roquin endonucleolytic repression, activation‑coupled mRNA stabilisation, stress‑granule–mediated translational silencing and activation‑relevant alternative splicing; each has direct evidence in T cells (Glasmacher et al., 2010; Vogel et al., 2013; Jeltsch et al., 2014; Mino et al., 2015; Moore et al., 2018; Ansa‑Addo et al., 2020; Oberdoerffer et al., 2008; Wu et al., 2006; Techasintana et al., 2015; Karginov et al., 2019; Piecyk et al., 2000).

Approach (in silico): We re‑analyse single‑cell RNA‑seq from melanoma tumours and quantify RBP activity with gene‑module scores and per‑gene overlays on UMAP manifolds, alongside a validated exhaustion score and pseudotime. We stratify by therapy (anti‑CTLA‑4, anti‑PD‑1, and anti‑CTLA‑4+PD‑1) and timepoint (Pre/Post), and compare responders vs. non‑responders. Within CD8+ T cells we define Tex states (CTL, pre‑Tex, Tex‑intermediate, terminal Tex) using exhaustion‑score quartiles and perform per‑patient statistics (Kruskal–Wallis across states; pairwise Mann–Whitney with BH‑FDR correction) to identify RBPs whose expression tracks Tex progression and response.

Dataset choice: We use the Sade‑Feldman melanoma cohort because it is large, clinically annotated (Pre/Post, therapy, response) and widely adopted, providing a representative TIL landscape for benchmarking (Sade‑Feldman et al., 2021). Importantly, by operating at the single‑cell transcriptomic level specifically in melanoma, our analysis minimises confounding from the heterogeneous tumour genomic landscape, enabling a focused test of RBP control of CD8+ Tex programmes under ICB.


Below is a concise “Methods”:

Input & pre-processing
- Matrix: 16,291 cells × 55,737 genes; 32 patients; 3 therapies; timepoints Pre/Post.
- Per cohort (therapy × timepoint), we filter genes detected in ≥3 cells; then compute highly variable genes (HVGs) for dimensionality reduction.

Linear embedding
- PCA on scaled HVGs → 50 PCs. These PCs summarize transcriptome variation and suppress noise.
- KNN graph on PCA space (30 PCs) using Euclidean distances; edges reflect local similarity.
	- Mathematically: build a neighborhood graph G=(V,E)G=(V,E) with adjacency AA, then normalize to connectivities (weights), typically symmetrized & renormalized (Scanpy pp.neighbors).

Nonlinear visualization
- UMAP on the KNN graph: UMAP optimizes a low‑dimensional embedding that preserves local topology by minimizing the cross-entropy between a fuzzy simplicial set in high‑D and low‑D. It’s purely for visualization (not statistics).
Exhaustion & RBP module scores
- Gene set scoring (Scanpy tl.score_genes) computes a per-cell z‑scored signature; we used (i) an exhaustion score and (ii) an RBP module (mean across present RBPs).

Diffusion Maps & DPT (pseudotime)
- We compute Diffusion Maps using nDC=15nDC=15 components.
	- Let PP be the Markov transition matrix from the KNN graph (row-normalized connectivities).
	- Eigen-decompose P=ΨΛΦ⊤P=ΨΛΦ⊤. The diffusion components (DCs) are columns of ΨΨ, ordered by eigenvalues λ1≥λ2≥⋯λ1≥λ2≥⋯. They parametrize low‑dimensional manifold geometry by random‑walk diffusion.
- DPT (Diffusion Pseudotime) measures a diffusion distance from a chosen root cell: roughly DPT(i)=∑k=2Kλkt(ψk(i)−ψk(root))2DPT(i)=∑k=2Kλkt(ψk(i)−ψk(root))2 for some kernel time tt (conceptually). The idea: if CTL→Tex progression is manifold-like and connected, DPT should increase monotonically from stem/TCF7-high through T‑CM/T‑EFF to intermediate/terminal Tex.
- Root selection should be a biologically early CD8 subset (TCF7/IL7R-high). If the graph is disconnected, DPT can become undefined or behave oddly (e.g., inf means or repeated 1.0/–1.0 eigenvalues).

Aggregation & statistics
- Patient‑level means (CD8-only options): we aggregate expression by patient (and later by CD8 stages) and compute mean ± SD across patients within each group (therapy × timepoint × response).
- Barplots compare Pre vs Post, Responder vs Non‑responder, and within CD8 stages.
