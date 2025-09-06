
# Why these 11 RNA‑binding proteins (RBPs) for CTL → Tex in melanoma?  
Target RBPs: ZFP36 (TTP), ZFP36L1, ZFP36L2, ZC3H12A (Regnase‑1), RC3H1 (Roquin‑1), RC3H2 (Roquin‑2), PCBP1, ELAVL1 (HuR), HNRNPLL, TIA1 and TIA1B (an isoform annotation of TIA1 in this dataset).

---

## 1) Rationale in one page

- Checkpoint therapy modulates chronic T‑cell activation. Tex biology is driven by sustained antigen and inflammatory signals. That program is controlled post‑transcriptionally by RBPs that degrade, stabilize or repress translation of immune mRNAs.  
- Our 11 RBPs cover the four levers of post‑transcriptional control that repeatedly emerge in T cells:  
  - ARE‑mediated decay: ZFP36 family (ZFP36/L1/L2) attenuate T‑cell activation and effector output by directly binding cytokine and signaling mRNAs (IL2, IFNG, TNF; CD28/ICOS/CTLA4, NFAT/NF‑κB components). Loss of ZFP36 accelerates activation and promotes an exhaustion‑like program downstream. (Moore et al., 2018; Petkau et al., 2022)  
  -   Regnase–Roquin axis:   ZC3H12A/Regnase‑1 (endoribonuclease) and RC3H1/2/Roquin‑1/2 (ROQ‑domain RBPs) cooperate to repress costimulatory/checkpoint mRNAs (e.g., ICOS,   OX40) and other inflammatory targets, shaping T‑cell fate. Disrupting their interaction induces autoimmunity and   enhances antitumor responses  , linking this pathway directly to immunotherapy. (Glasmacher  et al. , 2010; Vogel  et al. , 2013; Jeltsch  et al. , 2014; Mino  et al. , 2015; Behrens  et al. , 2021)  
  -   Splicing gatekeeper:     HNRNPLL   drives the CD45 (PTPRC) isoform switch essential for naïve→memory/effector transitions and activation thresholds in T cells. (Oberdoerffer  et al. , 2008; Wu  et al. , 2006)  
  -   Stability/translation & stress:  ELAVL1/HuR stabilizes activation‑induced immune transcripts and promotes T‑cell survival; TIA1/TIA1B enforces stress‑granule–mediated translational silencing of inflammatory mRNAs (e.g., TNF). (Techasintana  et al., 2015; Karginov  et al., 2019; Piecyk et al., 2000)  
  -   Intracellular checkpoint:  PCBP1 safeguards effector identity and restrains Treg‑skewing; T‑cell PCBP1 loss elevates PD‑1/TIGIT/VISTA and blunts antitumor immunity—directly relevant to melanoma immunotherapy. (Ansa‑Addo  et al., 2020)

These RBPs potentially sit at decisive nodes of T‑cell activation -> restraint and exhaustion. Under PD‑1/CTLA‑4 blockade, they are well‑positioned to influence whether CD8+ CTLs remain functional or progress toward Tex.

---

## 2) Rational for RBP (concise)

### ZFP36 family — ZFP36, ZFP36L1, ZFP36L2
- Directly   restrain activation and effector functions; CLIP and ribosome profiling in T cells show broad targeting of cytokines, co‑receptors and signaling components. Zfp36 KO accelerates activation and increases effector cytokines in vivo. (Moore  et al., 2018)  
-   Zfp36/Zfp36l1   set the   timing   and   potency   of CD8 effector differentiation—double perturbation speeds early effector fate. (Petkau  et al., 2022)

### Regnase‑1 & Roquin — ZC3H12A + RC3H1/RC3H2
- Roquin‑1/2   repress    ICOS  and  OX40  mRNAs, limiting Tfh‑skewing and chronic activation. (Glasmacher  et al. , 2010; Vogel  et al., 2013)  
- Regnase‑1 and Roquin   co‑regulate a common 3′UTR stem‑loop   in inflammatory mRNAs with   spatiotemporally distinct mechanisms. (Mino  et al. , 2015)  
- Antigen receptor signaling via   MALT1 cleaves   both proteins, releasing cooperative repression and promoting TH17 differentiation—a paradigm for chronic stimulation. (Jeltsch  et al., 2014)  
-   Disrupting Regnase–Roquin interaction   causes autoimmunity  and   enhances antitumor responses  , linking this axis to tumor immunity. (Behrens  et al. , 2021)

### HNRNPLL
- Master   splicing regulator   of CD45 (PTPRC) in T cells; required for proper activation thresholds and memory development. (Oberdoerffer  et al. , 2008; Wu  et al. , 2006)

### PCBP1
-   Intracellular immune checkpoint   in T cells; prevents Teff→Treg reprogramming, restrains inhibitory receptors (PD‑1/TIGIT/VISTA) and supports antitumor immunity. (Ansa‑Addo  et al. , 2020)

### ELAVL1 (HuR)
-   Stabilizes   activation‑induced transcripts and promotes survival/activation in T cells; strongly upregulated upon stimulation. (Techasintana  et al. , 2015; Karginov  et al. , 2019)

### TIA1 / TIA1B
-   Translational silencing   factor; nucleates stress granules to triage inflammatory mRNAs (e.g.,  TNF), tuning the magnitude/duration of responses. (Piecyk  et al. , 2000)

---

## 3) Relevance to melanoma & checkpoint therapy

- Chronic antigen exposure and inflammatory circuits in melanoma select for Tex‑like programs. Each RBP above controls core Tex drivers—cytokine bursts, co‑stimulatory/checkpoint receptors, TCR signaling, stress‑adaptive translation or splicing—precisely where anti‑PD‑1/anti‑CTLA‑4 act.  
- Engineering studies show that releasing Regnase–Roquin repression enhances antitumor function of T cells, underscoring the therapeutic leverage of this axis in checkpoint settings. (Behrens et al., 2021)

---

## 4) How this maps to the dataset

- In the melanoma single‑cell dataset (Sade‑Feldman et al.), our UMAPs show clear exhaustion gradients and diffuse RBP‑module signal across anti‑PD1 Pre/Post cohorts, consistent with the roles above (see umap_anti‑PD1_Pre_context.pdf, umap_anti‑PD1_Post_context.pdf).  
- Visualizing all 11 RBPs across therapies and Responder vs Non‑responder strata will let us test whether increased ARE‑decay (ZFP36 family) or Regnase–Roquin repression associates with Tex progression and clinical response.

---

## References 

- Ansa‑Addo, E.A., Huang, H., Riesenberg, B.  et al.  (2020) ‘RNA binding protein PCBP1 is an intracellular immune checkpoint for shaping T cell responses in cancer immunity’,  Science Advances , 6(22): eaaz3865. doi:10.1126/sciadv.aaz3865.

- Behrens, G., Hohn, C., de Jonge, L.S.  et al.  (2021) ‘Disrupting Roquin‑1 interaction with Regnase‑1 induces autoimmunity and enhances antitumor responses’,  EMBO Reports , 22(4): e51505. doi:10.15252/embr.202051505.

- Glasmacher, E., Hoefig, K.P., Vogel, K.U.  et al.  (2010) ‘Roquin binds inducible costimulator mRNA and effectors of mRNA decay to induce microRNA‑independent post‑transcriptional repression’,  Nature Immunology , 11(3), 281–289. doi:10.1038/ni.1830.

- Jeltsch, K.M., Hu, D., Brenner, S.  et al.  (2014) ‘Cleavage of roquin and regnase‑1 by the paracaspase MALT1 releases their cooperatively repressed targets to promote TH17 differentiation’,  Nature Immunology , 15(11), 1079–1089. doi:10.1038/ni.3008.

- Karginov, F.V., Musheev, M., Mobley, J.L.  et al.  (2019) ‘HuR controls apoptosis and activation response without compromising RNA homeostasis in T lymphocytes’,  Journal of Immunology , 203(7), 1836–1847. doi:10.4049/jimmunol.1900527.

- Mino, T., Murakawa, Y., Fukao, A.  et al.  (2015) ‘Regnase‑1 and Roquin regulate a common element in inflammatory mRNAs by spatiotemporally distinct mechanisms’,  Cell , 161(5), 1058–1073. doi:10.1016/j.cell.2015.04.029.

- Moore, M.J., Zhang, C., Gantman, E.C.  et al.  (2018) ‘ZFP36 RNA‑binding proteins restrain T cell activation and anti‑viral immunity’,  eLife , 7, e33057. doi:10.7554/eLife.33057.

- Oberdoerffer, S., Moita, L.F., Neems, D.  et al.  (2008) ‘Regulated alternative splicing controls  PTPRC  isoform expression and antigen receptor signaling in T cells’,  Nature Immunology , 9(7), 774–783. doi:10.1038/ni.1623.

- Petkau, G., Salem, R.M., Yates, A.J.  et al.  (2022) ‘The timing of differentiation and potency of CD8 effector T cells is set by Zfp36 and Zfp36l1’,  Nature Communications , 13, 5322. doi:10.1038/s41467-022-33012-0.

- Piecyk, M., Wax, S., Beck, A.R.P.  et al.  (2000) ‘TIA‑1 is a translational silencer that selectively regulates the expression of TNF‑α’,  Journal of Cell Biology , 151(1), 1–12. doi:10.1083/jcb.151.1.1.

- Techasintana, P., Chang, C., Chen, Y.  et al.  (2015) ‘Transcriptomic‑wide discovery of direct and indirect HuR targets in activated CD4 T cells’,  PLOS ONE , 10(5): e0129321. doi:10.1371/journal.pone.0129321.

- Vogel, K.U., Edelmann, S.L., Jeltsch, K.M.  et al.  (2013) ‘Roquin paralogs 1 and 2 redundantly repress the  Icos  and  Ox40  costimulator mRNAs and control follicular helper T cell differentiation’,  Immunity , 38(4), 655–668. doi:10.1016/j.immuni.2012.12.004.

- Wu, Z., Jia, X., de la Cruz, L.  et al.  (2006) ‘The RNA-binding protein Hnrpll regulates CD45 alternative splicing and modulates T cell receptor signaling’,  Journal of Biological Chemistry , 281(6), 3715–3724. doi:10.1074/jbc.M511792200.

- Reviews and mechanistic context: Rehage, N.  et al.  (2018)  Binding of NUFIP2 to Roquin promotes recognition and regulation of ICOS mRNA ,  Nature Communications  9: 299; Schaefer, J.S. (2015)  Roquin – a multifunctional regulator of immune homeostasis ,  International Journal of Molecular Sciences , 16(10), 24366–24384.

