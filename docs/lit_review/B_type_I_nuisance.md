# Topic B — Inflated type I error and nuisance structure in community-level tests relating traits/phylogeny to environment

Compiled 2026-10-20 for the Discussion of the PCPS / PCC simulation manuscript (MEE).

**Verification protocol.** Every DOI below was resolved on Crossref (`api.crossref.org/works/<DOI>`); authors, year, title, journal, volume(issue) and pages are copied from the Crossref record. "Source of claim" says what I actually read: **FT** = full text read from the local PDF in `references/`; **ABS** = abstract read (via Consensus, Wiley Scholar Gateway, Undermind or OpenAlex `abstract_inverted_index`). Nothing is attributed to a paper beyond what that source says. Statements marked *(interpretation)* are my inference, not the authors' claims.

**Search tools used.** Consensus (3 queries), Scholar Gateway semanticSearch (1 query; output mostly duplicated Duarte et al. 2018), Wiley Scholar Gateway full text (1 query; surfaced ter Braak et al. 2018), Undermind `search_papers` (1 query, after `get_orientation`), OpenAlex (3 keyword searches + abstract retrieval). None failed. Undermind also returned a bioRxiv preprint (Zelený, Helsen & Lee 2021, intraspecific CWM) and a 2026 MEE paper (Ortega et al., mean–variance scaling in ITV); I did not verify or read them, so they are not listed.

---

## 1. Verified references

### Core PCPS / PFW / PCC framework

**1. Pillar, V. D., & Duarte, L. d. S. (2010).** A framework for metacommunity analysis of phylogenetic structure. *Ecology Letters*, 13(5), 587–596. https://doi.org/10.1111/j.1461-0248.2010.01456.x
- Source: metadata only (Crossref; found via Undermind). Not read here.
- Relevance: original definition of phylogenetic fuzzy weighting and matrix P; cite for the method, not for type I error claims.

**2. Duarte, L. D. S., Debastiani, V. J., Freitas, A. V. L., & Pillar, V. D. (2016).** Dissecting phylogenetic fuzzy weighting: theory and application in metacommunity phylogenetics. *Methods in Ecology and Evolution*, 7(8), 937–946. https://doi.org/10.1111/2041-210x.12547
- Source: FT + Appendix S2 (local .docx).
- Claims: when neither phylogeny nor environment structured species distributions, site shuffle and taxa shuffle had correct type I error for Mantel, ADONIS and GLM on PCPS axes, but the largest-P criterion *underestimated* type I error (Table 1–2; e.g. Mantel largest-P rejection 0.016–0.023, GLM PCPS1 0.018–0.026). When only phylogeny structured distributions, taxa shuffle was conservative. GLM had acceptable power only for PCPS1, which they interpret as PCPS detecting "phylogenetic patterns related to more basal nodes". In Appendix S2 species tolerances are drawn as "random normally distributed values with standard deviation 10 and mean tolerance" (μtol = 5 or 10), i.e. niche breadth was independent of phylogeny; maximum abundance (h) was drawn uniform 0–30, also independent of phylogeny.
- Relevance: direct baseline. Our σ-signal scenarios are something their design did not include; our conservative behaviour of PERMANOVA/Mantel/GLM PCPS1 is consistent with their largest-P results.

**3. Peres-Neto, P. R., & Kembel, S. W. (2015).** Phylogenetic gradient analysis: environmental drivers of phylogenetic variation across ecological communities. *Plant Ecology*, 216(5), 709–724. https://doi.org/10.1007/s11258-014-0405-0
- Source: FT.
- Claims: row- (site) and column- (species) permutation tests each "produce elevated type I errors" for the fourth-corner and the PCC framework when used alone (citing Dray & Legendre 2008; Peres-Neto et al. 2012); "the only tests with a valid type I error rates are those based on the criteria of choosing the greatest P value between the row- and column-based permutation tests". When only phylogeny structures distributions, the column-based test is invalid and the row-based test is valid.
- Relevance: defines the PCC and the max test as used in our study. Their simulations did not vary niche breadth structure on the tree *(my reading of their simulation description; check their Appendix A before stating this in the manuscript)*.

**4. Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012).** Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. *Ecology*, 93(sp8). https://doi.org/10.1890/11-0494.1  (Crossref returned no page range; check before citing.)
- Source: ABS (OpenAlex).
- Claims: proposes an analytical framework linking niche properties based on traits and phylogeny to environmental and spatial variation; simulations "show that the procedure is robust under a variety of scenarios".
- Relevance: origin of the community simulation model (Gaussian response, tolerance, maximum abundance) that Duarte et al. (2016, 2018) and we build on; cite for the simulation lineage.

**5. Pavoine, S. (2016).** A guide through a family of phylogenetic dissimilarity measures among sites. *Oikos*, 125(12), 1719–1732. https://doi.org/10.1111/oik.03262
- Source: FT (abstract and Methods).
- Claims: compares PD-dissimilarity indices with permutation tests and null models along simulated gradients; results for a bat dataset depend on whether abundance is considered and on the index used. Type I error scenarios used pure-birth trees and species abundances independent of phylogeny and environment.
- Relevance: supports the point that abundance weighting changes the outcome of phylobetadiversity tests (our smaller inflation with abundance data), and that earlier type I evaluations used null scenarios without phylogenetically structured nuisance parameters.

### CWM / fourth-corner type I error debate

**6. Dray, S., & Legendre, P. (2008).** Testing the species traits–environment relationships: the fourth-corner problem revisited. *Ecology*, 89(12), 3400–3412. https://doi.org/10.1890/08-0349.1
- Source: ABS (OpenAlex).
- Claims: evaluates existing permutation models for the fourth-corner by simulation and proposes a two-step testing procedure.
- Relevance: origin of the combined row/column testing logic.

**7. ter Braak, C. J. F., Cormont, A., & Dray, S. (2012).** Improved testing of species traits–environment relationships in the fourth-corner problem. *Ecology*, 93(7), 1525–1526. https://doi.org/10.1890/12-0126.1
- Source: ABS (OpenAlex).
- Claims: none of the five Dray & Legendre (2008) permutation methods attained correct type I error in all cases; reporting the maximum of the individual P values (max test) "controls the type I error in all cases"; applicable to RLQ.
- Relevance: formal justification of the max test we use for PCPS/PCC.

**8. Peres-Neto, P. R., Dray, S., & ter Braak, C. J. F. (2017).** Linking trait variation to the environment: critical issues with community-weighted mean correlation resolved by the fourth-corner approach. *Ecography*, 40(7), 806–816. https://doi.org/10.1111/ecog.02302
- Source: FT.
- Claims: CWM and SNC correlations have "inconsistent parameter estimation and unacceptable significance rates when only the environment or only traits are structuring species distributions"; parametric and single permutation tests of CWM/SNC give inflated type I error; the row–column (max) test is unbiased and more powerful. In their simulations niche breadth for each species is `s` × U(0,1) and maximum abundance U(0.3, 1), both independent of traits and of phylogeny; `s` was varied (1, 5, 10) only to change correlation strength.
- Relevance: the central reference for "one-sided structure produces spurious significance". Our result adds that structure in a nuisance species parameter (σ) can do the same.

**9. Zelený, D. (2018).** Which results of the standard test for community-weighted mean approach are too optimistic? *Journal of Vegetation Science*, 29(6), 953–966. https://doi.org/10.1111/jvs.12688
- Source: FT.
- Claims: three hypothesis categories (A: assumes trait–composition link; B: assumes composition–environment link; C: assumes neither). The standard (row-based) test is inflated for B and C; Category C "requires the use of the max test as the only way to control for correct Type I error rate". Inflation of the standard test for C is negatively related to compositional β-diversity and positively related to the strength of the composition–environment link and to sample size.
- Relevance: (i) provides the vocabulary to state which null hypothesis each of our tests addresses; (ii) the positive relation between inflation and sample size parallels our largest inflation at N = 200 *(interpretation: Zelený's "sample size" is number of plots, not species pool; check our definition of N before drawing the parallel)*.

**10. Zelený, D., & Schaffers, A. P. (2012).** Too good to be true: pitfalls of using mean Ellenberg indicator values in vegetation analyses. *Journal of Vegetation Science*, 23(3), 419–431. https://doi.org/10.1111/j.1654-1103.2011.01366.x
- Source: ABS (OpenAlex).
- Claims: community means inherit compositional similarity, producing too-high explained variance and too many significant results; proposes the modified (column-based) permutation test.
- Relevance: earliest demonstration that weighted means carry compositional structure; analogous to matrix P inheriting co-occurrence structure.

**11. Hawkins, B. A., Leroy, B., Rodríguez, M. Á., Singer, A., Vilela, B., Villalobos, F., Wang, X., & Zelený, D. (2017).** Structural bias in aggregated species-level variables driven by repeated species co-occurrences: a pervasive problem in community and assemblage data. *Journal of Biogeography*, 44(6), 1199–1211. https://doi.org/10.1111/jbi.12953
- Source: ABS (OpenAlex).
- Claims: averages of random species attributes ("random intrinsic variables") explain moderate-to-large variance in richness; explained variance increases with species co-occurrence; assemblage-aggregated variables "cannot be taken at face value".
- Relevance: general mechanism behind spurious community-level signal; PCPS scores are aggregated species-level variables.

**12. ter Braak, C. J. F., Peres-Neto, P. R., & Dray, S. (2017).** A critical issue in model-based inference for studying trait-based community assembly and a solution. *PeerJ*, 5, e2885. https://doi.org/10.7717/peerj.2885
- Source: ABS (Consensus).
- Claims: with additional random variation in species responses to the environment, site-based resampling tests with GLM have type I error up to 90% at nominal 5%; it is "a problem of ignoring a random effect", not omitted-variable bias; the GLM-based pmax test controlled type I error but was slightly inflated with both missing traits and missing environmental variables.
- Relevance: closest conceptual analogue to our σ result. Species-specific variation in response (here, niche breadth) that is not part of the tested hypothesis inflates site-based tests.

**13. ter Braak, C. J. F., Peres-Neto, P. R., & Dray, S. (2018).** Simple parametric tests for trait–environment association. *Journal of Vegetation Science*, 29(5), 801–811. https://doi.org/10.1111/jvs.12666
- Source: ABS (Wiley Scholar Gateway full-text passage).
- Claims: the ZS-modified test and the fourth-corner max test "are conservative in different scenarios"; a parametric max test combining CWM and SNC regressions controls type I error with equal or higher power; a weighted version was inflated.
- Note: the Wiley citation line lists O. Wildi as fourth author; Crossref and Undermind list three authors. I used Crossref. Check the published PDF.
- Relevance: documented conservativeness of max tests supports our finding that PERMANOVA, Mantel and GLM PCPS1 become more conservative.

**14. ter Braak, C. J. F. (2019).** New robust weighted averaging- and model-based methods for assessing trait–environment relationships. *Methods in Ecology and Evolution*, 10(11), 1962–1971. https://doi.org/10.1111/2041-210x.13278
- Source: ABS (Consensus).
- Claims: because species vary greatly in frequency and abundance, regressions weighted by Hill's N2 combined in the max test best controlled type I error; the existing multilevel model showed bias and type I inflation, while MLM3 (with site-related random effects) did not; extensions to phylogeny are possible.
- Relevance: heterogeneity among species (frequency/abundance) is a recognised source of invalid tests; points to weighting and model-based remedies.

**15. Miller, J. E. D., Damschen, E. I., & Ives, A. R. (2019).** Functional traits and community composition: a comparison among community-weighted means, weighted correlations, and multilevel models. *Methods in Ecology and Evolution*, 10(3), 415–425. https://doi.org/10.1111/2041-210x.13119
- Source: ABS (Consensus).
- Claims: CWM regression had inflated type I error in some scenarios; weighted correlations with permutation had good type I control but low power; when species abundances are correlated among sites, a conditional parametric LRT had correct type I error but lower power.
- Relevance: power–type I trade-off; independent confirmation that CWM-type aggregation is unsafe.

**16. Lepš, J., & de Bello, F. (2023).** Differences in trait–environment relationships: implications for community weighted means tests. *Journal of Ecology*, 111(11), 2328–2341. https://doi.org/10.1111/1365-2745.14172
- Source: ABS (Consensus).
- Claims: CWM tests and species-level tests address different null hypotheses; the type I error of classic CWM tests "should not be a priori considered inflated" when the question is about community-level trait structure.
- Relevance: counterpoint to cite when defining what "type I error" means for GLM on PCPS axes (community-level vs. clade-level hypothesis).

### Phylogenetic non-independence, dispersion and null-model mis-specification

**17. Duarte, L. D. S., Debastiani, V. J., Carlucci, M. B., & Diniz-Filho, J. A. F. (2018).** Analyzing community-weighted trait means across environmental gradients: should phylogeny stay or should it go? *Ecology*, 99(2), 385–398. https://doi.org/10.1002/ecy.2081
- Source: FT.
- Claims: "even under moderate phylogenetic signal in both the trait used to define CWM and species distribution across communities, OLS models relating CWM variation to environmental gradients lead to inflated type I error"; phylogeny-corrected OLS controls type I error but may lose power. Simulation tolerances again drawn as random normal (SD 10, mean μtol = 10), i.e. independent of phylogeny.
- Relevance: phylogenetic signal as a source of false positives in PFW-related analyses; our study moves the signal from the trait to a niche parameter (σ).

**18. Li, D., & Ives, A. R. (2017).** The statistical need to include phylogeny in trait-based analyses of community composition. *Methods in Ecology and Evolution*, 8(10), 1192–1199. https://doi.org/10.1111/2041-210x.12767
- Source: ABS (Consensus).
- Claims: LMMs had "unacceptably high type I error rates when there was phylogenetic signal in both traits and the residual variation in species distributions" (e.g. from an unmeasured trait with phylogenetic signal); PLMMs fixed this.
- Relevance: phylogenetically structured σ acts like a latent trait with signal in residual species distributions *(interpretation)*; PLMM is a candidate remedy.

**19. Braga, J., ter Braak, C. J. F., Thuiller, W., & Dray, S. (2018).** Integrating spatial and phylogenetic information in the fourth-corner analysis to test trait–environment relationships. *Ecology*, 99(12), 2667–2674. https://doi.org/10.1002/ecy.2530
- Source: ABS (Consensus).
- Claims: spatial and phylogenetic autocorrelation can inflate type I error of the fourth-corner double-permutation test; Moran spectral randomization preserving phylogenetic and spatial structure "provides total control over type I error rates".
- Relevance: directly suggests a fix for our taxa shuffle: constrained (phylogeny-preserving) species randomization.

**20. Hardy, O. J. (2008).** Testing the spatial phylogenetic structure of local communities: statistical performances of different null models and test statistics on a locally neutral community. *Journal of Ecology*, 96(5), 914–926. https://doi.org/10.1111/j.1365-2745.2008.01421.x
- Source: ABS (Consensus).
- Claims: tests become liberal "when the randomization algorithm breaks down a structure in the original data set unrelated to the null hypothesis"; e.g. when species abundances are non-randomly distributed on the phylogeny, null models that preserve this structure (permuting species with similar abundances among tips) perform better.
- Relevance: the clearest general statement of our mechanism. Taxa shuffle breaks the phylogenetic structure of σ, which is unrelated to the environment–phylogeny H0.

**21. Tsang, T. P. N., & Cadotte, M. W. (2025).** Species overlap and phylogenetic relatedness result in community statistical non-independence (and what to do about it). *Ecology Letters*, 28(12), e70267. https://doi.org/10.1111/ele.70267
- Source: FT.
- Claims: shared species and shared evolutionary history among communities produce "compositional autocorrelation"; naïve linear regression has higher type I error and underestimates uncertainty, while mixed models with a community-level covariance matrix (from phylogenetic covariance upscaled by abundance) have lower type I error and similar or higher power.
- Relevance: recent and general framing; supports the argument that site-level PCPS scores are not independent observations and points to a model-based alternative to permutation.

**22. Anderson, M. J., & Walsh, D. C. I. (2013).** PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: what null hypothesis are you testing? *Ecological Monographs*, 83(4), 557–574. https://doi.org/10.1890/12-2010.1
- Source: ABS (Consensus).
- Claims: the Mantel test (and ANOSIM) are very sensitive to heterogeneous dispersions; PERMANOVA is largely unaffected in balanced designs; in unbalanced designs all tests can be too liberal or too conservative depending on which group has greater dispersion.
- Relevance: heterogeneous niche breadth produces heterogeneous dispersion in composition. Note that their designs are group-based, whereas ours use a continuous gradient, so the transfer is analogical *(interpretation)*.

**23. Warton, D. I., Wright, S. T., & Wang, Y. (2012).** Distance-based multivariate analyses confound location and dispersion effects. *Methods in Ecology and Evolution*, 3(1), 89–101. https://doi.org/10.1111/j.2041-210x.2011.00127.x
- Source: ABS (OpenAlex).
- Claims: distance-based analyses make implicit mean–variance assumptions out of step with real count data, which confounds location and dispersion effects; transformation does not solve this; multivariate GLMs do.
- Relevance: relevant to the presence/absence vs. abundance contrast and to why distance-based tests (PERMANOVA, Mantel on DP) respond differently from axis-based tests.

---

## 2. Candidate Discussion arguments

**A1. Earlier evaluations left nuisance species parameters unstructured.** Duarte et al. (2016, 2018; FT + Appendix S2), Peres-Neto et al. (2017; FT) and Pavoine (2016; FT) drew niche breadth and maximum abundance independently of phylogeny, so their type I error estimates hold only for a null in which these parameters carry no phylogenetic signal. Our design relaxes that assumption. Suggested phrasing: "Previous evaluations of PFW and fourth-corner–type tests generated species tolerances independently of phylogeny (Duarte et al., 2016, 2018; Peres-Neto et al., 2017), so their type I error estimates apply to a null in which nuisance niche parameters are phylogenetically unstructured."

**A2. Mechanism: taxa shuffle breaks structure unrelated to H0.** When σ is conserved, clades differ in how widely they spread along the gradient, which is a phylogenetic pattern in composition that is not an environment–phylogeny association in optima. Taxa shuffle randomizes tips and destroys this pattern, so the observed statistic is compared with a null distribution that is too narrow, making the column test liberal. The max test then fails whenever the site test also rejects. Support: Hardy (2008; ABS) gives the general principle; Braga et al. (2018; ABS) show phylogenetic autocorrelation inflating the double-permutation test; ter Braak et al. (2017 PeerJ; ABS) show that ignored species-level random variation in environmental responses inflates tests. The mechanism itself is *(interpretation)* and should be checked against our own diagnostics, for example the null distributions of the statistic under taxa shuffle with σ-signal vs. without.

**A3. Why PCPS2 and PCC but not PCPS1, PERMANOVA or Mantel.** Duarte et al. (2016; FT) show PCPS1 summarises basal-node gradients and carries most of the environmental power. A σ-driven pattern, in which clades with broad vs. narrow niches differ in occupancy breadth rather than position, may load on a subsequent orthogonal axis (PCPS2). Omnibus distance-based tests on DP average over all axes, which dilutes it. *(Interpretation. It needs a check: correlate the PCPS2 loadings with clade-level mean σ in the simulations.)* PCC weights species by their phylogenetic covariance (Peres-Neto & Kembel 2015; FT), which may make it similarly sensitive. Anderson & Walsh (2013; ABS) and Warton et al. (2012; ABS) explain why distance-based statistics respond differently to dispersion-type structure than to location-type structure.

**A4. Conservativeness of the max test is expected.** Duarte et al. (2016; FT) already found the largest-P criterion below nominal α under the full null. ter Braak et al. (2018; ABS) report that fourth-corner max tests "are conservative in different scenarios". Our more conservative PERMANOVA/Mantel/PCPS1 under σ-signal is in line with this. The practical cost is lost power, not false positives.

**A5. Sample size and inflation.** Zelený (2018; FT) shows that inflation of the standard CWM test increases with sample size and with the strength of the composition–environment link. Our peak inflation at N = 200 (largest species pool) fits the broader point that more data give more power to detect nuisance structure as well as real signal. *(Caveat: Zelený's sample size is the number of plots; ours is pool size. Present this as an analogy only.)*

**A6. Which null hypothesis is tested.** Using Zelený's (2018; FT) Categories A–C and Lepš & de Bello's (2023; ABS) distinction between community-level and species-level hypotheses, we can state that GLM on a PCPS axis with the max test targets a Category-C-like hypothesis (no assumed link of either phylogeny or environment to composition). Under that hypothesis σ-signal does create a real phylogeny–composition link, so the column test rejects correctly, and only the max test's site-side component protects against false conclusions. This makes the inflation a failure of the joint null, not of either component alone *(interpretation; consistent with Peres-Neto & Kembel 2015, FT, on what each permutation tests)*.

**A7. Abundance data reduce inflation.** Pavoine (2016; FT) shows that conclusions from phylobetadiversity indices depend on whether abundance is used. ter Braak (2019; ABS) shows that species heterogeneity in frequency and abundance affects test validity and that N2 weighting helps. A plausible reading: with abundances, the contribution of tail occurrences, which is where σ acts, is down-weighted relative to occurrences near the optimum *(interpretation)*. This could also explain why signal in maximum abundance had no effect, since h scales abundance without changing where a species occurs. *(Interpretation; I found no paper testing phylogenetic signal in maximum abundance.)*

**A8. Aggregated community variables inherit co-occurrence structure.** Matrix P and its PCPS scores are aggregates of species-level phylogenetic information. Zelený & Schaffers (2012; ABS), Hawkins et al. (2017; ABS) and Tsang & Cadotte (2025; FT) show that such aggregates carry compositional autocorrelation that standard tests ignore. This places our result in a general class of problems and is not a peculiarity of PFW.

**A9. Remedies to recommend or discuss.** (i) Phylogeny-constrained species randomization, e.g. Moran spectral randomization (Braga et al. 2018; ABS) or permuting species within classes of similar niche breadth, following Hardy's (2008; ABS) abundance-class swaps. (ii) Model-based tests with phylogenetic random effects: PLMM (Li & Ives 2017; ABS), MLM3 (ter Braak 2019; ABS), mixed models with compositional covariance (Tsang & Cadotte 2025; FT). (iii) Phylogenetic correction of the response (Duarte et al. 2018; FT). (iv) Report PCPS axes beyond the first only with caution, and check σ-signal empirically when niche-breadth estimates are available.

**A10. Tree imbalance had no effect on type I error.** I found no paper in this topic that tests how tree shape (Colless, γ) affects type I error of these community-level tests. Pavoine (2016; FT) used only pure-birth trees, and Duarte et al. (2016) simulated a different tree per replicate without manipulating shape. This can be stated as a gap: "[CITAR: study evaluating tree shape effects on community phylogenetic tests — not found in Topic B searches; see Topic A]".

---

## 3. Items NOT verified / to check

- Peres-Neto et al. (2012) page range (Crossref returned none).
- ter Braak et al. (2018) author list (Wiley citation line adds Wildi; Crossref does not).
- Zelený, Helsen & Lee (2021, bioRxiv) and Ortega et al. (2026, MEE): surfaced by Undermind; neither verified nor read.
- Whether Peres-Neto & Kembel (2015) Appendix A varied niche breadth relative to phylogeny: the Appendix (.doc) is in `references/appendices/` but I did not read it.
