# Topic A — Tree shape (γ, imbalance, size) and the performance of phylogenetic methods

Literature review for the Discussion (MEE manuscript, PCPS/PCC simulation study).
Compiled 2026-10-20.

**Verification protocol.** Every DOI below was resolved on Crossref (`api.crossref.org/works/<DOI>`).
Authors, year, title, journal, volume(issue) and pages are copied from the Crossref record. Pages
for Feng et al. 2017 came from OpenAlex because Crossref lists none. The "Source of claim" line
says what I actually read: abstract (Crossref/OpenAlex/Consensus/Undermind), full-text passage
(Wiley Scholar Gateway), or full PDF read through Undermind `read_pdfs`. Where I only inferred a
link to our results, the note is marked *[our interpretation]*.

**Tools used:** Consensus (5 searches), Undermind (orientation, 2 semantic searches, metadata
lookup, PDF reading of 4 papers), Wiley Scholar Gateway full text (1 search), Scholar Gateway
semantic search (2 searches; the output was too large to display inline, so I parsed the saved
JSON), OpenAlex (abstracts and biblio), Crossref (DOI verification). No tool failed.

**Our results (for reference).** Colless imbalance did not affect type I error. Very negative γ
reduced the power of the GLM on PCPS axis 1 (by about 0.40) and of PCC (by about 0.25), but did
not reduce PERMANOVA power. Low γ improved clade recovery by PCPS. Three empirical anuran regional
pools had γ between −2.5 and −3.9.

---

## 1. Verified references

### A1. Phylogenetic signal tests and tree size/shape

**1. Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K., & Thuiller, W. (2012).** How to measure and test phylogenetic signal. *Methods in Ecology and Evolution*, 3(4), 743–756. https://doi.org/10.1111/j.2041-210x.2012.00196.x
- *Source of claim:* abstract (OpenAlex).
- They simulated BM and other trait models on trees of different size, resolution and branch-length availability. Abouheif's C_mean and Pagel's λ outperformed Moran's I and Blomberg's K under BM. "Sample size influenced not only the uncertainty but also the expected values of most indices." Polytomies and missing branch lengths had "only negligible impacts". The guidelines also take into account index "sensitivities to different topologies of phylogenies".
- *Relevance:* This is the standard precedent for the claim that the expected value and power of phylogeny-based statistics depend on tree properties and not only on the trait process. That fits our finding that tree shape (γ) changes the power of the PCPS/PCC tests. Their size result is also a reason to report the species-pool size used in our design.

**2. Blomberg, S. P., Garland, T., & Ives, A. R. (2003).** Testing for phylogenetic signal in comparative data: behavioral traits are more labile. *Evolution*, 57(4), 717–745. https://doi.org/10.1111/j.0014-3820.2003.tb00285.x
- *Source of claim:* abstract (Consensus).
- The randomization test had correct type I error and "good power (0.8) for simulated datasets with 20 or more species". K "allows valid comparisons of the amount of phylogenetic signal across traits and trees". They introduced the ACDC transformation, in which the rate of evolution accelerates or decelerates. The ACDC parameter was significant for 40% of empirical traits, which suggests that branch-length structure (deep versus tippy) matters for fit.
- *Relevance:* This supports our type I error result. Permutation-based tests keep nominal α whatever the topology, and power is the quantity that depends on the tree. ACDC/early burst is the trait-level analogue of a low-γ tree, because in both most of the variance sits on deep branches *[our interpretation]*.

**3. Freckleton, R. P., Harvey, P. H., & Pagel, M. (2002).** Phylogenetic analysis and comparative data: a test and review of evidence. *The American Naturalist*, 160(6), 712–726. https://doi.org/10.1086/343873
- *Source of claim:* abstract (Consensus).
- λ is a powerful index with low type I error and is robust to incomplete phylogenetic information. In a review of 26 datasets, the cases where signal could be neither supported nor rejected "were found in small phylogenies and probably reflect a lack of power".
- *Relevance:* This is general support for the point that power, not type I error, is the property most sensitive to the size and structure of the tree.

### A2. Phylogenetic eigenvectors (PVR) and their critiques

**4. Diniz-Filho, J. A. F., de Sant'Ana, C. E. R., & Bini, L. M. (1998).** An eigenvector method for estimating phylogenetic inertia. *Evolution*, 52(5), 1247–1262. https://doi.org/10.1111/j.1558-5646.1998.tb02006.x
- *Source of claim:* abstract (Consensus).
- The paper introduces PVR, which is a PCoA of the phylogenetic distance matrix followed by regression on the eigenvectors. In 100 random cladograms, "the number of eigenvalues extracted by the broken-stick model can be fully explained by cladogram topology". Including more eigenvectors inflated type I error.
- *Relevance:* This is the most direct precedent for our mechanism. The eigenstructure of a phylogenetic distance matrix, and therefore what a PCPS axis captures, is determined by tree shape. In a low-γ tree a few deep splits dominate the eigenvalue spectrum. This plausibly concentrates PCPS axis 1 on one deep dichotomy, which would improve clade recovery but could misalign axis 1 with an environmentally sorted clade that is not the deepest split, and so reduce GLM power on axis 1 *[our interpretation; needs checking against our eigenvalue outputs]*.

**5. Rohlf, F. J. (2001).** Comparative methods for the analysis of continuous variables: geometric interpretations. *Evolution*, 55(11), 2143–2160. https://doi.org/10.1111/j.0014-3820.2001.tb00731.x
- *Source of claim:* abstract (OpenAlex).
- PIC is a special case of PGLS. "The importance of the eigenvalues and eigenvectors of the connection matrix, W, for the distribution of ρ is discussed". The phylogenetic autocorrelation method "has several problems that limit its usefulness in comparative studies" and cannot be made equivalent to PGLS.
- *Relevance:* This gives the geometric argument that the behaviour of eigenvector- and autocorrelation-based statistics depends on the spectrum of the phylogenetic matrix, and that spectrum is a property of tree shape. It supports interpreting the γ effect on PCPS/PCC as a consequence of the eigenstructure of matrix P.

**6. Freckleton, R. P., Cooper, N., & Jetz, W. (2011).** Comparative methods as a statistical fix: the dangers of ignoring an evolutionary model. *The American Naturalist*, 178(1), E10–E17. https://doi.org/10.1086/660272
- *Source of claim:* abstract (Consensus).
- In simulations, GLS was more efficient than PVR. It gave lower variance in parameter estimates, lower residual phylogenetic signal and lower type I error. The authors argue that the results "can very likely be generalized to eigenvector methods that control for space and both space and phylogeny".
- *Relevance:* PCPS axes are phylogenetic eigenvectors weighted by community composition, so this critique applies by analogy. It helps explain why a test on a single eigenvector axis (GLM on PCPS1) is more fragile than a test on the whole matrix (PERMANOVA on P). Note that Freckleton et al. focus on type I error, whereas our concern is power.

**7. Diniz-Filho, J. A. F., Bini, L. M., Rangel, T. F., Morales-Castilla, I., Olalla-Tárraga, M. Á., Rodríguez, M. Á., & Hawkins, B. A. (2012).** On the selection of phylogenetic eigenvectors for ecological analyses. *Ecography*, 35(3), 239–249. https://doi.org/10.1111/j.1600-0587.2011.06949.x
- *Source of claim:* abstract (Consensus).
- "Using a relatively small number of eigenvectors may be insufficient to control phylogenetic autocorrelation, leading to flawed conclusions". Selecting eigenvectors by minimizing residual autocorrelation performed best.
- *Relevance:* This supports the argument that analysing only PCPS axis 1 can miss the relevant phylogenetic signal when tree shape spreads it across several axes. It also suggests that our recommendation should favour multivariate tests (PERMANOVA on P) or several axes over a GLM on axis 1 alone.

### A3. Tree shape and community-phylogenetic metrics

**8. Mazel, F., Davies, T. J., Gallien, L., Renaud, J., Groussin, M., Münkemüller, T., & Thuiller, W. (2016).** Influence of tree shape and evolutionary time-scale on phylogenetic diversity metrics. *Ecography*, 39(10), 913–920. https://doi.org/10.1111/ecog.01694 (published online 2015)
- *Source of claim:* full PDF (Undermind `read_pdfs`) and abstract.
- They compared "tippy" and "stemmy" toy trees and pure-birth pools of 20–200 species, and quantified Colless imbalance and γ. MPD "is highly sensitive to deep branching structure", whereas PD and MNTD carry similar, shallower information. Branch-length distribution changed whether a metric registered clustering or overdispersion. Pool tree shape (Colless, γ) "did not influence the discrepancies between ses.PDFaith and ses.MPD".
- *Relevance:* This paper used exactly our two tree-shape descriptors, Colless and γ. Its main lesson is that metrics weight different phylogenetic depths. Matrix P/PCPS, like MPD, is dominated by deep structure, which may explain why low γ (stemmy, with deep splits) helps PCPS recover clades. The null effect of Colless on metric discrepancies parallels our null effect of imbalance on type I error.

**9. Molina-Venegas, R., & Roquet, C. (2014).** Directional biases in phylogenetic structure quantification: a Mediterranean case study. *Ecography*, 37(6), 572–580. https://doi.org/10.1111/j.1600-0587.2013.00442.x
- *Source of claim:* full PDF (Undermind `read_pdfs`) and abstract.
- "Tree shape influenced the accuracy loss derived from the lack of phylogenetic resolution, particularly for NRI and stdPD". In simulations that varied stemminess, the NRI bias changed with γ, while NTI was less affected (Fig. 3). Biases were generally towards type II errors.
- *Relevance:* This is a direct empirical precedent for γ modulating the performance of a basal (MPD-type) community metric. The γ effect fell on type II error (power), as in our results, and not on type I error.

**10. Kembel, S. W. (2009).** Disentangling niche and neutral influences on community assembly: assessing the performance of community phylogenetic structure tests. *Ecology Letters*, 12(9), 949–960. https://doi.org/10.1111/j.1461-0248.2009.01354.x
- *Source of claim:* abstract (Undermind).
- Metacommunity simulations with trait evolution and assembly showed that all tests were sensitive to phylogenetic signal in traits and abundances. Metrics most sensitive to niche processes were also more likely to detect non-random structure under dispersal assembly.
- *Relevance:* This is the closest simulation design to ours (a metacommunity combined with trait evolution). It supports our framing that test performance depends jointly on the assembly process and on how the trait maps onto the tree. Kembel did not vary tree shape, so our γ/Colless factorial extends his design.

**11. Vamosi, S. M., Heard, S. B., Vamosi, J. C., & Webb, C. O. (2009).** Emerging patterns in the comparative analysis of phylogenetic community structure. *Molecular Ecology*, 18(4), 572–592. https://doi.org/10.1111/j.1365-294x.2008.04001.x
- *Source of claim:* abstract (Consensus). Box 1 was seen via Scholar Gateway but not read in full.
- 23 of 39 studies (59%) found phylogenetic clustering. The authors "discuss the relationship between metrics of phylogenetic clustering and tree balance" and call for "multivariate descriptors" of community phylogenetic structure.
- *Relevance:* This is an early statement that tree balance can affect clustering metrics, and it motivates multivariate approaches such as PCPS. **I did not read the tree-balance section itself**, so cite it only for the general point.

### A4. γ statistic and the empirical prevalence of slowdowns

**12. Pybus, O. G., & Harvey, P. H. (2000).** Testing macro-evolutionary models using incomplete molecular phylogenies. *Proceedings of the Royal Society of London. Series B: Biological Sciences*, 267(1459), 2267–2272. https://doi.org/10.1098/rspb.2000.1278
- *Source of claim:* abstract (Consensus).
- The paper develops tests of constant speciation/extinction that are robust to incomplete sampling; this is the origin of the γ statistic and the MCCR test. When applied to published trees, the tests suggested that "in some cases, speciation rates have decreased through time".
- *Relevance:* This is the methods citation for γ.

**13. McPeek, M. A. (2008).** The ecological dynamics of clade diversification and community assembly. *The American Naturalist*, 172(6), E270–E284. https://doi.org/10.1086/593137
- *Source of claim:* abstract (OpenAlex). Summary statistics come from a secondary full-text passage in Machac et al. 2013 (ref. 15).
- In 245 chordate, arthropod, mollusk and plant phylogenies, "a majority of clades … showed rapid lineage accumulation early with a slowing more recently". A metacommunity model in which new species are ecologically differentiated reproduces this pattern. Machac et al. (2013) report the mean ± SD γ of McPeek's dataset as −0.737 ± 3.102 (I read this in Machac et al.'s full text, not in McPeek).
- *Relevance:* This establishes that negative γ is the empirical norm, so our low-γ scenarios are realistic and not extreme. It also links slowdowns to metacommunity processes, which is conceptually close to our framework. Our anuran pools (γ = −2.5 to −3.9) fall in the strongly negative tail of this distribution *[our comparison]*.

**14. Phillimore, A. B., & Price, T. D. (2008).** Density-dependent cladogenesis in birds. *PLoS Biology*, 6(3), e71. https://doi.org/10.1371/journal.pbio.0060071
- *Source of claim:* full PDF (Undermind `read_pdfs`) and abstract.
- "A third of all clades had a γ ≤ −1.645". "Thirteen of the 23 clades (57%) that include more than 20 species show significant slowdowns". Clade size correlated with γ (r = −0.58). Larger clades are expected to show stronger slowdowns even under constant rates, but the observed frequency exceeded that expectation.
- *Relevance:* This shows that significantly negative γ is common in species-rich clades. It also matters for our design because γ scales with tree size, so pool size and γ should not be confounded in our simulations or in our interpretation of the empirical pools.

**15. Machac, A., Storch, D., & Wiens, J. J. (2013).** Ecological causes of decelerating diversification in carnivoran mammals. *Evolution*, 67(8), 2423–2433. https://doi.org/10.1111/evo.12126
- *Source of claim:* full-text passages (Wiley Scholar Gateway) and abstract.
- Slowdowns were strong in carnivorans (γ = −0.935 ± 0.851) and comparable to McPeek's (2008) dataset. The effect of clade diversity on γ "may arise as a statistical artifact (i.e., greater clade diversity increases the ability of the gamma statistic to refute constant diversification)".
- *Relevance:* This is a secondary source for McPeek's γ distribution, and it adds a further caveat about γ and clade size.

**16. Cusimano, N., & Renner, S. S. (2010).** Slowdowns in diversification rates from real phylogenies may not be real. *Systematic Biology*, 59(4), 458–464. https://doi.org/10.1093/sysbio/syq032
- *Source of claim:* abstract (Consensus).
- "Oversampling deep nodes … strongly biases inferences toward downturns", which happens when sampling drops below about 80%. The authors suggest not inferring downturns from trees with less than 80% species sampling.
- *Relevance:* This is an important caveat for our empirical γ values (−2.5 to −3.9). Regional pools pruned from a large phylogeny (e.g., a global amphibian tree) are not incomplete clades in Cusimano and Renner's sense, but sampling and pruning still shape γ. We should state that our γ values describe the shape of the pool tree as used in the analysis, without implying a diversification process *[our interpretation]*.

**17. Fordyce, J. A. (2010).** Interpreting the γ statistic in phylogenetic diversification rate studies: a rate decrease does not necessarily indicate an early burst. *PLoS ONE*, 5(7), e11781. https://doi.org/10.1371/journal.pone.0011781
- *Source of claim:* abstract (Consensus).
- "The gamma statistic is extraordinarily sensitive to recent diversification rates, and does not necessarily detect early bursts". Its power was greater for recent rate decreases than for early bursts.
- *Relevance:* Negative γ mainly reflects few recent splits, that is, long terminal branches and short-to-moderate internal structure near the tips. For PCPS this means that species in a pool with very negative γ are phylogenetically nearly equidistant within clades and well separated between them. This plausibly sharpens the clade signal recovered by PCPS. At the same time it reduces the effective number of independent phylogenetic units available to a univariate GLM on axis 1 *[our interpretation]*.

**18. Liow, L. H., Quental, T. B., & Marshall, C. R. (2010).** When can decreasing diversification rates be detected with molecular phylogenies and the fossil record? *Systematic Biology*, 59(6), 646–659. https://doi.org/10.1093/sysbio/syq052
- *Source of claim:* abstract (Consensus).
- "A significantly negative γ value for a molecular phylogeny indicates a decreasing rate of diversification", but non-significant or positive γ can arise under several processes.
- *Relevance:* This is optional support for treating γ as a descriptor of tree shape rather than of process.

### A5. Anuran diversification

**19. Roelants, K., Gower, D. J., Wilkinson, M., Loader, S. P., Biju, S. D., Guillaume, K., Moriau, L., & Bossuyt, F. (2007).** Global patterns of diversification in the history of modern amphibians. *Proceedings of the National Academy of Sciences*, 104(3), 887–892. https://doi.org/10.1073/pnas.0608378104
- *Source of claim:* abstract (Consensus).
- Amphibian diversification shows "several episodes of accelerated" diversification rather than gradual accumulation. "Approximately 86% of modern frog species … descended from only five ancestral lineages that produced major radiations in the late Cretaceous and early Tertiary".
- *Relevance:* Concentrated, rapid radiations at moderate depth followed by many species are a scenario that tends to produce clade-structured trees with negative γ *[our interpretation; Roelants et al. do not report γ]*. This makes our empirical γ values plausible for anurans.

**20. Feng, Y.-J., Blackburn, D. C., Liang, D., Hillis, D. M., Wake, D. B., Cannatella, D. C., & Zhang, P. (2017).** Phylogenomics reveals rapid, simultaneous diversification of three major clades of Gondwanan frogs at the Cretaceous–Paleogene boundary. *Proceedings of the National Academy of Sciences*, 114(29), E5864–E5870. https://doi.org/10.1073/pnas.1704632114
- *Source of claim:* abstract (Consensus).
- Three clades (Hyloidea, Microhylidae, Natatanura), which together hold about 88% of extant anurans, "simultaneously underwent rapid diversification" at the K–Pg boundary.
- *Relevance:* Simultaneous early bursts in the major anuran lineages provide a phylogenetic basis for strongly negative γ in Neotropical anuran pools dominated by Hyloidea *[our interpretation; no γ reported in the abstract]*.

**21. Jiménez-Ortega, D., Valente, L., Dugo-Cota, Á., Rabosky, D. L., Vilà, C., & Gonzalez-Voyer, A. (2023).** Diversification dynamics in Caribbean rain frogs (*Eleutherodactylus*) are uncoupled from the anuran community and consistent with adaptive radiation. *Proceedings of the Royal Society B: Biological Sciences*, 290(1990), 20222171. https://doi.org/10.1098/rspb.2022.2171
- *Source of claim:* abstract (Consensus).
- In a 197-species Caribbean anuran assemblage, diversification followed "rapid diversification early in the radiation which gradually decreases towards the present", and equilibrium (diversity-dependent) dynamics influenced richness of all Caribbean anurans.
- *Relevance:* This is a direct example of a diversification slowdown in a *regional anuran assemblage*, the same unit as our species pools. It is the best anuran-specific support we found. γ itself is not quoted in the abstract.

### A6. Empirical tree imbalance

**22. Blum, M. G. B., & François, O. (2006).** Which random processes describe the tree of life? A large-scale study of phylogenetic tree imbalance. *Systematic Biology*, 55(4), 685–691. https://doi.org/10.1080/10635150600889625
- *Source of claim:* abstract/first paragraph (OpenAlex). The full text was not read.
- On a large database of published phylogenies, "the Yule model is inadequate to describe phylogenetic tree data", and many trees are "consistent with a simple branch split model first considered by Aldous in 1996".
- *Relevance:* Real trees are more imbalanced than Yule trees. Our imbalance gradient should therefore be described relative to Yule, and the empirical pools placed on it. Aldous's model corresponds to the β-splitting model with β = −1, which is more imbalanced than Yule (β = 0); this is my background knowledge and is not stated in the text I read, so check it in the full paper before citing the β value.

**23. Heard, S. B. (1992).** Patterns in tree balance among cladistic, phenetic, and randomly generated phylogenetic trees. *Evolution*, 46(6), 1818–1826. https://doi.org/10.1111/j.1558-5646.1992.tb01171.x
- *Source of claim:* abstract (Consensus).
- Across 208 literature trees, "imbalance decreases with tree size", and the degree of imbalance "is inconsistent with the uniform-rate speciation model".
- *Relevance:* This supports Blum & François (2006). Because imbalance depends on tree size, the Colless index should be normalized when tree sizes are compared.

---

## 2. Candidate arguments for the Discussion

1. **Type I error is robust to topology; power is not.** Permutation and likelihood tests of phylogenetic structure keep nominal α across topologies, and their power depends on tree size and structure (Blomberg et al. 2003; Freckleton et al. 2002; Münkemüller et al. 2012). Molina-Venegas & Roquet (2014) found that γ-dependent bias fell on type II errors. Our null effect of Colless imbalance on type I error, together with the γ effect on power, fits this general pattern.

2. **The γ effect works through the eigenstructure of the phylogenetic matrix.** The number and importance of phylogenetic eigenvectors is set by tree topology (Diniz-Filho et al. 1998), and the behaviour of autocorrelation/eigenvector statistics depends on the eigenvalues of the phylogenetic matrix (Rohlf 2001). In very negative-γ trees, a few deep splits dominate the spectrum, so PCPS axis 1 is captured by the deepest dichotomy. This (i) improves clade recovery, but (ii) reduces GLM power on axis 1 whenever the environmentally sorted clade is not the one that dominates axis 1. PERMANOVA uses the whole P matrix and is not affected. *[Mechanism is our interpretation; check it against the eigenvalue and axis-loading output.]*

3. **Single-axis tests on phylogenetic eigenvectors are fragile.** Critiques of PVR (Freckleton et al. 2011) and work on eigenvector selection (Diniz-Filho et al. 2012) show that a small number of eigenvectors can misrepresent phylogenetic structure. This supports recommending multivariate tests on P (PERMANOVA), or several PCPS axes, over a GLM on PCPS1 alone, especially for clades with strongly negative γ.

4. **Different methods weight different phylogenetic depths.** Mazel et al. (2016) showed that MPD-type metrics respond to deep branching while PD and MNTD track shallow structure, and that γ and Colless describe this tippy-versus-stemmy contrast. Matrix P and PCPS behave like deep-weighted metrics. This explains why stemmy (low-γ) trees favour PCPS clade recovery, and it predicts that methods weighting shallow structure (e.g., phylofactorization on tips) would respond differently. *[The prediction is our interpretation.]*

5. **Low γ is the empirical norm, especially in species-rich clades, and anurans are no exception.** Most of the 245 clades in McPeek (2008) showed slowdowns. A third of 45 bird clades, and 57% of large ones, had significantly negative γ (Phillimore & Price 2008). Anuran diversification was concentrated in a few rapid radiations (Roelants et al. 2007; Feng et al. 2017), and a regional anuran assemblage showed early-burst, diversity-dependent dynamics (Jiménez-Ortega et al. 2023). Our three anuran pools (γ = −2.5 to −3.9) are therefore representative, and the power loss we document for GLM/PCC under low γ is likely to matter in real amphibian datasets.

6. **Caveats on interpreting empirical γ.** γ is sensitive to sampling (Cusimano & Renner 2010), to recent branching (Fordyce 2010) and to clade size (Phillimore & Price 2008; Machac et al. 2013). We use γ as a descriptor of pool-tree shape, not as evidence for a diversification process. For the anuran pools, we should report sampling completeness relative to the source phylogeny.

7. **Imbalance: empirical trees are more imbalanced than Yule.** Published trees are more imbalanced than Yule expectations (Heard 1992; Blum & François 2006). Our range of Colless values should be placed relative to Yule and to the empirical pools. Because the methods were insensitive to imbalance across that range, the conclusions should hold for typical empirical imbalance.

## 3. Gaps / not verified

- I found **no published γ values for anuran regional assemblages** beyond the qualitative early-burst result of Jiménez-Ortega et al. (2023). The full text of that paper may report γ, but I did not read it.
- I did not read the tree-balance discussion in Vamosi et al. (2009) or the full text of Blum & François (2006).
- I found no study that tested PCPS, PCC or evoPCA performance across γ or imbalance gradients. This appears to be a gap that our study fills (based on these searches, not an exhaustive review).
- Billaud et al. (bioRxiv 2019, DOI 10.1101/586420) reports diversity decline in Archaeobatrachia but no γ. It was not included, and the published version was not verified.
- LeapSpace (Scopus) report 1 (2026-10-20, query on tree shape × phylogenetic signal): only Brock et al. 2011 (Syst. Biol. 60: 410–419, doi:10.1093/sysbio/syr007; abstract read on OpenAlex) added to the Discussion.
