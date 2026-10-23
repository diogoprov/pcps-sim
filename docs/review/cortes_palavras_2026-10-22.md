# Cortes para voltar aos limites de palavras (proposta, 2026-10-22)

Nada foi aplicado. Os cortes estão num script (no scratchpad da sessão) que roda
contra o `methods_results.qmd` atual e confere que cada trecho aparece uma única vez.

## Contagem

Regra: texto do Abstract até o fim das Referências, sem YAML, comentários, linhas
de tabela e legendas de figura, e sem a AI disclosure statement. Na versão de
7.996, essa regra dava 7.969, ou seja, ela conta cerca de 27 palavras a menos que
a contagem anterior. Os números abaixo já estão corrigidos para a escala antiga.

| | Atual | Com os cortes | Limite |
|---|---:|---:|---:|
| Texto principal (com refs, sem legendas/tabela) | ≈ 8.454 | **≈ 7.965** | 8.000 |
| Abstract (itens 1–4) | 367 | **346** | 350 |

O excesso não veio só do ponto 1. A referência de Drobniak et al. (2026), que
entrou com a declaração AIdIT depois da contagem de 7.996, acrescentou cerca de 100
palavras às referências. Se a declaração de IA (285 palavras) contar para o limite,
estes cortes não bastam. Nesse caso, é melhor perguntar ao editorial office se
declarações entram na contagem.

## Lógica

1. **Trechos recém-acrescentados (cerca de −250):** a definição das hipóteses
   nulas sai de *Scenarios and estimands* e se junta à frase do max test em
   *Methods compared*, onde já havia uma frase sobre as duas permutações, que
   antes ficava duplicada. Também saem a frase repetida nos Resultados e cerca de
   110 palavras da Discussão.
2. **Texto antigo, cortes de baixo custo (cerca de −220, incluindo as referências):** frases de apoio cuja
   citação só aparece nelas, de modo que a referência também sai (Diniz-Filho et
   al. 2012, Molina-Venegas & Roquet 2014, Ai & Ellwood 2022). Somam-se a elas uma
   repetição da Introdução (Feng et al. 2017 sobre a radiação dos anuros) e uma frase
   sobre pools aninhados (Gerhold et al. 2018 continua citado na Introdução).
3. **Abstract:** o item 3 fica sem "although the optima had no phylogenetic
   signal", porque "only niche breadth" já diz isso, e o item 4 fica sem a
   condição repetida.

Nenhum número de resultado sai do texto, exceto as medianas posteriores de poder da
RDA PCPS nos Resultados de S2σ, que estão no Apêndice S2 e não sustentam a
conclusão daquele parágrafo.

Pontos para você julgar:
- **D1** perde a ressalva de que "o nível baixo de γ era mais extremo com N = 200".
  Ela já está nos Métodos (Tree shapes), mas, se preferir manter, são +14 palavras.
- **D1** também perde a frase "Unfortunately, our diagnostic runs cannot tell…".
  O mecanismo continua com hedge ("seems to").
- **D8** tira Ai & Ellwood (2022). A frase seguinte ("Accordingly, more power was
  retained…") continua apoiada em Thompson et al. (2020). Se quiser manter Ai &
  Ellwood, troque esse corte pela frase do deslocamento 0,5 nos Resultados (−20).
- **D4** condensa as três evidências sobre competição numa frase. Revise se o
  resumo de cada estudo está fiel.

## Trechos

### A1 abstract item 3 (31 → 20 palavras, -11)

**Atual:**

> rejected the null hypothesis in up to 0.36 of the replicates with low γ and 200 species, and the PCC in up to 0.22, although the optima had no phylogenetic signal.

**Proposto:**

> rejected in up to 0.36 of replicates with low γ and 200 species, and the PCC in up to 0.22.

### A3 abstract item 3 (PERMANOVA) (8 → 5 palavras, -3)

**Atual:**

> whereas PERMANOVA was less sensitive to tree shape

**Proposto:**

> but not that of PERMANOVA

### A4 abstract item 3 (clades) (5 → 4 palavras, -1)

**Atual:**

> of the replicates on average

**Proposto:**

> of replicates on average

### A2 abstract item 4 (20 → 14 palavras, -6)

**Atual:**

> not reading a significant second PCPS axis or PCC as evidence of conserved optima when niche breadth may be conserved,

**Proposto:**

> not reading a significant second PCPS axis or PCC as evidence of conserved optima,

### M1 scenarios (270 → 156 palavras, -114)

**Atual:**

> The eight trait scenarios were the same in both kinds of simulation (@tbl-scenarios). Every test was evaluated with a site shuffle, whose null hypothesis is that $\mathbf{P}$ is independent of $E$, and a taxa shuffle, whose null hypothesis is that species are exchangeable among the tips of the tree, that is, that no attribute of their response to the gradient ($u_j$, $\sigma_j$ or $h_j$) carries phylogenetic signal. The max test controls the type I error when at least one of the two holds (Peres-Neto & Kembel, 2015; Duarte et al., 2016). S0 and S1 used a gradient independent of the one that generated the data, so the site-shuffle null hypothesis was true, and S2 used the generating gradient with optima and the other attributes without phylogenetic signal, so the taxa-shuffle null hypothesis was true; the rejection rate in these scenarios is therefore the type I error. In S2σ and S2h, the optima had no phylogenetic signal, but niche breadth or maximum abundance did. Thus, both null hypotheses were false, and a rejection is not a type I error of the max test. However, a significant test is usually read as evidence that lineages track the gradient through their optima, and this reading would be wrong in these two scenarios, since no permutation of the tips can separate signal in the optima from signal in the other attributes. We therefore report the rejection rate in S2σ and S2h as the probability that a test suggests conserved optima when only niche breadth or maximum abundance is conserved, and not as a type I error. The rejection rate in S3a–S3σ is the power.

**Proposto:**

> The eight trait scenarios were the same in both kinds of simulation (@tbl-scenarios). S0 and S1 used a gradient independent of the one that generated the data, and S2 used the generating gradient with no phylogenetic signal in any attribute of the species, so the null hypothesis of the max test was true and the rejection rate is the type I error. In S2σ and S2h, only niche breadth or maximum abundance carried phylogenetic signal. Both null hypotheses were then false, but a significant test would still be read as evidence that lineages track the gradient through their optima, since no permutation of the tips can separate signal in the optima from signal in the other attributes. We therefore report the rejection rate in these two scenarios as the probability that a test suggests conserved optima when they are not conserved, and not as a type I error. The rejection rate in S3a–S3σ is the power.

### M2 methods compared (35 → 66 palavras, +31)

**Atual:**

> Every test was evaluated with 999 permutations under two null models, the site shuffle, which permutes the sites, and the taxa shuffle, which permutes the tips of the tree while keeping its topology. We report

**Proposto:**

> Every test was evaluated with 999 permutations under two null models: the site shuffle permutes the sites, under the null hypothesis that $\mathbf{P}$ is independent of $E$, and the taxa shuffle permutes the tips of the tree while keeping its topology, under the null hypothesis that species are exchangeable on the tree, that is, that none of $u_j$, $\sigma_j$ and $h_j$ carries phylogenetic signal. We report

### R1 results S2σ (38 → 0 palavras, -38)

**Atual:**

> Since the taxa-shuffle null hypothesis was false in both scenarios, these rates are not type I errors of the max test, but they show how often each test would suggest conserved optima when only niche breadth was conserved.

**Proposto:**

*(removido)*

### R2 results RDA (47 → 23 palavras, -24)

**Atual:**

> The RDA PCPS rejected in 0.004–0.021 of the replicates in S2σ and lost much less power with low γ than GLM PCPS1 (posterior medians 0.65 with low γ against 0.76 in Yule trees, against 0.34 and 0.74 for GLM PCPS1, at N = 100 in S3b), whereas

**Proposto:**

> The RDA PCPS rejected in 0.004–0.021 of the replicates in S2σ and lost much less power with low γ than GLM PCPS1, whereas

### D1 discussion S2σ (375 → 264 palavras, -111)

**Atual:**

> Phylogenetic signal in niche breadth increased the rejection rate of only the second PCPS axis and the PCC, even though the optima had no phylogenetic signal. Previous evaluations could not detect this, because their simulations drew niche breadth independently of the phylogeny (Duarte et al., 2016, 2018; Peres-Neto et al., 2017). Strictly, these rejections are not type I errors of the max test, since the taxa shuffle assumes that species are exchangeable on the tree, and conserved σ violates this assumption as much as conserved optima do. Permutation tests reject whenever the randomization breaks down a structure of the data, whether or not this structure is the one of interest (Hardy, 2008), and no permutation of the tips can break the link between the phylogeny and the optima while keeping the link with niche breadth. The practical consequence is the same as that of an inflated type I error, because a significant test is usually read as evidence that lineages track the gradient through their optima. When σ was conserved, the second PCPS axis became more correlated with the gradient itself, whereas the curvature of the response, which tied the community-weighted mean of σ to the square of the gradient equally with and without signal, did not seem to drive the rejections (Appendix S2). This pattern seems to load on the axis orthogonal to the basal contrast of the first PCPS axis, which would explain why the increase was largest in trees with low γ and with 200 species, although the low level of γ was also more extreme at N = 200. Unfortunately, our diagnostic runs cannot tell how clades of broad-niche species come to covary linearly with the gradient within a replicate. Rotating the null axes by Procrustes removed this increase but made the test conservative and less powerful, whereas the RDA on both axes stayed below the nominal level and kept most of its power (Appendix S2). Thus, a significant second PCPS axis or PCC should not be taken as evidence of conserved optima when niche breadth may be phylogenetically conserved, and constrained randomizations that preserve phylogenetic structure (Braga et al., 2018) or model-based tests with phylogenetic random effects (Li & Ives, 2017; Tsang & Cadotte, 2025) seem to be promising alternatives.

**Proposto:**

> Phylogenetic signal in niche breadth increased the rejection rate of only the second PCPS axis and the PCC, even though the optima had no phylogenetic signal, which previous evaluations could not detect because they drew niche breadth independently of the phylogeny (Duarte et al., 2016, 2018; Peres-Neto et al., 2017). Strictly, these rejections are not type I errors, since conserved σ violates the exchangeability of species assumed by the taxa shuffle as much as conserved optima do, and permutation tests reject whenever the randomization breaks down a structure of the data, whether or not it is the one of interest (Hardy, 2008). However, the practical consequence is the same, because a significant test is usually read as evidence that lineages track the gradient through their optima. With conserved σ, the second PCPS axis became more correlated with the gradient itself, whereas its correlation with the square of the gradient did not increase (Appendix S2). This pattern seems to load on the axis orthogonal to the basal contrast of the first PCPS axis, which would explain why the increase was largest with low γ and 200 species. Rotating the null axes by Procrustes removed this increase but made the test conservative, whereas the RDA on both axes stayed below the nominal level (Appendix S2). Thus, a significant second PCPS axis or PCC should not be taken as evidence of conserved optima when niche breadth may be phylogenetically conserved, and constrained randomizations (Braga et al., 2018) or model-based tests with phylogenetic random effects (Li & Ives, 2017; Tsang & Cadotte, 2025) seem to be promising alternatives.

### D2 discussion γ (Diniz-Filho 2012) (25 → 0 palavras, -25)

**Atual:**

> This parallels the finding that a small number of phylogenetic eigenvectors may be insufficient to control phylogenetic autocorrelation in ecological analyses (Diniz-Filho et al., 2012).

**Proposto:**

*(removido)*

### D3 discussion imbalance (Molina-Venegas) (48 → 20 palavras, -28)

**Atual:**

> Similarly, the imbalance of the pool tree did not change how much community phylogenetic metrics disagreed (Mazel et al., 2016), and the bias of the net relatedness index caused by poor phylogenetic resolution depended on tree shape and consisted mostly of type II errors (Molina-Venegas & Roquet, 2014).

**Proposto:**

> Similarly, the imbalance of the pool tree did not change how much community phylogenetic metrics disagreed (Mazel et al., 2016).

### D4 discussion competition (50 → 33 palavras, -17)

**Atual:**

> Indeed, competition was weakly related to relatedness in vascular plants (Cahill et al., 2008), average fitness differences, but not stabilising niche differences, showed phylogenetic structure in 18 annual plants (Godoy et al., 2014), and competition can sometimes exclude less related taxa even when traits are conserved (Mayfield & Levine, 2010).

**Proposto:**

> Indeed, relatedness was a weak or inconsistent predictor of competitive outcomes in plants (Cahill et al., 2008; Godoy et al., 2014), and competition can also exclude less related taxa (Mayfield & Levine, 2010).

### D5 discussion pools (Feng) (19 → 3 palavras, -16)

**Atual:**

> , and the major lineages of frogs diversified rapidly at the Cretaceous–Paleogene boundary (Feng et al., 2017), but γ

**Proposto:**

> , but γ

### D6 discussion pools (Gerhold) (33 → 6 palavras, -27)

**Atual:**

> Comparing phylogenetic structure across nested species pools has been proposed as a way to detect the imprint of deep-past processes on local communities (Gerhold et al., 2018), but γ also tends to decrease

**Proposto:**

> However, γ also tends to decrease

### D7 recommendation (iii) (18 → 22 palavras, +4)

**Atual:**

> (iii) interpret the second PCPS axis and the PCC with caution when niche breadth may be phylogenetically conserved;

**Proposto:**

> (iii) not take a significant second PCPS axis or PCC as evidence of conserved optima when niche breadth may be phylogenetically conserved;

### D8 discussion mass effects (Ai & Ellwood) (22 → 4 palavras, -18)

**Atual:**

> Thompson et al. (2020), and a spatially implicit model failed to predict spatially explicit metacommunities at high dispersal (Ai & Ellwood, 2022).

**Proposto:**

> Thompson et al. (2020).

### ref Diniz-Filho, J. A. F., Bini, L (36 → 0 palavras, -36)

**Atual:**

> Diniz-Filho, J. A. F., Bini, L. M., Rangel, T. F., Morales-Castilla, I., Olalla-Tárraga, M. Á., Rodríguez, M. Á., & Hawkins, B. A. (2012). On the selection of phylogenetic eigenvectors for ecological analyses. *Ecography*, 35, 239–249. https://doi.org/10.1111/j.1600-0587.2011.06949.x

**Proposto:**

*(removido)*

### ref Molina-Venegas, R., & Roquet,  (20 → 0 palavras, -20)

**Atual:**

> Molina-Venegas, R., & Roquet, C. (2014). Directional biases in phylogenetic structure quantification: a Mediterranean case study. *Ecography*, 37, 572–580. https://doi.org/10.1111/j.1600-0587.2013.00442.x

**Proposto:**

*(removido)*

### ref Ai, D., & Ellwood, M. D. F. (2 (29 → 0 palavras, -29)

**Atual:**

> Ai, D., & Ellwood, M. D. F. (2022). A spatially implicit model fails to predict the structure of spatially explicit metacommunities under high dispersal. *Ecological Modelling*, 474, 110151. https://doi.org/10.1016/j.ecolmodel.2022.110151

**Proposto:**

*(removido)*

