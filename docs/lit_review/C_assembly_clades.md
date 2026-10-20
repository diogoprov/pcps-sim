# Revisão C — Processos de metacomunidade, detecção de estrutura ambiental/filogenética e identificação de clados

Data: 2026-10-20. Escopo: Discussion do manuscrito MEE (PCPS/PCC sob formatos de árvore e dinâmica de Thompson et al. 2020; recuperação de clado por PCPS, PCC, evoPCA e phylofactorization).

## Método de verificação

- Buscas: Wiley Scholar Gateway (`semanticSearch`, `search_wiley_fulltext`), Consensus (`search`, 2 chamadas), Undermind (`get_orientation`, `list_workspaces`, `search_papers` no workspace "Cap 1 - Marcos - Ecology"), OpenAlex (`/works?search=`). Na primeira tentativa o Undermind falhou porque não passei um workspace válido; depois funcionou.
- **Todos os DOIs abaixo foram conferidos na API do Crossref** (`api.crossref.org/works/<DOI>`). Autores, ano, título, periódico, volume e páginas foram copiados da resposta do Crossref.
- "Fonte da afirmação" indica o que eu de fato li: **abstract** (via OpenAlex ou Consensus), **full text** (PDF local em `references/`) ou **trecho de full text** (passagens do Wiley Scholar Gateway).
- Correções às pistas do briefing:
  - Leibold et al. 2010 saiu na **Ecology Letters**, não na PNAS.
  - A generalização da phylofactorization (Washburne et al. 2019) saiu na **Ecological Monographs**, não na MEE.
  - Não encontrei um artigo de Ovaskainen et al. 2019 chamado "How are species interactions structured in species-rich communities?" com comunidades virtuais. O artigo de Ovaskainen de 2019 que corresponde ao tema se chama "What can observational data reveal about metacommunity processes?" (Ecography). Uso este.
  - No lugar de Pontarp et al. 2019 (TREE, sobre o gradiente latitudinal; pouco pertinente aqui) uso Jang, Brännström & Pontarp 2022, que trata diretamente de dispersão, gradiente e forma da árvore.
  - Pigot & Etienne 2015 é um modelo nulo dinâmico (especiação alopátrica, colonização, extinção). Não é um teste de métodos com modelos de metacomunidade baseados em processos; entra só na lista suplementar.

---

## (1) Referências principais (20)

### A. Teoria de metacomunidades e mass effects

**1. Leibold et al. 2004**
Leibold, M. A., Holyoak, M., Mouquet, N., Amarasekare, P., Chase, J. M., Hoopes, M. F., Holt, R. D., Shurin, J. B., Law, R., Tilman, D., Loreau, M., & Gonzalez, A. (2004). The metacommunity concept: a framework for multi-scale community ecology. *Ecology Letters*, 7(7), 601–613. https://doi.org/10.1111/j.1461-0248.2004.00608.x
- Fonte: abstract.
- Relevância: define os quatro paradigmas (patch dynamics, species sorting, mass effects, neutral). É a referência-padrão para enquadrar os cenários SS vs ME. O abstract não traz previsões quantitativas sobre a relação espécie–ambiente; o "mass effects enfraquecem o filtro ambiental" deve ser citado a partir de Mouquet & Loreau 2003 e Thompson et al. 2020.

**2. Mouquet & Loreau 2003**
Mouquet, N., & Loreau, M. (2003). Community patterns in source-sink metacommunities. *The American Naturalist*, 162(5), 544–557. https://doi.org/10.1086/378857
- Fonte: abstract.
- Afirmação: a diversidade beta e gama "decline as dispersal increases because of increased homogenization of the metacommunity"; a diversidade alfa é máxima em dispersão intermediária.
- Relevância: base teórica direta para a perda de poder sob mass effects. A homogeneização reduz a variação composicional entre sítios que PERMANOVA, db-RDA e PCC precisam para detectar o gradiente.

**3. Cottenie 2005**
Cottenie, K. (2005). Integrating environmental and spatial processes in ecological community dynamics. *Ecology Letters*, 8(11), 1175–1182. https://doi.org/10.1111/j.1461-0248.2005.00820.x
- Fonte: abstract.
- Afirmação: em 158 conjuntos de dados, cerca de 50% da variação composicional foi explicada por ambiente e espaço. A maioria foi classificada como species sorting, seguida da combinação SS + mass effects.
- Relevância: justifica que os dois cenários simulados são empiricamente frequentes. A combinação SS + ME (o nosso cenário ME) é comum na natureza, então a perda de poder que observamos tem consequência prática.

**4. Thompson et al. 2020**
Thompson, P. L., Guzman, L. M., De Meester, L., Horváth, Z., Ptacnik, R., Vanschoenwinkel, B., Viana, D. S., & Chase, J. M. (2020). A process-based metacommunity framework linking local and regional scale community ecology. *Ecology Letters*, 23(9), 1314–1329. https://doi.org/10.1111/ele.13568
- Fonte: full text (PDF local).
- Afirmações (do texto):
  1. Taxas altas de dispersão facilitam dinâmica fonte–sumidouro, "thus promoting the spatial homogenisation of metacommunities".
  2. Em dispersão alta, "the metacommunity homogenises due to source-sink dynamics, spatial and temporal β-richness erodes".
  3. Por ser espacialmente explícito e com ambiente autocorrelacionado, no modelo "mass effects tend to only allow individual species to dominate subregions of the metacommunity", e os autores sugerem que "it is unlikely that dispersal can lead to full homogenisation".
- Relevância: é o simulador que usamos. O ponto (3) ajuda a explicar por que o GLM sobre PCPS2 manteve poder. Se ME produz dominância regional em vez de homogeneização total, um gradiente filogenético amplo ainda pode sobreviver num eixo de ordenação. Essa ligação é interpretação nossa.

### B. Simulações que avaliam métodos de inferência de montagem

**5. Ovaskainen, Rybicki & Abrego 2019**
Ovaskainen, O., Rybicki, J., & Abrego, N. (2019). What can observational data reveal about metacommunity processes? *Ecography*, 42(11), 1877–1886. https://doi.org/10.1111/ecog.04444
- Fonte: abstract.
- Afirmações:
  - Um modelo agent-based gera os quatro paradigmas, e um "virtual ecologist" analisa os dados com 18 métricas (beta diversidade, partição de variação, JSDM).
  - As comunidades mais previsíveis foram as de nicho em paisagens estáticas heterogêneas, "such as metacommunities following the species sorting paradigm or the mass effects paradigm".
  - Os JSDMs deram a maior resolução, e as abordagens carregam informação complementar.
- Relevância: é um contraponto útil. Para eles, ME não é um cenário de baixa previsibilidade, enquanto no nosso estudo ME derrubou o poder da maioria dos testes. Uma reconciliação possível (interpretação nossa): eles medem previsibilidade e capacidade de classificar paradigmas, não o poder de um teste de associação composição–ambiente com um gradiente filogenético plantado.

**6. Mason & Pavoine 2013**
Mason, N. W. H., & Pavoine, S. (2013). Does trait conservatism guarantee that indicators of phylogenetic community structure will reveal niche-based assembly processes along stress gradients? *Journal of Vegetation Science*, 24(5), 820–833. https://doi.org/10.1111/jvs.12033
- Fonte: abstract.
- Afirmações:
  - Com mass effects e priority effects incorporados, "none of the PD indices gave high power to detect assembly processes … even when trait variation was strongly conserved".
  - O conservatismo de traços no pool "was poorly maintained at the meta-community scale".
- Relevância: é o paralelo mais direto para os nossos resultados sob ME. Esta simulação independente também encontrou perda de poder de índices filogenéticos α e β quando mass effects atuam. Eles propõem um mecanismo (perda do sinal filogenético na escala da metacomunidade) que podemos discutir.

**7. Kembel 2009**
Kembel, S. W. (2009). Disentangling niche and neutral influences on community assembly: assessing the performance of community phylogenetic structure tests. *Ecology Letters*, 12(9), 949–960. https://doi.org/10.1111/j.1461-0248.2009.01354.x
- Fonte: abstract.
- Afirmações:
  - Numa simulação de metacomunidade, todos os testes foram sensíveis ao sinal filogenético nas abundâncias e nos traços.
  - Os métodos mais sensíveis a processos de nicho "were also more likely to find non-random patterns … under dispersal assembly".
  - As interpretações "should be modified to account for the influence of dispersal".
- Relevância: é um precedente de que a dispersão confunde testes filogenéticos de comunidade. Sustenta a nossa recomendação de reportar o cenário de dispersão e de não ler perda de poder como ausência de filtro.

**8. Kraft et al. 2007**
Kraft, N. J. B., Cornwell, W. K., Webb, C. O., & Ackerly, D. D. (2007). Trait evolution, community assembly, and the phylogenetic structure of ecological communities. *The American Naturalist*, 170(2), 271–283. https://doi.org/10.1086/519400
- Fonte: abstract.
- Afirmações:
  - Nearest-relative tests têm mais poder para competição; testes de parentesco total funcionam melhor para filtragem de habitat.
  - O tamanho da comunidade local em relação ao pool afeta fortemente o poder: "power increases with larger pool sizes for communities created by filtering but decreases for communities created by competition".
- Relevância:
  - Ajuda a interpretar por que a competição filogenética, que remove parentes próximos, mudou o poder dos testes baseados em composição.
  - Mostra que o desempenho depende de qual escala filogenética o processo afeta (nós terminais vs profundos).

**9. Münkemüller et al. 2020**
Münkemüller, T., Gallien, L., Pollock, L. J., Barros, C., Carboni, M., Chalmandrier, L., Mazel, F., Mokany, K., Roquet, C., Smyčka, J., Talluto, M. V., & Thuiller, W. (2020). Dos and don'ts when inferring assembly rules from diversity patterns. *Global Ecology and Biogeography*, 29(7), 1212–1229. https://doi.org/10.1111/geb.13098
- Fonte: abstract + trecho de full text (Wiley Scholar Gateway, seção "Solutions", sobre JSDMs e covariância residual).
- Afirmações:
  - Identifica 10 armadilhas na inferência a partir de padrões de α-diversidade, a maioria vinda de "an overly simplistic view of potential processes".
  - Recomenda métodos que abranjam escalas espaciais, ambientais e filogenéticas e a integração de modelagem mecanicista.
- Relevância:
  - Justifica o uso de um simulador baseado em processos (Thompson et al. 2020) em vez de só modelos estatísticos de resposta.
  - Apoia a nossa conclusão de que os testes devem ser interpretados à luz do regime de dispersão.

**10. Gerhold et al. 2015**
Gerhold, P., Cahill, J. F., Winter, M., Bartish, I. V., & Prinzing, A. (2015). Phylogenetic patterns are not proxies of community assembly mechanisms (they are far better). *Functional Ecology*, 29(5), 600–614. https://doi.org/10.1111/1365-2435.12425
- Fonte: abstract.
- Afirmações: lista pressupostos fracamente sustentados do uso de padrões filogenéticos como proxy, entre eles:
  - "(iii) trait similarity causes enhanced competition";
  - "(iv) competition causes species exclusion";
  - "(v) communities are at equilibrium";
  - "(vi) … the relative balance of the two can be thus quantified by a single parameter".
- Relevância:
  - O nosso cenário de competição filogeneticamente estruturada impõe exatamente (iii)–(iv), então é um caso favorável.
  - Mesmo assim, competição e filtragem interagiram de forma não aditiva com mass effects: a competição restaurou parte do poder e, ao mesmo tempo, alterou a árvore analisada. Isso ilustra (vi).

**11. Jang, Brännström & Pontarp 2022**
Jang, Y.-T., Brännström, Å., & Pontarp, M. (2022). The interactive effects of environmental gradient and dispersal shape spatial phylogenetic patterns. *Frontiers in Ecology and Evolution*, 10, 1037980. https://doi.org/10.3389/fevo.2022.1037980
- Fonte: abstract.
- Afirmações: num modelo eco-evolutivo,
  - gradiente íngreme com baixa dispersão produz "asymmetric phylogenies" e beta-diversidade filogenética alta;
  - dispersão alta com gradiente suave produz "symmetric phylogenies" e beta filogenética baixa.
- Relevância:
  - Dá suporte independente à ideia de que dispersão e forma da árvore não são independentes.
  - Nós mostramos o lado "ecológico" (a competição podou a árvore analisada e reduziu o Colless em ~75% em relação ao pool); eles mostram o lado "evolutivo".
  - Ressalva: os mecanismos diferem (diversificação vs filtragem/poda), então a analogia deve ser apresentada como tal.

### C. Competição e parentesco

**12. Mayfield & Levine 2010**
Mayfield, M. M., & Levine, J. M. (2010). Opposing effects of competitive exclusion on the phylogenetic structure of communities. *Ecology Letters*, 13(9), 1085–1093. https://doi.org/10.1111/j.1461-0248.2010.01509.x
- Fonte: abstract.
- Afirmação: dois tipos de diferença entre espécies (de nicho e de fitness) têm efeitos opostos, e "competition can sometimes eliminate more different and less related taxa, even when the traits … are phylogenetically conserved".
- Relevância:
  - O nosso modelo implementa só o caso em que a competição exclui parentes próximos (poda de terminais).
  - A Discussion deve dizer que o resultado não generaliza para competição dominada por diferenças de fitness, em que a competição produziria agrupamento filogenético e possivelmente outro efeito no formato da árvore.

**13. Cahill et al. 2008**
Cahill, J. F., Kembel, S. W., Lamb, E. G., & Keddy, P. A. (2008). Does phylogenetic relatedness influence the strength of competition among vascular plants? *Perspectives in Plant Ecology, Evolution and Systematics*, 10(1), 41–50. https://doi.org/10.1016/j.ppees.2007.10.001
- Fonte: abstract (via Consensus; o OpenAlex não tinha abstract).
- Afirmação: em 142 espécies de plantas vasculares, "the competition-relatedness relationship was weak and not significant".
- Relevância: o cenário de competição filogenética é um caso extremo plausível, não o caso típico. Os efeitos de poda que relatamos são, portanto, um limite superior.

**14. Godoy, Kraft & Levine 2014**
Godoy, O., Kraft, N. J. B., & Levine, J. M. (2014). Phylogenetic relatedness and the determinants of competitive outcomes. *Ecology Letters*, 17(7), 836–844. https://doi.org/10.1111/ele.12289
- Fonte: abstract (via Consensus).
- Afirmações:
  - Em 18 anuais da Califórnia, "Stabilising niche differences were unrelated to phylogenetic distance, while species' average fitness showed phylogenetic structure".
  - A coexistência não teve relação com a filogenia.
- Relevância: é evidência empírica de que o pressuposto do nosso cenário (competição mais forte entre parentes) nem sempre vale. Complementa Mayfield & Levine 2010.

### D. Métodos para identificar clados associados ao ambiente

**15. Washburne et al. 2017**
Washburne, A. D., Silverman, J. D., Leff, J. W., Bennett, D. J., Darcy, J. L., Mukherjee, S., Fierer, N., & David, L. A. (2017). Phylogenetic factorization of compositional data yields lineage-level associations in microbiome datasets. *PeerJ*, 5, e2969. https://doi.org/10.7717/peerj.2969
- Fonte: abstract.
- Afirmação: a phylofactorization é ao mesmo tempo redução de dimensionalidade, ordenação e "an inferential tool for identifying edges in the phylogeny along which putative functional ecological traits may have arisen".
- Relevância:
  - Como o método busca diretamente o ramo (edge) que separa clados com respostas distintas, faz sentido que a recuperação do clado plantado (0,44–0,88) não dependa de γ: o método testa cada partição da árvore, não a projeção global da estrutura filogenética.
  - Essa explicação é interpretação nossa.

**16. Washburne et al. 2019**
Washburne, A. D., Silverman, J. D., Morton, J. T., Becker, D. J., Crowley, D., Mukherjee, S., David, L. A., & Plowright, R. K. (2019). Phylofactorization: a graph partitioning algorithm to identify phylogenetic scales of ecological data. *Ecological Monographs*, 89(2), e01353. https://doi.org/10.1002/ecm.1353
- Fonte: abstract (via Consensus).
- Afirmação: generaliza o método para "any community ecological data" (presença–ausência e contagens, via GLM/GAM de família exponencial), e um dos exemplos confirma que "a large clade of Acidobacteria thrive in neutral soils".
- Relevância: é a referência metodológica correta para aplicar phylofactorization a dados de presença/abundância de metacomunidades fora de microbiomas, como fizemos.

**17. Pavoine 2016**
Pavoine, S. (2016). A guide through a family of phylogenetic dissimilarity measures among sites. *Oikos*, 125(12), 1719–1732. https://doi.org/10.1111/oik.03262
- Fonte: full text (PDF local).
- Afirmações (do texto):
  - evoPCA (Hellinger/Chord) aplica PCA centrada às abundâncias de "evolutionary units" (ramos ponderados pelo comprimento).
  - Uma propriedade importante de evoCA/evoNSCA/tb-evoPCA "is to allow for direct biplots with both evolutionary units and sites displayed; the dissimilarities among sites are then directly explained by identifying the evolutionary units that drove those differences".
  - No exemplo com morcegos, nós específicos (p. ex. o nó 21, Stenodermatinae) aparecem associados aos tipos de habitat.
- Relevância: explica mecanicamente por que o evoPCA recupera o clado. A unidade de análise é o ramo, então o ramo-tronco do clado plantado vira diretamente uma variável na ordenação.

**18. Borregaard et al. 2014**
Borregaard, M. K., Rahbek, C., Fjeldså, J., Parra, J. L., Whittaker, R. J., & Graham, C. H. (2014). Node-based analysis of species distributions. *Methods in Ecology and Evolution*, 5(11), 1225–1235. https://doi.org/10.1111/2041-210X.12283
- Fonte: abstract.
- Afirmação: SOS e GND percorrem "each node in the phylogeny and compares the distributions of descendant clades to a null model", com um null model balanceado que é "independent of phylogeny size".
- Relevância: é uma alternativa node-based que não testamos. A ideia de testar nó a nó é a mesma da phylofactorization e coerente com a superioridade dos métodos "por ramo" que encontramos. Candidata à seção de perspectivas.

**19. Leibold, Economo & Peres-Neto 2010**
Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity phylogenetics: separating the roles of environmental filters and historical biogeography. *Ecology Letters*, 13(10), 1290–1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
- Fonte: abstract.
- Afirmação: no método proposto, copépodes calanoides (menos vágeis) refletem legados biogeográficos, enquanto dafnídeos (mais vágeis) refletem adaptação à química da água e à morfologia do lago. A estrutura depende de "factors that affect the colonization rate".
- Relevância:
  - Mostra que a taxa de dispersão do grupo muda qual processo a estrutura filogenética da metacomunidade revela.
  - É o paralelo empírico do contraste SS vs ME. É também o antecessor conceitual da PCC (Peres-Neto & Kembel 2015).

### E. Anfíbios / anuros

**20. Duarte et al. 2014**
Duarte, L. D. S., Both, C., Debastiani, V. J., Carlucci, M. B., Gonçalves, L. O., Cappelatti, L., Seger, G. D. S., Bastazini, V. A. G., Brum, F. T., Salengue, E. V., & Bernardo-Silva, J. S. (2014). Climate effects on amphibian distributions depend on phylogenetic resolution and the biogeographical history of taxa. *Global Ecology and Biogeography*, 23(2), 213–222. https://doi.org/10.1111/geb.12089
- Fonte: abstract.
- Afirmações:
  - Com fuzzy weighting e partição de variação em 262 ecorregiões, a associação entre clima e nós filogenéticos profundos "showed a strong historical legacy", enquanto a distribuição de gêneros (nós rasos) "was mostly associated with climate".
  - Hyloidea aparece associada à origem gondwânica e a temperaturas médias altas.
- Relevância:
  - Aplica PCPS a anfíbios e mostra que o sinal ambiental depende da profundidade do nó.
  - Liga-se à nossa dependência de γ da PCPS: a distribuição de profundidade dos nós (γ) controla quais clados dominam os primeiros PCPS (interpretação nossa).

---

## Referências suplementares (verificadas no Crossref; uso opcional)

- **Leibold, M. A., & Chase, J. M. (2018).** *Metacommunity Ecology* (Monographs in Population Biology, Vol. 59). Princeton University Press. https://doi.org/10.1515/9781400889068
  - Crossref: "Metacommunity Ecology, Volume 59", tipo edited-book. Série e editora confirmadas pelo título da resenha na Quarterly Review of Biology (doi:10.1086/702366).
  - **Não li o livro**: só cite como referência geral, sem atribuir afirmações específicas.
- **Pigot, A. L., & Etienne, R. S. (2015).** A new dynamic null model for phylogenetic community structure. *Ecology Letters*, 18(2), 153–163. https://doi.org/10.1111/ele.12395
  - Abstract: a especiação alopátrica produz sobredispersão filogenética em relação ao null model clássico.
  - Útil se a Discussion tratar de padrões que um processo histórico gera sem competição.
- **Peres-Neto, P. R., Leibold, M. A., & Dray, S. (2012).** Assessing the effects of spatial contingency and environmental filtering on metacommunity phylogenetics. *Ecology*, 93(sp8). https://doi.org/10.1890/11-0494.1
  - O Crossref não retornou páginas; conferir antes de submeter.
  - Abstract: um framework analítico ligando traços/filogenia a variação ambiental e espacial, avaliado com simulações.
- **Duarte, L. D. S., Prieto, P. V., & Pillar, V. D. (2012).** Assessing spatial and environmental drivers of phylogenetic structure in Brazilian *Araucaria* forests. *Ecography*, 35(10), 952–960. https://doi.org/10.1111/j.1600-0587.2011.07193.x
  - Abstract: com PCPS + db-RDA, "the influence of explanatory factors on the phylogenetic structure was concentrated in few nodes" (samambaias arborescentes vs plantas com sementes; malvids). É o precedente empírico de PCPS usada para apontar nós.
- **Parra, J. L., McGuire, J. A., & Graham, C. H. (2010).** Incorporating clade identity in analyses of phylogenetic community structure: an example with hummingbirds. *The American Naturalist*, 176(5), 573–587. https://doi.org/10.1086/656619
  - Abstract: um framework "to identify the species and clades responsible for the inferred patterns of phylogenetic structure" e sua distribuição ao longo de gradientes. É o antecessor do node-based de Borregaard et al. 2014.
- **Hernández-Ordóñez, O., Santos, B. A., Pyron, R. A., Arroyo-Rodríguez, V., Urbina-Cardona, J. N., Martínez-Ramos, M., Parra-Olea, G., & Reynoso, V. H. (2019).** Species sorting and mass effect along forest succession: Evidence from taxonomic, functional, and phylogenetic diversity of amphibian communities. *Ecology and Evolution*, 9(9), 5206–5218. https://doi.org/10.1002/ece3.5110
  - Abstract: anfíbios estruturados principalmente por species sorting, mas 10 espécies presentes em todos os estágios sucessionais e beta baixa, o que "supports a mass effect".
  - É um exemplo empírico de SS + ME em anfíbios, com dimensão filogenética.
- **Leão-Pires, T. A., Luiz, A. M., & Sawaya, R. J. (2024).** Ecological constraints and trait conservatism drive functional and phylogenetic structure of amphibian larvae communities in the Brazilian Atlantic Forest. *Austral Ecology*, 49(9), e13580. https://doi.org/10.1111/aec.13580
  - Abstract: 19 de 33 comunidades de girinos com agrupamento filogenético significativo; diversidade de traços concentrada perto da raiz, indicando conservatismo.
  - Dá suporte à premissa de nicho filogeneticamente conservado em girinos.
- **Leibold, M. A., Rudolph, F. J., Blanchet, F. G., De Meester, L., Gravel, D., Hartig, F., Peres-Neto, P., Shoemaker, L., & Chase, J. M. (2022).** The internal structure of metacommunities. *Oikos*, 2022(1), oik.08618. https://doi.org/10.1111/oik.08618
  - Abstract: usa um modelo processual de colonização–extinção com JSDM para mostrar que espécies e sítios contribuem de forma desigual para a estrutura da metacomunidade.
  - Útil para argumentar que um teste global (PERMANOVA) pode perder um sinal concentrado num clado.
- **Provete, D. B. (2018).** Canopy cover and structural complexity affect the phylogenetic composition of a pond-dwelling Neotropical anuran metacommunity. *bioRxiv*. https://doi.org/10.1101/252478
  - Preprint do próprio autor; não li o conteúdo. Decidir se cabe citar.

---

## (2) Argumentos candidatos para a Discussion (rascunho em inglês)

**A1. Mass effects erode the compositional signal that all tests rely on.**
Dispersal-driven homogenization is a theoretical expectation of source–sink metacommunities: β and γ diversity decline as dispersal increases (Mouquet & Loreau 2003; Thompson et al. 2020). The loss of power we observed under mass effects is therefore expected. PERMANOVA lost 0.62 in power and PCC power nearly vanished. These losses reflect a reduction in the among-site compositional variation that every test partitions, not a failure specific to phylogenetic weighting. Mason & Pavoine (2013) likewise found that no phylogenetic diversity index had high power under mass effects, even with strongly conserved traits, and Kembel (2009) warned that dispersal alters the interpretation of community phylogenetic tests.
*Caveat:* Ovaskainen et al. (2019) classified mass-effects metacommunities among the most predictable ones. We should state that predictability (or classifying paradigms) differs from the power to detect a planted species/clade–environment association.

**A2. Why GLM on PCPS2 kept its power.**
Thompson et al. (2020) note that, in spatially explicit landscapes with autocorrelated environment, mass effects let individual species dominate subregions rather than homogenizing the whole metacommunity. A broad phylogenetic gradient can thus persist on a secondary ordination axis even when the total compositional signal is diluted.
*This is our interpretation*; it should be checked against the variance explained by PCPS1 vs PCPS2 under ME.

**A3. Phylogenetically structured competition is a favourable, not a typical, scenario.**
Our competition scenario assumes that relatedness increases competitive impact and leads to exclusion, the very assumptions Gerhold et al. (2015) flag as weakly supported. Empirical tests often find weak competition–relatedness relationships (Cahill et al. 2008) or find phylogenetic structure in fitness rather than niche differences (Godoy et al. 2014). In that case competition may exclude distant rather than close relatives (Mayfield & Levine 2010). The partial recovery of power under mass effects plus competition should therefore be read as an upper bound, conditional on a competition kernel that prunes close relatives.

**A4. Competition changes the tree being analysed, which confounds tree shape with process.**
Pruning close relatives reduced the Colless index of the analysed tree by ~75% relative to the pool. Hence, in realized communities, the "tree shape" factor is itself an outcome of assembly. This parallels the finding of Jang et al. (2022) that dispersal and gradient steepness jointly determine phylogenetic asymmetry, although their mechanism is diversification and ours is ecological filtering. Kraft et al. (2007) showed that the power of community phylogenetic tests depends on local community size relative to the pool, and that it moves in opposite directions under filtering and competition. Studies comparing methods across tree shapes should therefore report the shape of the realized (analysed) tree, not only that of the pool.

**A5. Branch-wise methods recover clades more reliably than global eigenvector methods.**
evoPCA and phylofactorization recovered the planted clade in 0.44–0.88 of replicates, regardless of γ. Both treat individual branches or edges as units:
- evoPCA ordinates "evolutionary units" and identifies those that drive differences among sites (Pavoine 2016);
- phylofactorization explicitly searches for the edge that best separates clades with contrasting responses (Washburne et al. 2017, 2019).

The node-by-node logic of the SOS/GND framework (Borregaard et al. 2014) and the clade-identity approach of Parra et al. (2010) belong to the same family and could be added to future comparisons.

**A6. PCPS's γ-dependence reflects which nodes dominate the leading eigenvectors.**
PCPS summarises phylogeny-weighted composition along a few axes. Empirical PCPS studies report that explanatory power concentrates in a few deep nodes (Duarte et al. 2012), and that deep and shallow nodes respond to different drivers (Duarte et al. 2014, amphibians of the Americas). Changes in node-depth distribution (γ) will change whether the planted clade's stem node is among those weighted by the leading PCPS axes, which is consistent with its 0.12–0.42 recovery and strong γ-dependence.
*This is our mechanistic interpretation, not a claim from those papers.* PCC (0.08–0.12) was not designed to locate clades: it tests a global phylogeny–composition–environment association (Leibold et al. 2010 is its conceptual antecedent).

**A7. Practical recommendation.**
Following Münkemüller et al. (2020), inference about assembly should account for the dispersal regime and use methods matched to the phylogenetic scale of the question:
- global tests (PCPS/PERMANOVA/PCC) to ask *whether* phylogeny-weighted composition tracks the environment;
- branch-wise methods (evoPCA, phylofactorization, node-based) to ask *which* clades do so.

Empirical anuran metacommunities show both conserved niches (Leão-Pires et al. 2024, tadpoles) and evidence for combined species sorting and mass effects (Hernández-Ordóñez et al. 2019), so both limitations apply to the taxa that motivate this work.

---

## Lacunas e pendências

- Páginas de Peres-Neto et al. 2012 (Ecology 93(sp8)) ausentes no Crossref; conferir no PDF.
- Leibold & Chase 2018 não foi lido; só cite de forma genérica.
- Não encontrei, nesta busca, outro estudo que compare PCPS/PCC com evoPCA/phylofactorization na recuperação de clados. A ausência sustenta a novidade do estudo, mas uma busca negativa não prova inexistência. Sugiro um deep search no Undermind antes de afirmar "first".
- Não verifiquei nenhuma referência sobre sinal filogenético em metacomunidades de girinos além de Leão-Pires et al. 2024 e do preprint do autor; se necessário, use [CITAR: phylogenetic signal tadpole metacommunity].
- LeapSpace (Scopus) report 2 (2026-10-20, query on PFW/PCPS): mostly empirical applications of PCPS (e.g. Brum et al. 2013; Cavalheri et al. 2015; Duarte et al. 2014 PLoS ONE); none added to the Discussion; candidates for the Introduction, to be verified when used.
