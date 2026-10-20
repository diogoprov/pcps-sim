# Memória do projeto — Simulações PCPS (capítulo de doutorado)

Última atualização: 2026-10-16

## Resumo

Capítulo da tese de doutorado (PPG Ecologia e Evolução, UFG) sobre as propriedades
do PCPS (Principal Coordinates of Phylogenetic Structure; Duarte 2011). Pergunta
original: como o formato da árvore (stemminess/γ de Pybus, desbalanço/Colless), o
tamanho do pool regional e o padrão de montagem da metacomunidade afetam o
autovalor do PCPS1. Está engavetado desde 2016; em 2026 foi reavaliado para ver se
há sobreposição com Duarte et al. 2016 (MEE).

Autores na última versão: D. B. Provete, V. J. Debastiani (contribuição igual),
L. D. S. Duarte, A. S. Melo.

## Onde estão os arquivos

**Desde 2026-10-16, o projeto vive no repositório git `~/Github/pcps-sim/`.** O
artigo é de autoria única de D. B. Provete; Leandro e Vanderlei não são coautores.

- `docs/`: esta memória e o `PLANO_SIMULACAO.md`.
- `references/`: PDFs e apêndices; **fora do git** (copyright). Inclui o draft de
  dez/2014 de Duarte et al. (`DuarteEtAl_DraftMEE_11Dez2014.pdf`).
- `legacy/`: material de 2012–2016, movido de `~/Downloads/Capitulo 1 - Simulações PCPS/`.
  - `legacy/manuscripts/`: `Provete et al. simulações PCPS MAR_2016.doc` é a
    **versão mais recente**. As versões anteriores (`Properties_*`, com comentários
    de LSD e ASM) também estão aqui.
  - `legacy/R/`, `legacy/figures/`, `legacy/notes/`: scripts e figuras antigos.
  - `legacy/dropbox/`: **cópia** dos arquivos pequenos de
    `~/Library/CloudStorage/Dropbox/Artigo simulação PCPS/` (scripts, `mean_sd_*.txt`,
    `Desenho Simulações PCPS.PDF`).
- **Ficaram no Dropbox**, sem cópia no repositório: os workspaces `res_*_workspace_1..5`
  (cerca de 4,6 GB), `dados_reduzidos.txt`, os gráficos grandes, `pcps_vs_nri.doc`
  (outro projeto) e `Decomposição phylobetadiversity/`.
- **Box** (`PPG Eco-Evo UFG/Capitulo 1 - Simulações PCPS/`): os originais seguem lá,
  intocados. Não baixar nada do Box.

## Desenho das simulações antigas (2013–2014)

- `sim.com()`: para cada comunidade, sorteia uma espécie de referência e filtra as
  candidatas pela similaridade filogenética (agrupamento: Sim = 0.7; sobredispersão:
  Sim2 = 0.85; filtro = 0.95). Abundâncias lognormais. 60 comunidades × 10 spp.
- 6 tipos de metacomunidade (Com_1–6): agrupada, aleatória, sobredispersa, mista
  (20/20/20), gradiente agrupamento→aleatório, gradiente sobredispersão→aleatório.
- Experimento de formato da árvore: `TreeSim::sim.bd.taxa`, 100 tips,
  λ/μ = 0.2/0, 0.35/0.15, 0.55/0.15, 0.75/0.15.
- Experimento de tamanho do pool: Yule com 20, 40, 80 e 160 tips.
- Matriz P via `SYNCSA::matrix.p`; PCoA. Nota: o script usa distância euclidiana,
  mas o texto diz Bray-Curtis com raiz quadrada. Verificar qual foi de fato usada.
- 10 árvores × 10 jobs × 500 réplicas por bloco, em 5 blocos. Rodado num MacPro do
  V. Pillar com o pacote `multicore`, que está obsoleto.
- Análise: ANCOVA de λ1 ~ Colless + γ + tipo de metacomunidade, com Hedges' d.

Resultados principais (MAR/2016): a estrutura da metacomunidade pesa mais que o
formato da árvore (R² = 0.20). Agrupada: d = 1.225. γ: d = 0.125. Colless:
d = −0.076. Tamanho do pool: d = −0.390.

### Problemas conhecidos do manuscrito antigo

- Contradição interna: os Resultados dizem que o γ simulado é "highly negative", mas
  os Métodos relatam γ entre 0.047 e 1.571.
- O espaço de árvores não é realista; o próprio texto admite isso.
- Não há gradiente ambiental, então não dá para medir poder nem erro tipo I.
- O autovalor do PCPS1 tem pouco valor inferencial como variável resposta.
- Hedges' d calculado a partir de correlação parcial e dependente do nível de
  referência: frágil.
- O experimento #3 do desenho (sinal filogenético BM/OU) nunca foi feito.

## Sobreposição com Duarte et al. 2016

Duarte, Debastiani, Freitas & Pillar 2016. Dissecting phylogenetic fuzzy weighting.
*Methods Ecol Evol* 7: 937–946. doi:10.1111/2041-210X.12547. **Lidos (2026-10-16):
o PDF publicado** (`~/Downloads/Methods Ecol Evol - 2016 - Duarte - Dissecting ...pdf`)
**e o Appendix S2** (`~/Downloads/mee312547-sup-0008-appendixs2.docx`).

- O que fazem: revisão teórica (matriz Q, self- e cross-belonging) e simulações no
  estilo de Peres-Neto et al. 2012 (curvas unimodais em E; ótimo dado por um traço F
  simulado com Grafen's ρ = 1 ou 0.0001). São 4 cenários: E+filogenia, só filogenia,
  só E e nenhum. Pool de 50, 100 ou 200 spp.; μtol = 5 ou 10; 100 comunidades;
  incidência; 1000 réplicas; gradiente de ruído (sd 5–150). Testes: Mantel, ADONIS
  e GLM sobre os PCPS (a versão publicada não tem RDA), com site shuffle, taxa
  shuffle e largest P. Comparam com COMDIST, COMDISTNT, Rao e UniFrac (só pool de
  100). A Tabela 2 traz o autovalor relativo do PCPS1, entre 0.31 e 0.57.
- **Árvores (Appendix S2):** `geiger::sim.bdtree`, Yule puro (b = 0.1, d = 0), uma
  árvore nova por réplica e mesmo processo gerador. O formato da árvore varia só
  estocasticamente em torno de Yule; γ e Colless não foram controlados nem
  reportados. Ramos de comprimento 0 foram trocados por 0.01. Só dois níveis de
  sinal (ρ de Grafen).
- Resultado-chave: o GLM tem poder só no PCPS1. No PCPS2, a rejeição pelo largest P
  fica em cerca de 0.01. O ADONIS em D_P é o mais robusto.
- **O draft ≠ o publicado.** A frase do draft de 2014 sobre topologia ("we suspect
  ... influenced by the topology ... deserves attention in future studies", na
  Discussão) **foi removida** da versão publicada. Não dá para citá-la; o argumento
  da lacuna tem que ser nosso.
- **Gancho:** a Introdução publicada afirma que a PFW "incorporates ... phylogenetic
  unbalance" e que isso a diferencia de outros métodos, mas nenhuma simulação varia
  o desbalanço. Essa propriedade central nunca foi testada.
- Sobreposição parcial: tamanho do pool e sinal filogenético.
- **Não sobrepõe: o efeito do formato da árvore.**

## Busca bibliográfica ampla (2026-10-16)

Feita com OpenAlex (citações diretas), Crossref, Consensus e web. **Limitação:** o
OpenAlex tem menos citações indexadas que o Google Scholar, e teses podem ter
escapado. O status de verificação de cada referência está entre parênteses.

**Citações de PCPS/PFW (349 trabalhos únicos varridos):** os citantes de Duarte 2016
(69), Duarte 2011 Oikos (107), Pillar & Duarte 2010 (213), Debastiani & Duarte 2014
e Duarte et al. 2018 Ecology. **Nenhum avalia PCPS/PFW variando Colless, γ ou o
tamanho do pool.** Os mais próximos:
- Debastiani, Bastazini & Pillar 2021 *Oecologia* 196: 633–647,
  doi:10.1007/s00442-021-04972-1 (OpenAlex): efeito da baixa resolução da árvore
  (politomias) sobre K, PD, MPD e MNTD. Não inclui PCPS.
- Debastiani & Duarte 2017 *Evol Biol* 44: 135–143, doi:10.1007/s11692-016-9396-1
  (metadados verificados): teste de Mantel para sinal filogenético.

**Crítica à CWM aplicada à matriz P.** CORREÇÃO (2026-10-16): Duarte et al. 2016
(site shuffle, taxa shuffle e largest P; Tabela 1) e Peres-Neto & Kembel 2015
(permutação de linhas, colunas e linhas/colunas) já aplicam o remédio padrão (max
test). A novidade está nas condições em que esse remédio pode falhar para P: formato
da árvore, sinal filogenético em σ e h sob H0, e testes por eixo. Ver
`PLANO_SIMULACAO.md`, seção 0. Literatura sobre a crítica:
- Peres-Neto, Dray & ter Braak 2017 *Ecography*, doi:10.1111/ecog.02302 (OpenAlex)
- Zelený 2018 *J Veg Sci*, doi:10.1111/jvs.12688 (OpenAlex)
- Duarte, Debastiani, Carlucci & Diniz-Filho 2018 *Ecology* 99: 385–398,
  doi:10.1002/ecy.2081 (OpenAlex): os autores do PCPS admitem a inflação do erro
  tipo I em CWM causada pelo sinal filogenético.
- Tsang & Cadotte 2025 *Ecol Lett* 28: e70267, doi:10.1111/ele.70267 (OpenAlex):
  "compositional autocorrelation" infla o erro tipo I em regressões entre
  comunidades; propõem mixed models.

**Formato da árvore × métricas (nenhum estudo sobre ordenações de phylobetadiversity):**
- Mazel et al. 2016 *Ecography* 39: 913–920, doi:10.1111/ecog.01694 (Crossref):
  efeito de formato e tamanho da árvore em PD, MPD e MNTD (só alfa). É o
  precedente mais próximo.
- Molina-Venegas & Rodríguez 2014 *Ecography*, doi:10.1111/j.1600-0587.2013.00442.x
  (Crossref): stemminess modula o viés de NRI e PD quando a resolução é baixa.
- Diniz-Filho et al. 1998 *Evolution* 52: 1247–1262,
  doi:10.1111/j.1558-5646.1998.tb02006.x (Crossref): no PVR, o número de autovalores
  retidos é "fully explained by cladogram topology". **É a base teórica para a
  pergunta autovalor × formato da árvore.**
- Pavoine 2016 *Oikos* 125: 1719–1732, doi:10.1111/oik.03262 (OpenAlex): compara
  ordenações de PD-dissimilarity ao longo de gradientes, sem variar a árvore. **É a
  referência publicada do evoPCA** (falta confirmar o título).
- **Correções:** Leprieur et al. 2012 **não** trata de tree balance (removido).
  Diniz-Filho et al. 2012 sobre seleção de autovetores saiu na *Ecography* 35:
  239–249, e não na *Evol Biol*.

**Faixas empíricas e ferramentas:**
- γ: McPeek 2008 *Am Nat* (maioria das árvores com γ < 0).
- Desbalanço: Blum & François 2006 *Syst Biol* (beta-splitting de Aldous; β ≈ −1
  **não verificado** no texto completo); Khurana et al. 2023 *Syst Biol*.
- Não há γ nem Colless publicados para anuros. Plano: calcular nas 100 árvores de
  Portik et al. 2023 *MPE*, doi:10.1016/j.ympev.2023.107907.
- Simulação: `treestats` (Janzen et al. 2024 *MPE*, desbalanço controlado),
  modelo β-α-η (Maliet et al. 2018 *Syst Biol*), `apTreeshape`, `DDD`, `TreeSim`,
  `treebalance`, `poweRbal`.
- Pendente: o pacote PCPS (Debastiani & Duarte 2014 *Front Biogeogr* 6) tem dois
  DOIs no Crossref (10.21425/f56322943 e 10.21425/f5fbg22943); confirmar qual é o
  canônico.

**Conclusão (2026-10-16): não houve garfo, mas o desenho antigo está datado.**

## Plano de retomada (proposto, ainda não aprovado)

Comparar três ordinações de composição filogenética sob formato de árvore
controlado:

1. **PCPS** (pacote `PCPS` no CRAN, de V. Debastiani)
2. **Phylogenetic gradient analysis**: Peres-Neto & Kembel 2015, *Plant Ecology*
   216: 709–724, doi:10.1007/s11258-014-0405-0. Li o Appendix A (salvo em
   `~/Downloads/PeresNeto_Kembel_2015_AppendixA.doc`): 100 sítios × 50 spp., uma árvore
   só, sinal via δ de Pagel. **Não há código público** (nem no CRAN, nem no GitHub,
   nem no suplemento); seria preciso reimplementar a partir das equações ou pedir o
   código aos autores.
3. **evoPCA** (Chord/Hellinger): Pavoine 2016, `adiv::evopcachord` e
   `adiv::evopcahellinger`. Pavoine, S. 2016. A guide through a family of
   phylogenetic dissimilarity measures among sites. *Oikos* 125: 1719–1732,
   doi:10.1111/oik.03262. Verificado no OpenAlex e no manual do adiv.

Esboço:

- Árvores com γ e Colless variados de forma independente e em faixa realista:
  diversity-dependence ou early burst para γ negativo; heterogeneidade de taxas ou
  `apTreeshape::rtreeshape` para o desbalanço; árvores empíricas de anuros como
  referência.
- Metacomunidades pelo protocolo Dray & Legendre / Peres-Neto et al. 2012, com sinal
  filogenético no ótimo de nicho, mais o filtro por espécie de referência do
  `sim.com`. Considerar `mcomsimr` ou `ecolottery`.
- Métricas: (i) erro tipo I e poder; (ii) recuperação do clado associado ao
  gradiente, a principal novidade; (iii) concentração de variação nos primeiros
  eixos.
- Pergunta central: o desempenho relativo dos métodos depende do formato da árvore?
- Alvo: MEE ou Ecography. Vanderlei é coautor natural.

**Comparação entre métodos (busca de 2026-10-16):** não encontrei nenhum estudo que
compare PCPS, PGA e evoPCA no mesmo desenho, que varie γ ou Colless ao avaliar métodos
filogenéticos de comunidade × ambiente, ou que use a recuperação do clado como
métrica. Varri os citantes de Peres-Neto & Kembel 2015 (21) e de Pavoine 2016 (~50).
Outros trabalhos próximos: Peres-Neto, Leibold & Dray 2012 *Ecology*,
doi:10.1890/11-0494.1; Root & Nelson 2011 *J Veg Sci*, doi:10.1111/j.1654-1103.2011.01320.x;
Washburne et al. 2017 *PeerJ*, doi:10.7717/peerj.2969 (phylofactorization, que
identifica clados associados a variáveis e seria um bom 4º método); Davies 2021
*Ecol Lett*, doi:10.1111/ele.13682 ("Ecophylogenetics redux"). Não li o texto completo
de Pavoine 2016, Tucker 2017 e Peres-Neto 2012 para descartar que algum deles tenha
variado o formato da árvore.

**A conferir:** a tese de D. B. Provete (2015, UFG; http://repositorio.bc.ufg.br/tede/handle/tede/5747)
está indexada no OpenAlex. Se o capítulo estiver nela, a tese é publicação prévia
(não muda a novidade, mas precisa ser declarada).

## Próximos passos

- [x] Passo 1 (2026-10-16): lidos Duarte et al. 2018 e o resumo de ter Braak et al. 2017 → o item 2 continua aberto (ver plano, seção 0).
- [x] O `mcomsimr` (Thompson et al. 2020) foi avaliado → desenho em duas camadas proposto (plano, seção 3b); aguarda o OK do Diogo.
- [x] Diogo aprovou o desenho em duas camadas: A = estática (Peres-Neto et al. 2012), B = `mcomsimr`.
- [x] Passo 2 (2026-10-16): Peres-Neto & Kembel 2015 lido inteiro → `docs/ESPEC_PCC.md`.
  O `picante` NÃO implementa PGA/PCC, apesar do que diz o artigo (repositório
  conferido). Dois achados para a Introdução:
  (i) p. 717: "tree shape ... Future studies will be important to assess the
  performance ... for different tree shapes", uma chamada publicada para o nosso eixo;
  (ii) p. 721: críticas à PCPS (teste só por linhas, correção de Lingoes, falta de
  ponderação → "very unstable regression estimates", citando simulações não
  publicadas de Peres-Neto & Dray). O objetivo 1 pode testar essa crítica.
- [~] Passo 3 (2026-10-16): `renv` iniciado (R 4.6.1). Funções da camada A em `R/`
  (`pfw.R`, `pcc.R`, `simulate_static.R`, `validation.R`); 7 testes unitários passam
  (`matrix_p` = `SYNCSA::matrix.p`; F e R² = `vegan::adonis2`; Mantel = `vegan::mantel`).
  **Validação rápida (200 réplicas/célula, `analysis/01_validate_layerA.R`):**
  - Duarte et al. 2016, Tabela 1: as taxas de rejeição batem (poder do ADONIS
    0,75/0,90/0,955 contra 0,766/0,871/0,931; env_only com site = 1 e taxa ≈ 0,05;
    Mantel phylo_only com largest P ≈ 0). As médias do Mantel batem
    (0,831/0,887/0,905 contra 0,840/0,884/0,906).
    **Discrepância:** o R² do ADONIS publicado só é reproduzido com Bray–Curtis
    **sem** raiz (0,830 contra 0,819 publicado; com raiz, 0,505), embora os Métodos
    digam "square-rooted". O Mantel bate com a raiz. **Confirmado no Appendix S3**
    (`references/appendices/Duarte_etal_2016_AppendixS3.doc`, 2026-10-17):
    `matrix.p.adonis` chama `adonis(p.matrix ~ envir)` com o `method = "bray"`
    padrão, sem raiz; só `matrix.p.mantel` tem `squareroot = TRUE`. O taxa shuffle
    deles é `picante::taxaShuffle` na matriz de distâncias, equivalente ao nosso.
    O site shuffle do ADONIS são as permutações do próprio `adonis`. A opção
    `sqrt_adonis` fica (padrão TRUE no estudo novo; FALSE para reproduzir Duarte).
    A conferir com 1000 réplicas: taxa shuffle do ADONIS em phylo_only, ~0,045
    contra 0,018–0,027 publicado.
  - Peres-Neto & Kembel 2015, Tabela 1 (PCC): bate com −½·PD sem quadrado, como
    impresso (poder 0,95/0,94 contra 0,945/0,947; env_only 0,055/0,065 contra
    0,054/0,052). Com o quadrado, o poder cai (0,925/0,89). Phylo_only col =
    0,23/0,195 contra 0,267/0,276; depende da árvore única deles.
  - Nota: o `PCPS::pcps` atual (1.0.9) **não** usa Lingoes (descarta autovalores
    negativos). A crítica (5) de Peres-Neto & Kembel pode não valer para a versão
    atual.
- [x] **Validação completa, 1000 réplicas (2026-10-17; 53 min com 16 workers).**
  Scripts: `analysis/01_validate_layerA.R 1000 16` e
  `analysis/02_compare_published.R 1000`. Alvos transcritos em `data/published/`.
  Critério: a taxa publicada cai no IC95% binomial da nossa (n = 1000).
  - **Duarte et al. 2016: validado.** 68 de 72 células (taxas site, taxa, largest P)
    dentro do IC; com cerca de 60 taxas distintas a 95%, ~3 fora seriam esperadas
    só por acaso. Fora: ADONIS both N = 100 (0,892 contra 0,871), ADONIS env_only
    N = 200 site (0,996 contra 0,999), ADONIS phylo_only N = 200 largest P (0,021
    contra 0,013), Mantel env_only N = 200 taxa (0,038 contra 0,059). As médias de
    ρ e R² batem até a 2ª casa em todas as células (ADONIS com BC sem raiz). A
    pendência do taxa shuffle em phylo_only (~0,045 com 200 réplicas) era ruído:
    com 1000, 0,019/0,017/0,026 contra 0,023/0,027/0,018.
  - **Peres-Neto & Kembel 2015 (PCC): validado nas células que não dependem da
    árvore.** neither e env_only batem (exceto neither μtol = 10 row/col: 0,009
    contra 0,004, na borda). O poder (0,925/0,916 sem quadrado; 0,948/0,926 com
    quadrado, contra 0,945/0,947) e o col em phylo_only (0,247/0,216 contra
    0,267/0,276) ficam abaixo, com uma árvore fixa diferente da deles. Com 1000
    réplicas, PD com ou sem quadrado não se distinguem (a vantagem do sem quadrado
    com 200 réplicas era ruído). Mantemos −½·PD como impresso.
  - **Dispersão entre árvores** (`analysis/03_pnk_tree_spread.R 20 200 16`, 2,5 min; formato em `analysis/04_pnk_tree_shape.R`):
    20 árvores Yule (`ape::rphylo(50, 1, 0)`), 200 réplicas cada. Poder do PCC
    entre árvores: 0,71–0,925 (μtol = 5) e 0,74–0,915 (μtol = 10), média ~0,85.
    Col em phylo_only: 0,135–0,275. Os valores publicados ficam no extremo superior,
    compatível com uma árvore única "favorável".
  - **Resultado preliminar para H1 (exploratório: 20 árvores, sem inferência).**
    Spearman entre árvores: poder do PCC × γ de Pybus = 0,90 (μtol = 5) e 0,87
    (μtol = 10); × Colless = −0,53/−0,56. O col em phylo_only também cresce com γ
    (0,69/0,73). Interpretação minha: árvores com nós mais perto das pontas (γ > 0)
    dão autovetores filogenéticos mais "informativos" para o gradiente. O formato da
    árvore já muda o poder do PCC em ~0,2, mesmo dentro de Yule. Precisa ser
    testado no desenho γ × Colless e para PFW/PCPS.
- [x] `renv.lock` criado (R 4.6.1, 101 pacotes; `force = TRUE` porque phangorn
  sugere BiocManager/BiocVersion, que não estão instalados). `.renvignore` exclui
  legacy/, references/, data/ e logs/.

- [ ] **Plano v0.1 em `PLANO_SIMULACAO.md` (2026-10-16): aguarda decisões do Diogo (seção 9) e PDFs (seção 10).** O eixo principal é a validade da inferência sobre P; a comparação de métodos é o objetivo 2. O PCC será reimplementado a partir do PDF e do Appendix A (decisão do Diogo).

- [x] Busca bibliográfica ampla (2026-10-16): lacuna confirmada (ver acima).
- [x] Ler o PDF publicado e o Appendix S2 de Duarte et al. 2016.
- [x] Código de Peres-Neto & Kembel 2015: não há código público.
- [ ] Conferir o texto completo de Pavoine 2016, Tucker et al. 2017 e Peres-Neto et
      al. 2012 quanto à variação do formato da árvore.
- [ ] Decidir se a phylofactorization (Washburne et al. 2017) entra como 4º método.
- [ ] Incluir no desenho o erro tipo I da matriz P sob sinal filogenético (crítica
      de Peres-Neto et al. 2017, Zelený 2018 e Tsang & Cadotte 2025).
- [ ] Escrever o plano de simulação detalhado e confirmar com o Diogo antes de
      programar.
- [ ] Reorganizar como research compendium (`data/`, `R/`, `analysis/`,
      `manuscript/`) com `renv`.

## Referências citadas nesta memória (status de verificação)

| Referência | Status |
|---|---|
| Duarte et al. 2016 MEE 7: 937–946, doi:10.1111/2041-210X.12547 | PDF publicado e Appendix S2 lidos |
| Peres-Neto & Kembel 2015 Plant Ecol 216: 709–724, doi:10.1007/s11258-014-0405-0 | metadados verificados; Appendix A lido |
| Pavoine 2016 Oikos 125: 1719–1732, doi:10.1111/oik.03262 | metadados verificados (OpenAlex e adiv); texto não lido |
| Demais referências da seção "Busca bibliográfica ampla" | metadados verificados (OpenAlex/Crossref), status indicado em cada item; textos não lidos |
| Peres-Neto et al. 2012; Dray & Legendre 2008; Duarte 2011 | citados nos manuscritos; não verifiquei |
| Leprieur et al. 2012 PLoS ONE 7: e42760 | o resumo não trata de tree balance; irrelevante para o formato da árvore |

## Passo 4: árvores (2026-10-17)

- [x] Gerador com γ e Colless independentes (`R/trees.R`). Detalhes em
  PLANO_SIMULACAO.md, seção 2b. Testes em `tests/testthat/test-trees.R`; 18
  expectativas passam no total.
- [x] Referência empírica (`analysis/05_tree_shape_empirical.R`). Árvore de Portik,
  Streicher & Wiens 2023 *MPE* (doi:10.1016/j.ympev.2023.107907; DOI verificado no
  OpenAlex). Cópia em `data/raw/` com SHA-256, vinda do compêndio dos girinos; é
  uma árvore só, não o posterior.
  - Anura tem 5.242 pontas.
  - Clados de 40–250 pontas: γ q10/q50/q90 = −5,6/−3,8/−1,5 (~50 spp.) até
    −8,9/−6,2/−4,3 (~200 spp.). Colless normalizado pelo Yule: q50 de 0,9 a 1,4;
    q90 de 2,4 a 3,0.
  - Subconjuntos aleatórios de espécies dão γ ainda mais negativo (efeito da
    amostragem incompleta).
- [x] Grade 3×3 × N (`analysis/06_tree_grid.R 200 12`, ~1 min com 12 workers).
  São 27 células com 200 árvores cada, em
  `data/derived/trees/tree_grid_k200.rds` (fora do git; regenerar com a seed
  20261017). As médias batem com os alvos e a correlação γ × Colless dentro das
  células fica ≈ 0.
- Rodado em paralelo com o Hmsc do japan-microcosm-bacteria (4 núcleos); a pedido
  do Diogo, usamos 12 workers.
- [ ] Obter o posterior de árvores datadas de Portik et al. 2023, para a incerteza
  dos quantis. **Onde:** projeto OSF `qgcf2` (https://osf.io/qgcf2/). A árvore
  usada no Zenodo da Jessyca (doi:10.5281/zenodo.20297485,
  `antipredator_v2.qmd` l. 212–214) é o arquivo `64ea07f5ac48a005e340256b`
  (download por `https://files.osf.io/v1/resources/qgcf2/providers/osfstorage/<id>?action=download`).
  O OSF fica fora do ar até segunda-feira, 2026-10-19 (status.cos.io). Listar os
  arquivos do projeto quando voltar e conferir se o SHA-256 da árvore única bate
  com o de `data/raw/`.
- [ ] Pools regionais (ex.: espécies do Brasil, do banco de girinos) como
  referência alternativa aos clados.
- [ ] Piloto: 3 células (Yule×Yule, γ low × Colless yule, γ yule × Colless high),
  com PFW (Mantel/ADONIS) e PCC nos cenários de H0. Confirmar o tempo antes.

## Piloto do formato de árvore (2026-10-17)

`analysis/07_pilot_tree_shape.R 200 16 999`: **45 min com 16 workers**, ou 12
CPU-h para 7.200 metacomunidades, ~6 s de CPU por réplica com PFW + PCC a 999
permutações. A estimativa de 25–30 min estava errada: o PCC a 999 permutações
custa mais que o PFW em N = 200.
- **Desenho:** 3 células (Yule×Yule, γ low × Colless yule, γ yule × Colless high),
  N = 50/100/200, 4 cenários, 200 réplicas e uma árvore nova por réplica.
- **Simulação:** ótimos por BM nos comprimentos reais (não Grafen); ADONIS e Mantel
  sobre sqrt(BC). Saída em `data/derived/pilot/pilot_reps200_perm999.csv`.
- **Precisão:** com 200 réplicas, o IC95% em torno de 0,05 é ~0,025–0,08. Tudo
  abaixo é exploratório.

**sqrt(BC) e propriedade euclidiana (checagem numérica, 60 matrizes P):** BC puro
deu autovalores negativos na PCoA em 100% dos casos (até 4% do maior); sqrt(BC),
em 0% (mínimo −1e−16). Decisão do Diogo: sqrt(BC) em todos os testes do estudo
novo.

Resultados:
- **Erro tipo I do componente ambiental** (env_only, taxa shuffle; neither): no
  geral perto de 0,05 nas três células e nos três métodos. Os valores de 0,08–0,095
  (ADONIS env_only γ low N = 100 e Colless high N = 200; PCC env_only Colless high
  N = 200) estão na borda do IC; conferir no estudo completo.
- **PFW no cenário phylo_only (só filogenia):** o taxa shuffle do ADONIS sobe para
  0,115 (γ low, N = 100), 0,08 (γ low, N = 200) e 0,09 (Yule, N = 200). Duarte
  (Grafen, BC sem raiz) tinha 0,018–0,027. O largest P fica ≤ 0,055. O Mantel
  continua baixo (taxa ≤ 0,06; largest P ≤ 0,01).
  Hipótese minha: a mudança vem de BM nos ramos reais + sqrt. Separar os dois
  efeitos.
- **PCC no cenário phylo_only:** o teste de colunas rejeita 0,13–0,39 e cresce com
  N, como em PNK (0,27). O row/col fica em 0,035–0,075. O valor de 0,075 (Colless
  high, N = 200) é possível inflação do max test; conferir.
- **Poder ("both")** — o γ baixo derruba o poder do PCC:

  | N | PCC: Yule | PCC: γ low | PFW: Yule | PFW: γ low |
  |---|---|---|---|---|
  | 50 | 0,865 | 0,325 | Mantel 0,705 | Mantel 0,49 |
  | 100 | 0,99 | 0,73 | ADONIS 0,77 | ADONIS 0,88 |
  | 200 | 1 | 0,94 | — | — |

  - PFW: efeito pequeno e de sinal instável. Com N = 50, cai no Mantel (tabela);
    com N = 100, o ADONIS sobe.
  - O Colless alto quase não muda o poder.
  - **Isso inverte a ordem dos métodos sob γ realista (N = 50: PFW > PCC)** →
    evidência preliminar para H4. Confirma a correlação poder × γ vista em
    `analysis/04`.
- **Estatísticas de árvore** (`analysis/08_tree_stats.R`, 40 s; 98 colunas por
  árvore, em `data/derived/trees/tree_stats_*.csv`): no cenário "both", o −log10(p)
  do PCC correlaciona com γ (ρ ≈ 0,45–0,49), treeness (0,43–0,50), comprimento
  médio dos ramos internos (−0,40 a −0,50) e MPD (−0,41 a −0,50). Com o PFW, |ρ| ≤
  0,19. Isso está confundido com o desenho de só 3 células; é exploratório.
- **Custo do estudo completo:** ~430 mil metacomunidades × ~6 s ≈ 700 CPU-h com
  PCC a 999 permutações. Opções: PCC com 499 permutações, ou só max test por
  réplica. Dimensionar a AWS a partir disso.

Próximos passos sugeridos:
1. Separar o efeito de BM nos ramos reais × Grafen e de sqrt × BC no taxa shuffle
   do ADONIS em phylo_only.
2. Adicionar o GLM sobre eixos PCPS (`pcps.sig`) e os cenários de nuisance (S2σ,
   S2h).
3. Pipeline `targets`, com checkpoint por réplica, para a AWS.

## GLM-PCPS, camada A completa e pipeline (2026-10-17)

- Decisão do Diogo: sqrt(BC) já é definitivo. Não vamos separar o efeito Grafen ×
  sqrt; o estudo completo mede o taxa shuffle em phylo_only com 1000 réplicas.
- **GLM sobre os eixos PCPS** (`R/pcps_glm.R`):
  - O F observado confere com `PCPS::pcps.sig` + `FUN.GLM`.
  - Com um eixo só, o Procrustes simétrico do `pcps.sig` não altera F (testado),
    por isso foi omitido.
  - `test_dp_all()` usa as mesmas permutações para ADONIS, Mantel e GLM.
    RSpectra calcula só os 2 primeiros autovetores.
  - Custo: PFW + GLM juntos levam 2,0/1,7/3,6 s por réplica (N = 50/100/200),
    contra 4,2/3,5/6,7 s separados. O PCC soma 0,6/1,5/4,2 s.
  - **O `PCPS` remove espécies ausentes de todos os sítios** antes de montar Q.
    Adotamos o mesmo no estudo novo: a árvore é podada, como com dados reais.
- **Camada A** (`R/simulate_layerA.R`): os 8 cenários do plano, λ de Pagel, σ e h
  com sinal por ranks (mantendo a marginal), incidência e abundância.
- **Pipeline** (`_targets.R`, `R/pipeline.R`, `analysis/run_pipeline.R <perfil>
  <workers>`):
  - Perfis: smoke, bench, pilot_glm e full; cada um tem seu próprio store,
    `_targets_<perfil>`.
  - Cada branch cruza uma célula de árvores com uma linha (cenário, dado, chunk).
  - Réplica r usa a árvore r da célula, num desenho pareado entre cenários.
  - O perfil smoke rodou de ponta a ponta, e a reexecução não tem nada
    desatualizado.
- **AWS** (`Dockerfile`, `aws/`, instruções em `aws/README.md`): instância spot
  única; o store vai para o S3 a cada 10 min e no aviso de interrupção, e a
  instância se encerra no fim.
  - **O build do Docker ainda não foi testado:** o Docker Desktop não estava
    rodando. Testar localmente, ou com o perfil bench numa c7a.2xlarge, antes do
    full.
- **Custo total previsto do projeto** (estimativa minha): objetivo 1 ≈ R$ 100–300.
  Mais camada B (mcomsimr; ainda sem medida, mais cara por réplica), objetivo 2
  (evoPCA e phylofactorization, ainda sem medida), bench/pilotos na AWS e
  reexecuções para revisores. Total provável de R$ 300–1.000.
- [ ] Mini-piloto do GLM (perfil pilot_glm: 3 células, 8 cenários, incidência e
  abundância, 100 réplicas; ~35–40 min com 16 workers) **agendado para quando os
  jobs do japan-microcosm (PIDs 5814 e 7391) terminarem.**

## AWS (2026-10-17)

- Conta ACCOUNT_ID. O CLI usa o perfil `pcps` (SSO do IAM Identity Center em
  us-east-2, usuário `diogo-cli`, permission set `iam-claude-code`). As máquinas
  rodam em **us-east-1**. O login SSO expira; para renovar, rodar
  `aws sso login --profile pcps`.
- O permission set **não** tem `servicequotas:GetServiceQuota` nem
  `ssm:GetParameters`. Por isso a AMI é resolvida com `ec2 describe-images`. A cota
  spot de vCPU ainda não foi verificada; o Diogo vai conferir no console ou
  acrescentar `ServiceQuotasReadOnlyAccess`.
- Recursos criados com aprovação do Diogo:
  - bucket `pcps-sim-ACCOUNT_ID`: privado, AES256, logs expiram em 90 dias,
    multipart incompleto é abortado em 7 dias;
  - role e instance profile `pcps-sim-ec2`: só o EC2 os assume, com acesso apenas
    a esse bucket.
- Preço spot em 2026-10-17 (us-east-1): c7a.48xlarge a US$ 3,27–4,07/h e
  c7a.2xlarge a US$ 0,18–0,24/h. O custo estimado foi aprovado pelo Diogo.
- **Bench** lançado em 2026-10-17: `i-00194c3103455d319` (c7a.2xlarge spot,
  8 workers, limite de 2 h), código no commit 8645914.
- **Bench na AWS (2026-10-17, 3ª tentativa, `i-07bbca48ab3ec0e6b`): OK.**
  - As tentativas 1 e 2 falharam no build: o `.dockerignore` excluía o
    `_targets.R`, e faltavam bibliotecas de sistema.
  - O rocker 4.6.1 não pegou binários; os 131 pacotes foram compilados em 240 s,
    com 8 vCPU.
  - Calibração das árvores: 5,8 min. Pipeline completo: 5,6 min de wall.
  - **Tempo por réplica (PFW + GLM + PCC, 999 permutações) numa c7a.2xlarge:
    6,0/7,65/14,3 s** para N = 50/100/200. No Mac M5 local, 2,6/3,1/7,9 s: a vCPU
    da AWS é ~2× mais lenta, e não 1,4× como eu tinha suposto.
  - **Estudo completo:** 432 mil réplicas × 9,3 s ≈ 1.120 vCPU-h. Numa
    c7a.48xlarge, ~6 h; spot a US$ 2,6–3,3/h dá ≈ US$ 16–20 (R$ 85–110). Com 32
    vCPU (a cota atual), ~35 h pelo mesmo custo.
  - Cota spot atual: 32 vCPU nas duas regiões; pedidos de 192 abertos em us-east-1
    e us-east-2.
  - Cotas de 192 vCPU spot aprovadas nas duas regiões (2026-10-17).
- **Mini-piloto GLM na AWS (2026-10-17, `i-0fed4d6774cc56ace`, c7a.48xlarge spot,
  us-east-2c, 190 workers, commit 17d41f5): OK.**
  - 3 células (Yule; γ baixo; Colless alto) × N 50/100/200 × 8 cenários ×
    incidência/abundância × 100 réplicas = 14.400 réplicas; ADONIS, Mantel,
    GLM PCPS1 e PCPS2 (sqrt-BC, 999 permutações), sem PCC. Nenhum p ausente.
  - Wall: 15,4 min de pipeline, ~21 min de instância (≈ US$ 1). CPU: 29,5 h,
    ≈ 7,4 s por réplica.
  - Resultados em `data/derived/pilot_glm/results.csv` (fora do git) e no S3.
  - **Erro tipo I do teste max** (α = 0,05):
    - S0 e S1: conservador em todos (0,006–0,026).
    - S2 e S2h: nominal em todos (0,049–0,054).
    - S2sigma: ADONIS, Mantel e GLM PCPS1 conservadores (~0,027), mas **GLM
      PCPS2 inflado: 0,103**. Cresce com N (0,062/0,110/0,138) e chega a 0,152
      na célula de γ baixo. O taxa shuffle no PCPS2 não controla o sinal
      filogenético em σ (amplitude de nicho).
  - **Poder do teste max** (S3a/S3b/S3sigma): ADONIS 0,51/0,73/0,67; Mantel
    0,50/0,66/0,57; GLM PCPS1 0,35/0,45/0,35; GLM PCPS2 ~0,01–0,03.
    - γ baixo derruba o poder do GLM PCPS1 em S3b (0,23 vs 0,56 em Yule); no
      ADONIS o efeito é pequeno (0,73 vs 0,77).
    - Abundância tem menos poder que incidência, sobretudo no GLM (0,28 vs 0,62 em
      S3b).
  - Espécies descartadas por ausência: em média 0,5/0,9/1,7 por N; nenhum sítio
    vazio.
- **Diagnóstico do S2sigma (2026-10-17, `analysis/09_diag_s2sigma.R`, local).**
  S2 (controle) e S2sigma; células Yule e γ baixo; N = 100/200; incidência; 100
  réplicas; 499 permutações; 5,7 min com 6 workers. Resultados em
  `data/derived/diag_s2sigma/`.
  - Taxa de rejeição do taxa shuffle:
    - GLM no PCPS2, um eixo por vez (como no pipeline): S2 0,04–0,08; S2sigma
      0,09–0,30 (0,30 em γ baixo com N = 200). A inflação se reproduz.
    - Com Procrustes simétrico dos eixos nulos 1–2 sobre os observados
      (`pcps.sig` com `choices = 1:2`): 0,00–0,02 nos dois cenários. Remove a
      inflação, mas fica muito conservador também em S2.
    - RDA nos eixos 1–2: S2 0,02–0,05; S2sigma 0,01–0,02.
    - GLM quadrático (env + env²) no PCPS2: 0,03–0,12. Não resolve.
  - CWM σ × env²: |r| ≈ 0,83–0,89 **nos dois cenários**. É um efeito de borda do
    gradiente, que existe com ou sem sinal em σ; não distingue os cenários.
  - Em S2sigma, o PCPS1 fica mais correlacionado com a CWM σ (0,63–0,83, contra
    0,44–0,61 em S2), e a variância do eixo 2 não muda (≈ 0,14).
  - Interpretação (minha): o clado de nicho amplo entra nos eixos principais e
    os gira. No taxa shuffle essa estrutura some, então o PCPS2 nulo não é o
    mesmo eixo que o PCPS2 observado. O problema é a **identidade dos eixos**
    no teste eixo a eixo sem rotação, e não uma associação real com o ambiente.
- **Histórico do git reescrito (2026-10-17)** com `git filter-repo`, antes de
  tornar o repositório público. Saíram de todos os commits: `legacy/`
  (manuscritos e dados não publicados, material de terceiros), a árvore
  `data/raw/portik2023_timetree.tre` (licença não verificada) e o número da conta
  da AWS (`ACCOUNT_ID`). Os dois arquivos continuam locais e estão no
  `.gitignore`; `aws/pack_code.sh` acrescenta a árvore ao tarball. **Os hashes de
  commit citados acima são anteriores à reescrita.** Backup do histórico antigo:
  `~/Github/pcps-sim_pre-filter_20261017.bundle`.
- Licenças: MIT (código) e CC BY 4.0 (texto, figuras, resultados); `CITATION.cff`.
- **Testes extras nos eixos PCPS (2026-10-17).** `test_dp_all(extra_axes = TRUE)`
  acrescenta `glm_pcps1_proc`, `glm_pcps2_proc` (GLM depois do Procrustes simétrico
  dos eixos nulos 1–2 sobre os observados, como `pcps.sig` com `choices = 1:2`),
  `rda_pcps` (RDA nos eixos 1–2) e as colunas `eig_share1/2` (fração da variação
  em cada PCPS). Os testes extras não consomem números aleatórios, então os
  originais não mudam. O custo marginal é desprezível. Perfil `axes`: igual ao
  `full`, sem PCC, com os extras; ≈ 890 vCPU-h, ≈ 4,7 h numa c7a.48xlarge.
- **Pasta Leandro-Vanderlei no Dropbox (inspecionada em 2026-10-17).**
  - Os dois scripts, o PDF do desenho e o `Testes PCPS.ppt` são idênticos às
    cópias em `legacy/`.
  - `simulações pcps .doc` (dez/2012): rascunho do artigo de forma da árvore, com
    Debastiani como coautor de mesma contribuição, Duarte e Melo. O objetivo era o
    efeito de Colless e γ na **fração da variação do PCPS1**. Resultado reportado:
    só γ tinha efeito (negativo), com F(2,10), ou seja apenas 13 árvores. Por isso
    o `eig_share` foi incluído no perfil `axes`.
  - `Artigo com aplicações ecológicas.doc`: esqueleto vazio (Debastiani, Provete,
    Duarte, Melo).
  - `pcps_vs_nri.doc`, `plot_nri_pcps.ppt`, `Result NRI.txt`: outro projeto (NRI
    vs PCPS1 em 12 metacomunidades de girinos da literatura, 167 comunidades).
    Pode servir de referência empírica para a composição das metacomunidades.
  - Os três PDFs `Grafico_*phylo_5_10_10_4_6_500.pdf` (~7 MB cada) não foram
    abertos.
- **Camada B, primeira versão (2026-10-17, `R/simulate_layerB.R`).**
  - Laço próprio, escrito a partir das equações de Thompson et al. 2020, sem
    copiar código: nem o Julia (`plthompson/Meta_com_framework`) nem o mcomsimr
    têm licença.
    - O mcomsimr também tem `rbind` dentro do laço, um bug com `max_r` vetorial e
      um ambiente padrão sem gradiente espacial.
  - Modelo: crescimento Gaussiano em E, com σ = DP (σ = √2·σ_Thompson);
    Beverton–Holt com α[j,k] = efeito de k sobre j; Poisson; emigração binomial;
    redistribuição multinomial ∝ `disp_mat %*% emigrantes`, uniforme por padrão
    (sítios permutáveis); ambiente estático U(0, 100); semeadura Poisson(0,5) nos
    passos 10–100; só o estado final é guardado.
  - Traits como na camada A: u por Pagel λ; σ e r_max (no papel de h, U(2, 8))
    com ou sem sinal, impostos por ranks.
  - Competição: `none`, `random` U(0, inter_max)·0,05, e `phylo`, com os mesmos
    valores atribuídos por rank de similaridade filogenética.
  - Velocidade: 1–2 s por metacomunidade (100 sítios × 100 espécies × 1000
    passos) no M5. Riqueza e número de espécies já estabilizam em 500 passos.
  - Equilíbrio de uma espécie sozinha: 79,5 (esperado 80). 17 testes passam.
  - **Problema aberto: a competição extingue muitas espécies.** Para N = 100,
    espécies restantes (média de 3 árvores, λ = 1):
    - inter_max = 0,05: 61–77
    - inter_max = 0,1: 45–66
    - inter_max = 0,2: 29–59
    - inter_max = 0,5 (Thompson): 15–40
    - sem competição: 89–99
    Com dispersão alta (a = 0,5) perdem-se mais espécies. O N efetivo deixa de
    ser o tamanho da árvore. Decisão pendente com o Diogo.
- **Revisão de código (2026-10-17, skill code-review em duas frentes, 4b69fdc..6d941a8)
  e correções.**
  - Traits e finalização compartilhados entre as camadas A e B
    (`sim_species_traits`, `finish_metacommunity`), `pcoa_axes` compartilhado, o
    site shuffle de `test_dp_all` reusa `stats_dp`, o perfil `axes` é derivado de
    `full` com `modifyList`, `stopifnot` quando `extra_axes` não tem pfw + glm, e
    a dispersão uniforme é calculada sem `disp_mat`.
  - Verificado: camada A, `test_dp_all` (com e sem extras) e `test_pcps_glm` dão
    resultados **idênticos** aos de antes com a mesma semente. A camada B só muda
    pela normalização dos pesos da multinomial: com os pesos antigos, fica idêntica.
  - Ponto não verificado: a ordem dentro do passo (Poisson, depois a dispersão,
    como no mcomsimr). No artigo, a dispersão pode estar dentro da Poisson.
    Precisa do PDF de Thompson et al. 2020, que não está em `references/`.
- **Decisões da camada B (Diogo, 2026-10-17):** competição fraca (inter_max
  0,05–0,1) e **calibrar o tamanho da árvore** para que o N depois da exclusão
  fique ≈ 50/100/200 (caminho 2). Zenodo 10.5281/zenodo.3833035 (código de
  Thompson, v2.0): licença "Other (Open)", sem descrição; não é uma licença de
  reuso definida.
- **Eficiência (2026-10-17).** Profiling (Rprof): na dinâmica da camada B, 59% do
  tempo vai no `%*%`; no `test_dp_all`, 32% vai no `vegan::vegdist`. O BLAS é o
  de referência no Mac e também no Docker: o rocker/r-ver usa `libblas-dev`
  (`install_R_source.sh`).
  - `d_p`: Bray–Curtis como manhattan / (soma das linhas), 2,5× mais rápido que o
    `vegdist`. Diferença ≤ 4e-12 nas estatísticas, p-valores idênticos.
  - `mc_dynamics`: depois da última semeadura, as espécies extintas saem das
    contas. Resultado **idêntico** com a mesma semente e ~35% mais rápido
    (pool 150, competição phylo).
  - Pendente: OpenBLAS no Docker (`libopenblas-dev` +
    `OPENBLAS_NUM_THREADS=1`), a testar num bench depois da rodada `axes`.
- **Calibração da camada B (2026-10-17, 23:33–23:55, 16 workers, 10 réplicas,
  commit 7df9a16).** 66/72 jobs convergiram (±5%). Pool necessário:
  sorting/none = N; mass_effects/none ≈ 1,07 N; sorting/phylo ≈ 1,3 N (50),
  1,7–1,9 N (100), 2,6–3 N (200); mass_effects/phylo ≈ 2–2,4 N (50), 4,2–5 N
  (100; 4/6 sem convergir, realizado 88–97) e, para N = 200, pools de 2.348–2.857
  pontas dão só 148–161 espécies: **teto de coexistência regional** (dispersão
  alta homogeneiza os sítios; mais pool = mais parentes próximos competindo;
  interpretação não testada). Pools tão podados também afastam a árvore
  analisada da Yule.
- **Decisão (Diogo, 2026-10-18): opção 1** — retirar só a célula
  mass_effects/phylo × N = 200 (desenho desbalanceado em uma célula; a
  inviabilidade é resultado). No piloto, regra geral: cai a célula com N
  realizado a mais de 15% do alvo (hoje, exatamente essas 8 linhas); as N = 100
  não convergidas ficam, com o N realizado registrado por réplica.
- **Piloto da camada B (2026-10-18, 00:10–00:59, 16 workers, commit 61cb68f).**
  100 réplicas × 88 células (11 arquétipo × N, 8 cenários), nperm 499, árvores
  Yule; 8.800 réplicas, nenhuma falha, nenhum sítio vazio. Arquivo:
  `data/derived/pilot_layerB/pilot_layerB_reps100_perm499.csv`. EP Monte Carlo de
  uma taxa de 0,05 com 300 réplicas ≈ 0,013.
  - Erro tipo I (max test): ≤ nominal em S0, S1, S2 e S2h, exceto Mantel em S2h
    (0,07; ~1,8 EP). **S2sigma só infla em mass_effects/none**: ADONIS 0,09,
    GLM PCPS1 0,083, RDA 0,083, Mantel 0,07. sorting/phylo é muito conservador
    em S2sigma/S2h (≤ 0,017).
  - A inflação do GLM PCPS2 em S2sigma da camada A (0,10) **não se repete**
    (0,037–0,07). O Procrustes continua mais conservador.
  - O site shuffle sozinho infla muito em S2sigma/S2h (até 1,0); o taxa shuffle
    controla, e o max test herda isso.
  - Poder (S3a/S3b/S3sigma): sorting 0,49–0,91 (Mantel o maior), mass effects
    0,11–0,36 (ADONIS o maior). Sobe com N (S3b ADONIS 0,46 → 0,73 de 50 para
    200). GLM PCPS2 quase sem poder (≤ 0,19).
  - **Confundimento:** a exclusão por competição filogenética poda parentes
    próximos e torna a árvore analisada fortemente negativa em γ (γ médio −0,6 a
    −3,5 com phylo; ≈ 0 sem competição). Forma da árvore e arquétipo estão
    confundidos na camada B; as células de forma precisam controlar a árvore
    *depois* da poda, ou registrar γ realizado como covariável.
- **Rodada full na AWS (2026-10-17 ~19:05 → 2026-10-18 03:17, ~8 h; i-0fd983111cb9e7873).**
  432.000 réplicas completas (8 cenários × 9 células de forma × N 50/100/200 ×
  abundância/incidência; ADONIS, Mantel, GLM PCPS1/2, PCC; 5 testes). Arquivo:
  `data/derived/full/results.csv` (fora do git) e no S3. O token SSO expirou
  ~03:00 e o monitor ficou cego; a instância encerrada some do describe
  ("None"), e o loop esperava "terminated" — monitor novo trata os dois casos.
  - Erro tipo I (max test): S2 e S2h nominais em todos os testes (0,046–0,055);
    S0/S1 conservadores (0,004–0,033).
  - **S2sigma: GLM PCPS2 0,12** (abund. 0,09, incid. 0,15), **cresce com N**
    (0,07 / 0,11 / 0,18) e com γ baixo (0,15 em g_low vs 0,075 em Yule).
    **PCC também infla em S2sigma** (0,08; incidência 0,11; N 200 0,10), mais
    em árvores Yule (0,09–0,096). ADONIS, Mantel e GLM PCPS1 ficam ≤ 0,04.
  - Poder: ADONIS 0,50/0,73/0,67 (S3a/b/σ), Mantel 0,46/0,64/0,55, PCC
    0,40/0,81/0,84, GLM PCPS1 0,28/0,38/0,28, PCPS2 ≈ 0,01. Incidência > abundância
    para ADONIS/Mantel/GLM. γ baixo derruba o GLM PCPS1 (S3b 0,21 vs 0,58 Yule) e
    o PCC (0,68 vs 0,94); ADONIS quase não muda (0,71–0,75). Colless sem efeito
    visível.
- **Rodada axes lançada 2026-10-18 08:58** (i-056c0789ae71860bc, commit 861e862),
  ~5 h.
- **Rodada axes (2026-10-18 08:58 → 15:36, ~6,6 h; i-056c0789ae71860bc).**
  432.000 réplicas completas (1000 por célula; 7 testes: ADONIS, Mantel, GLM
  PCPS1/2, GLM PCPS1/2 com Procrustes, RDA nos eixos 1–2). Arquivo:
  `data/derived/axes/results.csv` (fora do git) e no S3.
  - Repete o full: GLM PCPS2 em S2sigma 0,116 (até 0,21 com incidência e
    N 200; ~0,14 com γ baixo, ~0,075 em Yule).
  - **Procrustes:** elimina a inflação (PCPS2 0,002 em S2sigma), mas é
    **conservador demais em toda parte** (PCPS2 0,002 em S2/S2h, PCPS1 0,024 em
    S2) e custa poder (S3b PCPS1 0,29 vs 0,38 sem Procrustes).
  - **RDA nos eixos 1–2:** controla o erro tipo I (S2 0,050, S2h 0,049, S2sigma
    0,024, S0/S1 ≤ 0,011) e tem poder intermediário (S3a/b/σ 0,41/0,57/0,41),
    abaixo do ADONIS (0,50/0,73/0,67) e acima do GLM PCPS1. Também perde
    poder com γ baixo (S3b 0,51 vs 0,66 em Yule).
  - **eig_share1 depende de γ** (S0: 0,29 em γ baixo vs 0,36 em Yule; S3b: 0,37
    vs 0,52); eig_share2 ≈ 0,14–0,15 em todas as células. Reproduz, na grade
    completa, o efeito de γ no autovalor da PCPS1 do rascunho de 2012
    (Debastiani & Provete) e explica a perda de poder do GLM PCPS1 com γ baixo.
  - Leitura: entre os testes nos eixos, a RDA é a alternativa válida; o GLM por
    eixo isolado não é recomendável para o PCPS2. ADONIS segue o mais poderoso e
    o mais insensível à forma da árvore.
- **OSF de Portik et al. 2023 conferido (2026-10-18).** Componente "Final Trees"
  (osf.io/qgcf2, dentro de osf.io/3cbsh), **licença GPL-3.0** (copyright Daniel
  Portik 2023). `TreePL-Rooted_Anura_bestTree.tre` (osf.io/download/akupm) tem o
  **mesmo SHA-256** de `data/raw/portik2023_timetree.tre`. Não há posterior
  bayesiano: a incerteza vem de `TreePL-Rooted_Anura_bootstraps.tre` (22,9 MB,
  bootstraps datados com TreePL) e dos ICs dos nós. Pendente (confirmar com o
  Diogo): baixar os bootstraps e recalcular os quantis empíricos de γ e Colless
  com a incerteza. `data/raw/README.md` atualizado (sai o `[CITAR: ...]`).
- **BLAS no Docker (2026-10-18, `analysis/12_bench_blas.R`, Docker arm64 no Mac).**
  Correção: o rocker/r-ver 4.6.1 **já liga o R ao OpenBLAS pthread** (arm64 e
  amd64; `libopenblas0-pthread`), não à BLAS de referência como eu tinha dito.
  Sem limite de threads, cada worker do crew dispara threads do OpenBLAS que
  competem com os outros workers. Com `OPENBLAS_NUM_THREADS=1`, por réplica:
  dinâmica B 1,27 → 1,18 s, `test_dp_all` (extras, 999) 6,94 → 3,14 s,
  `test_pcc` (999) 4,08 → 1,24 s; estatísticas e p-valores **idênticos**.
  Dockerfile agora fixa `OPENBLAS_NUM_THREADS=1` e `OMP_NUM_THREADS=1`. Na AWS,
  com 190 workers, o ganho deve ser igual ou maior (não medido lá). O R nativo do
  Mac usa a BLAS de referência (libRblas).
- **Opção B para o γ na camada B (Diogo, 2026-10-18).** A célula da grade define o
  pool: a árvore do pool, com o tamanho calibrado, é sorteada com os alvos de γ
  e Colless da célula para o N alvo, com (β, c) **recalibrados para cada tamanho
  de pool** e a mesma rejeição por tolerância da camada A (sem rejeição pela
  árvore podada). A revisão de especificação mostrou que usar os (β, c) de N num
  pool maior dava γ −17,7 em vez de −5,6 com 450 pontas, recriando o
  confundimento com o arquétipo. Sem competição (pool ≈ N), o pool equivale a uma
  árvore da camada A. Pools independentes entre cenários (sem o pareamento da
  camada A). Cada réplica
  registra `gamma_pool`, `colless_pool`, `gamma` e `colless` da árvore podada e
  `n_realized`; a análise vai modelar a rejeição em função do γ realizado
  (ex.: rejeição ~ s(γ) + arquétipo + N + cenário). A pergunta causal sobre γ
  fica com a camada A. Implementado em `R/pipeline.R` (`design_table_B`,
  `run_chunk_B`), perfis `layerB` (5 células: 4 quinas + centro; 4 arquétipos;
  8 cenários; N 50/100/200 sem mass_effects/phylo × 200; incidência; 1000
  réplicas; PFW, GLM, PCC e extras) e `smoke_B`. Tamanhos de pool agora em
  `data/calibrated/layerB_pool_sizes.csv` (rastreado e copiado para a imagem).
  O aviso de "não convergido" virou `message()` em `feasible_pools()`: as células
  ficam de propósito.
- **Revista-alvo: Methods in Ecology and Evolution** (Diogo, 2026-10-18). Conferir
  as normas de figura da MEE antes da submissão (não verificadas).
- **Figura 1** (2026-10-18): esquema em `manuscript/figures/fig1_schematic.*`,
  gerado por `analysis/13_fig1_schematic.R`. Aprovada pelo Diogo. Pendências de
  desenho: painel 5 é ilustrativo; "999 permutações" e "×1000" valem para full e
  axes.
- **Rodada layerB na AWS** (2026-10-18, 18:43, i-0c37dab3d146dd037, commit
  e6e5188): 8.800 chunks de 50 réplicas, ~490 s por chunk. Às 21:47, 254 dos
  ~3.700 chunks prontos tinham falhado com "subscript out of bounds":
  `sim_tree_cell()` não achava árvore de pool dentro da tolerância (±0,5 em γ e
  Colless) nas 50 tentativas padrão. Taxa de aceitação: ~8–12% nas piores
  combinações célula × n_pool (g_yule__c_high, g_low__c_high; 2000 árvores
  cada), 10–21% nas outras. Com 50 tentativas, ~0,5% das réplicas falham, e
  basta uma para derrubar um chunk de 50.
  Corrigido com `max_tries = 5000` em `run_chunk_B()`, com
  `tar_cue(depend = FALSE)` nos `results` da camada B. Assim, só os chunks que
  falharam rodam de novo quando a instância for relançada sobre o store do S3,
  e os prontos ficam (testado num pipeline de brinquedo). Os chunks que deram
  certo não ficam viesados: aceitar na 1ª ou na 500ª tentativa tira árvores da
  mesma distribuição condicional. O cue precisa sair antes de qualquer mudança
  que deva refazer todos os chunks: ele ignora mudanças em `profile` e nas
  funções chamadas. **TODO:** reverter o cue num commit logo após a reexecução.
  Passos da reexecução: esperar a instância atual terminar (ela sincroniza o
  store), `aws/pack_code.sh` (a instância lê o código do S3, não do git) e
  `aws/launch_spot.sh ... layerB`.
- **Reexecução e fim da layerB (2026-10-19).** A 1ª instância terminou à 01:17
  com 816 de 8.800 chunks com erro; a reexecução (i-004f36ee7b9f1d21f, commit
  7f1ed73, 01:18 → 02:11) refez só esses. **440.000 réplicas completas**, 8
  testes, nenhum p ausente. `data/derived/layerB/results.csv` (fora do git) e
  no S3. ~508 s por chunk. Custo estimado das duas: ~7,5 h × US$ 2,6–3,2 ≈
  US$ 20–24. O `tar_cue(depend = FALSE)` foi revertido: os branches da 1ª
  instância ficaram com o hash antigo de `run_chunk_B()`, então um `tar_make`
  novo nesse store refaria quase tudo (o results.csv já é o final).
  - Pool: γ e Colless da árvore do pool batem com os alvos em todos os
    arquétipos (opção B funciona). Competição filogenética: γ podado mais
    negativo (pool Yule → −1,9/−2,4) e Colless podado ~45% menor com sorting
    (0,65–0,75 vs 1,17–1,33) e ~75% menor com mass effects (0,27–0,36).
  - Erro tipo I (max): tudo ≤ 0,075 exceto, na borda, PCC em S2h sem competição
    (0,076–0,078) e GLM PCPS2 em S2sigma com sorting/none (0,073; 0,096 com
    N = 200; repete a camada A, mais fraco). sorting/phylo: ADONIS e RDA muito
    conservadores em S2/S2h/S2sigma (0,013–0,023). A inflação do ADONIS em
    S2sigma/mass effects do piloto não se confirma (0,057).
  - Poder S3b: sorting ADONIS 0,74–0,81, Mantel 0,82–0,84, PCC 0,82–0,84;
    mass effects ADONIS 0,19–0,27, **PCC ~0,02 (sem poder)**. γ baixo do pool
    derruba GLM PCPS1 (0,76 → 0,34, sorting/none), PCC (0,95 → 0,69) e RDA
    (0,79 → 0,63); ADONIS insensível com sorting. Repete a camada A.
  - Tabular pelo γ realizado confunde com N; falta o modelo com covariáveis.
- **Resumo do projeto em `docs/resumo_projeto.qmd`** (2026-10-19, pedido do
  Diogo): estado, tabelas das três rodadas e próximos passos. HTML renderizado
  fica fora do git.
- **GAM da camada B pelo γ realizado** (2026-10-19, `analysis/14_gam_layerB.R`,
  `R/rejection_gam.R`): `mgcv::gam` (REML) binomial, rejeição (p_max ≤ 0,05)
  ~ s(γ) + s(Colless) da árvore analisada, **um modelo por teste × cenário ×
  arquétipo × N** (704 modelos de 5000 réplicas, ~30 s com 8 núcleos). Versões
  descartadas nas revisões: (1) por teste × cenário com arquétipo e N aditivos
  e média marginal (misturava sorting e mass effects e extrapolava γ entre N);
  (2) curvas sobre todo o intervalo de γ (interpolavam nos buracos entre as
  células do pool). Curvas condicionais à mediana do Colless, IC 95% pontual,
  trechos a > 0,25 de qualquer γ observado ficam NA (metade da tolerância do
  pool). Saídas em
  `data/derived/layerB_gam/` (fora do git).
  - O γ realizado se agrupa nas células do pool (r dentro do estrato
    0,83–0,9997; > 0,99 sem competição), então o contraste entre extremos das
    curvas repete o contraste entre células. S3b sorting/none N = 100: GLM PCPS1
    0,39 → 0,84, PCC 0,70 → 0,99, RDA 0,67 → 0,84; ADONIS e Mantel quase
    planos; com N = 200 ADONIS 0,97 → 0,90 e Mantel 0,98 → 0,95.
  - **Efeito da poda** (`pruning_gamma_slope`: glm por estrato com célula +
    γ do pool + (γ − γ do pool) + Colless; cada réplica sorteia sua árvore do
    pool dentro da tolerância, então a célula sozinha não fixa o pool; versão
    anterior só com a célula foi descartada na 3ª revisão). Só com
    sorting/phylo a mudança de γ pela poda tem efeito próprio: 35 de 72 estratos
    de poder com IC > 0 e nenhum < 0 (16, 15 e 4 com N = 50, 100, 200). PCC até
    0,46 logit/γ, Mantel até 0,29, ADONIS/GLM PCPS1/RDA 0,10–0,22. Sem
    competição a poda quase não muda γ (DP 0,1–0,3); mass/phylo muda muito
    (DP ~1,2) mas só 2 de 48 estratos com efeito. Erro tipo I: 27 de 440 com
    IC excluindo zero, nos dois sentidos (~22 ao acaso).
  - Erro tipo I nas curvas: IC inferior > 0,05 em algum ponto em 41 estratos,
    > 0,075 em 6: GLM PCPS2 em S2sigma sorting N = 200 (sobe com γ baixo, até
    0,13) e PCC em S2h N = 100–200 (até 0,135).
  - Colless com p < 0,01 em 67 de 704 modelos (efeitos pequenos). edf máx. 4,7
    com k = 6 (k' = 5), apertado em poucos estratos; o γ tem só ~3 grupos.
  - Limites: aditivo em γ e Colless (sem `te()`); IC pontuais, sem correção
    para 704 modelos; em S0/S1 a rejeição também cai com γ baixo (testes mais
    conservadores), então parte da perda de poder pode ser calibração.
    Exploratório; o modelo final segue em `brms`.
- **Decisões do Diogo sobre as pendências** (2026-10-19):
  - Bootstraps de Portik: sim. Baixados os 100 bootstraps datados (TreePL) em
    `data/raw/portik2023_timetree_bootstraps.tre` (fora do git; SHA-256 no
    README de `data/raw/`). A árvore única bate byte a byte com o
    `TreePL-Rooted_Anura_bestTree.tre` do OSF. Falta usar os bootstraps para a
    incerteza dos quantis empíricos de γ e Colless.
  - Objetivo 2: **evoPCA e phylofactorization entram**. Diogo pediu mais
    detalhes sobre a métrica de recuperação de clados.
  - **Sem segundo simulador** (revisto no mesmo dia). A camada B já é o
    framework de Thompson et al. 2020, que o Diogo considera suficiente e
    padrão na área (ele cita artigos recentes de Peres-Neto que o usam; não
    verificados aqui).
  - **Recuperação de clados: clado plantado.** Sortear um nó com 10–30% das
    pontas, deslocar o ótimo dessas espécies ao longo de E; profundidade do nó
    como fator (H3). Escore por nó = |correlação| entre pertença ao clado e os
    escores das espécies no eixo mais ligado a E (PCPS, evoPCA, PCC);
    phylofactor dá a aresta. Métricas: acerto do nó, posto do nó plantado e
    Jaccard entre o clado recuperado e o verdadeiro.
  - Pools regionais empíricos: sim, com os dados de literatura do Diogo. Cópia
    local em `~/Downloads/Capitulo X - NRI vs. PCPS/` (projeto NRI vs PCPS):
    ~15 metacomunidades de girinos com composição e árvore datada. As árvores
    têm só 9–34 espécies (γ −2,7 a 1,4; Colless normalizado −0,4 a 2,5), bem
    abaixo do N = 50–200 do desenho. **Decisão do Diogo: opção (b)** — juntar
    as listas por região, atualizar os nomes para a taxonomia de Portik e podar
    a árvore de Portik para formar pools regionais de 50–200 espécies.
  - Figura 1: o painel 5 fica ilustrativo. Rótulos ajustados ("1000 time
    steps"; "×1000 replicates per cell") e o ícone da RDA virou um triplot
    (sítios, espécies como cruzes e a seta de E).
  - **Autoria (correção):** a ideia da simulação foi do Diogo (o PDF com
    anotações à mão é dele). V. J. Debastiani ajudou na implementação e nos
    detalhes, mas o desenho atual é outro. Agradecimentos: V. J. Debastiani,
    L. D. S. Duarte, A. S. Melo (orientador do doutorado), e a banca: M. V.
    Cianciaruso, T. Rangel e T. Quental (conferir nomes completos e iniciais).
- **Figura por célula da camada B** (2026-10-19, `analysis/15_fig_layerB_cells.R`):
  taxa de rejeição por célula do pool (IC de Wilson) no γ realizado mediano,
  barra horizontal = quantis 2,5–97,5% do γ realizado, forma = Colless da
  célula. Uma figura por N em `data/derived/layerB_gam/cell_rates_n*.png`.
  Substitui as curvas do GAM para leitura (sem interpolar entre células).
- **Thompson et al. 2020 lido (2026-10-19; PDF em ~/Downloads).** A ordem da
  camada B confere com as Eqs. 1, 2 e 4: Ñ(t+1) = Poisson(N r / (1 + Σ α N)),
  depois emigração binomial a partir de Ñ(t+1) com probabilidade a_i, depois
  imigração com probabilidade ∝ Σ_{y≠x} E_iy (kernel exponencial no artigo;
  uniforme aqui, por desenho). Semeadura Poisson(0,5) a cada 10 passos nos 100
  primeiros também confere. Diferenças de desenho já documentadas: ambiente
  estático U(0, 100) em vez de autocorrelado no espaço e no tempo, kernel
  uniforme e sítios sem coordenadas.
- **Objetivo 2: escores por método e piloto** (2026-10-19; `R/clade_methods.R`,
  `analysis/16_pilot_clade.R`). Incidência. Cada método de ordenação escolhe,
  entre os 4 primeiros eixos, o mais correlacionado com E:
  - PCPS: correlação de cada espécie (coluna de P) com o eixo (igual a
    `PCPS::pcps()$correlations`, testado).
  - PCC: valores do autovetor filogenético cuja coluna PCC mais se
    correlaciona com E. A 1ª versão combinava todos os autovetores pelos β, o
    que cancela a árvore (vira ranking das espécies por E); corrigido na
    revisão.
  - evoPCA (`adiv::evopcahellinger`, abundance = FALSE com incidência e TRUE
    com abundância): correlação das
    colunas de espécie com o eixo. As colunas de ramo saturam com incidência
    (o ramo do clado está presente onde qualquer descendente está) e os
    `co` escalam com o comprimento do ramo; rejeitadas.
  - phylofactor binomial: objetivo do `gpf(algorithm = "phylo")` com
    cbind(S, F) ~ phylo * env (pvDeviance) calculado para todo nó, com dados
    agregados por sítio × lado da partição; o melhor nó bate com o 1º fator
    do `gpf` (teste). `phylofactor::matrix.to.phyloframe` quebra no R ≥ 4
    (`class(x) != "matrix"`); a tabela longa é montada à mão.
  - Piloto (3 células × N 50/100/200 × profundidade old/young × shift 1/2 ×
    50 réplicas = 1800; 12 núcleos): **1,2 min**, 0,13/0,20/0,47 s por réplica
    para N = 50/100/200. Acerto do nó (shift 2, N = 100, old/young): evoPCA
    0,79/0,73, phylofactor 0,76/0,73, PCPS 0,37/0,12, PCC 0,11/0,11. O PCPS
    acerta mais clados antigos (coerente com a H3); o PCC, limitado a 4
    autovetores, quase não isola clados de 10–30% das pontas. Idade relativa
    média: old 0,67, young 0,49.
  - Pendências: cronometrar cada método; validar o phylofactor contra o `gpf`
    em mais sementes.
- **Objetivo 2: rodada completa local** (decisões do Diogo, 2026-10-19: PCC só
  com os 4 primeiros autovetores; `PhyloFactor` composicional para abundância;
  deslocamentos 0,5/1/2; rodar local às 21h em vez da AWS).
  `analysis/17_run_clade.R`: 9 células × N 50/100/200 × old/young × 3
  deslocamentos × 1000 réplicas com incidência, mais abundância com N = 100:
  216 blocos, 864 mil linhas (4 métodos). Fase 1 sorteia 1000 árvores **novas**
  por célula e N com os (β, c) calibrados da camada A. Desenho pareado: a
  réplica r usa a árvore r e uma semente por (árvore, profundidade, r), então
  clado plantado, ótimos e gradiente são os mesmos entre deslocamentos e tipos
  de dado. Arquivos com `n_reps` no nome, gravação atômica, retomável. Phylofactor
  composicional (`phylofactor_ilr_node_scores`): variância explicada da
  balance ILR de log(x + 0,5) regredida em E, igual ao `choice = "var"`; o
  melhor nó bate com o `PhyloFactor` (teste). Bug achado na revisão: nó
  verdadeiro sem escore (glm sem convergir) ficava com posto 1; agora NA.
- **Objetivo 2: rodada completa feita** (2026-10-19, local, 16 núcleos,
  18h53 → 20h05: árvores 1,7 min, blocos 1h09; ~0,29 s por réplica).
  `data/derived/clade/results.csv` (fora do git): 864.000 linhas (216 blocos ×
  1000 réplicas × 4 métodos), nenhum NA em acerto, posto ou Jaccard. O log tem
  ≥ 50 avisos (provavelmente glm sem convergir em algum nó; nenhum nó plantado
  ficou sem posto).
  - Acerto do nó plantado, incidência, shift 2, antigo/recente: N = 100 —
    evoPCA 0,79/0,66, phylofactor 0,79/0,67, PCPS 0,38/0,15, PCC 0,11/0,09;
    N = 200 — evoPCA 0,88/0,80, phylofactor 0,88/0,82, PCPS 0,42/0,16. Com
    shift 0,5 todos ficam em 0,05–0,21.
  - Abundância (N = 100, shift 2): evoPCA 0,78/0,65, phylofactor
    composicional 0,68/0,56, PCPS 0,32/0,10, PCC 0,09/0,08.
  - **Formato da árvore** (incidência, N = 100, shift 2): o PCPS depende muito
    do γ: clados antigos 0,55 (γ baixo) → 0,41 (médio) → 0,19 (Yule);
    recentes 0,28 → 0,14 → 0,03. evoPCA e phylofactor quase não mudam com γ
    (0,78–0,80 antigos; 0,65–0,69 recentes). Colless pesa pouco (PCPS
    0,22 alto → 0,31 Yule). Leitura preliminar: H3 sustentada para o PCPS
    (γ muito negativo favorece a recuperação de clados basais); H4 sustentada
    (a distância entre PCPS e os demais muda com γ). Falta o modelo formal.
- **Modelos formais (brms) feitos** (2026-10-19; `R/formal_models.R`,
  `analysis/18_formal_models.R`, `analysis/19_shape_contrasts.R`; saídas em
  `data/derived/models/`, fora do git). Binomial agregado por condição (1000
  réplicas): erro tipo I (1350 condições), poder (810) e clados (864). Forma da
  árvore como **fatores** γ (Yule/médio/baixo) × Colless (Yule/médio/alto) com
  interação, porque cada célula tem um só alvo por eixo e a camada A não guardou
  o γ de cada árvore (o PLANO pedia contínuo; decisão minha, a confirmar com o
  Diogo). Tipo I e poder: `test * scenario * (g + cl) + test * (g:cl + n_f +
  data_type)`; clados: `method * (g * depth + cl + g:cl + shift_f + n_f +
  data_type)`. **Beta-binomial** (a 1ª versão com efeito aleatório por condição
  batia no treedepth máximo em 99% das transições). Priors: normal(0, 1,5) nos
  b; intercepto normal(−2,9, 1,5) no tipo I e normal(0, 1,5) nos outros; φ
  gamma(2, 0,01). 4 cadeias × 2000 iterações, cmdstanr.
  - Convergência: R-hat máx. 1,00; ESS bulk ≥ 1237, tail ≥ 1816; 0 divergências;
    0 transições no treedepth máximo. PPC: 97% das condições no IC preditivo de
    95% e 58–66% no de 50% (levemente conservador). Dispersão quasibinomial:
    3,0 (tipo I), 13,8 (poder), 5,0 (clados).
  - Contrastes (N = 100, incidência; mediana [IC 95%]; γ baixo − Yule no Colless
    Yule, Colless alto − Yule no γ Yule). O modelo impõe os mesmos efeitos de
    forma (escala logit) para todo N e tipo de dado; os contrastes em N = 100
    não testam se o padrão vale nos outros N.
    - Tipo I, γ baixo: GLM PCPS2 em S2sigma +0,078 [0,062; 0,095] (0,08 →
      0,16); PCC em S2sigma −0,021 [−0,035; −0,007]; em S1, GLM PCPS2 −0,008 e
      PCC −0,013; em S0, 4 de 5 testes ficam mais conservadores (−0,011 a
      −0,019; no PCC o IC inclui 0). Colless alto: |Δ| ≤ 0,005, ICs com 0.
    - Poder, γ baixo: GLM PCPS1 −0,41 a −0,46; PCC −0,20 a −0,36; Mantel −0,09
      a −0,18; ADONIS ≈ 0 em S3b/S3sigma (−0,09 em S3a); GLM PCPS2 −0,02 (poder
      ≈ 0,01–0,04). Colless alto: −0,01 a −0,08. Ordem dos métodos muda: em S3a
      o GLM PCPS1 cai de 2º para 4º; em S3b o PCC cai de 1º para 3º.
    - Clados, γ baixo: PCPS +0,31 (antigos) e +0,28 (recentes); evoPCA e
      phylofactor ≈ 0; PCC −0,12 (antigos) e +0,08 (recentes). Colless alto:
      −0,02 a −0,07. evoPCA e phylofactor sempre no topo.
  - Leitura das hipóteses (revisada na revisão de spec):
    - H1: sem evidência nesta resolução (Colless não muda o tipo I; ICs ±0,01).
      A H1 fala do taxa shuffle; a resposta modelada é o p_max.
    - H2: não testada como formulada (falta o contraste de cenário S2sigma/S2h
      vs S0/S1). As taxas em Yule sugerem inflação só no GLM PCPS2 (0,084) e
      no PCC (0,098) em S2sigma; o efeito de γ vai em sentidos opostos nos dois.
    - H3: parcial. Com γ baixo o PCPS recupera mais clados basais, mas também
      mais recentes; a concentração da variação no PCPS1 não entrou. Não há
      nível γ > Yule, então a metade "γ positivo" não pode ser testada.
    - H4: sustentada para poder (a ordem muda em S3a e S3b) e, nos clados, só
      entre PCPS e PCC.
  - Limites e próximos passos do modelo: γ e Colless como fatores (o PLANO
    pedia contínuos; a confirmar com o Diogo); sem interações forma × N, forma
    × tipo de dado nem g:cl × cenário; PPC só in-sample (cobertura otimista; a
    dispersão de 13,8 no poder sugere estrutura faltando, provavelmente forma ×
    N). Falta um PPC agrupado por N × célula e o contraste de cenário para H2.
- **Modelos formais, 2ª versão** (2026-10-19, pedido do Diogo: interações
  forma × N e forma × tipo de dado, contraste de cenário para H2, PPC
  agrupado). Diogo aceitou γ e Colless como fatores. Fórmulas: camada A
  `test * scenario * (g + cl) + test * (g:cl + n_f + data_type) +
  test:(g + cl):(n_f + data_type)`; clados `method * (g * depth + cl + g:cl +
  shift_f + n_f + data_type) + method:(g + cl):(n_f + data_type + shift_f)`.
  Beta-binomial, mesmos priors. Convergência: R-hat ≤ 1,005, ESS bulk ≥ 914,
  tail ≥ 1712, 0 divergências, 0 treedepth. PPC in-sample: 95% → 0,97–0,98;
  agrupado por teste/método × N × célula: 95% → 1,00, 50% → 0,94–0,99
  (conservador). Dispersão quasibinomial do poder segue alta (13), então resta
  estrutura não modelada (provavelmente cenário × N ou tipo de dado × cenário).
  Contrastes em `data/derived/models/contrasts.csv` (`analysis/19_shape_contrasts.R`).
  - **H1:** |tipo I − 0,05| com Colless alto − Yule: nenhum dos 75 contrastes
    (teste × cenário H0 × N) tem IC que exclua 0; máx. |mediana| 0,011. Sem
    evidência (não "sem efeito"): só Colless alto vs Yule, e a parte
    comparativa da H1 (taxa shuffle vs métodos só de distância) não foi testada.
  - Contrastes só com incidência; forma × tipo de dado foi ajustada mas não
    reportada. Nos clados, abundância só existe com N = 100, então efeitos de
    tipo de dado em outros N seriam extrapolação. O PPC agrupado (95% → 1,00)
    quase não discrimina.
  - **H2** (S2sigma ou S2h − S2, no Colless Yule; N = 50/100/200): sinal em σ
    infla o tipo I só no GLM PCPS2 (Yule +0,03; γ baixo +0,09/+0,11/+0,16) e
    no PCC (Yule +0,045; γ baixo +0,02; os ICs de Yule e γ baixo se sobrepõem
    no PCC); ADONIS, GLM PCPS1 e Mantel ficam mais conservadores (−0,01 a
    −0,025). Sinal em h: nenhum efeito. **H2 rejeitada como enunciada** ("em
    todos os métodos"); o efeito de σ fica restrito a GLM PCPS2 e PCC. O
    modelo não tem cenário × N, então a variação com N desse contraste vem só
    da escala logit; e o contraste só foi feito no Colless Yule.
  - **Poder, γ baixo − Yule** (N = 50/100/200): GLM PCPS1 −0,32 a −0,46 (piora
    com N); PCC −0,07 a −0,48 (S3b −0,48/−0,23/−0,10 e S3sigma −0,43/−0,19/
    −0,07 melhoram com N; S3a piora, −0,20 → −0,35);
    Mantel −0,01 a −0,23 (melhora com N); ADONIS −0,17 a +0,08 (com N = 200 o
    γ baixo até ajuda). A forma × N importa: os contrastes mudam bastante com N.
  - **Clados, γ baixo − Yule:** PCPS +0,23 a +0,45 (cresce com N; antigos >
    recentes); PCC −0,10 a −0,14 (antigos) e +0,07 a +0,09 (recentes); evoPCA
    e phylofactor ≈ 0 em todos os N.
- **Decisões do Diogo para o manuscrito (2026-10-19):** um único artigo com as
  camadas A e B; a **camada B entra no texto principal como 2º grande
  objetivo** (desempenho dos métodos sob processos de montagem da
  metacomunidade, que os artigos de benchmark e os que propuseram os testes
  não avaliaram), com modelo formal próprio. Pools regionais: juntar as
  metacomunidades empíricas do Diogo em conjuntos maiores e harmonizar a
  nomenclatura (talvez com o pacote AmphiNom). Manuscrito: Métodos e
  Resultados em .docx com figuras no texto, números de linha e de página, na
  voz do Diogo, seguindo DiRenzo et al. 2023 (MEE 14: 203–217) e Williams et al.
  2024 (MEE 15: 1926–1939); resumo dos guias em `docs/GUIAS_RELATO_MEE.md`
  (Williams só pelo resumo: o PDF não pôde ser baixado; o Diogo precisa
  fornecê-lo). Log de uso de IA em `docs/AI_INTERACTION_LOG.md`.
- **Modelo formal da camada B** (2026-10-19, `analysis/20_formal_models_layerB.R`;
  saídas em `data/derived/models/B_<type_I|power>/`). Beta-binomial agregado,
  6 testes (sem as variantes Procrustes, como nas tabelas do resumo), célula do
  pool como fator de 5 níveis (as 4 quinas e o centro). 1ª versão
  (`test * archetype * (scenario + cell) + test * n_f`) descartada na revisão:
  os desvios do PPC se concentravam em N = 50 e 200 e a dispersão do poder era
  10,7. Versão final: `test * archetype * (scenario + cell) + test * n_f *
  (archetype + scenario + cell)`; dispersão 2,2 (tipo I) e 3,9 (poder); desvios
  do PPC espalhados entre os N. Tipo I com max_treedepth 12 (com 10, 38 de 4000
  transições batiam no limite). Convergência: R-hat ≤ 1,01, ESS bulk ≥ 835,
  tail ≥ 1596, 0 divergências, 0 treedepth. PPC agrupado: 95% → 1,00.
  Atenção: `file_refit = "on_change"` ignorava mudanças em `control` e
  recarregou um ajuste antigo sem avisar; agora `file_refit = "always"`.
  Confusão de desenho: a célula é a forma do pool; sob competição filogenética
  a poda muda a árvore analisada, então o contraste de competição mistura o
  processo com a forma. Mass effects com competição não tem N = 200.
  - Poder (S3b, N = 100), arquétipo − sorting sem competição (Yule / γ baixo;
    1ª versão do modelo, a conferir com a final): mass effects sem competição
    derruba todos (PCC −0,95/−0,75; Mantel −0,76/−0,79; ADONIS −0,63/−0,72; GLM
    PCPS1 −0,55/−0,26; RDA −0,58/−0,55), exceto o GLM PCPS2 (+0,09/+0,08, poder
    de base quase nulo); mass effects com competição: PCC −0,93, Mantel −0,55 a
    −0,63. Competição filogenética com sorting (modelo final): −0,12 a +0,11
    nos cenários de poder (S3b: −0,08 a +0,10).
  - Tipo I: no modelo final, nenhuma condição com N = 100 tem IC inferior >
    0,075 (máx. mediana 0,088). As taxas brutas passam do limite quase só com
    N = 200: PCC em S2h (mass effects/Yule 0,143; sorting/Yule 0,111) e GLM
    PCPS2 em S2sigma com sorting (0,104–0,129).
- **Guias da MEE** (docs/GUIAS_RELATO_MEE.md, v0.2 com o PDF completo de
  Williams et al. 2024). Ajustes propostos ao fluxo: MCSE em todas as taxas e
  justificativa de K = 1000 e das 999 permutações; arquivo com DOI (Zenodo);
  estudo de caso com as metacomunidades empíricas (item 9; a decidir);
  alvo estatístico, versões, sementes, ausência de pré-registro e desvios do
  PLANO nos Métodos; corrigir o comentário de `rejection_rates()` (diz que
  calcula IC binomial e não calcula).
- **Pools regionais (definidos em 2026-10-19)** (`analysis/21_regional_pools.R`;
  dados locais em `data/raw/nri_pcps/`, fora do git: compilados da literatura
  e com dados não publicados cedidos ao Diogo; saídas em `data/derived/pools/`).
  Decisões do Diogo: três pools aninhados; estudo de caso com as metacomunidades
  empíricas (pode virar um box); não revisar planilhas; excluir indeterminadas;
  localidades Crivelari = Paraná, Vasconcelos e Motta = noroeste de SP (as demais
  inferidas, marcadas em `datasets.csv`); confirmados Vitreorana uranoscopa,
  Elachistocleis ovalis (Mainara) e Ololygon trapicheiroi (Mainara "Stra"); Both
  tem L. latrans. Fontes dos nomes: .nex de cada metacomunidade; Both pelos nomes
  completos de `Both_tree.txt` (Lma = Limnomedusa macroglossa, Loc = L. latrans);
  Mainara por uma tabela de códigos inferida dos clados de `Mainara_phylo.txt` e
  da fauna do litoral SP/RJ (18 nomes "kept: inferred (Mainara)"). Regras:
  indeterminadas (sp., aff., cf., gr., "Xsp") e códigos ambíguos excluídos; nome
  que o ASW não conhece trocado pelo nome do ASW mais próximo (≤ 2 edições,
  único; pode aplicar também a sinonímia do ASW, por exemplo A. callipygius → A.
  albosignatus, L. chaquensis → L. macrosternum); casamento com Portik pela
  sinonímia do ASW. Fora da árvore de Portik: Crossodactylus dispar, Scinax
  caldarum, Ololygon hiemalis. **123 espécies.** Pools (γ e Colless
  normalizado pelo Yule na árvore de Portik; entre parênteses a variação entre os
  100 bootstraps, mediana [2,5%; 97,5%]; não é IC da estimativa pontual, porque
  os bootstraps diferem em topologia e têm raiz mais profunda):
  - Mata Atlântica (SE + S), n = 60: γ −2,51 (−2,21 [−2,63; −1,81]), Colless 0,67
    (0,86 [0,65; 1,10]);
  - Mata Atlântica + interior de SP, n = 81: γ −3,16 (−2,81 [−3,29; −2,38]),
    Colless 0,67 (0,83 [0,58; 1,14]);
  - América do Sul, n = 123: γ −3,90 (−3,60 [−4,13; −3,19]), Colless 1,99 (2,05
    [1,12; 2,26]).
  Em relação à grade (alvos com N = 50: γ −5,58/−3,79/0, Colless 0/0,94/2,39;
  N = 100: γ −6,54/−4,94/0, Colless 0/1,31/2,74): γ dos três fica entre os
  níveis Yule e médio; Colless entre Yule e médio nos dois primeiros e entre
  médio e alto na América do Sul.
- **Estudo de caso e manuscrito (2026-10-19/20).** `analysis/23_case_study.R`
  (funções de nomes movidas para `R/taxonomy.R`): testes e métodos de clado em
  Bocaina (Provete et al. 2014, Hydrobiologia 734: 69–79), sul do Brasil (Both,
  Melo, Cechin & Hartz 2011, Acta Oecologica 37: 140–145; confirmado pelo Diogo)
  e noroeste de SP ("vitao" = Prado & Rossa-Feres 2014, Austral Ecology 39:
  540–547), incidência, profundidade como gradiente e dossel na Bocaina (o
  gradiente do preprint Provete 2018, bioRxiv 10.1101/252478). Nenhum teste
  rejeita com profundidade (p max ≥ 0,11); dossel na Bocaina: site shuffle
  p ≈ 0,05, taxa shuffle p = 0,23–0,56. γ: −1,59, −1,53, +1,39. Clados: na
  Bocaina PCC, evoPCA e phylofactor escolhem o mesmo clado (Dendropsophus +
  Scinax); no NW de SP PCC e evoPCA, os microhilídeos.
  Manuscrito (`manuscript/methods_results.qmd` → .docx com linhas e páginas
  numeradas; referência `manuscript/reference.docx` gerada por
  `make_reference_docx.py`): Métodos e Resultados na voz do Diogo, com H1–H4,
  versões, ICs, posto/Jaccard e Box 1; Material Suplementar
  (`manuscript/supporting_information.qmd`; Tabelas S1–S3 de
  `analysis/24_si_tables.R`); figuras de `analysis/22_manuscript_figures.R`.
  Atenção registrada no texto: o GLM PCPS2 em S2σ chega a 0,36 de tipo I bruto
  com N = 200 e γ baixo (o modelo, sem cenário × N, subestima).
- **Modelos formais, 3ª versão da camada A (2026-10-20).** Após a revisão de
  estilo (subestimação do tipo I bruto com N = 200), entrou `test:scenario:n_f`
  nos modelos de tipo I e poder; max_treedepth 12 em todos os modelos (com 10,
  uma transição do tipo I batia no limite). Convergência: R̂ ≤ 1,01, ESS bulk
  ≥ 908, tail ≥ 1.477, 0 divergências, 0 treedepth; dispersão quasibinomial
  2,1 (tipo I) e 8,5 (poder, antes 13). Novos contrastes: GLM PCPS2 em S2σ − S2
  com γ baixo +0,05/+0,10/+0,21 (N = 50/100/200; antes +0,09/+0,11/+0,16); PCC
  +0,05 em Yule e +0,02 com γ baixo (N = 100); poder do PCC com γ baixo −0,25
  (N = 100), −0,48 (N = 50), −0,06 (N = 200). Texto, figuras e Tabelas S1–S3
  atualizados.
- **Normas da MEE** (PDF salvo pelo Diogo em `references/MEE_author_guidelines.pdf`,
  fora do git): artigo padrão com 7.000–8.000 palavras incluindo referências;
  espaço duplo, coluna única, linhas e páginas numeradas; página de título em
  arquivo separado (autores, agradecimentos, contribuições, Data Availability);
  resumo numerado 1–4 (≤ 350 palavras) e declaração de dados/código para
  revisão logo abaixo dele; até oito palavras-chave em ordem alfabética;
  figuras e tabelas podem ficar no texto; citações "(Jones, 1998)", com DOI em
  todas as referências; "Fig. 1", "Figs 2–4", "Table 1"; P para valores de
  probabilidade; números de um a nove por extenso. As instruções técnicas de
  figura (formato, resolução) ficam num link que o PDF não traz.

## Reestruturação do manuscrito (2026-10-20)
- Pedido do autor após ler o docx: equações nativas do Word (math LaTeX no qmd → OMML via pandoc), subscritos em math, Methods/Results enxutos, detalhes movidos para Appendix S1 (métodos estendidos) e S2 (resultados estendidos) em manuscript/supporting_information.qmd, que tem lista de referências própria. Results sem subtítulos; Box 1 mantido.
- Figuras: Figs 2, 3 e 5 em YlGnBu com rótulo branco em célula escura; Fig. 4 e Box 1 em Okabe–Ito; ícone RDA da Fig. 1 com quatro setinhas; Box 1 com ponto = mediana bootstrap e losango = árvore pontual. Rótulos σ nos painéis.
- Novo no S1, tirado do código: σ = √2 σ_T de Thompson, imigração multinomial, só o estado final.

## Submissão: anonimização, Zenodo e Discussion (2026-10-20)
- MEE é duplo-cega (confirmado no PDF de normas). Manuscrito e SI sem autor; links do GitHub fora do texto; folha de rosto separada (title_page.qmd) com DOIs do Zenodo (conceito 10.5281/zenodo.23500104; v0.1.0 10.5281/zenodo.23500105).
- `tools/make_review_bundle.sh` gera review_bundle/pcps-sim-review.zip (gitignored), anonimizado e checado; inclui README_REVIEW.md (tools/review_bundle_README.md). Exclui notas internas, fontes do manuscrito e os scripts 21 e 23 (nomeiam provedores de dados). Depósito anônimo no Zenodo a fazer pelo Diogo, no modelo do registro 23438248 (criador "Anonymous, Anonymous").
- PENDENTE (decisão do Diogo): analysis/21_regional_pools.R no repositório público (desde 2026-10-19) lista provedores, localidades e a tabela código→espécie de Mainara (dados possivelmente não publicados).
- Figuras de submissão em manuscript/figures/submission/ (PDF + PNG 600 dpi). Fig. 1 a 180 mm tem textos de ~3 pt: precisa de novo layout.
- Discussion redigida (≈1.650 palavras) com revisão de literatura em docs/lit_review/A–C. Total do texto ≈ 8.500 palavras sem Introduction (limite 7.000–8.000 com referências): cortes necessários.
- Mainara = Jordani, M. X., Melo, L. S. O., Queiroz, C. S., Rossa-Feres, D. C., & Garey, M. V. (2017). Tadpole community structure in lentic and lotic habitats: richness and diversity in the Atlantic Rainforest lowland. Herpetological Journal, 27, 299–306 (sem DOI; PDF na BHS). Publicado; as 20 espécies da tabela de códigos e as coordenadas (23°21'36"S, 44°50'52"W, Picinguaba, Ubatuba) conferem com o artigo (2026-10-20). A tabela ficou em data/raw/nri_pcps/mainara_codes.csv (fora do git).
- (2026-10-20, tarde) Box 1 e os métodos do estudo de caso foram para o Appendix S3 (figura agora é Fig. S1; submission/Figure_S1). Discussion enxuta (~1.300 palavras); mais detalhes de Methods no S1. Total ≈ 8.150 palavras com abstract, legendas, tabela e referências. Introduction redigida (objetivos e H1–H4 saíram de Methods). O depósito anônimo do Zenodo (10.5281/zenodo.23500661) foi feito antes dessa mudança; atualizar com nova versão antes de submeter.
