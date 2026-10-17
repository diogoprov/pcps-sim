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
