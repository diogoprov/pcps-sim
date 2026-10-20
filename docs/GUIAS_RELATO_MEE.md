# Guias da MEE para relatar estudos de simulação, e como aplicá-los ao pcps-sim

Versão 0.1, 2026-10-19. Documento interno de trabalho, em português. As citações
literais ficam em inglês, como no original.

## 0. Fontes e o que foi realmente lido

| Guia | Fonte usada | Cobertura |
|---|---|---|
| DiRenzo, Hanks & Miller | PDF publicado (versão *Early View*, pp. 1–15, baixado pela UFMS em 06/12/2022) em `~/Documents/Ecologia-Converted2.Data/PDF/4146382399/DiRenzo-2022.pdf`; texto extraído com Ghostscript | **texto completo lido** (Resumo, seções 1–5, Box 1, Box 2, Tabela 1, Data availability) |
| Williams et al. | Só o **resumo**, pela API do DOAJ (registro do artigo); metadados no Crossref | **itens do checklist não recuperados** |

Tentativas que falharam (registro para não repetir):

- **Wiley Scholar Gateway:** 6 buscas (~105 passagens) e nenhuma veio da *Methods
  in Ecology and Evolution* (ISSN 2041-210X). Pelo visto, a MEE não está no corpus.
  O resultado salvo antes (`…1792461979646.txt`) também não tem passagens de
  nenhum dos dois guias.
- **`references/mee_2041-210X.14030.pdf` e `references/mee_2041-210X.14415.pdf`
  não são PDFs.** São páginas HTML do Cloudflare ("Just a moment...", 5820 bytes).
  Convém apagar e trocar pelos PDFs verdadeiros.
- **Wiley `pdfdirect` via curl:** deu 403. **Undermind:** "PDF not available".
  **Chrome:** a extensão não estava conectada. **EcoEvoRxiv:** não achei preprint.
- O Williams-2024.pdf que existe localmente é **outro artigo** (Williams, Zipkin &
  Brodie 2024, *Nat Commun*).
- O Williams et al. 2024 é de acesso aberto (CC BY-NC-ND, segundo o Unpaywall).
  Para completar a seção 2, basta baixar o PDF pelo navegador em
  https://doi.org/10.1111/2041-210X.14415 e salvar em `references/`.

---

## 1. DiRenzo, Hanks & Miller (2023)

**Referência (Crossref):** DiRenzo, G. V., Hanks, E., & Miller, D. A. W. (2023).
A practical guide to understanding and validating complex models using data
simulations. *Methods in Ecology and Evolution*, 14(1), 203–217.
https://doi.org/10.1111/2041-210X.14030. Publicado online em 18/11/2022; o
fascículo impresso é de janeiro de 2023. O PDF lido é o *Early View* (paginação
1–15), por isso as páginas citadas abaixo são as do *Early View*. Tipo:
*Perspective*. Licença: CC BY.

**Natureza do guia.** Não é um checklist de relato numerado. É uma taxonomia de
seis tipos de simulação, cada um com *Objective*, *Simulation settings*, *Model
fitting* e *Example*, mais um algoritmo geral (Box 1) e recomendações gerais. Os
itens abaixo foram derivados por mim dessa estrutura. A numeração é minha.

### 1.1 Itens/recomendações

| # | Item | Citação literal (inglês) | Onde entra |
|---|---|---|---|
| D1 | Declarar o objetivo e as estatísticas de interesse | "It is important to clearly identify the goal of any simulation study" (p. 3) | Métodos (início da seção de simulação) |
| D2 | Descrever os três passos: (1) simular dados com um *data generating model* M; (2) ajustar o modelo/algoritmo A; (3) resumir por Monte Carlo | "(1) simulate one to many unique datasets using a data generating model (referred to as M)" (p. 3; Box 1) | Métodos |
| D3 | Especificar as *simulation settings*: modelo M, parâmetros p, variáveis independentes X e o número K de réplicas | "it is common for simulation studies to consider a range of simulation settings" (Box 1) | Métodos (tabela de cenários) |
| D4 | Definir as saídas/estatísticas O do algoritmo (estimativa, largura ou cobertura do IC, p-valor, estatística de teste) | "may be an estimate of a parameter, the width or coverage of a credible interval or a p value" (Box 1) | Métodos |
| D5 | Calcular o **erro padrão de Monte Carlo** (MCSE; Koehler et al. 2009) e usá-lo para decidir se K basta. IC aproximado: S ± 1,96·MCSE | "This provides a straightforward way to assess whether or not more simulations are needed." (Box 1) | Métodos (justificativa de K); Resultados (IC das taxas) |
| D6 | Checagem visual de K: gráfico do MCSE para K crescente, para mostrar convergência | Box 1: "A good visual check of the effect of replication size K" | Supporting Information |
| D7 | **Basic validation simulation**: ajustar o modelo aos dados, simular um conjunto com as estimativas e reajustar, para checar identificabilidade, viabilidade e erros de código | "determine whether the model, fitting algorithm and code can generate realistic parameter estimates" (p. 5) | Métodos/SI. Os autores pedem isso no mínimo, sempre que houver modelo novo (p. 13) |
| D8 | **Determining statistical properties**: um *full simulation study* com centenas a milhares de réplicas | "100s to 1000s of simulations are performed" (p. 9) | Métodos |
| D9 | Métricas de desempenho definidas explicitamente. *Accuracy* (erro absoluto médio, MSE), *precision* (largura do IC), *bias* (média das estimativas − verdade), *coverage* (proporção de ICs que contêm o verdadeiro) | Coverage: "how often does the true parameter value fall within the range of the 95% CI?" (p. 9) | Métodos (definições); Resultados |
| D10 | **Goodness-of-fit** por simulação (PPC, DHARMa, p-valor de Monte Carlo) | "determine if the statistical model applied in the analysis can generate the observed data" (p. 10) | Métodos; Resultados/SI |
| D11 | **Simulation-based study design / power**: variar p, X e n; dá para usar *space-filling designs* (hipercubo latino) | "may include the model parameters (p), the explanatory variables (X) or the sample sizes" (p. 11) | Métodos |
| D12 | Definição operacional de poder (no exemplo: proporção de IC 95% que exclui zero, com o sinal correto) | "did not overlap zero and showed the same sign" (p. 11) | Métodos (definição) |
| D13 | **Assessing statistical robustness**: simular sob *misspecification* e comparar com um *baseline* em que M coincide com o modelo de análise | "compared to data from a data generating model that matches the statistical model as a baseline" (p. 12) | Métodos (desenho de cenários); Resultados |
| D14 | **Comparing the efficacy of different approaches**: o mesmo conjunto simulado é analisado por vários modelos; comparar por MSE, MSPE e cobertura | "the data are analysed using multiple statistical models" (p. 12) | Métodos; Resultados |
| D15 | Relatar requisitos e tempo computacionais (Tabela 1, em "Determining statistical properties") | Tabela 1: "What are the computational requirements/time of running the model?" | Métodos/SI |
| D16 | Publicar o código que simula e ajusta, para reprodutibilidade (o próprio artigo põe o código em Supplement S1 e numa *USGS software release* com DOI) | "the inclusion of model code increases reproducibility and transparency" (p. 13) | Data availability; SI |

**Contexto empírico do guia (do próprio artigo):** em 50 artigos da MEE de 2021,
18 de 23 que propunham método novo usaram simulação, mas só 5 de 23 mostraram
recuperação de parâmetros conhecidos e 9 de 23 mostraram propriedades
estatísticas (acurácia, precisão, viés, cobertura) (p. 2; Tabela S2).

**Limite de aplicação.** O guia trata de validar **modelos de estimação**
(hierárquicos, bayesianos). O pcps-sim avalia **testes de hipótese** por
permutação. Viés e cobertura de parâmetros (D9) só se aplicam diretamente aos
modelos `brms` de segunda etapa, não aos testes. Para os testes, as estatísticas
de interesse (D4) são a taxa de rejeição (erro tipo I e poder). Essa adaptação é
interpretação minha.

---

## 2. Williams et al. (2024)

**Referência (Crossref):** Williams, C., Yang, Y., Lagisz, M., Morrison, K.,
Ricolfi, L., Warton, D. I., & Nakagawa, S. (2024). Transparent reporting items
for simulation studies evaluating statistical methods: Foundations for
reproducibility and reliability. *Methods in Ecology and Evolution*, 15(11),
1926–1939. https://doi.org/10.1111/2041-210X.14415. Online em 30/09/2024.

### 2.1 O que o resumo permite afirmar (verificado, DOAJ)

| # | Ponto | Citação literal | Onde entra |
|---|---|---|---|
| W-a | São propostos **11 itens de relato**, que estendem a orientação de Morris et al. (2019, *Stat Med* 38: 2074) | "we propose 11 reporting items for statistical simulation studies" | Métodos inteiro |
| W-b | Os itens cobrem **três estágios** | "These items span across three stages: planning, coding and analysis." | Métodos (planejamento, código) e Resultados (análise) |
| W-c | O artigo esclarece a terminologia dos componentes estatísticos e os propósitos dos estudos de simulação | "We also clarify the terminology related to statistical components" | Métodos (vocabulário) |
| W-d | Num levantamento de 100 artigos de ecologia e evolução, poucos relataram incerteza de Monte Carlo | "only a small proportion of articles reported a Monte Carlo uncertainty (17%; 17 out of 98)" | Resultados (MCSE nas taxas) |
| W-e | Muitos artigos não deram o código | "32% (32 out of 100) articles did not provide code" | Data/code availability |
| W-f | Propõe **registro** (*registration*) de estudos de simulação | "we propose the registration of statistical simulation studies" | Métodos (declaração de registro ou de desvios do plano) |
| W-g | Recomenda ferramentas de ciência aberta para a reprodutibilidade | "we discuss the benefits of open science tools" | Data availability; SI |

### 2.2 Os 11 itens

**Não recuperado.** Não tive acesso ao texto completo. Por isso não listo os
nomes dos 11 itens, a distribuição deles entre os três estágios, a redação de
cada item, nem o que o artigo diz sobre *seeds*, versões de software, número de
réplicas ou apresentação de resultados. Não inferi esses itens a partir de
Morris et al. 2019 nem de outras fontes. A estrutura ADEMP (*aims,
data-generating mechanisms, estimands, methods, performance measures*) é
atribuída a Morris et al. 2019. **Não verifiquei**, nestas fontes, se Williams et
al. organizam os itens por ADEMP.

---

## 3. Como aplicar ao pcps-sim

Fontes do projeto: `docs/PLANO_SIMULACAO.md` (v0.2), `docs/resumo_projeto.qmd`
(2026-10-19), `docs/MEMORIA_PROJETO.md`, `_targets.R`, `R/pipeline.R`,
`analysis/17_run_clade.R`, `analysis/18_formal_models.R` e `renv.lock`.
Legenda: **OK** = o projeto já tem e é só relatar; **PARCIAL** = existe, mas
falta algo; **LACUNA** = não existe ou não está documentado.

| Item | O que o projeto tem | Estado | O que fazer no manuscrito |
|---|---|---|---|
| D1, W-a (objetivo) | Objetivo 1 (validade da inferência sobre P/PCPS: tipo I e poder) e objetivo 2 (comparação de métodos, recuperação de clados); hipóteses H1–H4 (PLANO §1) | OK | Abrir a seção de Métodos com os objetivos e as hipóteses. Dizer que H1–H4 foram formuladas antes das rodadas (datas no git) |
| D2, D3 (*data-generating mechanism*) | Duas camadas. **A:** estática, Peres-Neto, Leibold & Dray 2012, 100 sítios, E ~ U(0,100), resposta unimodal + Poisson. **B:** dinâmica de Thompson et al. 2020, reimplementada, com 4 arquétipos (*sorting*/*mass effects* × sem competição/competição filogenética). 8 cenários (S0, S1, S2, S2sigma, S2h, S3a, S3b, S3sigma); grade 3 × 3 γ × Colless; N = 50/100/200; incidência e abundância | OK | Tabela de cenários no texto. Fórmulas completas e o gerador de árvores (beta-splitting + taxa k^c, tolerância ±0,5) no SI. Declarar que o gerador é escolha de desenho própria, não um modelo publicado (PLANO §2b) |
| D3 (fatores fixos e excluídos) | Espaço/dispersão fora da camada A; célula *mass effects* + competição + N = 200 retirada por inviabilidade | PARCIAL | Relatar as exclusões e o motivo (teto de coexistência) nos Métodos, e o desvio do plano original (biblioteca de vários geradores → gerador de dois controles) |
| D4 (estatísticas de interesse) | p-valor do teste max (largest P de site e taxa shuffle); taxa de rejeição a α = 0,05; para clados: acerto, posto do nó e Jaccard (PLANO §5); fração de variação no PCPS1–2 | PARCIAL | O resumo mostra só o acerto. Relatar também o posto e o Jaccard, ou justificar por que ficaram de fora. Dizer explicitamente que a resposta é o teste max e não o taxa shuffle isolado (ressalva já registrada para H1/H2) |
| D5, W-d (MCSE) | O PLANO §5 e o resumo dizem que o IC95% de Monte Carlo em 0,05 com 1000 réplicas é ±0,014 (confere: 1,96·√(0,05·0,95/1000) ≈ 0,0135) | **LACUNA** nos Resultados | Pôr o MCSE (ou IC de Monte Carlo) em **todas** as taxas das tabelas. No poder, o MCSE é maior: no pior caso, p = 0,5, IC ≈ ±0,031 com 1000 réplicas (cálculo meu). Nos contrastes do `brms` já há IC posterior; dizer que eles incorporam a variação de Monte Carlo |
| D5, D6 (justificar K) | K = 1000 réplicas por condição (camadas A e B, clados); a justificativa implícita é o ±0,014 em 0,05 e a comparação com Duarte et al. 2016 (1000 réplicas) | PARCIAL | Uma frase nos Métodos justificando K pela precisão desejada na taxa nominal e pela comparabilidade. No SI, opcional: o gráfico do MCSE contra K (D6) |
| Número de permutações | 999 permutações (`nperm = 999L` em `R/pipeline.R`, camadas A e B); 499 em alguns pilotos | **LACUNA** (justificativa) | Justificar as 999 e lembrar que o p-valor por permutação tem resolução de 1/1000. Este ponto não vem dos dois guias; é recomendação minha |
| D7 (*basic validation*) | Validação contra Duarte et al. 2016 (68 de 72 células no IC; Tabela 1) e contra Peres-Neto & Kembel 2015 (PCC reimplementado; células independentes da árvore); implementações próprias conferidas com `SYNCSA`, `vegan` e `PCPS`; testes unitários (`tests/testthat`) | OK | Relatar nos Métodos (resumo) e no SI (tabelas completas), incluindo a discrepância do Bray–Curtis sem raiz no código de Duarte et al. (Appendix S3) |
| D8, D9 (métricas e definições) | Tipo I, poder e critério de Bradley [0,025; 0,075]. A referência de Bradley 1978 está como `[CITAR: Bradley 1978, não verificado]` | PARCIAL | Definir cada métrica nos Métodos. **Verificar Bradley 1978** antes de citar. Viés e cobertura (D9) só se aplicam aos modelos de segunda etapa |
| D10 (adequação do modelo de 2ª etapa) | `brms` beta-binomial, 4 cadeias × 2000 iterações; priors explícitos (normal(0; 1,5); intercepto normal(−2,9; 1,5) no tipo I; φ ~ gamma(2; 0,01)); R-hat ≤ 1,005, ESS bulk ≥ 914, 0 divergências; PPC in-sample (97–98% das condições no IC95%); sobredispersão residual no poder (13) | PARCIAL | Relatar priors, diagnósticos e PPC no SI. Declarar como limitação que o PPC é in-sample e pouco discriminante, que a sobredispersão do poder segue alta e que falta o termo cenário × N |
| D11 (desenho e poder) | Fatores N, cenário, forma da árvore e tipo de dado em grade fatorial completa (não *space-filling*) | OK | Dizer que o desenho é fatorial completo. Na camada B, que só 5 células da grade foram usadas (4 cantos + centro) |
| D12 (definição de poder) | Taxa de rejeição a α = 0,05 nos cenários S3a/S3b/S3sigma | OK | Definir o poder como taxa de rejeição sob H1 |
| D13 (*misspecification* e *baseline*) | S2 (sem sinal em σ/h) é o *baseline* de S2sigma e S2h; os contrastes S2sigma − S2 e S2h − S2 estão no modelo formal | OK | Apresentar S2sigma/S2h explicitamente como cenários de *nuisance structure* contra o *baseline* S2 |
| D14 (comparação de métodos nos mesmos dados) | Desenho pareado: todos os testes rodam na mesma metacomunidade; nos clados, o mesmo clado, ótimos e gradiente entre deslocamentos e tipos de dado | PARCIAL | Relatar o pareamento. **Lacuna:** evoPCA e phylofactorization estão só no objetivo 2 (clados), não no tipo I nem no poder; o resumo marca isso como pendente. Ou incluir, ou delimitar o objetivo 2 |
| D15 (custo computacional) | Tabela de rodadas e durações (resumo, "Infraestrutura e custo"); AWS c7a.48xlarge, 190 workers, Docker rocker 4.6.1, OpenBLAS com 1 thread | OK | Pôr no SI: tempo por rodada, hardware e paralelização |
| *Seeds* (planejamento/código; item não recuperado de Williams) | `targets` com `seed = 20261018` e uma semente por branch (`_targets.R`); clados: `base_seed = 20261019` com sementes por (árvore, profundidade, réplica) (`analysis/17_run_clade.R`); `brms` `seed = 20261019`; pilotos com `furrr` com sementes fixas | OK no código, **LACUNA** no texto | Declarar a estratégia de sementes no SI e no README do compêndio |
| Versões de software | `renv.lock`: R 4.6.1, 240 pacotes (por exemplo brms 2.23.0, rstan 2.32.7, vegan 2.7-6, PCPS 1.0.9, SYNCSA 1.3.5, ape 5.8-1, adiv 2.2.1, phylofactor 0.0.1, mgcv 1.9-4, targets 1.12.0, crew 1.3.3); Dockerfile | OK no repositório, **LACUNA** no texto | Citar o R e os pacotes principais com versão nos Métodos; apontar o `renv.lock` e a imagem Docker. Confirmar o DOI canônico do PCPS (pendente na memória) |
| W-e, D16 (código e dados) | GitHub público (MIT/CC BY 4.0; `CITATION.cff` v0.1.0). Os `results.csv` (milhões de linhas) e os modelos ficam **fora do git**, só no S3 privado. A árvore de Portik et al. 2023 fica fora do git (licença) | **LACUNA** | Arquivar uma versão congelada no Zenodo (DOI) com o código e os resultados brutos por réplica, ou ao menos as taxas agregadas por condição, para os Resultados serem reprodutíveis sem rodar de novo. Declarar como obter a árvore de Portik (OSF; SHA-256) |
| W-f (registro e desvios) | Não há pré-registro formal. Há um plano datado e versionado (`PLANO_SIMULACAO.md` v0.2, 2026-10-16) com desvios registrados (§2b; opção B do γ na camada B; GAM exploratório → `brms` final) | PARCIAL | Declarar que o estudo não foi registrado e relatar os desvios do plano, com o histórico do git. Atenção: o histórico foi reescrito em 2026-10-17 (`git filter-repo`), então os hashes anteriores não valem |
| Relato completo, sem seleção | Todos os testes e cenários estão nas tabelas do resumo, inclusive os desfavoráveis (GLM PCPS2, Procrustes conservador, H2 rejeitada) | OK | Manter tudo no SI: tabelas completas por célula × N × tipo de dado × teste, inclusive as falhas de controle do erro tipo I |
| Falhas e réplicas perdidas | `error = "continue"` no `targets`; jobs de calibração não convergidos (66 de 72); clados "sem NA" | PARCIAL | Relatar o número de réplicas com falha ou NA por rodada (inclusive zero), e os critérios de exclusão. Recomendação minha; o item correspondente de Williams não foi recuperado |
| Análises exploratórias × confirmatórias | GAM da camada B, explicitamente exploratório (704 modelos, sem correção para multiplicidade); `brms` como modelo final | PARCIAL | Separar no texto. A camada B ainda não tem modelo formal (`analysis/20_formal_models_layerB.R` não está commitado) |

### 3.1 Principais lacunas, por prioridade

1. **MCSE em todas as taxas relatadas** (D5; Williams: 17 de 98 artigos o
   fizeram), em especial no poder, onde o IC de Monte Carlo chega a ±0,03.
2. **Arquivo permanente de código e resultados** (Zenodo, com DOI). Hoje os
   resultados ficam só num S3 privado.
3. **Justificar K = 1000 e as 999 permutações** nos Métodos.
4. **Sementes e versões de software** declaradas no texto ou no SI. No código,
   já estão resolvidas.
5. **Ausência de registro e desvios do plano**, declarados explicitamente.
6. **Escopo do objetivo 2:** evoPCA e phylofactor só na recuperação de clados.
7. **Referências pendentes de verificação:** Bradley 1978; DOI do pacote PCPS.
8. **Obter o PDF de Williams et al. 2024** e completar a seção 2.2, para mapear
   os 11 itens um a um.
