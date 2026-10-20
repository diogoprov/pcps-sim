# Guias da MEE para relatar estudos de simulação, e como aplicá-los ao pcps-sim

Versão 0.2, 2026-10-19 (seção 2 completada com o PDF de Williams et al. 2024). Documento interno de trabalho, em português. As citações
literais ficam em inglês, como no original.

## 0. Fontes e o que foi realmente lido

| Guia | Fonte usada | Cobertura |
|---|---|---|
| DiRenzo, Hanks & Miller | PDF publicado (versão *Early View*, pp. 1–15, baixado pela UFMS em 06/12/2022) em `~/Documents/Ecologia-Converted2.Data/PDF/4146382399/DiRenzo-2022.pdf`; texto extraído com Ghostscript | **texto completo lido** (Resumo, seções 1–5, Box 1, Box 2, Tabela 1, Data availability) |
| Williams et al. | PDF publicado (versão final, fascículo 15(11), pp. 1926–1939) em `references/Williams_2024_MEE_transparent_reporting.pdf`; texto extraído com `pdftools::pdf_text`; a Figura 2 (imagem, com os 11 itens) foi lida pela página renderizada (p. 1932) | **texto completo lido** (Resumo, seções 1–5, Box 1, Figura 1, Figura 2). Supporting Information (S1, Tabelas S1–S2, Figuras S1–S6) **não lida** |

Tentativas que falharam antes de o PDF de Williams et al. ser obtido (registro histórico; o problema foi resolvido com o PDF salvo em `references/`):

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
- O Williams et al. 2024 é de acesso aberto (CC BY-NC-ND; confirmado na
  primeira página do PDF).

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

Tipo: *Perspective*. Licença: CC BY-NC-ND. Páginas citadas abaixo são as do
fascículo (1926–1939). Os autores citam DiRenzo et al. como "2022" (data online);
a referência final é DiRenzo, Hanks & Miller (2023), *MEE* 14(1): 203–217.

### 2.1 Pontos gerais (verificados no texto completo)

| # | Ponto | Citação literal | Onde (no artigo) |
|---|---|---|---|
| W-a | São **11 itens de relato**, que estendem Morris et al. (2019, *Stat Med* 38: 2074) | "we propose 11 reporting items for statistical simulation studies" | Resumo |
| W-b | Os itens cobrem **três estágios**; o relato deixa de ser um estágio separado (como em Morris et al.) e passa a atravessar os três | "These items span across three stages: planning, coding and analysis." (Resumo); "our proposed items integrate reporting components within each main stage" (p. 1930) | Resumo; §3.1; Figura 1 |
| W-c | O estágio de planejamento **segue o ADEMP** (Morris et al. 2019). Verificado: os itens 1–5 são A, D, E, M, P | "Planning (following ADEMP)" (Figura 2) | Figura 2 |
| W-d | Itens **novos** em relação a Morris et al.: 6 (*Implementation*), 8 (*Code for performance measure estimates*) e 9 (*Case Study*), marcados com asterisco | "we introduce three steps: 'Implementations', 'Code for performance measure estimates' and 'Case Study'" (p. 1930) | §3.1; Figura 2 |
| W-e | Os itens 1–5 vão no texto principal, não no SI | "These five core steps are typically reported in the Introduction or the Methods section of a manuscript and not in the supporting materials or information." (p. 1930) | §3.1.1 |
| W-f | Levantamento de 100 artigos (2017–2022): só 17% relataram incerteza de Monte Carlo e 1% justificou o número de réplicas | "only a small proportion of articles reported a Monte Carlo uncertainty (17%; 17 out of 98)" (Resumo); "only one study (1%; 1 out of 93) in our survey provided justification" (p. 1933) | Resumo; Box 1; §3.2 |
| W-g | 32% não deram código | "32% (32 out of 100) articles did not provide code" (Resumo) | Resumo; Box 1 |
| W-h | **Sementes:** relatar a semente aleatória, se definida | "if a random seed is defined in the code for data generation, reporting it enhances reproducibility" (p. 1931) | §3.1.2 |
| W-i | **Versões:** implementações e versões de pacotes podem mudar resultados | "default settings and versioning of software and packages across different software and programming languages can result in minor differences" (p. 1931) | §3.1.2 |
| W-j | **Relato completo:** relatar todos os resultados de desempenho planejados | "if only a selection of performance results is reported, the reader cannot fully interpret the findings" (p. 1933) | §3.1.3 |
| W-k | **MCSE** em gráfico ou tabela, no material suplementar | "we also strongly encourage the reporting of Monte Carlo standard errors (MCSE) in simulation studies either graphically or in tabular format in the supporting materials" (p. 1933) | §3.1.3 |
| W-l | **Justificativa do número de réplicas** a priori, não *post hoc* | "we propose that a more nuanced understanding of the simulation results variability can be achieved by providing background to the choice of the number of simulation replications" (p. 1933) | §3.1.3 |
| W-m | **Registro** do estudo de simulação (o termo é *registration*, não *pre-registration*); citam o modelo ADEMP de registro de Siepe et al. (2023); desvios são esperados e devem ser declarados; resultados de piloto podem entrar no registro | "These deviations and additions can be added and should be acknowledged." (p. 1934); "We suggest to report any pilot results within a registration" (p. 1934) | §4.2 |
| W-n | O registro é sobretudo importante ao comparar método novo com existentes (viés de relato seletivo) | "registration is particularly valuable when comparing a new method or procedure to existing ones" (p. 1934) | §4.2 |

### 2.2 Os 11 itens (Figura 2, p. 1932)

Texto da coluna *Reporting item* e da coluna *Expected reporting section*,
transcrito da Figura 2. Asterisco = item novo em relação a Morris et al. (2019).

| # | Estágio | Passo | *Reporting item* (literal) | *Expected reporting section* |
|---|---|---|---|---|
| W1 | Planning | (A) Aims | "Clearly define the objectives and hypotheses your simulation is designed to explore or test." | Introduction, Methods |
| W2 | Planning | (D) Data-generating mechanisms | "Describe the processes used to generate simulated data, referencing underlying models or theories." | Methods |
| W3 | Planning | (E) Estimand (statistical target) | "Identify and describe the specific statistical target or task that your simulation study is set up to estimate or test." | Methods |
| W4 | Planning | (M) Methods | "Outline the simulation methods employed, including experimental design and analysis techniques." | Introduction, Methods |
| W5 | Planning | (P) Performance measures | "List all metrics used to assess the method's performance, with a rationale for their selection." | Methods, Results |
| W6 | Coding | Implementation* | "Report all software and/or programming language used for simulations, including versions and relevant citations." | Methods, Results, Supplementary |
| W7 | Coding | Code for simulated data analysis | "Provide the code or a reference to where the code can be accessed for replicating the simulation study steps." | Methods, Results, Supplementary |
| W8 | Coding | Code for performance measure estimates* | "Make the code available or cite its location for calculating performance measure estimates." | Methods, Results, Supplementary |
| W9 | Analysis | Case Study* | "Summarise the real-world scenario that the simulation study exemplifies or is applicable to." | Results, Supplementary |
| W10 | Analysis | Estimates of performance | "Report the calculated performance metrics for the statistical methodologies evaluated (e.g. graphical or tabular)." | Results, Supplementary |
| W11 | Analysis | Monte Carlo standard errors (SEs) of performance measures | "Present the Monte Carlo standard errors or a justification of the number of simulation runs." | Results, Supplementary |

Observações:

- O artigo não tem uma categoria de *code repository* na coluna de seção. O
  código (W7, W8) é posto em "Methods, Results, Supplementary", e os exemplos da
  Figura 2 apontam para repositórios com DOI (Dryad, Appendix S1). Repositório
  com DOI é, portanto, a forma de cumprir W7–W8; a escolha do repositório é
  minha recomendação.
- W11 tem a forma "MCSE **ou** justificativa do número de réplicas". O texto
  (p. 1933) recomenda as duas coisas (W-k, W-l).
- Os itens são apresentados como protótipo, a ser refinado por consenso da
  comunidade (§5, p. 1934); não é ainda um *guideline* formal tipo CONSORT ou
  PRISMA.
- As 24 perguntas do levantamento (Tabela S2) estão no Supporting Information,
  que **não li**.

---

## 3. Como aplicar ao pcps-sim

Fontes do projeto: `docs/PLANO_SIMULACAO.md` (v0.2), `docs/resumo_projeto.qmd`
(2026-10-19), `docs/MEMORIA_PROJETO.md`, `_targets.R`, `R/pipeline.R`,
`analysis/17_run_clade.R`, `analysis/18_formal_models.R` e `renv.lock`.
Legenda: **OK** = o projeto já tem e é só relatar; **PARCIAL** = existe, mas
falta algo; **LACUNA** = não existe ou não está documentado.

### 3.1 Os 11 itens de Williams et al., um a um

| Item | Onde o artigo pede | O que o pcps-sim tem | Estado | O que fazer no manuscrito |
|---|---|---|---|---|
| W1 Aims | Introduction, Methods | Objetivo 1 (validade de ADONIS/Mantel/GLM-PCPS/RDA/PCC: tipo I e poder) e objetivo 2 (recuperação de clado plantado com PCPS/PCC/evoPCA/phylofactor); hipóteses H1–H4 (PLANO §1) | OK | Objetivos e hipóteses no fim da Introdução e no início dos Métodos |
| W2 Data-generating mechanisms | Methods | Camada A estática (Peres-Neto, Leibold & Dray 2012); camada B por processo (arquétipos de Thompson et al. 2020, reimplementados); 8 cenários; gerador de árvores com grade 3 × 3 γ × Colless calibrada em clados de anuros de Portik et al. 2023; N = 50/100/200 | OK | Citar os modelos de origem (o item pede "referencing underlying models or theories"). Justificar a escolha (Williams: 25% não justificaram). Declarar que o gerador de árvores é desenho próprio |
| W3 Estimand / statistical target | Methods | Objetivo 1: tarefa = teste de hipótese, alvo = H0 de ausência de relação filogenética-ambiental (por cenário). Objetivo 2: alvo = o nó do clado plantado. 2ª etapa: parâmetros dos modelos `brms` | PARCIAL | Nomear explicitamente a tarefa e o alvo de cada objetivo com o vocabulário do artigo (*task*/*target*). Dizer qual H0 cada teste avalia, e que a resposta é o teste max |
| W4 Methods | Introduction, Methods | Cinco testes de hipótese + quatro métodos de recuperação, desenho fatorial completo e pareado (mesma metacomunidade para todos os testes) | PARCIAL | Descrever o desenho (fatores, níveis, pareamento). Delimitar: evoPCA e phylofactor só entram no objetivo 2 |
| W5 Performance measures | Methods, Results | Tipo I, poder, critério de Bradley; acerto do nó, posto e Jaccard (clados) | PARCIAL | Listar todas as métricas **com justificativa** (o item pede "with a rationale for their selection"). Bradley 1978 continua **não verificado** |
| W6 Implementation | Methods, Results, Supplementary | `renv.lock` (R 4.6.1; brms 2.23.0, vegan 2.7-6, PCPS 1.0.9, SYNCSA 1.3.5, phylofactor 0.0.1, targets 1.12.0 etc.); Docker rocker 4.6.1; AWS c7a.48xlarge | OK no repositório, **LACUNA** no texto | Citar R e pacotes com versão e referência nos Métodos; tabela completa no SI. Confirmar o DOI do PCPS |
| W7 Code for simulated data analysis | Methods, Results, Supplementary | GitHub público com o pipeline `targets` (simulação + testes) e `renv`; sementes fixas | PARCIAL | Arquivar versão congelada com DOI (Zenodo) e citá-la em Methods e Data availability. Declarar as sementes (W-h) |
| W8 Code for performance measure estimates | Methods, Results, Supplementary | `R/validation.R` (`rejection_rates`), `R/clade_recovery.R`, `R/formal_models.R`, `analysis/18_formal_models.R`, `analysis/20_formal_models_layerB.R` | PARCIAL | Apontar esses arquivos explicitamente. **Obs.:** o comentário de `rejection_rates()` diz "with binomial 95% CI", mas a função só calcula a média; o IC não é calculado ali |
| W9 Case Study | Results, Supplementary | **Nenhuma aplicação a dados reais.** Há só a calibração empírica das árvores (Portik et al. 2023) e, como ideia, as ~15 metacomunidades de girinos do projeto NRI vs PCPS (árvores de 9–34 espécies, abaixo do N do desenho) | **LACUNA** | Decidir: (a) aplicar os testes a uma metacomunidade real como exemplo, ou (b) declarar que não há estudo de caso e descrever o cenário real a que a simulação se aplica (metacomunidades de anuros; faixa empírica de γ e Colless). O item só exige (b) ("Summarise the real-world scenario") |
| W10 Estimates of performance | Results, Supplementary | Tabelas de taxas por condição no resumo, inclusive resultados desfavoráveis; contrastes `brms` | OK | Figuras/tabelas no texto principal; **todas** as células × N × tipo de dado × teste no SI (W-j) |
| W11 MCSE ou justificativa de K | Results, Supplementary | O ±0,014 em 0,05 com 1000 réplicas está no PLANO §5, mas as tabelas não trazem MCSE por taxa | **LACUNA** | MCSE (ou IC de Monte Carlo) em todas as taxas, em tabela ou gráfico no SI (W-k), **e** justificativa a priori de K = 1000 nos Métodos (W-l) |
| Registro (W-m, §4.2; não é item numerado) | — | Sem registro. Plano datado (`PLANO_SIMULACAO.md` v0.2) com desvios anotados; pilotos documentados | **LACUNA** (declarar) | Dizer que o estudo não foi registrado; relatar os desvios do plano e os pilotos. Opcional: registrar agora o plano da parte ainda não rodada (modelo de Siepe et al. 2023, citado por Williams; **não verifiquei** Siepe et al. 2023 diretamente) |

### 3.2 Mapeamento detalhado (DiRenzo et al. e pontos gerais de Williams et al.)

| Item | O que o projeto tem | Estado | O que fazer no manuscrito |
|---|---|---|---|
| D1, W1 (objetivo) | Objetivo 1 (validade da inferência sobre P/PCPS: tipo I e poder) e objetivo 2 (comparação de métodos, recuperação de clados); hipóteses H1–H4 (PLANO §1) | OK | Abrir a seção de Métodos com os objetivos e as hipóteses. Dizer que H1–H4 foram formuladas antes das rodadas (datas no git) |
| D2, D3, W2 (*data-generating mechanism*) | Duas camadas. **A:** estática, Peres-Neto, Leibold & Dray 2012, 100 sítios, E ~ U(0,100), resposta unimodal + Poisson. **B:** dinâmica de Thompson et al. 2020, reimplementada, com 4 arquétipos (*sorting*/*mass effects* × sem competição/competição filogenética). 8 cenários (S0, S1, S2, S2sigma, S2h, S3a, S3b, S3sigma); grade 3 × 3 γ × Colless; N = 50/100/200; incidência e abundância | OK | Tabela de cenários no texto. Fórmulas completas e o gerador de árvores (beta-splitting + taxa k^c, tolerância ±0,5) no SI. Declarar que o gerador é escolha de desenho própria, não um modelo publicado (PLANO §2b) |
| D3 (fatores fixos e excluídos) | Espaço/dispersão fora da camada A; célula *mass effects* + competição + N = 200 retirada por inviabilidade | PARCIAL | Relatar as exclusões e o motivo (teto de coexistência) nos Métodos, e o desvio do plano original (biblioteca de vários geradores → gerador de dois controles) |
| D4 (estatísticas de interesse) | p-valor do teste max (largest P de site e taxa shuffle); taxa de rejeição a α = 0,05; para clados: acerto, posto do nó e Jaccard (PLANO §5); fração de variação no PCPS1–2 | PARCIAL | O resumo mostra só o acerto. Relatar também o posto e o Jaccard, ou justificar por que ficaram de fora. Dizer explicitamente que a resposta é o teste max e não o taxa shuffle isolado (ressalva já registrada para H1/H2) |
| D5, W11, W-k (MCSE) | O PLANO §5 e o resumo dizem que o IC95% de Monte Carlo em 0,05 com 1000 réplicas é ±0,014 (confere: 1,96·√(0,05·0,95/1000) ≈ 0,0135) | **LACUNA** nos Resultados | Pôr o MCSE (ou IC de Monte Carlo) em **todas** as taxas das tabelas. No poder, o MCSE é maior: no pior caso, p = 0,5, IC ≈ ±0,031 com 1000 réplicas (cálculo meu). Nos contrastes do `brms` já há IC posterior; dizer que eles incorporam a variação de Monte Carlo |
| D5, D6, W11, W-l (justificar K) | K = 1000 réplicas por condição (camadas A e B, clados); a justificativa implícita é o ±0,014 em 0,05 e a comparação com Duarte et al. 2016 (1000 réplicas) | PARCIAL | Uma frase nos Métodos justificando K pela precisão desejada na taxa nominal e pela comparabilidade. No SI, opcional: o gráfico do MCSE contra K (D6) |
| Número de permutações | 999 permutações (`nperm = 999L` em `R/pipeline.R`, camadas A e B); 499 em alguns pilotos | **LACUNA** (justificativa) | Justificar as 999 e lembrar que o p-valor por permutação tem resolução de 1/1000. Este ponto não vem dos dois guias; é recomendação minha |
| D7 (*basic validation*) | Validação contra Duarte et al. 2016 (68 de 72 células no IC; Tabela 1) e contra Peres-Neto & Kembel 2015 (PCC reimplementado; células independentes da árvore); implementações próprias conferidas com `SYNCSA`, `vegan` e `PCPS`; testes unitários (`tests/testthat`) | OK | Relatar nos Métodos (resumo) e no SI (tabelas completas), incluindo a discrepância do Bray–Curtis sem raiz no código de Duarte et al. (Appendix S3) |
| D8, D9 (métricas e definições) | Tipo I, poder e critério de Bradley [0,025; 0,075]. A referência de Bradley 1978 está como `[CITAR: Bradley 1978, não verificado]` | PARCIAL | Definir cada métrica nos Métodos. **Verificar Bradley 1978** antes de citar. Viés e cobertura (D9) só se aplicam aos modelos de segunda etapa |
| D10 (adequação do modelo de 2ª etapa) | `brms` beta-binomial, 4 cadeias × 2000 iterações; priors explícitos (normal(0; 1,5); intercepto normal(−2,9; 1,5) no tipo I; φ ~ gamma(2; 0,01)); R-hat ≤ 1,005, ESS bulk ≥ 914, 0 divergências; PPC in-sample (97–98% das condições no IC95%); sobredispersão residual no poder (13) | PARCIAL | Relatar priors, diagnósticos e PPC no SI. Declarar como limitação que o PPC é in-sample e pouco discriminante, que a sobredispersão do poder segue alta e que falta o termo cenário × N |
| D11 (desenho e poder) | Fatores N, cenário, forma da árvore e tipo de dado em grade fatorial completa (não *space-filling*) | OK | Dizer que o desenho é fatorial completo. Na camada B, que só 5 células da grade foram usadas (4 cantos + centro) |
| D12 (definição de poder) | Taxa de rejeição a α = 0,05 nos cenários S3a/S3b/S3sigma | OK | Definir o poder como taxa de rejeição sob H1 |
| D13 (*misspecification* e *baseline*) | S2 (sem sinal em σ/h) é o *baseline* de S2sigma e S2h; os contrastes S2sigma − S2 e S2h − S2 estão no modelo formal | OK | Apresentar S2sigma/S2h explicitamente como cenários de *nuisance structure* contra o *baseline* S2 |
| D14 (comparação de métodos nos mesmos dados) | Desenho pareado: todos os testes rodam na mesma metacomunidade; nos clados, o mesmo clado, ótimos e gradiente entre deslocamentos e tipos de dado | PARCIAL | Relatar o pareamento. **Lacuna:** evoPCA e phylofactorization estão só no objetivo 2 (clados), não no tipo I nem no poder; o resumo marca isso como pendente. Ou incluir, ou delimitar o objetivo 2 |
| D15 (custo computacional) | Tabela de rodadas e durações (resumo, "Infraestrutura e custo"); AWS c7a.48xlarge, 190 workers, Docker rocker 4.6.1, OpenBLAS com 1 thread | OK | Pôr no SI: tempo por rodada, hardware e paralelização |
| *Seeds* (W-h, p. 1931; ligado a W7) | `targets` com `seed = 20261018` e uma semente por branch (`_targets.R`); clados: `base_seed = 20261019` com sementes por (árvore, profundidade, réplica) (`analysis/17_run_clade.R`); `brms` `seed = 20261019`; pilotos com `furrr` com sementes fixas | OK no código, **LACUNA** no texto | Declarar a estratégia de sementes no SI e no README do compêndio |
| Versões de software (W6, W-i) | `renv.lock`: R 4.6.1, 240 pacotes (por exemplo brms 2.23.0, rstan 2.32.7, vegan 2.7-6, PCPS 1.0.9, SYNCSA 1.3.5, ape 5.8-1, adiv 2.2.1, phylofactor 0.0.1, mgcv 1.9-4, targets 1.12.0, crew 1.3.3); Dockerfile | OK no repositório, **LACUNA** no texto | Citar o R e os pacotes principais com versão nos Métodos; apontar o `renv.lock` e a imagem Docker. Confirmar o DOI canônico do PCPS (pendente na memória) |
| W7, W8, W-g, D16 (código e dados) | GitHub público (MIT/CC BY 4.0; `CITATION.cff` v0.1.0). Os `results.csv` (milhões de linhas) e os modelos ficam **fora do git**, só no S3 privado. A árvore de Portik et al. 2023 fica fora do git (licença) | **LACUNA** | Arquivar uma versão congelada no Zenodo (DOI) com o código e os resultados brutos por réplica, ou ao menos as taxas agregadas por condição, para os Resultados serem reprodutíveis sem rodar de novo. Declarar como obter a árvore de Portik (OSF; SHA-256) |
| W-m (registro e desvios) | Não há pré-registro formal. Há um plano datado e versionado (`PLANO_SIMULACAO.md` v0.2, 2026-10-16) com desvios registrados (§2b; opção B do γ na camada B; GAM exploratório → `brms` final) | PARCIAL | Declarar que o estudo não foi registrado e relatar os desvios do plano, com o histórico do git. Atenção: o histórico foi reescrito em 2026-10-17 (`git filter-repo`), então os hashes anteriores não valem |
| Relato completo, sem seleção (W10, W-j) | Todos os testes e cenários estão nas tabelas do resumo, inclusive os desfavoráveis (GLM PCPS2, Procrustes conservador, H2 rejeitada) | OK | Manter tudo no SI: tabelas completas por célula × N × tipo de dado × teste, inclusive as falhas de controle do erro tipo I |
| Falhas e réplicas perdidas | `error = "continue"` no `targets`; jobs de calibração não convergidos (66 de 72); clados "sem NA" | PARCIAL | Relatar o número de réplicas com falha ou NA por rodada (inclusive zero), e os critérios de exclusão. Recomendação minha; Williams et al. não têm item específico sobre réplicas com falha |
| Análises exploratórias × confirmatórias | GAM da camada B, explicitamente exploratório (704 modelos, sem correção para multiplicidade); `brms` como modelo final | PARCIAL | Separar no texto. `analysis/20_formal_models_layerB.R` já está no git; conferir se o modelo formal da camada B foi ajustado e diagnosticado antes de relatar |

### 3.3 Principais lacunas, por prioridade

1. **MCSE em todas as taxas relatadas** (D5; W11; Williams: 17 de 98 artigos o
   fizeram), em especial no poder, onde o IC de Monte Carlo chega a ±0,03.
2. **Arquivo permanente de código e resultados** (W7, W8; Zenodo, com DOI). Hoje
   os resultados ficam só num S3 privado.
3. **Justificar K = 1000** (W11, W-l; só 1 de 93 artigos o fez) **e as 999
   permutações** nos Métodos.
4. **Sementes e versões de software** (W6, W-h, W-i) declaradas no texto ou no
   SI. No código, já estão resolvidas.
5. **Estudo de caso (W9):** não há aplicação a dados reais; decidir entre um
   exemplo empírico ou só descrever o cenário real a que a simulação se aplica.
6. **Ausência de registro e desvios do plano** (W-m), declarados explicitamente.
7. **Escopo do objetivo 2** (W4): evoPCA e phylofactor só na recuperação de clados.
8. **Justificativa das métricas** (W5) e alvo estatístico nomeado (W3).
9. **Referências pendentes de verificação:** Bradley 1978; DOI do pacote PCPS;
   Siepe et al. 2023 (citado por Williams, não lido).
