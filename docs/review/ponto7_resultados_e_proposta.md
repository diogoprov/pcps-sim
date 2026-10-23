# Ponto 7: resultados e proposta de texto (2026-10-22)

Script: `analysis/38_raw_contrasts.R` (commit 1884ccb), com saídas em
`data/derived/models/raw_contrasts.csv` e `model_discrepancy.csv`. O script recalcula
os 1.534 contrastes dos modelos a partir das taxas brutas, com intervalos de Monte
Carlo que respeitam o desenho:
- **não pareados** entre células do grid e entre arquétipos;
- **pareados por árvore** em S2σ/S2h × S2.

Os erros-padrão usam a estimativa de Jeffreys. Nada foi alterado no manuscrito.

## 1. O que a checagem mostrou

**As conclusões se mantêm.** Conferi uma a uma as afirmações do texto contra os
contrastes brutos:

| Afirmação no texto | Modelo | Bruto (IC de Monte Carlo 95%) | Situação |
|---|---|---|---|
| H1: desbalanço não altera o erro tipo I | máx. ±0,012; nenhum IC exclui 0 | máx. ±0,017; 4 de 150 excluem 0 (o acaso daria ~7,5) | mantém |
| γ baixo, GLM PCPS1, S3b, N = 100 | −0,40 | −0,42 (−0,46 a −0,38) | mantém |
| γ baixo, PCC | −0,25 | −0,27 (−0,30 a −0,24) | mantém |
| γ baixo, PERMANOVA "não mudou" | 0,01 (−0,02 a 0,05) | **+0,05 (0,02 a 0,08)** | o poder até aumentou: "did not decrease" |
| Desbalanço reduz o poder "no máximo 0,09" | 0,09 | **0,13** (Mantel, N = 50) | trocar pelo bruto |
| S2σ, GLM PCPS2, γ baixo, N = 50/100/200 | 0,05 / 0,10 / 0,21 | 0,04 / 0,11 / **0,30** | o modelo encolhe o maior efeito |
| S2σ, GLM PCPS2, Yule, N = 100 | 0,02 | 0,03 (0,01–0,05) | mantém |
| S2σ, PCC, Yule, N = 100 | 0,05 | **0,08** (0,06–0,10) | o modelo encolhe |
| S2σ, PCC, Yule, N = 200 (não citado) | 0,07 | **0,16** | o modelo encolhe |
| Mass effects, Yule pool, S3b | −0,55 a −0,93 | −0,56 a −0,97 | mantém |
| Competição: PERMANOVA perde 0,46 em vez de 0,62 | 0,46 / 0,62 | 0,46 / 0,57 | mantém |
| Clados: PCPS ganha com γ baixo; evoPCA e phylofactor não mudam | (16 de 16 concordam) | idem | mantém |

Em 1.307 dos 1.534 contrastes, os intervalos bruto e do modelo concordam quanto a
excluir zero e quanto ao sinal. Dos 227 que discordam, 178 são diferenças pequenas
de erro tipo I, sobretudo entre arquétipos da camada de processos, que o texto não
interpreta.

**O ponto fraco são as condições individuais, não os contrastes.**
- As medianas do beta-binomial se afastam da taxa bruta em mais de 0,02 em 699 das
  5.664 condições, com máximo de 0,13.
- Em 574 dessas condições, a mediana cai fora do IC bruto de 95%. Nos dados de
  poder isso acontece em cerca de um terço das condições (268 de 810); na
  recuperação de clados, em cerca de 18% (152 de 864).
- O modelo com efeito aleatório por condição reproduz as taxas brutas.
- Ou seja, a estrutura fixa do beta-binomial resume bem os efeitos de formato, mas
  não cada célula.

**Isso afeta as figuras.** As Figs. 2, 3 e 5 mostram em cada célula a mediana
posterior do beta-binomial (N = 100, incidência). Na Fig. 3, 45 células de poder
ficam fora do IC bruto, com diferença de até 0,09; o pior caso é o Mantel em S3a,
com 0,61 no modelo contra 0,52 bruto.

**Uma frase dos Métodos ficou enganosa:**
> The beta-binomial accounts for the dispersion of the rates, which was 2.1–8.5
> times larger than binomial in quasibinomial models with the same structure.

O próprio cabeçalho do `35_mean_structure.R` diz que essa "dispersão" é falta de
ajuste das médias, porque cada linha agregada é uma condição. Proponho cortar a
frase (C3).

## 2. Decisão pendente: as figuras

- **(A) Recomendado:** as Figs. 2, 3 e 5 passam a mostrar as **taxas brutas** em
  cada célula, e a Fig. 4 mostra a taxa bruta com o IC de Monte Carlo. As legendas
  trocam "posterior medians" por "raw rates". Os modelos continuam nos contrastes
  do texto. É coerente com T1 e é o que o Codex pediu.
- **(B)** Manter as medianas do modelo e só acrescentar a Tabela S2 de taxas brutas
  e a frase T2. Dá menos trabalho, mas deixa células com valores fora do IC bruto.

## 3. Texto proposto para o manuscrito

Com T1–T5 e as compensações C1–C3, o texto fica em **≈ 7.981 palavras**. Na opção
(A), as legendas das figuras mudam, mas legendas não entram na contagem.

### T1 Methods: raw rates as primary evidence (+13 palavras)

**Atual:**

> We summarize the results as posterior rates and contrasts between tree shapes or archetypes, with medians and 95% credible intervals (CI), and report the raw rates where they exceed the modelled ones.

**Proposto:**

> We treat the raw rates and their Monte Carlo contrasts, which respect the pairing of the design, as the primary evidence and the posterior medians and 95% credible intervals (CI) of the models as a synthesis, and report the raw contrasts where the two differ.

### T2 Results: model fit (+18 palavras)

**Atual:**

> The predicted rates were largely insensitive to the priors and to the lack of fit of the fixed structure (Appendix S2).

**Proposto:**

> The predicted rates were largely insensitive to the priors; the beta-binomial medians of single conditions departed from the raw rates by up to 0.13, but the model contrasts agreed with the raw ones in every conclusion below (Appendix S2).

### T3 Results: H1 (+8 palavras)

**Atual:**

> was within ±0.012, and no 95% CI excluded zero, with incidence or abundance data. The posterior probability of a change smaller than 0.025 exceeded 0.95 in 297 of 300 conditions, and the taxa shuffle alone was not affected either (Appendix S2).

**Proposto:**

> was within ±0.012 (±0.017 in raw rates), and no 95% CI excluded zero, with incidence or abundance data; four of 150 raw Monte Carlo intervals did, fewer than expected by chance. The posterior probability of a change smaller than 0.025 exceeded 0.95 in 297 of 300 conditions (Appendix S2).

### T4 Results: S2σ in raw paired contrasts (+6 palavras)

**Atual:**

> Compared with the same scenario without signal (S2), signal in σ (S2σ) increased the rejection rate of GLM PCPS2 by 0.02 (95% CI 0.01–0.03) in Yule trees at N = 100, and by 0.05 (0.04–0.06), 0.10 (0.08–0.12) and 0.21 (0.19–0.23) with low γ at N = 50, 100 and 200, with raw rates of up to 0.36. The PCC increased by 0.05 (0.03–0.06) in Yule trees and by 0.02 (0.01–0.03) with low γ at N = 100, with raw rates of up to 0.22.

**Proposto:**

> Compared with the same scenario without signal (S2), and paired by tree, signal in σ (S2σ) increased the rejection rate of GLM PCPS2 by 0.03 (95% Monte Carlo interval 0.01–0.05) in Yule trees at N = 100, and by 0.04 (0.01–0.06), 0.11 (0.09–0.14) and 0.30 (0.27–0.33) with low γ at N = 50, 100 and 200, with raw rates of up to 0.36. The PCC increased by 0.08 (0.06–0.10) in Yule trees and by 0.04 (0.02–0.06) with low γ at N = 100, with raw rates of up to 0.22.

### T5 Results: power, PERMANOVA and imbalance (+6 palavras)

**Atual:**

> whereas the power of PERMANOVA did not change (0.01, −0.02 to 0.05), and high imbalance reduced power by at most 0.09.

**Proposto:**

> whereas the power of PERMANOVA did not decrease (raw contrast 0.05, 95% Monte Carlo interval 0.02–0.08), and high imbalance reduced power by at most 0.13 (raw rates).

### C1 compensation: priorsense detail (-10 palavras)

**Atual:**

> for the coefficients and for the predicted rate of every condition.

**Proposto:**

> .

### C3 Methods: dispersion sentence (lack of fit of the means) (-23 palavras)

**Atual:**

> The beta-binomial accounts for the dispersion of the rates, which was 2.1–8.5 times larger than binomial in quasibinomial models with the same structure.

**Proposto:**

> 

### C2 compensation: performance measure sentence (-2 palavras)

**Atual:**

> The performance measure of each test was its rejection rate across the $K = 1{,}000$ replicates of a condition, and that of each clade-recovery method was its proportion of hits.

**Proposto:**

> The performance measures were the rejection rate of each test and the proportion of hits of each clade-recovery method across the $K = 1{,}000$ replicates of a condition.

## 4. Parágrafo proposto para o Apêndice S2 (novo subtítulo)

> ## Raw contrasts and fit of the formal models
>
> For every contrast of the formal models, we computed the raw contrast with a 95% Monte Carlo interval (script 38_raw_contrasts.R; file raw_contrasts.csv), unpaired between cells of the grid and between archetypes, whose replicates used different trees or communities, and paired by tree between S2σ or S2h and S2. Standard errors used the Jeffreys estimate of each rate, and the contrast on the distance from 0.05 used a parametric bootstrap, which is approximate near 0.05. The raw and model intervals agreed on whether and in which direction a contrast differed from zero in 1,307 of 1,534 contrasts, and 178 of the 227 disagreements were small differences in type I error, mostly between archetypes, which we do not interpret. However, the beta-binomial models shrank the largest increases under signal in σ (GLM PCPS2 at N = 200 with low γ: raw 0.30, 95% Monte Carlo interval 0.27–0.33, against a posterior median of 0.21; PCC at N = 200 in Yule trees: 0.16 against 0.07). Their medians for single conditions departed from the raw rates by more than 0.02 in 699 of 5,664 conditions, by up to 0.13, and fell outside the raw 95% interval in 574 of them, mostly in power and clade recovery (file model_discrepancy.csv), whereas a model with a random effect per condition reproduced the raw rates. Thus, the fixed structure of the beta-binomial models summarizes the effects of tree shape and assembly, but not every condition, and the raw rates of Table S2 should be preferred for single conditions.

(A Tabela S2 já traz todas as taxas brutas com o EP de Monte Carlo.)
