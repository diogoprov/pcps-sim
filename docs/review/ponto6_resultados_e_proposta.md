# Ponto 6: resultados e proposta de texto (2026-10-22)

O código e as saídas estão no commit ae0f1f3 (`analysis/39_clade_sensitivity.R`,
`data/derived/models/clade_sensitivity.csv` e `clade_sensitivity_ranking.csv`).
A rodada teve 24 condições × 1.000 réplicas e levou 35 min com 4 workers. Ela
**reproduziu exatamente a rodada principal** (acerto, rank e Jaccard, 72.000
linhas). Nada foi alterado no manuscrito.

## 1. Resultados

**Número de eixos: não muda o ranking.**

| Método | 2 eixos | 4 eixos (atual) | 8 eixos |
|---|---:|---:|---:|
| evoPCA, deslocamento 2 | 0,739 | 0,739 | 0,739 |
| PCPS, deslocamento 2 | 0,234 | 0,234 | 0,234 |
| PCC, deslocamento 2 | 0,069 | 0,106 | 0,118 |

(médias das 6 combinações de célula × profundidade)

- O evoPCA escolheu **sempre o eixo 1**: conferi à parte e o eixo mais
  correlacionado com E é sempre o primeiro. Por isso, 2, 4 e 8 eixos dão resultados
  idênticos.
- No PCPS, a diferença máxima foi de 0,001.
- O PCC perde 0,02–0,06 com 2 eixos e ganha até 0,02 com 8. Os intervalos pareados
  excluem zero, mas o PCC continua muito abaixo dos outros métodos.

**Objetivo da phylofactorization: o atual dá falsos positivos.**

| Cenário | Objetivo completo (atual) | Só interação com E |
|---|---:|---:|
| Controle de prevalência (falso positivo) | **0,069–0,116** | 0,001–0,004 |
| Deslocamento 0 (acaso) | 0,005–0,010 | 0,003–0,005 |
| Deslocamento 1 | 0,296 (média) | 0,322 (+0,010 a +0,042) |
| Deslocamento 2 | 0,754 (média) | 0,760 (−0,010 a +0,017) |

- No controle de prevalência, o objetivo completo "acha" o clado cerca de **dez
  vezes** mais que o acaso.
- O objetivo só de interação elimina esses falsos positivos e não perde nada:
  recupera tanto quanto o completo com deslocamento 2 e um pouco mais com
  deslocamento 1.
- Os outros métodos ficam no nível do acaso no controle: evoPCA 0,005, PCPS 0,026
  e PCC 0,029.

**Ranking entre métodos** (`clade_sensitivity_ranking.csv`):
- **Deslocamento 2:** evoPCA e phylofactorization são os dois melhores nas 36
  combinações de célula × profundidade × eixos × objetivo, e a phylofactorization
  fica em primeiro em 30 delas.
- **Deslocamento 1:** na célula de γ baixo, o PCPS passa a phylofactorization em 9
  combinações. Isso é coerente com o ganho do PCPS com γ baixo que já está no
  texto. O evoPCA fica em primeiro em todas as 36.

**O controle funcionou como planejado.** A prevalência foi de 0,32 dentro do clado
contra 0,47 fora, parecida com a diferença do deslocamento 2 (0,35 contra 0,47).
Nenhuma réplica perdeu o clado (n = 1.000 em todas as condições), ao contrário do
que um dos revisores temia.

**Consequência para o artigo.** A recomendação de usar evoPCA ou
phylofactorization se mantém, mas a phylofactorization precisa do objetivo
restrito à interação com o gradiente.

## 2. Texto proposto para o manuscrito

O texto fica em **≈ 7.982 palavras**, com cortes compensatórios C1–C4: duas frases
de apoio e duas repetições. Phillimore & Price (2008) e Mazel et al. (2016)
continuam citados em outros pontos.

### M1 Methods: sensitivity analysis (replaces the rank/Jaccard sentence, which moves to Appendix S1) (+17 palavras)

**Atual:**

> A method recovered the clade when its best-scoring clade was the planted one or its complement. We also recorded the rank of the planted node among all nodes and the Jaccard similarity between the best-scoring and the planted clades.

**Proposto:**

> A method recovered the clade when its best-scoring clade was the planted one or its complement. At N = 100, we also chose the axis among the first two or eight, restricted the objective of phylofactorization to the interaction between the split and $E$, and added a planted clade that differed only in prevalence (Appendix S1).

### R1 Results: rank/Jaccard clause moves to Appendix S2 (-10 palavras)

**Atual:**

> , and ranks and Jaccard similarities agreed with this ranking (Appendix S2)

**Proposto:**

> (Appendix S2)

### R2 Results: sensitivity (+59 palavras)

**Atual:**

> and H4 was supported for clade recovery through the change in rank between PCPS and PCC.

**Proposto:**

> and H4 was supported for clade recovery through the change in rank between PCPS and PCC. The number of axes did not change the recovery by evoPCA or PCPS and changed that of the PCC by at most 0.06. However, the default objective of phylofactorization recovered a clade that differed only in prevalence in 0.07–0.12 of the replicates, against 0.001–0.004 when restricted to the interaction with $E$, which recovered planted clades as often (Appendix S2).

### D1 Discussion: branch-wise methods (+15 palavras)

**Atual:**

> Thus, branch-wise methods seem to be the appropriate choice to identify which clades respond to a gradient.

**Proposto:**

> Thus, branch-wise methods seem to be the appropriate choice to identify which clades respond to a gradient, provided that the objective of phylofactorization does not reward differences in prevalence unrelated to it.

### D2 Recommendation (v) (+10 palavras)

**Atual:**

> (v) use evoPCA or phylofactorization, rather than PCPS axes, to identify the clades associated with a gradient.

**Proposto:**

> (v) use evoPCA, or phylofactorization with an objective restricted to the interaction with the gradient, rather than PCPS axes, to identify the clades associated with a gradient.

### C1 compensation: pool size (Phillimore & Price 2008 stays cited elsewhere) (-39 palavras)

**Atual:**

> Interestingly, γ became more negative as the pool widened from the Atlantic Forest to South America. However, γ also tends to decrease with the number of species (Phillimore & Price, 2008), so this trend may partly reflect pool size.

**Proposto:**

*(removido)*

### C2 compensation: Mazel et al. 2016 (stays cited in the Introduction) (-20 palavras)

**Atual:**

> Similarly, the imbalance of the pool tree did not change how much community phylogenetic metrics disagreed (Mazel et al., 2016).

**Proposto:**

*(removido)*

### C3 compensation: empirical reference (repeats the end of the Introduction) (-23 palavras)

**Atual:**

> To place real regional species pools on the grid of tree shapes and to illustrate the methods with real data, we computed γ and the normalized Colless index of three nested regional pools of anurans built from 14 tadpole metacommunities, and applied all tests and clade-scoring methods to three published tadpole metacommunities (Appendix S3).

**Proposto:**

> We computed γ and the normalized Colless index of three nested regional pools of anurans built from 14 tadpole metacommunities, and applied all methods to three published tadpole metacommunities (Appendix S3).

### C4 compensation: Introduction (Phillimore & Price 2008 stays cited in the Discussion) (-10 palavras)

**Atual:**

> in most clades of a large compilation of phylogenies (McPeek, 2008) and in many clades of birds (Phillimore & Price, 2008)

**Proposto:**

> in most clades of a large compilation of phylogenies (McPeek, 2008)

## 3. Texto proposto para o Apêndice S1 (no fim de "Node scores of the clade-recovery methods")

> We also recorded the rank of the planted node among all nodes and the Jaccard similarity between the best-scoring and the planted clades. In a sensitivity analysis (script 39_clade_sensitivity.R), we reanalysed the metacommunities of the main run at N = 100 with incidence data, in the Yule, low-γ and high-Colless cells and with shifts of 0, 1 and 2, using the same trees and seeds, so that the default settings reproduced the main run exactly. The three ordinations chose the axis most correlated with E among the first two, four or eight, and phylofactorization used either its default objective (the sequential deviance of the split and of its interaction with E) or the deviance of the interaction alone. As a control, we planted a clade without a shift whose maximum abundance was multiplied by 0.2, so that it differed from the rest in prevalence (0.32 against 0.47 on average, similar to the difference produced by a shift of 2, 0.35 against 0.47) but had no association with E; a hit in this control is a false recovery. All variants were computed on the same metacommunity, so their contrasts are paired by replicate.

## 4. Texto proposto para o Apêndice S2 (nova seção após "Clade recovery with abundance data")

> ## Sensitivity of clade recovery
>
> The number of axes did not change the recovery by evoPCA, whose first axis was always the one most correlated with E, and changed that of PCPS by at most 0.001 (file clade_sensitivity.csv). The PCC recovered fewer clades with two axes (−0.02 to −0.06 with a shift of 2; paired 95% Monte Carlo intervals excluded zero) and slightly more with eight (up to +0.02), but remained far below evoPCA and phylofactorization. In the prevalence control, the default objective of phylofactorization recovered the planted clade in 0.069–0.116 of the replicates, about ten times its rate without a planted effect (0.005–0.010), whereas the interaction-only objective recovered it in 0.001–0.004, and the other methods stayed at or below their chance levels (evoPCA 0.005 on average, PCPS 0.026 and PCC 0.029). With a planted shift, the interaction-only objective recovered the clade as often as the default one with a shift of 2 (difference −0.010 to +0.017) and slightly more often with a shift of 1 (+0.010 to +0.042). With a shift of 2, evoPCA and phylofactorization were the two best methods in all 36 combinations of cell, depth, number of axes and objective, and phylofactorization was first in 30 of them; with a shift of 1, evoPCA was first in all combinations, and PCPS overtook phylofactorization in the low-γ cell, as expected from its gain with low γ (file clade_sensitivity_ranking.csv). Thus, the ranking of the methods does not depend on these choices, but the default objective of phylofactorization also responds to differences in prevalence unrelated to the gradient.

## 5. Decisão pendente

Aplico tudo (texto principal, S1 e S2), re-renderizo e faço o commit?
