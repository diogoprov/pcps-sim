# Ponto 6: sensibilidade da recuperação de clados (plano, 2026-10-22)

Plano para aprovação. Nada roda sem o seu ok.

## Pergunta

A recomendação do abstract é "use evoPCA ou phylofactorization para identificar
clados". O Codex aponta duas escolhas do procedimento que podem mudar o ranking dos
métodos:

1. **Número de eixos:** as ordenações (PCPS, PCC e evoPCA) escolhem o eixo mais
   correlacionado com E entre os **quatro primeiros**.
2. **Objetivo binomial da phylofactorization:** ele soma o termo `phylo` (diferença
   de prevalência entre os lados da partição) ao termo `phylo:env`, de modo que
   **premia diferenças de prevalência sem relação com E**. O próprio comentário em
   `R/clade_methods.R` já reconhece isso.

## Desenho proposto

- **Variantes:** são calculadas sobre a **mesma metacomunidade**, então o desenho é
  pareado.
  - PCPS, PCC e evoPCA com 2, 4 e 8 eixos;
  - phylofactorization com o objetivo atual e com um objetivo que só usa a
    interação com E (`d2 − d3`). As duas versões saem dos mesmos ajustes, sem custo
    extra.
- **Cenários:**
  - deslocamento 0 (acaso);
  - deslocamento 1;
  - deslocamento 2;
  - **controle de prevalência (novo):** o clado plantado tem *h* × 0,2 e nenhum
    deslocamento dos ótimos. Isso cria uma diferença de prevalência entre o clado e
    o resto sem nenhuma associação com E. No piloto, a prevalência ficou em 0,35
    dentro do clado contra 0,50 fora, uma diferença de tamanho parecido com a que o
    deslocamento 2 produz (0,37 contra 0,46). Um "acerto" nesse cenário é um falso
    positivo.
- **Células:** Yule/Yule, γ baixo/Yule e Yule/desbalanço alto. São as três dos
  contrastes que estão no texto.
- **Profundidade:** clados basais e recentes.
- **N e dados:** N = 100, incidência (é onde está o objetivo binomial).
- **Réplicas:** 1.000 por condição. Elas reaproveitam as árvores e o esquema de
  sementes do `17_run_clade.R` (árvore, profundidade, réplica). Assim, nos cenários
  de deslocamento, a variante "4 eixos + objetivo atual" **reproduz exatamente** a
  rodada principal, o que serve de controle de que nada mudou.
- **Tamanho:** 3 células × 2 profundidades × 4 cenários = 24 condições, ou 24.000
  réplicas.

## Benchmark

O piloto rodou 5 réplicas × 2 células × 3 cenários × 2 profundidades, com as 12
variantes por réplica, em um único núcleo, com a máquina em carga 3 de 18 núcleos.

| | Segundos por réplica |
|---|---:|
| Total, as 12 variantes | 0,23 (máx. 1,6) |
| evoPCA (3 números de eixos) | 0,14 |
| phylofactorization (2 objetivos) | 0,06 |
| PCPS / PCC (3 números de eixos) | 0,01 / 0,007 |

**Estimativa:** 24.000 réplicas × 0,23 s ≈ **1,5 h de CPU**, ou seja, **cerca de 25
min com 4 workers**. Usaria no máximo 4 workers, para não competir com o run da
outra sessão. Uma extensão para N = 200 deve custar algo entre 3 e 4 vezes isso.
Proponho deixá-la de fora, a menos que o resultado com N = 100 mude o ranking.

Os resultados do piloto, só como indicação, com 10 réplicas por célula:
- **evoPCA:** acerta 0,9 com o deslocamento 2 nos três números de eixos e 0 no
  controle de prevalência.
- **Phylofactorization:**
  - com o objetivo atual, acerta 0,85 com o deslocamento 2 e 0,10 no controle;
  - com o objetivo só de interação, acerta 0,90 com o deslocamento 2 e 0 no
    controle.
- **PCPS:** acerta 0,25–0,30 nos três números de eixos.

Se esse padrão se confirmar, a recomendação continua de pé. Mas a
phylofactorization deveria usar o objetivo de interação, porque o atual dá falsos
positivos quando os clados diferem só em prevalência.

## Implementação

- **`analysis/39_clade_sensitivity.R`:**
  - usa `pcps/pcc/evopca_species_scores()` com `n_axes` (o parâmetro já existe);
  - acrescenta uma função para os dois objetivos da phylofactorization e o
    simulador do controle de prevalência;
  - grava chunks atômicos por condição, como o 17.
- **`R/`:** os dois acréscimos (`phylofactor_node_scores(objective = ...)` e
  `h_factor` em `sim_metacommunity_clade()`) entram com testes em `tests/testthat`.
  Os padrões atuais não mudam, então a rodada principal continua reproduzível.
- **Saídas:** `data/derived/clade_sensitivity/` e um resumo `clade_sensitivity.csv`
  em `data/derived/models/`, com taxas brutas e contrastes pareados de Monte Carlo
  entre as variantes (pareados por metacomunidade).
- Depois: revisão de código com os dois revisores, commit e proposta de texto para
  o Apêndice S2. No texto principal, só entra alguma frase se o ranking mudar ou se
  o objetivo da phylofactorization tiver que ser trocado na recomendação.

**Esforço:** cerca de 2–3 h de trabalho (código, testes, revisão, análise e texto),
mais ~25 min de máquina.

## Decisões para você

1. Rodar agora, antes do run da madrugada, ou amanhã?
2. Só N = 100 (recomendo), ou incluir N = 200?
3. *h* × 0,2 no controle de prevalência está bom? É o fator que, no piloto, dá uma
   diferença de prevalência parecida com a do deslocamento 2.
