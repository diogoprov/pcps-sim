# Fluxo para os PRs no treestats e no PCPS (2026-10-23)

Nesta pasta há duas contribuições prontas. Cada uma tem a função com o help em
roxygen, o arquivo `.Rd` gerado e um patch que se aplica limpo ao clone atual do
repositório do pacote.

| Pacote | O que entra | Arquivos | Testado |
|---|---|---|---|
| **treestats** (GPL-3, Thijs Janzen) | `sim_tree_shape(n, beta, c)`: árvore ultramétrica com desbalanço (β-splitting) e tempos de ramificação (taxa k^c) controlados separadamente | `treestats/R/sim_tree_shape.R`, `treestats/tests/testthat/test-sim_tree_shape.R`, `treestats/treestats-sim_tree_shape.patch` | 15 testes passam; o `.Rd` passa no `tools::checkRd` |
| **PCPS** (GPL-2, Vanderlei Debastiani) | `permut.torus()`, que gera as translações no toro (com simetrias), e o argumento `permutation.site` em `pcps.sig()`, `matrix.p.sig()` e `matrix.p.null()`, que substitui o site shuffle livre | `PCPS/R/permut.torus.R`, `PCPS/PCPS-permutation_site.patch`, `PCPS/check_permutation_site.R` | todas as checagens passam (sem o argumento, o resultado é idêntico ao atual; com ele, há uma estatística nula por permutação e o P é calculado sobre elas) |

## Antes de abrir qualquer PR

1. **Anonimato.** O manuscrito está em revisão duplo-cega. Um PR público na sua
   conta do GitHub, descrevendo exatamente os métodos do artigo, pode identificar
   você para um revisor que procure. Recomendo abrir os PRs **depois do aceite**,
   ou pelo menos não citar o manuscrito nem o repositório `pcps-sim` até lá. A
   decisão é sua.
2. **Licença.** O código no `pcps-sim` é MIT e o autor é você, então você pode
   contribuí-lo sob GPL-3 (treestats) e GPL-2 (PCPS).
3. **Uso de IA.** O código foi escrito com assistência do Claude, como está no log
   do projeto. Alguns mantenedores pedem que isso seja declarado. Os textos de PR
   abaixo trazem uma linha opcional para isso.
4. **Referência pendente no PCPS.** O help de `permut.torus()` está sem
   referência para o teste de translação no toro, porque eu não verifiquei
   nenhuma. Há um `TODO` no topo do arquivo. Acrescente uma referência que você
   tenha conferido, ou remova o `TODO` antes do PR.
5. **Autoria no PCPS.** As funções do PCPS têm uma linha `@author`. Acrescente a
   sua em `permut.torus.R` (`#' @author Seu Nome <email>`).

## treestats

O treestats é um pacote de estatísticas de árvores. Como um simulador pode estar
fora do escopo dele, **abra primeiro uma issue** perguntando ao mantenedor se
tem interesse.

**1. Issue** (rascunho, em https://github.com/thijsjanzen/treestats/issues):

> **Title:** Proposal: simulate trees with independently controlled imbalance and branching times
>
> Hi Thijs, I have been using treestats to compute the gamma statistic and the Colless index of simulated trees. To build a grid of tree shapes where both vary independently, I wrote a small simulator: the topology comes from Aldous' beta-splitting model and the internode intervals from a pure-birth process where each lineage splits at rate k^c (c = 0 is Yule, c < 0 gives negative gamma). Since gamma depends only on the intervals and Colless only on the topology, beta mostly moves Colless and c mostly moves gamma. It is pure R (ape and stats only), documented and tested. Would a function like `sim_tree_shape(n, beta = 0, c = 0)` fit in treestats? If so, I can open a pull request.

**2. Se ele aceitar**, no terminal:

```bash
gh repo fork thijsjanzen/treestats --clone
```

```bash
cd treestats && git checkout -b sim-tree-shape
```

```bash
git apply ~/Github/pcps-sim/contrib/treestats/treestats-sim_tree_shape.patch
```

O patch acrescenta `R/sim_tree_shape.R`, o teste, os quatro `.Rd` e
`export(sim_tree_shape)` no `NAMESPACE`. Se você preferir regenerar o help com
`devtools::document()`, desfaça a mudança que o roxygen faz em
`Config/roxygen2/version` no `DESCRIPTION`, porque ela não faz parte da
contribuição.

**3. Compilar e checar.** O treestats tem código C++ e, nesta máquina, a
compilação falhou no linker (`library 'emutls_w' not found`). Isso indica que
falta o gfortran do toolchain do CRAN para macOS. Instale-o (veja a página de
ferramentas do CRAN para macOS) e rode, no R, dentro da pasta do pacote:

```r
devtools::test(filter = "sim_tree_shape")
devtools::check()
```

**4. NEWS.md.** Acrescente uma linha no topo, no formato do arquivo, por exemplo:
`- Added sim_tree_shape(), which simulates trees with independently controlled topology (beta-splitting) and branching times (lineage-dependent pure birth).`

**5. Commit, push e PR:**

```bash
git add -A && git commit -m "Add sim_tree_shape(): trees with independent control of imbalance and branching times"
```

```bash
git push -u origin sim-tree-shape
```

```bash
gh pr create --repo thijsjanzen/treestats --title "Add sim_tree_shape()" --body-file ~/Github/pcps-sim/contrib/treestats/PR_BODY.md
```

Rascunho do corpo do PR (já salvo em `contrib/treestats/PR_BODY.md`; edite o número da issue):

> Adds `sim_tree_shape(n, beta = 0, c = 0)`, discussed in #<issue>. Topology from Aldous' beta-splitting model; internode intervals from a pure-birth process with per-lineage rate k^c, rescaled to root depth 1; branching times assigned to nodes in a uniform random linear extension of the tree order. Only uses ape and stats. Tests check the output (binary, ultrametric, root depth 1), argument checks, and that c moves gamma without moving the Yule-normalized Colless index and beta does the opposite. The documentation explains how to reach target values of gamma and Colless by calibration and rejection.
>
> *(Optional)* The code was written with the assistance of an AI model (Claude, Anthropic) and checked by me.

## PCPS

**O GitHub do PCPS está desatualizado.** O `master` é a versão 1.0.8, de 2021; o
CRAN tem a 1.0.9, de junho de 2026. No GitHub, `FUN.ADONIS` ainda chama
`vegan::adonis`, que o vegan atual já removeu. As três funções que o patch
altera são idênticas nas duas versões, então o patch se aplica de qualquer
forma. Ainda assim, convém **escrever antes ao Vanderlei**, por e-mail ou por
issue, perguntando se ele aceita PRs no GitHub e se pode subir a 1.0.9 para lá.
Você tem contato direto com ele, e ele está nos agradecimentos do artigo.

**1. Mensagem** (rascunho):

> Hi Vanderlei, in a simulation study with spatially autocorrelated sites I needed a restricted site shuffle in `pcps.sig()`: the unrestricted shuffle was strongly liberal (rejection rates of 0.18–0.70 under the null), whereas torus translations of a regular grid kept the type I error near the nominal level. I wrote a minimal extension: a new argument `permutation.site` in `pcps.sig()`, `matrix.p.sig()` and `matrix.p.null()` (a matrix of site permutations that replaces `SYNCSA::permut.vector()` for the site shuffle; P values are then computed over its rows), and a helper `permut.torus()` that returns all torus translations of a grid, with its reflections and rotations. With the default `permutation.site = NULL` the results are unchanged. The GitHub master is still at 1.0.8; could you push 1.0.9 so that I can open a pull request against it?

**2. Fork, branch e patch:**

```bash
gh repo fork vanderleidebastiani/PCPS --clone
```

```bash
cd PCPS && git checkout -b permutation-site
```

```bash
git apply ~/Github/pcps-sim/contrib/PCPS/PCPS-permutation_site.patch
```

O patch altera `R/matrix.p.null.R`, `R/pcps.sig.R` e `R/matrix.p.sig.R`,
acrescenta `R/permut.torus.R`, inclui `export(permut.torus)` no `NAMESPACE` e o
arquivo novo no `Collate` do `DESCRIPTION`, e atualiza `man/pcps.sig.Rd` e
`man/matrix.p.null.Rd`. Ele acrescenta também `man/permut.torus.Rd`.

**Não rode `devtools::document()` sem conferir o diff.** O roxygen atual
reformata o `NAMESPACE` inteiro e vários `.Rd` do PCPS (que foi gerado com o
roxygen 7.1). O patch já traz só as mudanças necessárias. Se você editar o
roxygen de `permut.torus.R` (para a referência ou o `@author`), regenere e
inclua no commit apenas `man/permut.torus.Rd`.

**3. Checar** (no R, na pasta do pacote):

```r
devtools::install()
source("~/Github/pcps-sim/contrib/PCPS/check_permutation_site.R")
devtools::check()
```

**4. Commit, push e PR:**

```bash
git add -A && git commit -m "Restricted site permutations: permutation.site argument and permut.torus()"
```

```bash
git push -u origin permutation-site
```

```bash
gh pr create --repo vanderleidebastiani/PCPS --title "Restricted site permutations (permutation.site, permut.torus)" --body-file ~/Github/pcps-sim/contrib/PCPS/PR_BODY.md
```

Rascunho do corpo do PR (já salvo em `contrib/PCPS/PR_BODY.md`):

> Adds an optional argument `permutation.site` to `pcps.sig()`, `matrix.p.sig()` and `matrix.p.null()`: a matrix with one permutation of the sampling units per row that replaces the unrestricted site shuffle (for instance to keep spatial autocorrelation). P values of the site shuffle are computed over its rows. Adds `permut.torus(nrow, ncol, symmetries, runs)`, which returns all torus translations of a regular grid (optionally with reflections and rotations, 8 × nrow × ncol − 1 permutations for a square grid). With the default `permutation.site = NULL`, results are identical to the current version (checked with the same seed). Invalid permutation matrices are refused.
>
> *(Optional)* The code was written with the assistance of an AI model (Claude, Anthropic) and checked by me.

## Depois do merge

- Se o treestats incorporar `sim_tree_shape()`, o `pcps-sim` pode passar a usá-la,
  numa versão futura, em vez de `R/trees.R`. Os resultados não mudam, porque a
  função é a mesma. Só não faça isso antes de publicar, para não mexer na
  reprodutibilidade do artigo.
- Cite as versões dos pacotes que incluírem as contribuições, se você as usar em
  trabalhos futuros.
