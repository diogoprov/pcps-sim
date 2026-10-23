# Ponto 1: estimandos de S2σ/S2h (proposta, 2026-10-22)

Proposta para revisão do Diogo. Nada foi editado no manuscrito. Nenhuma simulação
nova; os números novos vêm de `data/derived/diag_s2sigma/diag_s2sigma.csv`
(script `analysis/09_diag_s2sigma.R`, já no repositório).

## 0. Ajuste no meu próprio diagnóstico

Na conversa eu disse que a correlação da CWM de σ com E² mostrava uma associação
filogenética com o gradiente. Conferindo o CSV, isso precisa ser refinado:

| Cenário | Célula | N | rejeição GLM PCPS2 (linear) | rejeição GLM PCPS2 (quadrático) | \|r\| PCPS2 × E | \|r\| PCPS2 × E² | \|r\| CWM σ × E² | K de Blomberg de σ |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| S2  | γ baixo | 100 | 0,07 | 0,07 | 0,29 | 0,86 | 0,87 | 0,26 |
| S2  | γ baixo | 200 | 0,04 | 0,05 | 0,34 | 0,84 | 0,92 | 0,21 |
| S2σ | γ baixo | 100 | 0,16 | 0,03 | 0,44 | 0,80 | 0,87 | 0,88 |
| S2σ | γ baixo | 200 | 0,30 | 0,10 | 0,59 | 0,73 | 0,93 | 0,96 |
| S2  | Yule    | 100 | 0,08 | 0,05 | 0,32 | 0,40 | 0,88 | 0,06 |
| S2σ | Yule    | 100 | 0,11 | 0,10 | 0,40 | 0,31 | 0,88 | 0,64 |

(100 réplicas por linha, incidência, taxa shuffle, 499 permutações; medianas.)

- A CWM de σ × E² é **igual** em S2 e S2σ (0,87–0,93). Ela é um efeito de borda do
  modelo de resposta gaussiano e não tem nada de filogenético em si.
- O que muda em S2σ é a correlação **linear** do PCPS2 com E (0,29 → 0,44 e
  0,34 → 0,59 com γ baixo). O GLM quadrático rejeita *menos* (0,03–0,10), então as
  rejeições vêm do componente linear, não da curvatura.

Conclusão: o argumento do Codex continua de pé, mas o mecanismo não é "a curvatura
vira filogenética". Em cada réplica, a composição filogenética covaria
linearmente com E quando σ é conservado. O mecanismo exato fica em aberto e deve
ser escrito com hedge. O ponto lógico central não depende dele (seção 1).

## 1. Tabela de hipóteses por cenário (para o Diogo; não vai para o manuscrito)

Notação: $H_0^{site}$, composição ponderada pela filogenia ($\mathbf{P}$)
independente de E (site shuffle); $H_0^{taxa}$, espécies intercambiáveis nas pontas
da árvore, isto é, nenhum atributo da resposta ($u_j$, $\sigma_j$, $h_j$) tem sinal
filogenético (taxa shuffle); $H_0^{max} = H_0^{site} \cup H_0^{taxa}$ (o max test
controla o erro tipo I quando pelo menos uma é verdadeira); $H_0^{u}$, ótimos sem
sinal filogenético, com σ e h livres.

| Cenário | $H_0^{site}$ | $H_0^{taxa}$ | $H_0^{max}$ | $H_0^{u}$ | O que a taxa de rejeição estima |
|---|---|---|---|---|---|
| S0  | V | V | V | V | erro tipo I |
| S1  | V | F | V | F | erro tipo I |
| S2  | F | V | V | V | erro tipo I (componente filogenético) |
| S2σ | F | **F** | **F** | V | rejeição com σ conservado: não é erro tipo I de $H_0^{max}$; é a frequência com que o teste sugeriria ótimos conservados sem que eles o sejam |
| S2h | F | **F** | **F** | V | idem, com h conservado |
| S3a | F | F | F | F | poder |
| S3b | F | F | F | F | poder |
| S3σ | F | F | F | F | poder (com σ conservado) |

Ponto-chave: em S2σ/S2h as duas hipóteses componentes são falsas, então o max test
rejeitar não é erro tipo I **da hipótese que ele testa**. Só é "erro" em relação a
$H_0^{u}$, que nenhuma permutação de pontas consegue isolar: permutar pontas quebra
ao mesmo tempo o vínculo da filogenia com $u$ e com σ. É isso que o texto atual
ofusca quando escreve "inflated the type I error".

## 2. Tabela teste × cenário (para o Diogo; base do texto)

Todos os testes usam o mesmo par de permutações, então a hipótese nula é a mesma
(seção 1). O que muda é a estatística e, portanto, a qual parte da violação de
$H_0^{taxa}$ cada teste é sensível.

| Teste | Estatística | Alvo da estatística | S2σ (estático, N = 100–200, incidência) | S2h |
|---|---|---|---|---|
| PERMANOVA | pseudo-F de $D_P$ ~ E | associação multivariada de $\mathbf{P}$ com E | ficou mais conservador | sem mudança |
| Mantel | r($D_P$, \|ΔE\|) | idem, por pares de sítios | mais conservador | sem mudança |
| GLM PCPS1 | F do eixo 1 ~ E | contraste filogenético dominante | mais conservador | sem mudança |
| GLM PCPS2 | F do eixo 2 ~ E | segundo contraste (em parte um arco do eixo 1) | rejeição até 0,36 | sem mudança |
| RDA PCPS1–2 | F multivariado | os dois contrastes juntos, invariante a rotação | 0,004–0,021 | (só na 2ª rodada) |
| GLM PCPS com Procrustes | F após rotação dos eixos nulos | idem, eixo a eixo | 0,002–0,009 | (só na 2ª rodada) |
| PCC | var. do R² ajustado entre autovetores | vários contrastes filogenéticos | rejeição até 0,22 | > 0,075 só na camada de processos |

## 3. Texto proposto para o manuscrito (`methods_results.qmd`)

### 3.1. Scenarios and estimands (substitui o parágrafo da linha 73 e a tabela)

> The eight trait scenarios were the same in both kinds of simulation (@tbl-scenarios). Every test was evaluated with a site shuffle, whose null hypothesis is that $\mathbf{P}$ is independent of $E$, and a taxa shuffle, whose null hypothesis is that species are exchangeable among the tips of the tree, that is, that no attribute of their response to the gradient ($u_j$, $\sigma_j$ or $h_j$) carries phylogenetic signal. The max test controls the type I error when at least one of the two holds (Peres-Neto & Kembel, 2015; Duarte et al., 2016). S0 and S1 used a gradient independent of the one that generated the data, so the site-shuffle null hypothesis was true, and S2 used the generating gradient with optima and the other attributes without phylogenetic signal, so the taxa-shuffle null hypothesis was true; the rejection rate in these scenarios is therefore the type I error. In S2σ and S2h, the optima had no phylogenetic signal, but niche breadth or maximum abundance did. Thus, both null hypotheses were false, and a rejection is not a type I error of the max test. However, a significant test is usually read as evidence that lineages track the gradient through their optima, and this reading would be wrong in these two scenarios, since no permutation of the tips can separate signal in the optima from signal in the other attributes. We therefore report the rejection rate in S2σ and S2h as the probability that a test suggests conserved optima when only niche breadth or maximum abundance is conserved, and not as a type I error. The rejection rate in S3a–S3σ is the power. For aim (iii), the estimand was the probability that a method identifies a planted clade.

| Scenario | Optima (λ) | Signal in σ or h | Gradient analysed | Null hypothesis of the max test | Estimand |
|---|---|---|---|---|---|
| S0 | 0 | none | independent | true | type I error |
| S1 | 1 | none | independent | true | type I error |
| S2 | 0 | none | generating | true | type I error |
| S2σ | 0 | σ | generating | false | rejection with conserved σ only |
| S2h | 0 | h | generating | false | rejection with conserved h only |
| S3a | 0.5 | none | generating | false | power |
| S3b | 1 | none | generating | false | power |
| S3σ | 1 | σ | generating | false | power with conserved σ |

### 3.2. Abstract

Item 3, primeira e segunda frases (o resto do item fica igual):

> Tree imbalance did not affect the type I error of any test. When only niche breadth carried phylogenetic signal, the linear model on the second PCPS axis rejected the null hypothesis in up to 0.36 of the replicates with low γ and 200 species, and the PCC in up to 0.22, although the optima had no phylogenetic signal.

Item 4, a recomendação correspondente:

> We recommend reporting γ of the analysed tree, not relying on a single test, not reading a significant second PCPS axis or PCC as evidence of conserved optima when niche breadth may be conserved, considering abundance data for the PCC, and using evoPCA or phylofactorization to identify clades associated with a gradient.

### 3.3. Results (substitui o parágrafo da linha 118)

> Phylogenetic signal in niche breadth increased the rejection rate of only two tests (@fig-type1), even though the optima had no phylogenetic signal. Compared with the same scenario without signal (S2), signal in σ (S2σ) increased the rejection rate of GLM PCPS2 by 0.02 (95% CI 0.01–0.03) in Yule trees at N = 100, and by 0.05 (0.04–0.06), 0.10 (0.08–0.12) and 0.21 (0.19–0.23) with low γ at N = 50, 100 and 200, with raw rates of up to 0.36. The PCC increased by 0.05 (0.03–0.06) in Yule trees and by 0.02 (0.01–0.03) with low γ at N = 100, with raw rates of up to 0.22. In contrast, PERMANOVA, GLM PCPS1 and the Mantel test became more conservative under signal in σ (by 0.005–0.030), and signal in maximum abundance (S2h) did not change the rejection rate of any test. Since the taxa-shuffle null hypothesis was false in both scenarios, these rates are not type I errors of the max test, but they show how often each test would suggest conserved optima when only niche breadth was conserved. The RDA PCPS rejected in 0.004–0.021 of the replicates in S2σ and lost much less power with low γ than GLM PCPS1 (posterior medians 0.65 with low γ against 0.76 in Yule trees, against 0.34 and 0.74 for GLM PCPS1, at N = 100 in S3b), whereas rotating the null axes by Procrustes removed the increase of GLM PCPS2 but made it conservative (Appendix S2). As a result, H2 was rejected as stated.

Legenda da Fig. 2 (linha 120):

> Rejection rate of each test over the γ × Colless grid at N = 100 with incidence data (posterior medians of the beta-binomial model). Rows S0 to S2 are type I errors; in S2σ and S2h the optima had no phylogenetic signal, but niche breadth or maximum abundance did. Notice that rows S0 to S2h are almost uniform across the grid, whereas GLM PCPS2 doubles its rejection rate in the low-γ cells of S2σ.

Camada de processos (linha 126), só a frase sobre a hipótese nula:

> Under the null hypotheses, no lower CI bound exceeded 0.075 at N = 100, but raw rates of GLM PCPS2 in S2σ and of the PCC in S2h did at N = 100 and 200 (Appendix S2).

(Aqui o plural "null hypotheses" basta se a seção 3.1 já estiver clara. Alternativa:
"In S0–S2h, no lower CI bound…".)

### 3.4. Discussion (substitui o parágrafo da linha 142)

> Phylogenetic signal in niche breadth increased the rejection rate of only the second PCPS axis and the PCC, even though the optima had no phylogenetic signal. Previous evaluations could not detect this, because their simulations drew niche breadth independently of the phylogeny (Duarte et al., 2016, 2018; Peres-Neto et al., 2017). Strictly, these rejections are not type I errors of the max test, since the taxa shuffle assumes that species are exchangeable on the tree, and conserved σ violates this assumption as much as conserved optima do. Permutation tests reject whenever the randomization breaks down a structure of the data, whether or not this structure is the one of interest (Hardy, 2008), and no permutation of the tips can break the link between the phylogeny and the optima while keeping the link with niche breadth. The practical consequence is the same as that of an inflated type I error, because a significant test is usually read as evidence that lineages track the gradient through their optima. When σ was conserved, the second PCPS axis became more correlated with the gradient itself, whereas the curvature of the response, which tied the community-weighted mean of σ to the square of the gradient equally with and without signal, did not seem to drive the rejections (Appendix S2). This pattern seems to load on the axis orthogonal to the basal contrast of the first PCPS axis, which would explain why the increase was largest in trees with low γ and with 200 species, although the low level of γ was also more extreme at N = 200. Unfortunately, our diagnostic runs cannot tell how clades of broad-niche species come to covary linearly with the gradient within a replicate. Rotating the null axes by Procrustes removed this increase but made the test conservative and less powerful, whereas the RDA on both axes stayed below the nominal level and kept most of its power (Appendix S2). Thus, a significant second PCPS axis or PCC should not be taken as evidence of conserved optima when niche breadth may be phylogenetically conserved, and constrained randomizations that preserve phylogenetic structure (Braga et al., 2018) or model-based tests with phylogenetic random effects (Li & Ives, 2017; Tsang & Cadotte, 2025) seem to be promising alternatives.

Saiu a frase sobre ter Braak et al. (2017), que tratava de erro tipo I em testes
por sítio na relação traço–ambiente; ela reforçava o enquadramento antigo. Se
quiser mantê-la como analogia, ela pode entrar depois da frase de Hardy (2008).

### 3.5. H2 na Introdução (linha 45)

O texto atual diz que o sinal em σ ou h "would inflate the type I error of all tests".
Duas opções:

- **(a) manter H2 como foi formulada** (é a hipótese que foi de fato proposta antes
  dos resultados) e deixar a seção 3.1 explicar por que o resultado não é erro
  tipo I. É a opção mais honesta do ponto de vista do registro.
- **(b) reescrever:** "(ii) phylogenetic signal in niche breadth or maximum
  abundance would increase the rejection rate of all tests even without signal in
  the optima, with a magnitude that depends on tree shape (H2)".

Recomendo (a).

## 4. Outros arquivos que precisam acompanhar (não redigidos ainda)

- **Cover letter** (`cover_letter.qmd:19`): "phylogenetic signal in niche breadth
  inflated the type I error of…" → "…led the test on the second PCPS axis and the
  PCC to reject the null hypothesis often even though the optima had no
  phylogenetic signal".
- **SI** (`supporting_information.qmd`): linha 81 (subtítulo "Type I error" engloba
  S2σ; trocar "inflation" por "increase"); linha 105 (acrescentar a correlação
  linear do PCPS2 com E e a comparação da CWM de σ × E² entre S2 e S2σ, com os
  números da seção 0); linha 109 ("where the type I error was inflated" → "where
  the rejection rate increased"); linhas 91 e 133 (mesma troca de vocabulário).
- **Código e figuras:** `scenarios_layerA()` rotula S2σ/S2h como `type_I_phylo`,
  e a Fig. 2 e os modelos formais agrupam por esse rótulo. Proponho **não** mexer
  nos modelos (eles continuam estimando taxas de rejeição) e só trocar a legenda e,
  se houver, o título do painel ("Type I error" → "Rejection rate under true or
  partly true null hypotheses", ou separar S2σ/S2h numa faixa própria). Preciso
  olhar `22_manuscript_figures.R` antes de propor algo concreto.
- **Contagem de palavras:** a seção 3.1 cresce cerca de 110 palavras e a 3.4 cerca de
  80; o texto está perto do limite de 8.000. Se for preciso, dá para compensar
  encurtando a 3.4 (a frase "Unfortunately…" pode sair).

## 5. Referências

Nenhuma referência nova. Todas as citadas já estão no manuscrito: Hardy 2008,
Braga et al. 2018, Li & Ives 2017, Tsang & Cadotte 2025, Peres-Neto & Kembel 2015,
Duarte et al. 2016, 2018, Peres-Neto et al. 2017.
