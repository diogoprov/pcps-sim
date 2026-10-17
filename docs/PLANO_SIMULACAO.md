# Plano de simulação — inferência sobre composição ponderada pela filogenia

Versão 0.2 — 2026-10-16 — **enquadramento aprovado pelo Diogo; nada foi programado ainda**

Autoria: D. B. Provete (autor único). Documento interno de trabalho, em português. O texto do manuscrito será em inglês.
O contexto e as referências verificadas estão em `MEMORIA_PROJETO.md`.

---

## 0. Correção de premissa (ler antes)

A crítica à CWM, na forma mais simples, **já foi tratada para a matriz P**. Eu disse
o contrário antes, e estava errado.

- A crítica (Peres-Neto, Dray & ter Braak 2017; Zelený 2018) diz que testar CWM ~ E
  só com permutação de sítios infla o erro tipo I quando E estrutura a composição
  mas o atributo não se relaciona com E. A solução proposta é um teste de linhas e
  colunas com o maior P (max test). [Conferir no texto de Peres-Neto et al. 2017:
  só os metadados foram verificados.]
- A matriz P = QW é, formalmente, uma CWM: cada fuzzy set (linha de Q) funciona
  como um "atributo filogenético".
- **Duarte et al. 2016 já fizeram esse teste** (Tabela 1, lida no PDF publicado).
  No cenário "only environment", o site shuffle rejeita 99,4–100% das vezes, o taxa
  shuffle fica perto de 0,05 e eles adotam o largest P, citando Ter Braak, Cormont &
  Dray 2012.
- **Peres-Neto & Kembel 2015 também fizeram** (Appendix A, lido): permutação por
  linha, por coluna e linha/coluna. Pelo que lembro, a phylogenetic community
  composition deles também é uma CWM de autovetores filogenéticos. **Não verifiquei
  isso; preciso do PDF.**

**Consequência.** O eixo principal não pode ser "a crítica à CWM vale para P". Tem
que ser: **em que condições o remédio padrão (taxa shuffle + max test) falha para P?**
Três condições ficaram de fora dos dois estudos e são específicas da PFW:

1. **Formato da árvore.** Q depende da topologia via self-belonging, e Duarte et al.
   apresentam isso como uma virtude do método. O taxa shuffle permuta espécies entre
   pontas com self-belonging muito desigual quando a árvore é desbalanceada. Ninguém
   avaliou se o teste continua válido nesse caso. Os dois estudos usaram só árvores
   Yule ou uma única árvore.
2. **Estrutura filogenética "incômoda" sob H0.** Nos dois estudos, sob o H0 da
   hipótese 2 ("E afeta a composição, mas não via filogenia"), só o ótimo de nicho é
   aleatório na árvore. A amplitude de nicho (σ) e a abundância máxima (h) são
   sorteadas sem sinal. Na natureza, prevalência e amplitude de nicho costumam ter
   sinal. Se σ ou h têm sinal, P carrega estrutura filogenética mesmo que os ótimos
   não tenham. Resta saber se isso gera falsos positivos no teste do componente
   filogenético. É o análogo direto do problema da CWM em um nível acima: o max test
   protege contra "atributo sem relação com E", mas não necessariamente contra
   "estrutura filogenética em outra coisa que não o ótimo".
3. **Testes por eixo (PCPS + GLM, `pcps.sig`).** Duarte et al. só testaram o
   procedimento com o ajuste procrustiano do null em árvores Yule. A validade e a
   interpretação por clado dependem do espectro de autovalores. O espectro, por sua
   vez, depende da topologia (cf. Diniz-Filho et al. 1998 para o PVR).

**Checado em 2026-10-16: o item 2 continua aberto, e agora tem dois precedentes
que o sustentam.**

- **ter Braak, Peres-Neto & Dray 2017**, *PeerJ* 5: e2885, doi:10.7717/peerj.2885.
  Só li o resumo (via OpenAlex); a página da PeerJ deu 403. O p-max test controla
  o erro tipo I nos modelos usuais, mas "still gave slightly inflated error in more
  complex models that included both missing (but important) traits and missing
  (but important) environmental variables". Trabalham com atributos, sem filogenia.
  → σ e h com sinal são justamente "atributos ausentes mas importantes" e
  estruturados pela filogenia.
- **Duarte, Debastiani, Carlucci & Diniz-Filho 2018**, *Ecology* 99: 385–398. Lido
  em `references/`. Quando o atributo e a distribuição das espécies têm sinal
  filogenético moderado ou alto, o OLS CWM ~ E infla o erro tipo I **mesmo usando
  trait shuffle**. A correção que propõem (PGLS ou autovetores; NTDR) devolve o
  erro tipo I ao nominal, mas com perda de poder. Limites do estudo: só CWM de
  atributos, não a matriz P; só árvores Yule (`sim.bdtree`, λ = 0.1); μtol = 10 com
  σ sem sinal. Eles mesmos dizem que o fourth corner "possibly suffer[s] from the
  same shortcomings", **sem testar**.
- **Conclusão:** ninguém testou se o taxa shuffle ou o max test sobre P/PCPS falha
  quando a filogenia estrutura outros atributos da espécie que não o ótimo, nem
  como o formato da árvore modula isso. O enquadramento da seção 0 se mantém.
  Duarte et al. 2018 é o precedente mais próximo e precisa ser citado com destaque.
- **Para o desenho:** incluir a correção filogenética de Duarte et al. 2018
  (autovetores) como um método extra no objetivo 1, para ver se ela resolve o
  problema também para P.

---

## 1. Objetivos e hipóteses

**Objetivo 1 (eixo principal): validade da inferência sobre P/PCPS.** Estimar o erro
tipo I dos testes sobre D_P (ADONIS, Mantel) e sobre eixos PCPS (GLM) com site
shuffle, taxa shuffle e largest P, em função de:
(a) formato da árvore (γ, desbalanço);
(b) sinal filogenético em σ e h sob H0;
(c) tamanho do pool;
(d) incidência × abundância.

**Objetivo 2: comparação de métodos.** Comparar, nas mesmas metacomunidades, PFW/PCPS,
PCC (Peres-Neto & Kembel 2015), evoPCA (Pavoine 2016) e phylofactorization (Washburne
et al. 2017) quanto a:
(i) erro tipo I;
(ii) poder;
(iii) recuperação dos clados associados ao gradiente.

**Hipóteses de trabalho.** São inferências minhas, ainda não testadas.

- H1. O desbalanço da árvore afasta o erro tipo I do taxa shuffle sobre D_P do
  nominal mais do que o dos métodos baseados só em distâncias. Motivo: o
  self-belonging é heterogêneo entre pontas.
- H2. Sinal filogenético em σ ou h infla o erro tipo I do componente filogenético
  (taxa shuffle e largest P) em todos os métodos. A magnitude depende do formato da
  árvore.
- H3. Com γ muito negativo (nós basais longos), o PCPS concentra variação no PCPS1 e
  recupera bem os clados basais, mas mal os recentes. Com γ positivo, o espectro se
  achata e a recuperação piora.
- H4. A ordem dos métodos em poder e recuperação de clados muda com o formato da
  árvore.

---

## 2. Árvores

**Objetivo:** variar γ e desbalanço de forma (aproximadamente) independente, em
faixas realistas.

1. **Biblioteca.** Gerar cerca de 20 mil árvores ultramétricas com vários
   geradores, normalizar a profundidade da raiz para 1 e calcular γ (Pybus & Harvey
   2000) e o Colless normalizado (`treestats`, Janzen et al. 2024; ou `treebalance`).
   Geradores:
   - Yule e birth–death com taxas constantes: γ ≈ 0 ou positivo; desbalanço baixo.
   - diversity-dependence (`DDD`) e rate-shift (`TreeSim::sim.rateshift.taxa`):
     γ negativo.
   - Heterogeneidade de taxas entre clados: desbalanço alto. Opções: modelo β-α-η
     de Maliet et al. 2018 (**verificar se há implementação em R**); ou topologia
     por beta-splitting (`apTreeshape`) com tempos de ramificação atribuídos depois.
2. **Grade.** Uma grade 3 × 3 (γ baixo/médio/alto × desbalanço baixo/médio/alto),
   com limites definidos pelos quantis das árvores empíricas. As árvores são
   sorteadas da biblioteca por amostragem estratificada até cada célula ter
   árvores suficientes. As células sem árvores suficientes ficam registradas: elas
   indicam onde a correlação natural entre γ e desbalanço impede o controle.
3. **Referência empírica.** Calcular γ e Colless em subárvores regionais (50–200
   spp.) das 100 árvores datadas de anuros de Portik et al. 2023, para mostrar que a
   grade cobre o espaço real.
4. **Tamanho do pool (N):** 50, 100 e 200 spp., como em Duarte et al. 2016, para comparação direta.

## 3. Metacomunidades

O núcleo é o protocolo de Peres-Neto, Leibold & Dray 2012, usado tanto por Duarte
et al. 2016 quanto por Peres-Neto & Kembel 2015. Assim os resultados ficam
comparáveis com os deles.

- 100 sítios; E ~ U(0, 100).
- Resposta unimodal: L_ij = h_j · exp(−(E_i − u_j)² / (2σ_j²)), com Poisson em
  seguida. Analisar incidência (como nos dois estudos) **e abundância** (decidido).
- Ótimo u_j: atributo simulado na árvore por BM com λ de Pagel ∈ {0, 0.5, 1},
  reescalado para [−1, 101]. Uso λ e não o ρ de Grafen porque λ separa melhor sinal
  de formato da árvore. A manipulação do ρ interage com a profundidade dos nós.
- Novidade: σ_j (média μtol = 10, dp 10, truncada em > 0) e h_j (U(0, 30)) podem
  ser sorteados sem sinal ou por BM na árvore (λ = 1), depois reescalados para as
  mesmas distribuições marginais.

**Cenários:**

| Código | u_j | σ, h | E usado na análise | Estima |
|---|---|---|---|---|
| S0 | λ = 0 | sem sinal | E0 aleatório | tipo I (H1 e H2) |
| S1 | λ = 1 | sem sinal | E0 | tipo I (H1) |
| S2 | λ = 0 | sem sinal | E original | tipo I do componente filogenético |
| **S2σ** | λ = 0 | **σ com sinal** | E original | **tipo I sob estrutura incômoda** |
| **S2h** | λ = 0 | **h com sinal** | E original | **idem** |
| S3a | λ = 0.5 | sem sinal | E original | poder |
| S3b | λ = 1 | sem sinal | E original | poder |
| S3σ | λ = 1 | σ com sinal | E original | poder sob estrutura incômoda |

Fica fora (pode virar extensão): espaço e dispersão. Peres-Neto, Leibold & Dray 2012
já tratam disso, e Tsang & Cadotte 2025 tratam da autocorrelação composicional.

### 3b. Segunda camada: simulação baseada em processo (`mcomsimr`)

Proposta do Diogo (2026-10-16): usar o `mcomsimr` (github.com/plthompson/mcomsimr),
que implementa o framework de Thompson et al. 2020 *Ecol Lett* 23: 1314–1329,
doi:10.1111/ele.13568. A citação vem do README; o artigo ainda não foi lido.

Conferi o código de `simulate_MC()`:
- Crescimento independente da densidade: Gaussiano em E, com
  `r = max_r · exp(−((z_j − E)/(2·breadth_j))²)`.
- Competição Beverton–Holt com uma matriz α (`N·r / (1 + N %*% int_mat)`) e Poisson.
- Dispersão com kernel exponencial ou `disp_mat` dado pelo usuário.
- Burn-in e séries temporais.

Pontos de entrada para a filogenia:
- `env_optima`: vetor fornecido pelo usuário → ótimos simulados por BM na árvore.
- `env_niche_breadth`: aceita um valor por espécie, segundo a documentação. Falta
  confirmar em `env_traits()`. Daí sai σ com sinal.
- `int_mat`: **competição estruturada pela filogenia**. Por exemplo,
  α_jk = f(similaridade filogenética), que é a hipótese competition–relatedness.
  É uma nova fonte de estrutura incômoda sob H0, mecanística e ancorada em teoria.
- `landscape`, `disp_mat` e `env.df`: espaço e dispersão controlados.
- **Não há h por espécie**, porque `max_r` é comum a todas. Seria preciso adaptar a
  função (fork) ou abrir mão de h com sinal nesta camada.

**Decisão proposta: um desenho em duas camadas.**
- **Camada A (estática; Peres-Neto et al. 2012):** validação (reproduzir Duarte et
  al. 2016 e Peres-Neto & Kembel 2015), comparabilidade, varredura completa da grade
  de árvores e custo baixo.
- **Camada B (`mcomsimr`):** cenários principais com realismo de processo e três
  arquétipos de Thompson et al. 2020:
  1. species sorting, com dispersão intermediária;
  2. mass effects, com dispersão alta;
  3. competição estruturada pela filogenia.

  Cobre um subconjunto das células de árvore: os 4 cantos e o centro da grade
  γ × desbalanço.

**Cuidados com a camada B:**
1. Com espaço e dispersão limitada, os sítios deixam de ser permutáveis, e o site
   shuffle fica inválido por autocorrelação espacial, independentemente da
   filogenia. Para não confundir as duas fontes de erro, o baseline usa sítios sem
   estrutura espacial (`disp_mat` uniforme e E sem autocorrelação). O espaço entra
   como fator separado, com teste espacialmente restrito (por exemplo, Moran
   spectral randomization, `adespatial::msr`). Esse fator pode ficar para um
   segundo artigo.
2. O pacote só existe no GitHub (não está no CRAN): fixar o commit no `renv`. O
   README avisa que as simulações do artigo foram feitas em Julia e que o pacote
   pode dar resultados diferentes.
3. Custo: cada metacomunidade roda cerca de 2200 passos de tempo (padrão:
   initialization 200, burn-in 800, timesteps 1200). Medir no piloto; dá para
   encurtar a série e usar só o estado final, ou a média dos últimos passos.

**Do desenho antigo (2013–2016):** o gerador `sim.com`, com espécie de referência e
filtro por similaridade, sai do desenho, porque não tem gradiente ambiental e não
permite estimar erro tipo I nem poder. Podemos citar o manuscrito antigo como
motivação.

## 4. Métodos e testes

| Método | Estatística | Implementação |
|---|---|---|
| PFW / D_P | ADONIS (R²) e Mantel sobre √Bray–Curtis de P | `SYNCSA::matrix.p`, `vegan::adonis2`; reescrever o laço de permutação para eficiência |
| PCPS por eixo | GLM do eixo k ~ E (k = 1–4) | `PCPS::pcps`, `pcps.sig`; checar o null procrustiano |
| PCC (Peres-Neto & Kembel 2015) | conforme o artigo | **reimplementar a partir do PDF e do Appendix A** |
| evoPCA | RDA dos eixos ~ E, ou F global | `adiv::evopcahellinger` e `adiv::evopcachord` |
| Phylofactorization | fatores (arestas) associados a E | pacote `phylofactor` (GitHub; **verificar** disponibilidade e manutenção) |

**Permutações:** site shuffle; taxa shuffle (permutação das pontas, topologia
fixa); largest P. A phylofactorization tem inferência própria; tratar à parte.

**Reimplementação do PCC** (especificação em `docs/ESPEC_PCC.md`; o `picante` não tem o código). Incluir também a versão β (MPD + db-RDA) e uma variante da PCPS com padronização e regressão ponderadas, para testar a crítica (7) de Peres-Neto & Kembel 2015. Precisamos do PDF de Peres-Neto & Kembel 2015. A
validação é reproduzir o desenho do Appendix A (100 sítios × 50 spp., μtol = 5 ou
10, 1 gradiente) e conferir se as taxas de rejeição batem com as publicadas, dentro
do erro de Monte Carlo. Se não baterem, aí sim vale escrever ao Pedro.

**Validação da PFW.** Reproduzir a Tabela 1 de Duarte et al. 2016 nas árvores Yule
(N = 50, 100, 200; μtol = 10). É o teste de sanidade de todo o pipeline.

## 5. Métricas

- **Erro tipo I e poder:** taxa de rejeição a α = 0,05. Com 1000 réplicas, o IC95%
  de Monte Carlo em torno de 0,05 é de cerca de ±0,014. Critério de aceitabilidade:
  intervalo liberal de Bradley, [0,025; 0,075].
- **Recuperação de clados.** O "efeito verdadeiro" de cada nó interno é a diferença
  padronizada entre as médias de u_j dos dois clados-filhos, ou a contribuição do nó
  ao BM simulado; isso ainda precisa ser decidido no piloto. Cada método produz um
  escore por nó:
  - PCPS: correlação entre a pertença ao clado e os escores das espécies no eixo
    associado a E;
  - PCC e evoPCA: loadings das espécies agregados por clado;
  - phylofactor: as arestas escolhidas.

  Métrica: correlação de Spearman entre efeito verdadeiro e escore, e
  precision/recall dos k nós de maior efeito.
- **Espectro:** autovalor relativo dos eixos 1 a 4. Isso recupera a pergunta
  original do capítulo, agora como variável descritiva e não como resposta
  principal.

## 6. Análise dos resultados da simulação

- Tratar γ e o desbalanço como **contínuos**, porque cada árvore tem seu valor.
  Modelo binomial da taxa de rejeição em função de:
  - spline de γ × desbalanço;
  - N, cenário, método, tipo de dado (incidência/abundância);
  - interações relevantes.
- Ajuste em `brms` com priors explícitos (por exemplo normal(0, 1.5) nos
  coeficientes na escala logit), diagnósticos (R-hat, ESS, divergências) e PPC.
  Se ficar pesado, usar um GAM binomial (`mgcv`) para explorar e reservar o `brms`
  para o modelo final.
- **Figura-chave:** superfícies da taxa de rejeição sobre o plano γ × desbalanço,
  por método.

## 7. Custo e execução

Estimativa por alto, que **o piloto vai confirmar**:
- Objetivo 1: 8 cenários × 9 células de árvore × 3 N × 2 tipos de dado × 1000
  réplicas, com 999 permutações de taxa shuffle. São cerca de 430 mil
  metacomunidades.
- O taxa shuffle refaz Q, P e D_P a cada permutação. Isso é barato para N ≤ 150,
  mas o ADONIS a cada permutação pesa.

Opções para reduzir o custo:
1. 499 permutações e 500 réplicas nos cenários de poder.
2. Objetivo 2 só com N = 100 (incidência e abundância).
3. **Execução na AWS** (não há cluster). Opções: EC2 com muitos núcleos
   (instâncias c7i/c7a) e `targets` + `crew`, ou AWS Batch com `crew.aws.batch`.
   Primeiro medir o custo no piloto local e depois dimensionar. Usar instâncias
   spot, com checkpoint por réplica (o `targets` já faz isso).

**Piloto (primeira etapa de programação):** 50 réplicas por cenário em 3 células
de árvore, para medir tempo por réplica e checar a reprodução de Duarte et al. 2016.

## 8. Estrutura do compendium

```
pcps-sim/
  R/            funções: árvores, simulação de metacomunidades, métodos, testes, métricas
  analysis/     _targets.R, scripts do piloto e das rodadas completas
  data/raw/     árvores de Portik et al. 2023 (imutáveis)
  data/derived/ biblioteca de árvores, resultados da simulação
  manuscript/   .qmd (inglês)
  renv.lock
```

## 9. Decisões (2026-10-16)

1. Enquadramento: o eixo principal é "condições em que o max test falha para P" —
   **aprovado**.
2. Computação: AWS.
3. N = 50, 100 e 200.
4. Incidência e abundância.
5. Autoria: Provete como autor único. O PCC será reimplementado a partir do PDF e
   do Appendix A, sem pedir o código ao Pedro.

Ainda pendente: checar se ter Braak et al. 2017 *PeerJ* e Duarte et al. 2018
*Ecology* já cobrem o item 2 da seção 0.

## 10. Leituras necessárias antes de programar

Os PDFs estão em `references/` (fora do git). Já baixados: Duarte et al. 2016 e o
Appendix S2, Peres-Neto & Kembel 2015 e o Appendix A, Peres-Neto et al. 2017,
Zelený 2018, Duarte et al. 2018, Pavoine 2016 e Tsang & Cadotte 2025. Faltam os
dois de acesso aberto (ter Braak et al. 2017; Washburne et al. 2017).
- Peres-Neto & Kembel 2015 *Plant Ecol* (PDF principal): **obrigatório** para
  reimplementar o PCC.
- Peres-Neto, Dray & ter Braak 2017 *Ecography*; Zelený 2018 *J Veg Sci*.
- ter Braak et al. 2017 *PeerJ* (acesso aberto); Duarte et al. 2018 *Ecology*.
- Pavoine 2016 *Oikos*; Washburne et al. 2017 *PeerJ* (acesso aberto).
- Tsang & Cadotte 2025 *Ecol Lett*.
