# Especificação: phylogenetic gradient analysis (Peres-Neto & Kembel 2015)

Versão 0.1 — 2026-10-16. Fonte: o PDF (`references/PeresNeto_Kembel_2015_PlantEcol.pdf`,
pp. 711–716) e o Appendix A (`references/appendices/`). Algumas equações saíram
corrompidas na extração do PDF; a reconstrução está marcada com **[R]** e precisa
ser conferida na figura do PDF.

## 0. Código disponível

O artigo (p. 716) diz: "Functions implementing both frameworks ... are available in
the picante R package". **Isso não se confirma.** Clonei `github.com/skembel/picante`
(último commit em 2026-08-23) e procurei por `peres`, `eigenvector`, `cmdscale` e
`weighted mean`: não há nenhuma função de PGA ou PCC. Vamos reimplementar.

## 1. Framework A: phylogenetic beta-diversity + db-RDA ("PGA-β")

1. PB (n × n) = MPD entre comunidades (= COMDIST), com incidência.
2. A = −½ PB², centrada duplamente: G = (I − 11ᵀ/n) A (I − 11ᵀ/n).
3. H = E (EᵀE)⁻¹ Eᵀ. pseudo-F = [tr(HGH)/p] / [tr((I−H)G(I−H))/(n−p−1)].
   O PDF diz (n − p); usar a forma de McArdle & Anderson 2001 e conferir **[R]**.
4. Na prática, equivale a `vegan::dbrda` ou `adonis2` sobre `picante::comdist`.
   Implementar direto com álgebra matricial, para ganhar velocidade nas
   permutações.

## 2. Framework B: phylogenetic community composition (PCC)

Entradas: L (n × k), de incidência ou abundância (a abundância é extensão nossa;
eles só usaram incidência); árvore ultramétrica; E (n × p).

1. **Autovetores filogenéticos.** PD = distância patrística (k × k). B = −½ PD.
   O PDF escreve PD **sem elevar ao quadrado**, ao contrário de PB. Isso equivale à
   PCoA de √PD, ou seja, aos autovetores da matriz de covariância BM duplamente
   centrada. Manter como está escrito e testar a alternativa (−½ PD²) na validação.
   Centrar duplamente e fazer a eigen-decomposição. PEV = k × (k − 1) autovetores
   com autovalor > 0. A escala dos autovetores é irrelevante depois do passo 2.
2. **Padronização ponderada.** w_j = soma da coluna j de L (prevalência; com
   abundância, a abundância total). Para cada autovetor:
   x̄_w = Σ w_j PEV_j / Σ w_j,  s²_w = Σ w_j (PEV_j − x̄_w)² / Σ w_j.
   O denominador é tr(W_k), sem −1. PEV_std = (PEV − x̄_w)/s_w.
3. **PCC** (n × (k − 1)) = diag(1/r_i) · L · PEV_std, com r_i = soma da linha i de
   L (riqueza). É a CWM dos autovetores padronizados.
4. **E padronizado** com pesos r_i: média e desvio ponderados, com denominador
   Σ r_i.
5. **Regressão ponderada.** PCĈ = E_std (E_stdᵀ W_n E_std)⁻¹ E_stdᵀ W_n PCC, com
   W_n = diag(r). Para um único preditor, a inclinação é igual à correlação do
   fourth corner.
6. **R² por autovetor** **[R]**. PL_i = L ∘ (1_n PEV_std,iᵀ): as presenças viram o
   valor padronizado da espécie. Como os autovetores têm média ponderada 0 e
   variância 1, tr(PL_iᵀPL_i) = Σ_j w_j = soma de L, igual para todos os autovetores.
   R²_i = tr(PCĈ_iᵀ W_n PCĈ_i) / tr(PL_iᵀ PL_i).
   R²adj,i = 1 − (n − 1)/(n − p − 1) · (1 − R²_i).
7. **Estatística de teste: variância dos R²adj,i** entre os k − 1 autovetores. Não
   é a média. Racional: autovetores iniciais devem explicar mais quando o nicho tem
   sinal filogenético. O denominador da variância não importa para um teste de
   permutação.
8. **R² reportado:** a média só dos R²adj,i positivos (recomendação deles, p. 715).
9. **Ordenação:** PCA de PCĈ, ponderada por r **[a decidir]**. Escores das espécies
   no PC1 = PEV_std · loadings do PC1. Esses escores alimentam a métrica de
   recuperação de clados (plano, seção 5).

## 3. Testes

- **Linhas:** permutar as linhas de L em relação a E e recalcular tudo. Os pesos
  W_n mudam de posição, então recalcular a padronização de E.
- **Colunas:** permutar as colunas de L (espécies em relação às pontas da árvore) e
  recalcular w, PEV_std e PCC.
- p = (#{null ≥ obs} + 1)/(nperm + 1); **p_max = max(p_linha, p_coluna)**.
- Eles usaram 99 permutações nas simulações e 999 nos dados empíricos.

## 4. Validação contra a Tabela 1 (100 sítios × 50 spp., 1 gradiente)

Desenho do Appendix A:
- **uma única árvore** Yule de 50 spp. (ape), fixa em todas as réplicas;
- traço do ótimo por `phytools::fastBM`, reescalado para [−1, 101]; E ~ U(0, 100);
- h ~ U(0, 30); σ ~ N(μtol, 10). σ negativo não é tratado no texto: usar |σ| ou
  truncar, e registrar a escolha **[ambíguo]**;
- Poisson → incidência; 1000 réplicas; 99 permutações.

Alvos (taxas de rejeição publicadas, PCC):

| Cenário | μtol | row/col | row | col |
|---|---|---|---|---|
| neither | 5 | 0.004 | 0.058 | 0.042 |
| neither | 10 | 0.004 | 0.043 | 0.050 |
| só ambiente | 5 | 0.054 | 1.000 | 0.054 |
| só ambiente | 10 | 0.052 | 1.000 | 0.052 |
| só filogenia | 5 | 0.030 | 0.043 | 0.267 |
| só filogenia | 10 | 0.040 | 0.050 | 0.276 |
| ambos (poder) | 5 | 0.945 | 1.000 | 0.945 |
| ambos (poder) | 10 | 0.947 | 1.000 | 0.947 |

PGA-β, poder row/col: 0.937 (μtol = 5) e 0.939 (μtol = 10).

**Critério:** cada taxa reproduzida deve cair no IC95% binomial da taxa publicada
(n = 1000). Exceção: o poder depende da árvore única, que não temos. Usar várias
árvores Yule e esperar só a ordem de grandeza.

## 5. Afirmações de Peres-Neto & Kembel 2015 que interessam ao nosso eixo 1

Todas da p. 721 (Discussão):
- "(2) ... RLQ and PCPS are usually tested using the row permutation approach".
- "(5) the PCPS approach is based on a Bray–Curtis dissimilarity ... decomposed via
  PCoA using the Lingoes correction ... this correction has been shown to reduce
  the statistical power".
- "(7) the PCPS approach also uses weights to build their site by phylogeny matrix,
  though neither this matrix is standardized using these weights nor the regression
  approach using Redundancy Analysis is weighted, leading to **very unstable
  regression estimates (Peres-Neto and Dray unpublished simulations)**".
  → É uma crítica à PFW/PCPS **nunca publicada com dados**. O nosso objetivo 1
  pode testá-la diretamente: PCPS com e sem padronização e regressão ponderadas.
- p. 717, **sobre o formato da árvore**: "it is likely that tree shape (e.g., tree
  symmetry) impacts our ability to link phylogeny, species distributions, and
  environmental variation. Future studies will be important to assess the
  performance of the proposed frameworks for different tree shapes, which is
  beyond the scope of the present study."
  → **É a chamada publicada que justifica o eixo do formato da árvore.** Ela
  substitui a frase do draft de Duarte et al., que não chegou à versão publicada.
