# Revisão crítica do código e plano de correção

- **Autoria deste documento:** Codex, assistente de inteligência artificial da OpenAI, baseado em GPT-6.
- **Data:** 22 de outubro de 2026 (America/Campo_Grande).
- **Solicitante:** Diogo B. Provete.
- **Origem:** crítica e plano apresentados pelo Codex nesta conversa, consolidados a pedido do solicitante.
- **Natureza:** revisão metodológica simulada na perspectiva de um revisor da *Methods in Ecology and Evolution*. Não constitui parecer da revista nem revisão por um especialista humano independente.
- **Escopo desta entrega:** documentação da avaliação e do plano; as correções propostas ainda não foram implementadas. A autoria deste documento não implica coautoria do artigo.

## 1. Parecer geral e verificações realizadas

**Recomendação simulada: revisão maior (major revision).** A simulação tem uma base computacional consistente e resultados potencialmente publicáveis, mas algumas conclusões centrais ainda não estão demonstradas pelo desenho atual. A maior preocupação é a definição da hipótese nula, seguida pela interpretação causal dos contrastes.

Foram examinados as funções de simulação, os testes, o pipeline, as análises, o manuscrito e seu suplemento. Toda a suíte de `tests/testthat` foi executada com as dependências locais do projeto e passou. Na camada estática completa, foram conferidas 2.160 condições de teste, cada uma com 1.000 réplicas: 2.160.000 linhas, sem valores ausentes em `p_max`. Isso corresponde a 432.000 metacomunidades avaliadas por cinco testes. As simulações completas não foram reexecutadas, nem os modelos bayesianos reajustados nesta revisão.

Pontos fortes: modularização clara; testes contra implementações de referência; sementes controladas; dependências registradas; inclusão do observado no cálculo dos P-valores de permutação; cenários nulos distintos; transparência sobre mudanças no desenho; análises de sensibilidade já realizadas. Não foi encontrada uma falha aritmética evidente que invalide toda a simulação. Passar nos testes computacionais, porém, não demonstra a validade das interpretações científicas.

As referências ao código abaixo são relativas à raiz do repositório e descrevem a versão inspecionada nesta data.

## 2. Oito pontos da crítica

### 2.1. A hipótese nula de S2σ/S2h precisa ser demonstrada

**Localização:** `R/simulate_layerA.R`, funções `scenarios_layerA()` e `sim_species_traits()`; seção “Scenarios and estimands” do manuscrito.

S2σ é classificado como nulo porque os ótimos não têm sinal filogenético. Entretanto, a largura do nicho tem sinal e entra diretamente na resposta ao ambiente:

`L_ij = h_j exp[-(E_i - u_j)^2 / (2 sigma_j^2)]`.

Assim, clados podem diferir sistematicamente na forma da resposta ambiental, mesmo sem ótimos filogeneticamente estruturados. Isso pode produzir associação entre composição filogenética e ambiente, especialmente com respostas não lineares, limites do gradiente e conversão em incidência. “Ausência de sinal nos ótimos” não implica “ausência de associação entre composição filogenética e ambiente”.

A permutação das pontas também destrói a estrutura de σ presente no gerador. O max test não garante validade se os testes componentes não forem válidos sob suas respectivas hipóteses nulas.

As rejeições destacadas foram confirmadas: 0,356 para PCPS2 em S2σ, N = 200, incidência, γ baixo e desbalanço Yule; e até 0,218 para PCC em S2σ. Esses números demonstram taxas de rejeição elevadas nesse cenário, mas não demonstram, isoladamente, inflação do erro tipo I sob uma hipótese nula apropriada.

**Correção exigida:** definir matematicamente o alvo de cada teste. Distinguir uma associação ambiental filogenética em sentido amplo de conservadorismo especificamente dos ótimos. Até resolver isso, descrever S2σ como avaliação da “taxa de rejeição sob conservadorismo da largura do nicho”. Examinar também S2h e cenários de competição filogenética. Esta é uma questão de validade do estimando, não um erro aritmético demonstrado.

### 2.2. O contraste de competição mistura efeitos

**Localização:** `R/simulate_layerB.R`, `competition_matrix()` e `archetypes_layerB()`.

A função implementa competição `none`, `random` e `phylo`, mas os arquétipos usam apenas `none` e `phylo`. O contraste muda simultaneamente presença de competição interespecífica, estrutura filogenética e seleção das espécies sobreviventes. A mudança da árvore analisada, já reconhecida no manuscrito, acrescenta outro caminho.

**Correção exigida:** incluir competição aleatória com a mesma distribuição de coeficientes. Separar ausência/presença de competição de organização aleatória/filogenética. Sem esse controle, interpretar o efeito total do tratamento, evitando atribuí-lo exclusivamente à estrutura filogenética.

### 2.3. A versão espacial muda também a estocasticidade da dispersão

**Localização:** `R/simulate_layerB.R`, ramo espacial de `mc_dynamics()`.

A versão não espacial usa redistribuição multinomial e conserva emigrantes por espécie. A espacial sorteia imigrantes por Poisson entre destinos. A aproximação está declarada, mas sua inocuidade não foi demonstrada.

Condicionalmente aos emigrantes, a variância do total de imigrantes é zero na redistribuição conservativa e igual ao total esperado de emigrantes no esquema Poisson. Probabilidades pequenas por destino não eliminam essa diferença nem recuperam as covariâncias entre destinos.

**Correção exigida:** comparar com redistribuição multinomial por origem, conservando indivíduos. Até lá, diferenças entre versões não devem ser atribuídas exclusivamente à espacialidade ou à limitação da dispersão.

### 2.4. Ótimos com marginal fixa não isolam completamente o efeito de análise

**Localização:** `R/simulate_layerA.R`, opção `optima = "uniform"`; `analysis/36_optima_uniform.R`; interpretação no suplemento.

Os valores uniformes continuam sendo atribuídos pelos postos de um caráter browniano gerado na árvore. Mudar os tempos de ramificação ainda pode mudar quais espécies recebem ótimos semelhantes. A marginal é fixa, mas a estrutura conjunta entre espécies não é.

**Correção exigida:** limitar a conclusão atual a “o efeito não é explicado apenas pela marginal dos ótimos”. Para separar os caminhos, manter comunidades e caracteres fixos e modificar somente os comprimentos dos ramos na análise, preservando a correspondência entre espécies.

### 2.5. O resultado sobre PCPS depende da variante inferencial

**Localização:** `R/pcps_glm.R`; `tests/testthat/test-pcps_glm.R`; resultados de Procrustes no suplemento.

O suplemento mostra que Procrustes remove a rejeição elevada de PCPS2 em S2σ, embora torne o procedimento conservador. A escolha da configuração é, portanto, central. A documentação de `PCPS::pcps.sig` descreve ajuste procrustiano das ordenações nulas.

**Correção exigida:** distinguir testes individuais sem rotação conjunta, testes com Procrustes nos eixos selecionados e RDA dos eixos. Verificar não apenas estatísticas observadas, mas também estatísticas nulas sob permutações controladas. Não generalizar o comportamento de uma configuração para “PCPS” sem qualificação.

### 2.6. Recuperação de clados compara procedimentos específicos de extração

**Localização:** `R/clade_methods.R`, funções de escores de PCC, evoPCA e phylofactorization.

PCC é restrito aos quatro primeiros autovetores, com escolha pela associação ao ambiente. EvoPCA usa correlações das colunas de espécies, em vez dos carregamentos dos ramos. Phylofactorization avalia diretamente contrastes dos nós. Essas escolhas podem ser defensáveis, mas fazem parte do procedimento comparado e influenciam o ranking.

O objetivo binomial de phylofactorization também recompensa diferenças de prevalência entre grupos, mesmo sem associação com E; o comentário do código reconhece esse comportamento.

**Correção exigida:** testar sensibilidade ao número de eixos e ao objetivo de phylofactorization baseado apenas na interação com E. Delimitar as recomendações ao procedimento completo. O controle com deslocamento zero é útil, mas não resolve sozinho a comparabilidade dos objetivos.

### 2.7. Modelos formais não recuperam automaticamente o pareamento

**Localização:** `R/formal_models.R`; `analysis/18_formal_models.R`, `20_formal_models_layerB.R`, `31_paired_contrasts.R` e `35_mean_structure.R`.

Árvores são compartilhadas entre cenários estáticos, e métodos são aplicados às mesmas comunidades. A agregação seguida de modelos beta-binomiais não representa integralmente essas dependências. A beta-binomial acomoda dispersão adicional, mas não recupera o pareamento perdido.

A análise pareada existente é um ponto forte e deve ser estendida aos contrastes principais. Na sensibilidade com efeito por condição, diferenças medianas pequenas convivem com diferenças máximas de aproximadamente 0,068 para taxas de tipo I estáticas e 0,117 para poder estático. As condições responsáveis por essas discrepâncias precisam ser examinadas.

**Correção exigida:** apresentar taxas brutas e contrastes com incerteza de Monte Carlo como evidência principal; usar modelos para síntese. Verificar as conclusões nas condições com discrepâncias relevantes, sem depender apenas de convergência MCMC ou cobertura preditiva global.

### 2.8. Desbalanço e dinâmica precisam de conclusões delimitadas

**Localização:** `analysis/30_h1_equivalence.R`; `data/derived/models/h1_equivalence.csv`; `mc_dynamics()`.

A equivalência atingiu probabilidade posterior acima de 0,95 em 297 de 300 contrastes, mas não em três; o mínimo foi 0,596. Isso apoia efeitos geralmente pequenos no desenho estudado, não ausência universal de efeito.

Também não foi identificada uma demonstração suficiente de robustez da escolha do estado final em 1.000 passos. Sem recolonização externa após o período de semeadura, extinções são absorventes; não se deve exigir um equilíbrio permanente não nulo. É necessário avaliar robustez ao horizonte temporal e à inicialização. Os arquétipos devem ser sustentados pelo comportamento realizado, incluindo evidências de fontes e sumidouros, e não apenas pelo valor do parâmetro de dispersão.

**Correção exigida:** qualificar a equivalência e testar alguns horizontes temporais e condições iniciais.

## 3. Plano de correção

### 3.1. Ordem de trabalho

Resolver primeiro os pontos **1, 5 e 7**, antes de novas rodadas grandes. Eles podem mudar a interpretação e a seleção das simulações adicionais. Depois executar controles e diagnósticos em subconjuntos, seguidos por rodadas confirmatórias. Não há necessidade demonstrada de descartar ou repetir todo o estudo.

### 3.2. Ponto 1 — hipóteses e estimandos

1. Preparar uma tabela por teste e cenário: pergunta científica, hipótese estatística, dependências do gerador e justificativa da permutação.
2. Demonstrar analítica ou numericamente como clados com diferentes larguras de nicho respondem ao gradiente.
3. Decidir se S2σ/S2h avaliam robustez a outros componentes do nicho ou erro tipo I sob um alvo explicitamente restrito aos ótimos.
4. Atualizar classificações, tabelas e interpretação. Não inventar uma permutação restrita sem demonstrar sua validade.

**Conclusão:** separação explícita entre erro tipo I, sensibilidade a componentes alternativos do nicho e poder. Desenvolver um novo teste condicional, caso necessário, constitui ampliação do escopo e não está incluído nos prazos abaixo.

### 3.3. Ponto 2 — controle de competição

Adicionar `random` ao desenho. Começar com N = 100, células Yule e γ baixo com desbalanço Yule, dispersão baixa/alta e cenários selecionados após o ponto 1. Compartilhar os valores dos coeficientes entre versões aleatória e filogenética, alterando sua atribuição aos pares de espécies.

Verificar riqueza realizada e separar efeito total do tratamento de análises com riqueza aproximadamente controlada. Não tratar automaticamente a riqueza sobrevivente como covariável externa: ela é consequência do tratamento.

**Conclusão:** contrastes identificáveis entre presença de competição e sua organização filogenética.

### 3.4. Ponto 3 — dispersão espacial conservativa

Implementar multinomiais por origem. Verificar conservação por espécie, frequências de destino e casos simples. Comparar Poisson e multinomial quanto a riqueza, abundância, ocupação e rejeição dos testes.

Se o efeito for desprezível para as conclusões, documentar a sensibilidade. Se for relevante, adotar a versão conservativa e repetir as condições espaciais necessárias.

**Conclusão:** impacto da aproximação quantificado e escolha do algoritmo justificada.

### 3.5. Ponto 4 — experimento de mecanismo de γ

Usar a mesma topologia, espécies correspondentes e duas cronologias, Yule e γ baixo. Cruzar árvore de geração dos caracteres e árvore usada na análise:

| Geração | Análise |
|---|---|
| Yule | Yule |
| Yule | γ baixo |
| γ baixo | Yule |
| γ baixo | γ baixo |

Cada comunidade gerada é analisada com as duas cronologias. O experimento é um diagnóstico de mecanismo: analisar com árvore diferente da geradora não representa uma avaliação usual com filogenia corretamente especificada.

**Conclusão:** separar o que muda por geração dos caracteres do que muda apenas na análise.

### 3.6. Ponto 5 — variantes de PCPS

Comparar configurações de eixo 1, eixo 2 e ambos os eixos com `PCPS::pcps.sig`. Verificar estatísticas observadas e, quando possível, as mesmas permutações. Padronizar nomes dos procedimentos e trazer Procrustes e RDA aos resultados principais.

**Conclusão:** equivalência computacional documentada e recomendações específicas por configuração.

### 3.7. Ponto 6 — recuperação de clados

Testar 2, 4 e 8 eixos, inicialmente com N = 100. Comparar o objetivo binomial atual com um objetivo que isole a interação com E. Manter deslocamento zero e acrescentar um controle com prevalências distintas entre clados, mas sem associação sistemática com E.

**Conclusão:** estabilidade do ranking avaliada e limites das regras de extração explícitos.

### 3.8. Ponto 7 — taxas brutas e contrastes

Reanalisar primeiro os resultados existentes. Calcular taxas, intervalos de Monte Carlo e contrastes por réplica, respeitando a dependência real: árvore compartilhada entre cenários estáticos e comunidade compartilhada entre métodos. Não presumir pareamento entre tratamentos dinâmicos apenas por terem o mesmo número de réplica.

Examinar as condições com maiores diferenças entre modelos e verificar os contrastes que sustentam as recomendações. Usar os modelos para síntese, sem substituir a evidência direta das taxas.

**Conclusão:** conclusões principais sustentadas por taxas brutas e contrastes que preservam o desenho.

### 3.9. Ponto 8 — tempo, inicialização e equivalência

Guardar estados e diagnósticos das mesmas trajetórias em 500, 1.000 e 2.000 passos, com mais de uma inicialização. Comparar riqueza, composição, abundância e desempenho inferencial. Registrar evidências de fontes e sumidouros para justificar mass effects.

Apresentar os três contrastes sem equivalência apoiada e substituir afirmações universais por conclusões restritas às condições estudadas.

**Conclusão:** robustez temporal documentada e alcance das inferências corretamente delimitado.

## 4. Estimativas de esforço e computação

Estimativas preliminares de esforço técnico, baseadas na inspeção do código e nos tempos registrados no projeto; não são um compromisso de prazo nem um benchmark das novas implementações.

| Atividade | Horas de trabalho |
|---|---:|
| 1. Hipóteses e estimandos | 8–14 |
| 2. Controle de competição | 4–7 |
| 3. Dispersão espacial conservativa | 5–8 |
| 4. Separação dos caminhos de γ | 6–10 |
| 5. Validação das variantes de PCPS | 4–7 |
| 6. Sensibilidade na recuperação de clados | 6–10 |
| 7. Contrastes e modelos | 6–10 |
| 8. Equivalência e robustez temporal | 5–9 |
| Subtotal dos oito pontos | **44–75** |
| Integração de texto, figuras, tabelas e verificação final | **6–10** |
| **Total** | **50–85** |

O total corresponde aproximadamente a **7–12 dias úteis** de trabalho concentrado. A integração já está incluída nas 50–85 horas; esta explicitação evita a ambiguidade da estimativa apresentada inicialmente na conversa.

Cronograma indicativo:

1. **Dias 1–2:** ponto 1, auditoria de PCPS e contrastes existentes.
2. **Dias 3–5:** implementação dos controles e pilotos dos pontos 2, 3, 4, 6 e 8.
3. **Dias 6–8:** rodadas confirmatórias selecionadas e análise.
4. **Dias 9–12:** integração, verificação das conclusões e atualização do manuscrito.

Começar com **50–100 réplicas por condição** para diagnóstico de comportamento e custo. Usar **1.000 réplicas nas condições confirmatórias** destinadas a estimar taxas finais. O piloto não substitui a precisão da rodada confirmatória.

Reservar preliminarmente **12–48 horas de computação em infraestrutura semelhante à AWS usada no projeto**, para os subconjuntos propostos. Esse intervalo precisa de benchmark, especialmente para dispersão conservativa e horizontes temporais maiores. Uma reexecução ampla da camada dinâmica pode ultrapassá-lo. Tempo de máquina não equivale a horas de trabalho técnico; parte pode ocorrer durante outras atividades.

A decisão sobre o ponto 1 é a maior incerteza do plano. Se a solução for sobretudo redefinir estimandos e reformular conclusões, grande parte dos resultados pode ser preservada. Se envolver desenvolver um novo teste condicional ou repetir extensamente a camada dinâmica, será necessário reestimar o escopo.

## 5. Avaliação editorial final

Preservar o estudo, mas resolver a definição das hipóteses antes da submissão. O resultado mais defensável atualmente é que a configuração dos testes, os tempos de ramificação e os processos de montagem alteram suas taxas de rejeição e recuperação de clados. A afirmação de erro tipo I elevado sob hipótese nula válida exige demonstração adicional.

## 6. Fontes metodológicas consultadas

- ter Braak: [Fourth-corner correlation is a score test statistic in a log-linear trait–environment model that is useful in permutation testing](https://doi.org/10.1007/s10651-017-0368-0). Fundamentação do esquema max e suas hipóteses.
- [Documentação do pacote PCPS](https://vanderleidebastiani.r-universe.dev/PCPS/doc/manual.html). Descrição do ajuste procrustiano em `pcps.sig`.
- Thompson et al. (2020): [A process-based metacommunity framework linking local and regional scale community ecology](https://doi.org/10.1111/ele.13568). Referência do modelo dinâmico e da interpretação dos processos.

As críticas sobre o repositório derivam da inspeção do código, do manuscrito e dos resultados locais. As hipóteses de mecanismo e as correções propostas não devem ser confundidas com resultados já demonstrados por novas simulações.
