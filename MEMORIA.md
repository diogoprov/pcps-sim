# MEMORIA.md — pcps-sim

Consolidação, em 23/10/2026, da memória do Claude Code que estava na pasta antiga do projeto
na Box (`~/.claude/projects/-Users-diogoborgesprovete-Library-CloudStorage-Box-Box-PPG-Eco-Evo-UFG-Capitulo-1---Simula--es-PCPS/memory/`).
Aquela pasta está desatualizada: **este repositório é o projeto**.

O histórico cronológico completo e o plano continuam em
[`docs/MEMORIA_PROJETO.md`](docs/MEMORIA_PROJETO.md) e
[`docs/PLANO_SIMULACAO.md`](docs/PLANO_SIMULACAO.md). Aqui ficam as decisões que valem e o
porquê delas.

## Origem e escopo (avaliação de 16/10/2026)

O capítulo (manuscrito de março/2016, coautores Debastiani, Duarte e Melo) **não se sobrepõe**
a Duarte et al. 2016 (*MEE*) no núcleo — o efeito do formato da árvore; só parcialmente em
tamanho do pool e sinal filogenético. O desenho antigo foi considerado datado.

Lido o PDF publicado e o Appendix S2 (árvores Yule puras, γ e Colless não controlados): a
frase do rascunho sobre topologia foi removida da versão publicada. Varredura no OpenAlex
(~420 citantes) confirmou a lacuna — ninguém testou PCPS ou outras ordenações filogenéticas
sob formato de árvore controlado.

Decisão do Diogo: o eixo principal passa a ser a **validade da inferência sobre P**, e a
comparação de métodos vira o objetivo 2. O `max test` já foi usado por Duarte et al. 2016 e
por Peres-Neto & Kembel 2015, então a novidade está nas **condições em que ele falha**. O PCC
será reimplementado a partir do PDF, sem pedir o código ao Pedro.

## Autoria (decisão do Diogo, 19/10/2026)

A ideia da simulação é do Diogo (o PDF com anotações à mão é dele); V. J. Debastiani ajudou na
implementação. **Debastiani e Duarte vão nos agradecimentos, não como coautores, e a ideia
original não deve ser atribuída a Debastiani.**

A versão de 2016 listava Debastiani como coautor de mesma contribuição — isso mudou. Autor
único hoje: Diogo B. Provete.

## Declaração de uso de IA

Todo enunciado em manuscrito ou preprint segue o formato **AIdIT v0.1** (Drobniak et al. 2026,
*Research Integrity and Peer Review* 11:57, doi:10.1186/s41073-026-00230-1), nunca um
parágrafo livre estilo COPE.

Montar a partir do template do app oficial (github.com/szymekdr/AI_reporting, `code/app.R`):
"We used the following AI engines in the present study: <nome (v.; data/hora)>; … Area(s) of
generative AI usage: …", mais as frases adicionais (verificado/aprovado; conformidade ética;
prompts sob pedido ou em Supplementary Materials) e um comentário livre se houver. Vai como
seção "AI disclosure statement" do texto principal, citando Drobniak et al. 2026.

**Pedir ao Diogo para confirmar motores, áreas, datas e frases adicionais** — não inventar. O
texto do template é usado literalmente; conferir a fonte em vez de adivinhar.

## Trabalho de código

Explicitar suposições e perguntar quando ambíguo; mínimo de código que resolve; mudanças
cirúrgicas (cada linha alterada ligada ao pedido); critério de verificação definido antes de
implementar. Depois de um bloco coerente de mudanças, rodar revisão em duas frentes —
aderência a padrões e aderência ao spec — contra o último commit revisado, e corrigir o que
aparecer **antes** de dizer que está pronto. Pedido explícito do Diogo em 17/10/2026.

## Ambiente e arquivos

- Não baixar do Box nem forçar sincronização: a cópia local foi para
  `~/Downloads/Capitulo 1 - Simulações PCPS/` e os PDFs avulsos para `~/Downloads/`.
  A sincronização do Box trava.
- A execução pesada vai para a AWS; confirmar o plano antes de rodar algo demorado.
- Registrar no fim de cada sessão uma entrada em
  [`docs/AI_INTERACTION_LOG.md`](docs/AI_INTERACTION_LOG.md) — é a base da declaração AIdIT.
