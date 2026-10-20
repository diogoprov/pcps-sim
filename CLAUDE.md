# pcps-sim — instruções do projeto

- **Leia `docs/MEMORIA_PROJETO.md` e `docs/PLANO_SIMULACAO.md` no início de cada sessão.**
  Eles registram o histórico, as decisões, as referências verificadas e os próximos passos.
- Autor único: Diogo B. Provete. Claude trabalha como pós-doc.
- **Faça commit a cada etapa concluída** (pedido do Diogo, 2026-10-16). Mensagens em
  português, curtas e descritivas.
- Mantenha `docs/MEMORIA_PROJETO.md` atualizado depois de decisões ou resultados novos.
- Ao fim de cada sessão, acrescente uma entrada datada em `docs/AI_INTERACTION_LOG.md`
  (em inglês; o que o assistente fez, o que o Diogo decidiu, correções feitas após
  revisão). Pedido do Diogo, 2026-10-19, para a declaração de uso de IA.
- `references/` (PDFs) fica fora do git por causa do copyright. Os dados brutos em
  `data/raw/` são imutáveis.
- `legacy/` guarda o material de 2012–2016, só para consulta. Não edite.
- Não baixe nada do Box nem force a sincronização do Dropbox. Os originais grandes
  (workspaces `.RData`) ficam no Dropbox.
- Stack: R (tidyverse, `|>`), `targets`, `renv`; `brms` com priors explícitos e
  diagnósticos; Quarto para o manuscrito, que será em inglês.
- A execução pesada vai para a AWS. Antes de rodar algo demorado, confirme o plano.
