#!/usr/bin/env bash
# Anonymized copy of the repository for double-anonymous peer review
# (Methods in Ecology and Evolution): git HEAD without the author's internal
# notes, title page and manuscript sources, with names, the GitHub account and
# the licence holder removed. Fails if an identifying string is left.
# Usage: tools/make_review_bundle.sh   (writes review_bundle/pcps-sim-review.zip)
set -euo pipefail

root="$(git rev-parse --show-toplevel)"
out_dir="$root/review_bundle"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
dest="$work/pcps-sim"
mkdir -p "$dest" "$out_dir"
git -C "$root" archive HEAD | tar -x -C "$dest"

# internal notes, AI log, citation metadata and manuscript sources (the
# manuscript is submitted separately); the figures and SI tables stay
rm -f "$dest"/CLAUDE.md "$dest"/CITATION.cff \
      "$dest"/docs/MEMORIA_PROJETO.md "$dest"/docs/AI_INTERACTION_LOG.md \
      "$dest"/docs/resumo_projeto.qmd "$dest"/docs/GUIAS_RELATO_MEE.md \
      "$dest"/manuscript/*.qmd "$dest"/manuscript/*.docx "$dest"/manuscript/*.py
rm -rf "$dest"/docs/lit_review "$dest"/tools
# internal working notes in Portuguese (first person, decision dates)
rm -f "$dest"/docs/PLANO_SIMULACAO.md "$dest"/docs/ESPEC_PCC.md
# scripts that read the third-party tadpole data (not distributed); they list
# the data providers and their localities. Their outputs are included.
rm -f "$dest"/analysis/21_regional_pools.R "$dest"/analysis/23_case_study.R
cp "$root/tools/review_bundle_README.md" "$dest/README_REVIEW.md"

# names in the remaining files (BSD sed, as on macOS)
scrub() {
  sed -i '' -E \
    -e 's/Copyright \(c\) ([0-9]{4}) Diogo B\. Provete/Copyright (c) \1 Anonymous (for peer review)/' \
    -e '/^Author: Diogo B\. Provete\.?$/d' \
    -e 's#https?://github\.com/diogoprov/pcps-sim#[repository URL withheld for review]#g' \
    -e "s/the author's //g; s/the author confirmed/it was confirmed/g; s/shared with the author/shared with the authors/g" \
    -e '/brazilian-tadpole-knowledge-gaps/d; s/Copied on ([0-9-]+) from compendium$/Copied on \1 from a local copy./' \
    -e 's/\(Diogo, ([0-9-]+):/(decision of \1:/g' \
    -e 's/Proposta do Diogo/Proposta/g; s/aprovado pelo Diogo/aprovado/g' \
    -e 's/Autoria: D\. B\. Provete \(autor único\)\. //; s/Autoria: Provete como autor único\./Autoria: autor único./' \
    "$1"
}
while IFS= read -r -d '' f; do scrub "$f"; done < <(find "$dest" -type f \
  \( -name '*.R' -o -name '*.md' -o -name '*.sh' -o -name 'LICENSE*' -o -name 'Dockerfile' \) -print0)

# self-citations stay (third person, as in the manuscript); anything else fails
identifying="diogo|diogoprov|dbprovete|borges|author's|brazilian-tadpole|mainara|crivelari|rogeio|decio|goias|/Users/"
left="$(grep -rIl -i -E "$identifying" "$dest" || true)"
if [[ -n "$left" ]]; then
  echo "Identifying strings left in:" >&2
  grep -rIn -i -E "$identifying" "$dest" >&2
  exit 1
fi
# a line is accepted only if every "Provete" on it is a citation of 2014 or 2018
grep -rIn -i 'provete' "$dest" | sed -E 's/Provete( et al\.)?,? \(?(2014|2018)//gI' |
  grep -i 'provete' >&2 && {
  echo "Unexpected 'Provete' above" >&2; exit 1; } || true

rm -f "$out_dir/pcps-sim-review.zip"
(cd "$work" && zip -qr "$out_dir/pcps-sim-review.zip" pcps-sim)
echo "Written: $out_dir/pcps-sim-review.zip ($(du -h "$out_dir/pcps-sim-review.zip" | cut -f1))"
