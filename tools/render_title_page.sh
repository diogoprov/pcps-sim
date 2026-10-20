#!/usr/bin/env bash
# Submission copy of the title page with the corresponding e-mail, which is
# kept out of the public repository. Writes manuscript/title_page_submission.docx
# (gitignored).
# Usage: tools/render_title_page.sh <e-mail>
set -euo pipefail
email="${1:?usage: tools/render_title_page.sh <e-mail>}"
root="$(git rev-parse --show-toplevel)"
cd "$root/manuscript"
sed "s/CORRESPONDING_EMAIL/$email/" title_page.qmd > title_page_submission.qmd
trap 'rm -f title_page_submission.qmd' EXIT
quarto render title_page_submission.qmd --to docx --quiet
echo "Written: $root/manuscript/title_page_submission.docx"
