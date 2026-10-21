#!/usr/bin/env python3
"""Submission copy of the cover letter with the reviewers' e-mails, which are
kept out of the public repository (manuscript/private/reviewer_emails.tsv:
surname<TAB>e-mail). Writes manuscript/cover_letter_submission.docx (gitignored).
Usage: tools/render_cover_letter.py
"""
import subprocess
from pathlib import Path

root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
ms = root / "manuscript"
text = (ms / "cover_letter.qmd").read_text()
for line in (ms / "private" / "reviewer_emails.tsv").read_text().splitlines():
    surname, email = line.split("\t")
    text = text.replace("REVIEWER_EMAIL_" + surname.replace(" ", "_"), email)
qmd = ms / "cover_letter_submission.qmd"
qmd.write_text(text)
try:
    subprocess.run(["quarto", "render", qmd.name, "--to", "docx", "--quiet"], cwd=ms, check=True)
finally:
    qmd.unlink()
print(f"Written: {ms / 'cover_letter_submission.docx'}")
