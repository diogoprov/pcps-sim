#!/usr/bin/env python3
"""Preprint: a single docx with the title page and the main text, not
anonymized (authors, affiliations, statements, and the public repository
instead of the anonymous archive for review). Built from
manuscript/methods_results.qmd and manuscript/title_page.qmd, so it can be
rebuilt after every revision. The e-mail is not kept in the repository: the
output (manuscript/preprint/, gitignored) is written with the e-mail given.
Usage: tools/build_preprint.py <e-mail>
"""
import re
import subprocess
import sys
from pathlib import Path

email = sys.argv[1]
root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
ms = root / "manuscript"
public = ("https://github.com/diogoprov/pcps-sim, archived at Zenodo "
          "(https://doi.org/10.5281/zenodo.23500104)")


def body(text):
    """Text after the YAML header, without HTML comments."""
    text = text.split("---", 2)[2]
    return re.sub(r"<!--.*?-->", "", text, flags=re.S).strip()


main = (ms / "methods_results.qmd").read_text()
title = re.search(r'^title: "(.*)"$', main, flags=re.M).group(1)
yaml = main.split("---", 2)[1]
yaml = re.sub(r'^subtitle: .*\n', 'subtitle: "Preprint"\n', yaml, flags=re.M)
main = body(main)
tp = body((ms / "title_page.qmd").read_text())

# title-page block (authors to word count) and the statements (# sections)
head, statements = tp.split("# Acknowledgements", 1)
head = head.replace("CORRESPONDING_EMAIL", email)
head = re.sub(r"^(Word count|Article type):.*$\n?", "", head, flags=re.M)
statements = "# Acknowledgements" + statements
statements = re.sub(r" For peer review, an anonymized copy is available at \S+\.", "", statements)

# public repository instead of the anonymous archive
main = re.sub(r"Data and code availability for review: (.*?) are available in an anonymized "
              r"archive for review \([^)]*\); the public repository and its Zenodo archive "
              r"will be cited after review\.",
              lambda m: f"Data and code availability: {m.group(1)} are available at {public}.", main)
main = re.sub(r"Code is available in an anonymized archive for review \([^)]*\),",
              f"Code is available at {public},", main)
assert "anonymized" not in main, "an anonymized-archive statement was left in the main text"

# statements go before the references
text, refs = main.split("# References (verified)", 1)
out = f"---{yaml}---\n\n{head.strip()}\n\n{text.strip()}\n\n{statements.strip()}\n\n# References{refs}"
qmd = ms / "preprint.qmd"
qmd.write_text(out)
try:
    subprocess.run(["quarto", "render", qmd.name, "--to", "docx", "--output-dir", "preprint"],
                   cwd=ms, check=True)
finally:
    qmd.unlink()
print(f"Written: {ms / 'preprint' / 'preprint.docx'} ({title})")
