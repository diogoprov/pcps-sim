"""Reference .docx for the manuscript: pandoc's default styles plus continuous
line numbers, page numbers centred in the footer and double line spacing.

Usage: python3 manuscript/make_reference_docx.py
(writes manuscript/reference.docx; needs quarto on the PATH)
"""
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

OUT = Path(__file__).with_name("reference.docx")

FOOTER = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:p>
    <w:pPr><w:jc w:val="center"/></w:pPr>
    <w:r><w:fldChar w:fldCharType="begin"/></w:r>
    <w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>
    <w:r><w:fldChar w:fldCharType="separate"/></w:r>
    <w:r><w:t>1</w:t></w:r>
    <w:r><w:fldChar w:fldCharType="end"/></w:r>
  </w:p>
</w:ftr>
"""


def main():
    tmp = Path(tempfile.mkdtemp())
    default = tmp / "default.docx"
    with open(default, "wb") as fh:
        subprocess.run(["quarto", "pandoc", "--print-default-data-file", "reference.docx"],
                       stdout=fh, check=True)
    with zipfile.ZipFile(default) as z:
        parts = {name: z.read(name) for name in z.namelist()}

    # footer part, its relationship and content type
    parts["word/footer1.xml"] = FOOTER.encode()
    rels = parts["word/_rels/document.xml.rels"].decode()
    rels = rels.replace(
        "</Relationships>",
        '<Relationship Id="rIdFooter1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
        'Target="footer1.xml"/></Relationships>')
    parts["word/_rels/document.xml.rels"] = rels.encode()
    types = parts["[Content_Types].xml"].decode()
    types = types.replace(
        "</Types>",
        '<Override PartName="/word/footer1.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/>'
        "</Types>")
    parts["[Content_Types].xml"] = types.encode()

    # section properties: footer reference and continuous line numbering
    doc = parts["word/document.xml"].decode()
    # (schema order: footer references first, line numbering after the margins)
    doc = re.sub(r"<w:sectPr([^>]*)>",
                 r'<w:sectPr\1><w:footerReference w:type="default" r:id="rIdFooter1"/>',
                 doc, count=1)
    # pandoc's default section has only footnotePr: add A4, 2.5 cm margins and
    # the line numbering after it
    doc = re.sub(r"(</w:footnotePr>)",
                 r'\1<w:pgSz w:w="11906" w:h="16838"/>'
                 r'<w:pgMar w:top="1418" w:right="1418" w:bottom="1418" w:left="1418" '
                 r'w:header="709" w:footer="709" w:gutter="0"/>'
                 r'<w:lnNumType w:countBy="1" w:restart="continuous"/>', doc, count=1)
    assert "lnNumType" in doc
    if 'xmlns:r=' not in doc[:2000]:
        doc = doc.replace(
            "<w:document ",
            '<w:document xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" ',
            1)
    parts["word/document.xml"] = doc.encode()

    # double spacing for every paragraph by default
    styles = parts["word/styles.xml"].decode()
    styles = re.sub(r"<w:pPrDefault>.*?</w:pPrDefault>",
                    '<w:pPrDefault><w:pPr><w:spacing w:after="0" w:line="480" '
                    'w:lineRule="auto"/></w:pPr></w:pPrDefault>',
                    styles, count=1, flags=re.S)
    parts["word/styles.xml"] = styles.encode()

    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in parts.items():
            z.writestr(name, data)
    shutil.rmtree(tmp)
    print(f"Written: {OUT}")


if __name__ == "__main__":
    main()
