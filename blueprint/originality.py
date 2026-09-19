"""Flag word sequences shared between a written chapter (src/chNN-*/*.md) and its source chapters.

    python blueprint/originality.py 06 [--n 8] [--extra C16]

Sources are the DAG nodes listed for the chapter in book.yaml (plus --extra nodes).
Text is extracted with pdftotext; scanned PDFs without a text layer are reported as
unchecked. Shared n-grams made only of mathematical boilerplate are ignored.
"""
import argparse
import re
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "regression"
BOOK_DIRS = {
    "R": "linear-models-in-statistics-rencher-schaalje-2e-chapters",
    "S": "linear-regression-analysis-seber-lee-2e-chapters",
    "C": "plane-answers-to-complex-questions-linear-models-christensen-5e-chapters",
    "N": "regression-analysis-theory-methods-and-applications-sen-srivastava-chapters",
    "A": "foundations-of-linear-and-generalized-linear-models-agresti-chapters",
    "F": "regression-models-methods-and-applications-2e-chapters",
}
COMMON = {"the", "of", "a", "is", "and", "to", "in", "that", "for", "be", "by", "with", "if", "then", "on", "as", "an", "are", "this", "it"}


def source_pdf(node: str) -> Path:
    book, rest = node[0], node[1:]
    d = SRC / BOOK_DIRS[book]
    if rest.isdigit():
        return next(d.glob(f"{int(rest):02d}-*.pdf"))
    return next(d.glob(f"*appendix-{rest.lower()}-*.pdf"))


def words(text: str):
    return re.findall(r"[a-z]+(?:'[a-z]+)?", text.lower())


def md_to_text(md: str) -> str:
    md = re.sub(r"```.*?```", " ", md, flags=re.S)
    md = re.sub(r"\\\(.*?\\\)", " MATH ", md, flags=re.S)
    md = re.sub(r"\\\[.*?\\\](\{#[^}]*\})?", " ", md, flags=re.S)
    md = re.sub(r"\\begin\{center\}.*?\\end\{center\}", " ", md, flags=re.S)
    md = re.sub(r"^:::.*$", " ", md, flags=re.M)
    md = re.sub(r"@[a-z]+-[A-Za-z0-9-]+", " ", md)
    md = re.sub(r"\]\([^)]*\)", "]", md)
    return md


def tex_to_text(tex: str) -> str:
    tex = re.sub(r"(?<!\\)%.*", "", tex)
    tex = re.sub(r"\\begin\{(equation|align|lstlisting)\*?\}.*?\\end\{\1\*?\}", " ", tex, flags=re.S)
    tex = re.sub(r"\\\[.*?\\\]", " ", tex, flags=re.S)
    tex = re.sub(r"\$[^$]*\$", " MATH ", tex)
    tex = re.sub(r"\\(cref|Cref|label|ref|eqref|pyfile|gen|input|inputgen|includegraphics|textcite|parencite|Textcite)(\[[^\]]*\])?\{[^}]*\}", " ", tex)
    tex = re.sub(r"\\[a-zA-Z]+\*?", " ", tex)
    return re.sub(r"[{}\[\]]", " ", tex)


def ngrams(ws, n):
    return {tuple(ws[i:i + n]) for i in range(len(ws) - n + 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("chapter")
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--extra", nargs="*", default=[])
    a = ap.parse_args()
    book = yaml.safe_load((ROOT / "blueprint/book.yaml").read_text())
    ch = next(c for c in book["chapters"] if c["num"] == int(a.chapter))
    files = sorted((ROOT / "src").glob(f"ch{int(a.chapter):02d}-*/*.md"))
    mine = words(md_to_text("\n".join(f.read_text() for f in files)))
    mine_ng = ngrams(mine, a.n)
    print(f"chapter {a.chapter}: {len(mine)} words of prose checked, n = {a.n}")
    nodes = [s[0] for s in ch["sources"]] + a.extra
    total = 0
    for node in nodes:
        pdf = source_pdf(node)
        txt = subprocess.run(["pdftotext", str(pdf), "-"], capture_output=True, text=True).stdout
        src = words(txt)
        if len(src) < 200:
            print(f"  {node}: UNCHECKED (no text layer: {pdf.name})")
            continue
        shared = [g for g in mine_ng & ngrams(src, a.n) if sum(w not in COMMON and w != "math" for w in g) >= 3]
        total += len(shared)
        print(f"  {node}: {len(src)} source words, {len(shared)} shared {a.n}-grams")
        for g in sorted(shared)[:40]:
            print("     · " + " ".join(g))
    print(f"total shared: {total}")


if __name__ == "__main__":
    main()
