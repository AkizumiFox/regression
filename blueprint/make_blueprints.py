"""Write authoring/blueprints/chNN.md for every chapter from blueprint/book.yaml.

A blueprint is the brief a writer (or drafting agent) starts from: the chapter's place in
the book, its planned sections, the source sections its coverage is checked against, its
page budget and its prerequisites from the DAG. Written chapters keep their blueprint as a
record. Regenerate after editing book.yaml.
"""
import json
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
BOOK_NAMES = {"R": "Rencher and Schaalje", "S": "Seber and Lee", "C": "Christensen",
              "N": "Sen and Srivastava", "A": "Agresti", "F": "Fahrmeir et al."}


def prerequisites(book, dag):
    cut = book.get("cut", {})
    uses = defaultdict(set)
    for ch in book["chapters"]:
        for node, _ in ch["sources"]:
            uses[node].add(ch["num"])
    first = {n: min(c) for n, c in uses.items()}
    req = defaultdict(set)
    for a, b in dag["hard"]:
        a, b = cut.get(a, a), cut.get(b, b)
        if a is None or b is None or a not in first:
            continue
        for c in uses.get(b, ()):
            if first[a] < c:
                req[c].add(first[a])
    return req


def main():
    book = yaml.safe_load((ROOT / "blueprint/book.yaml").read_text())
    dag = json.loads((ROOT / "blueprint/dag.json").read_text())
    written = {int(p.name[2:4]) for p in (ROOT / "src").glob("ch[0-9][0-9]-*") if p.is_dir()}
    req = prerequisites(book, dag)
    parts = {p["num"]: p for p in book["parts"]}
    titles = {c["num"]: c["title"] for c in book["chapters"]}
    out = ROOT / "authoring" / "blueprints"
    out.mkdir(parents=True, exist_ok=True)
    for ch in book["chapters"]:
        part = parts[ch["part"]]
        lines = [f"# Chapter {ch['num']}: {ch['title']}", "",
                 f"Directory: `src/ch{ch['num']:02d}-{ch['slug']}/`. "
                 f"State: {'written' if ch['num'] in written else 'not started'}.", "",
                 f"**Part {part['num']}, {part['title']}.** {part.get('drops') or ''} {part['model']}".replace("  ", " "), ""]
        if ch.get("new"):
            lines += ["**New chapter.** None of the six planning sources treats this topic in depth; it is original writing.", ""]
        if ch.get("note"):
            lines += [f"*{ch['note']}*", ""]
        lines += [f"**Page budget:** about {ch['pages']} source-book pages "
                  f"(about {round(ch['pages'] * 0.6)} A4 pages of this book's PDF).", ""]
        pre = sorted(req.get(ch["num"], []))
        lines += ["**Prerequisites (from the DAG):** " + (", ".join(f"Chapter {p} ({titles[p]})" for p in pre) or "none") + ".", ""]
        lines += ["## Planned sections", ""] + [f"{i}. {s}" for i, s in enumerate(ch["outline"], 1)] + [""]
        if ch["sources"]:
            lines += ["## Coverage is checked against", ""]
            for node, secs in ch["sources"]:
                n = dag["nodes"][node]
                where = "" if secs == "all" else f", {secs}"
                lines.append(f"- {BOOK_NAMES[n['book']]}, ch. {n['ch']} ({n['title']}){where}")
            lines.append("")
        (out / f"ch{ch['num']:02d}.md").write_text("\n".join(lines))
    print(f"wrote {len(book['chapters'])} blueprints ({len(written)} chapters written)")


if __name__ == "__main__":
    main()
