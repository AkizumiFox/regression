#!/usr/bin/env python3
"""List numbers in a chapter's prose that no script records, and recorded values the prose never quotes.

    python tools/check_prose_numbers.py 13 [--all]

tools/check_numbers.py guarantees that the scripts still compute the recorded values
(code/_generated/chNN/*.json); it cannot see whether the prose quotes them correctly. This
lists every decimal with three or more digits after the point in the chapter's Markdown
(outside code cells) that is not a recorded value of that chapter. Each one must be a hand
calculation shown in the text, a standard constant, or a typo. With --all it also lists
recorded decimals the prose never quotes (usually values used only in figures).
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    ch = f"ch{int(sys.argv[1]):02d}"
    recorded = {}
    for f in (ROOT / "code" / "_generated" / ch).glob("*.json"):
        recorded.update(json.loads(f.read_text()))
    clean = lambda v: v.replace("\\ensuremath{-}", "-").lstrip("-")
    values = {clean(v) for v in recorded.values()}
    unknown = []
    text_all = ""
    for page in sorted((ROOT / "src").glob(f"{ch}-*/*.md")):
        text = re.sub(r"```.*?```", "", page.read_text(), flags=re.S)
        text_all += text
        for lineno, line in enumerate(text.splitlines(), 1):
            for m in re.finditer(r"(?<![\w.])(\d+\.\d{3,})(?!\d)", line):
                if m.group(1) not in values:
                    unknown.append(f"{page.name}:{lineno}: {m.group(1)}   {line.strip()[:110]}")
    print(f"{ch}: {len(unknown)} prose decimals not recorded by a script")
    for u in unknown:
        print("  " + u)
    if "--all" in sys.argv:
        unquoted = sorted(k for k, v in recorded.items()
                          if re.fullmatch(r"\d+\.\d+", clean(v)) and clean(v) not in text_all)
        print(f"{len(unquoted)} recorded decimals not quoted in the prose:")
        for k in unquoted:
            print(f"  {k} = {recorded[k]}")


if __name__ == "__main__":
    main()
