#!/usr/bin/env python3
"""Run every script in code/chNN/ and confirm the numbers the book quotes are still what they compute,
and that every runnable cell in src/ still matches the script it was taken from.

Each script asserts its own mathematical claims and records the values the text quotes
(code/_generated/chNN/*.json). tools/numbers.json is the list of values the Markdown in
src/ was written with. A script that fails, or a quoted value that has changed, is an error:
update the prose (and then tools/numbers.json with --accept) if the change is intended.

    MPLBACKEND=Agg python tools/check_numbers.py [--accept] [--no-run]
"""
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RECORD = ROOT / "tools" / "numbers.json"


def check_cells():
    """Each {.python .run #cell-<script>-<region>} cell must consist of lines of code/chNN/<script>.py."""
    import re
    problems = 0
    cells = 0
    for md in sorted((ROOT / "src").glob("ch*/*.md")):
        chapter = md.parent.name[:4]
        scripts = {p.stem.replace("_", "-"): p for p in (ROOT / "code" / chapter).glob("*.py")}
        for m in re.finditer(r"```\{\.python \.run #cell-([^}\s]+)\}\n(.*?)```", md.read_text(), re.S):
            cells += 1
            cid, code = m.group(1), m.group(2)
            script = next((p for stem, p in sorted(scripts.items(), key=lambda kv: -len(kv[0])) if cid.startswith(stem)), None)
            if script is None:
                print(f"CELL {md.relative_to(ROOT)} #cell-{cid}: no matching script in code/{chapter}/")
                problems += 1
                continue
            source = {line.strip() for line in script.read_text().splitlines()}
            missing = [l for l in code.splitlines() if l.strip() and l.strip() not in source]
            if missing:
                print(f"CELL {md.relative_to(ROOT)} #cell-{cid}: lines not in {script.relative_to(ROOT)}: {missing[:3]}")
                problems += 1
    print(f"{cells - problems}/{cells} runnable cells match their scripts")
    return problems


def main():
    accept = "--accept" in sys.argv
    if check_cells():
        return 1
    if "--no-run" not in sys.argv:
        env = dict(os.environ, PYTHONPATH=str(ROOT / "code"), MPLBACKEND="Agg")
        failed = []
        for script in sorted((ROOT / "code").glob("ch*/*.py")):
            r = subprocess.run([sys.executable, str(script)], cwd=ROOT, env=env, capture_output=True, text=True)
            status = "ok" if r.returncode == 0 else "FAILED"
            print(f"{status:6s} {script.relative_to(ROOT)}")
            if r.returncode:
                failed.append(script)
                print(r.stderr[-2000:])
        if failed:
            print(f"{len(failed)} script(s) failed")
            return 1
    current = {}
    for f in sorted((ROOT / "code" / "_generated").glob("ch*/*.json")):
        current.update(json.loads(f.read_text()))
    if accept or not RECORD.exists():
        RECORD.write_text(json.dumps(current, indent=1, ensure_ascii=False, sort_keys=True) + "\n")
        print(f"recorded {len(current)} values in {RECORD.relative_to(ROOT)}")
        return 0
    recorded = json.loads(RECORD.read_text())
    changed = {k: (v, current.get(k)) for k, v in recorded.items() if current.get(k) != v}
    for k, (old, new) in changed.items():
        print(f"CHANGED {k}: text says {old!r}, script now gives {new!r}")
    print(f"{len(recorded) - len(changed)}/{len(recorded)} quoted values unchanged")
    return 1 if changed else 0


if __name__ == "__main__":
    sys.exit(main())
