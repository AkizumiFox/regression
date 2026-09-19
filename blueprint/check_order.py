"""Verify the chapter order in book.yaml against the prerequisite DAG.

A hard edge a -> b is violated when some chapter uses (part of) b before
any chapter has used a. Cut nodes are replaced by the node that subsumes
them. Documented exemptions are reported but do not fail the check.
Soft edges are reported as warnings only.
"""
import json
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).parent


def load():
    dag = json.loads((HERE / "dag.json").read_text())
    book = yaml.safe_load((HERE / "book.yaml").read_text())
    return dag, book


def main():
    dag, book = load()
    cut = book.get("cut", {})
    exempt = {tuple(e["edge"]): e["why"] for e in book.get("exemptions", [])}
    nodes = set(dag["nodes"])

    first_use = {}
    errors = []
    for ch in book["chapters"]:
        for node, _ in ch["sources"]:
            if node not in nodes:
                errors.append(f"ch.{ch['num']}: unknown DAG node {node}")
            first_use.setdefault(node, ch["num"])

    def resolve(n):
        while n in cut:
            n = cut[n]
            if n is None:
                return None
        return n

    used = set(first_use)
    unused = sorted(n for n in nodes if n not in used and n not in cut)
    for n in unused:
        errors.append(f"DAG node {n} ({dag['nodes'][n]['title']}) is neither used nor cut")

    violations, exempted, warnings = [], [], []
    for kind, edges in (("hard", dag["hard"]), ("soft", dag["soft"])):
        for a0, b0 in edges:
            a, b = resolve(a0), resolve(b0)
            if a is None or b is None or a == b:
                continue
            if a not in first_use or b not in first_use:
                continue
            # every use of b must come at or after the first use of a
            uses_b = [ch["num"] for ch in book["chapters"] for n, _ in ch["sources"] if n == b]
            early = [c for c in uses_b if c < first_use[a]]
            if not early:
                continue
            msg = f"{a0}->{b0}: {b} used in ch.{early} before {a} (first in ch.{first_use[a]})"
            if kind == "soft":
                warnings.append(msg)
            elif (a0, b0) in exempt:
                exempted.append(f"{msg}\n      exemption: {exempt[(a0, b0)]}")
            else:
                violations.append(msg)

    nums = [c["num"] for c in book["chapters"]]
    if nums != list(range(1, len(nums) + 1)):
        errors.append("chapter numbers are not consecutive")

    print(f"{len(book['chapters'])} chapters, {sum(c['pages'] for c in book['chapters'])} pp budget, "
          f"{len(used)} DAG nodes used, {len(cut)} cut")
    for title, items in (("ERRORS", errors), ("HARD VIOLATIONS", violations),
                         ("EXEMPTED", exempted), ("soft warnings", warnings)):
        print(f"\n{title}: {len(items)}")
        for m in items:
            print("  " + m)
    unused_exemptions = [e for e in exempt if not any(m.startswith(f"{e[0]}->{e[1]}:") for m in exempted)]
    for e in unused_exemptions:
        print(f"  note: exemption {e} is not needed")
    return 1 if errors or violations else 0


if __name__ == "__main__":
    sys.exit(main())
