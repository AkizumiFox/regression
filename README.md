# Regression: From Projection to Distribution

One source, Pandoc Markdown in `src/`, built by `build.py` into both the website and the PDFs
(the same engine as the linear-algebra book).

| Task | Command |
|---|---|
| Write with live preview | `./build.py serve` → http://127.0.0.1:8000 (rebuilds and reloads on save) |
| Build the HTML site | `./build.py html` |
| Build the book PDF | `./build.py book` (writes `_build/pdf/book.pdf`) |
| Everything | `./build.py all` (see `BUILD.md`) |
| Validate references and numbering | `./build.py check` |
| Run the book's code, check its numbers and cells | `MPLBACKEND=Agg .venv/bin/python tools/check_numbers.py` |
| Chapter order against the prerequisite DAG | `.venv/bin/python blueprint/check_order.py` |
| Originality against the source books | `.venv/bin/python blueprint/originality.py 06` |
| Publish | `./build.py deploy --push` → regression.akizumifox.com (see `DEPLOY.md`) |

First time: `python3 -m venv .venv && .venv/bin/pip install -r code/requirements.txt`.

## Reading the site

A toolbar at the top of each page folds the contents (left) and "On this page" (right) and
switches light/dark. Press `/` to search, ← and → to turn pages. Solutions start folded. Code
blocks are live Python cells: **Run** executes them in the browser (Pyodide). The sidebar links
to a list of all results and a dependency graph.

## Layout

- `src/`: the written chapters (1–25), one Markdown file per section, with the figures
  of each chapter beside its Markdown (`name.pdf` for print, `name.svg` for the web). `src/index.md`
  is the preface, how-to-use and notation page.
- `config/config.json`: chapters, parts, environments, deploy target. `latex/macros.tex`: the
  notation macros, used by both the PDF and the website.
- `code/chNN/*.py`: the scripts behind every figure and computed number. They assert their own
  claims; `regbook` writes figures into `src/` and the quoted values to `code/_generated/`.
  `tools/numbers.json` records the values the text quotes; `tools/check_numbers.py` fails if a
  script fails, a quoted value changes, or a runnable cell no longer matches its script.
- `blueprint/`: the plan. `book.yaml` (all 45 chapters: sources, outlines, budgets), `dag.json`
  (the prerequisite graph of the six source books), `check_order.py`, `originality.py`,
  `coverage/chNN.md` (source topics mapped to where the book covers them), `make_blueprints.py`.
- `authoring/`: `STYLE.md` (how to write a section), `NOTATION.md`, `STATUS.md`, and
  `blueprints/chNN.md` (the brief for every chapter, written or not).
- `build/`, `filters/`, `templates/`, `latex/*.sty`, `tests/`, `widgets/`: the engine, copied from
  the linear-algebra book (see "Another book" in `BUILD.md`).
- `regression/`: the six source books (PDFs, copyrighted, gitignored; used only for coverage
  checklists and originality checks, never published).
- `_build/` (output) and `site/` (deploy clone) are gitignored.

## History

The book was first written as a LaTeX edition (`book/`, compiled with LuaLaTeX). On 2026-09-19 it
was converted to this Markdown source with a one-off converter and the LaTeX edition was deleted.
