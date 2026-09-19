# Writing standard for the book

For anyone who writes a chapter: me, or a drafting agent. Read this, `authoring/NOTATION.md`
and the chapter's brief in `authoring/blueprints/`. Chapter 6 (`src/ch06-projections/`) is the
approved exemplar of register, density and rigour; read a few of its sections first.

**Hard rules**

- **No copying.** The source books in `regression/` are copyrighted. Use them only to build the
  coverage checklist (`blueprint/coverage/chNN.md`): which topics and results the chapter must
  contain. Prose, proofs, examples, data and exercises are our own. Run
  `blueprint/originality.py NN` (at n = 8 and n = 6) and rewrite any shared phrase that is not
  stock mathematics.
- **Every number comes from code.** A computed number, table entry or figure comes from a script
  in `code/chNN/` that asserts the claims it illustrates, records the quoted values with
  `regbook.Generated`, and saves figures with `regbook.figure_path`. After writing, run
  `tools/check_numbers.py --accept` once to record the quoted values, then plain
  `tools/check_numbers.py` from then on.
- **Every claim is true, every hypothesis stated.** Full proofs of every central result.
- **Proofs cite only earlier results.** Forward pointers in prose ("Chapter 11 uses this") are fine.
- **Book-wide conventions** in `NOTATION.md` (noncentrality without ½, nonnegative definite, `\T`,
  bold upright matrices and vectors: `\mathbf{a}` for a bold Latin letter without a macro,
  `\boldsymbol{\uptheta}` for a bold lowercase Greek letter, `\boldsymbol{\Sigma}` for capital
  Greek; never `\bm` or `\boldsymbol{\theta}`).

## Page skeleton

One Markdown file per section: `src/chNN-slug/MM-slug.md`, first line `# Title` (no math in
titles). `##` for subsections, `###` only inside Exercises. The chapter's `index.md` has the
title, a one-paragraph story, **What you need.** (earlier chapters and results, linked), and a
**Roadmap** listing the sections. The last page of each chapter is `NN-summary-and-notes.md`:
`## Summary` (a short list of the chapter's results, linked), `## Notes and sources` (the
source books and primary literature, honestly credited) and `## References`.

## Blocks

```markdown
::: {#thm-proj-fwl}
[Frisch–Waugh–Lovell]

Statement.
:::

::: {.proof}
Proof.
:::
```

- Numbered: `thm`, `lem`, `cor`, `prp` (one counter), `def` (its own counter). Unnumbered and
  titled: `exm` (examples), `exr` (exercises). Small blocks: `.proof`, `.remark`, `.idea`
  (key ideas), `.solution`, `.warning`, `.check`.
- Label ids: `thm-<chapterprefix>-<name>` (prefixes: `mat`, `rv`, `mvn`, `qf`, `proj`, …).
- Equations: `\[ … \]{#eq-name}`. A labelled display ends its paragraph (blank line after it).
- References: `@thm-name`, `@eq-name`, `@exm-name`. Sections: a Markdown link,
  `[Section 6.4](../ch06-projections/04-least-squares.html)`. Chapters not yet published are plain text.
- Labelled lists: `::: {.enumerate options="label=(\alph*)"}` around a numbered list.
- Figures: both formats, from the chapter directory:

  ```markdown
  ::: {when-format="html"}
  ![**Figure 6.6.1.** Caption.](fwl_added_variable.svg){#fig-proj-fwl width=100%}
  :::

  ::: {when-format="pdf"}
  ![Caption.](fwl_added_variable.pdf){width=100%}
  :::
  ```

  The web caption carries the number; in print LaTeX numbers the figure itself (the two agree
  section by section). Refer to it as `[Figure 6.6.1](#fig-proj-fwl)`. TikZ may be written inline
  as raw LaTeX, followed by a caption paragraph `[**Figure 6.5.1.** Caption.]{#fig-...}`.
- Code: a runnable cell whose lines are taken from the chapter's script (checked by
  `tools/check_numbers.py`). The first cell from a script on a page includes the imports and data
  loading it needs.

  ````markdown
  ```{.python .run #cell-fwl-fwl}
  import numpy as np
  ...
  ```
  ````

## Exercises

At the end of each section page: `## Exercises`, then `### A. Check your understanding`
(routine), `### B. Practice` (core, used later) and `### C. Going deeper` (hard). Each exercise
is `::: {#exr-...}` titled `[A1]`, `[B2]`, …; a worked solution follows it as `::: {.solution}`
(folded on the web). Give solutions to the core exercises that later chapters rely on.
