# Vendored sleek

The LaTeX preamble is based on **sleek**, François Rozet's style
(<https://github.com/francois-rozet/sleek>, v1.01, 2020-10-29). The copies here
came from a local checkout at `~/coding/packages/`, which is outside this
repository and is committed nowhere; referring to it by absolute path would
break the build for everyone else and in CI, so the two files the book uses are
vendored into the tree.

| file | what it gives the book |
|------|------------------------|
| `sleek.sty` | page geometry, running-head frame, paragraph and table typography, link colours, a few helper macros |
| `sleek-theorems.sty` | the `mdthm` key for `\declaretheorem`, and the `thicc` mdframed style the environments are framed in |

**Not vendored**, because the book uses none of it: `sleek-title.sty` (the book
has its own title page, and `\maketitle` is `\let` to `\relax`),
`sleek-code.sty`, `sleek-diagrams.sty` and `sleek-listings.sty` (a different
book's minted listings, flowchart TikZ styles and `listings` colour schemes,
with Chinese captions).

`latex/preamble.tex` and `latex/preamble-book.tex` load these two, then layer
the book on top: `latex/latex-template.sty` (macros, colours, helpers),
`latex/chapter-style.sty` (chapters, `\part` banners, running heads) and
`latex/theorem-envs.sty` (the environments `filters/theorems.lua` emits, framed
in the sleek style, coloured from the book's config).

## Changes to `sleek.sty`

Every removal is left in the file as a comment saying what went and why.

1. **`\ProvidesPackage{sleek}` → `\ProvidesPackage{latex/sleek/sleek}`**, so the
   name matches the path it is loaded by. `latex/latex-template.sty` and
   `latex/chapter-style.sty` already do this.
2. **Removed the whole `biblatex` block** (`backend=biber, style=numeric-comp,
   sorting=none, maxbibnames=99, backref=true`) and the two
   `\DefineBibliographyStrings` that configure it. Neither edition of this book
   has a bibliography, so it bought nothing but a `biber` dependency and a
   "Please (re)run Biber" warning on every build.
3. **Removed `\abs`, `\norm`, `\diff` and `\tr`.** The book defines three of the
   four itself, and those definitions stay the only ones:
   - `\abs`: `latex/latex-template.sty` declares it with
     `\DeclarePairedDelimiter`, which has a starred form sleek's plain
     `\left|...\right|` does not;
   - `\norm`: `latex/macros.tex` declares it with `\newcommand`, so sleek's
     would not have been a silent override but a hard "command already
     defined" error;
   - `\tr`: `latex/latex-template.sty` declares it with `\DeclareMathOperator`,
     which is `\newcommand` underneath — the same hard error;
   - `\diff`: the book neither defines nor uses it, so dropping it costs
     nothing and keeps this copy identical to the sibling book's.
   Keeping any of sleek's would have silently reflowed mathematics in all 45
   chapters.
4. **Removed `\renewcommand{\headrulewidth}{0pt}`.** Both sides `\renewcommand`
   it; the book's 0.4pt rule (in `latex/chapter-style.sty`) wins, and is now the
   only `\headrulewidth` in the preamble. The book has running heads with
   chapter titles, which want the rule; sleek has no chapters.
5. **Removed the `hyperref` load.** hyperref is loaded exactly once, by
   `templates/latex/template.tex`, after `latex/macros.tex` — it has to be,
   because it reclaims `\H` and the preamble restores the bold matrix **H**
   afterwards. Sleek's options (`pdfusetitle, colorlinks, linktoc=page,
   allcolors=blue`) reach that one load through `\PassOptionsToPackage` in
   `latex/preamble.tex`, so the sleek link style survives without a second load
   and without an option clash. `templates/latex/template.tex` asks for
   `colorlinks=true, linkcolor=blue, urlcolor=blue`, which is what
   `allcolors=blue` says, so no link changes colour.
6. **Changed the geometry** from `top=3.5cm,bottom=2.5cm,right=2.5cm,left=2.5cm`
   to a uniform `2.5cm`. Left and right are sleek's own values either way, so the
   measure — and with it every line break in the book — is identical; sleek's
   only real change was a 1cm deeper top margin, which costs about 5% more pages
   and buys nothing the 0.4pt head rule does not already give. This line is now
   the single place the page is decided: `latex/preamble.tex` loads `geometry`
   with no options at all.
7. **Removed `\DeclareOption{noheader}`** (it does `\newgeometry`, fighting the
   line above) **and `\DeclareOption{french}`** (it needs babel-french, siunitx
   and the bibliography strings, none of which this book loads). The `parindent`
   option is kept.
8. **Removed packages the book does not use**: `eurosym`, `esint`, `siunitx`,
   `float`, `multicol`, `csquotes`. Checked against the generated LaTeX: no
   `\euro`, no `\oiint`, no `\si`/`\SI`, no `[H]` float, no `multicols`, no
   `\enquote`.
9. **Dropped the `parfill` option from `parskip`** (`\RequirePackage{parskip}`,
   not `\RequirePackage[parfill]{parskip}`). `parfill` sets
   `\parfillskip=30pt plus 1fil`, i.e. 30pt of glue at the end of every
   paragraph that can stretch but not shrink. In a book this dense that is
   expensive: it adds 30pt to the natural width of every single-line paragraph,
   which turns lines that merely reached the margin into overfull boxes, and in
   at least one place makes the only feasible line break disappear so a formula
   runs into the margin. Measured in the sibling linear-algebra book, which
   runs the same engine and the same preamble, over its 262 section PDFs: 363
   overfull boxes before any of this, 357 with sleek and `parfill`, **290** with
   sleek without it. It was never switched on here. The rest of `parskip` is kept — `latex/chapter-style.sty` sets `\parskip`
   and `\parindent` afterwards anyway, so what the package still contributes is
   its list-spacing corrections.
10. **Added `\RequirePackage{lmodern}`** next to sleek's `\RequirePackage[T1]{fontenc}`.
   This is the one thing added rather than removed. T1 without `cm-super` (not
   installed here) falls back to bitmap Type 3 fonts; Latin Modern is the
   standard Type 1 answer and is metric-compatible with the Computer Modern the
   book used under OT1.

   Going from OT1 to T1 is the one change in this file the book had to answer.
   `latex/latex-template.sty` maps the caron letters the prose uses (Š in
   "Šidák", č, ě) to `\accent` plus the base letter, because the vector
   shorthands shadow `\v`. Slot 20 is the caron in OT1 but *guillemotright* in
   T1, so "Šidák" would have set as "»Sidák"; those nine
   `\DeclareUnicodeCharacter` lines now name slot 7, T1's caron. T1 also gives
   the accented letters real glyphs, so "Scheffé" now comes out of the PDF as
   one character instead of `e` plus a combining acute.

## Changes to `sleek-theorems.sty`

1. **`\ProvidesPackage`** renamed to the path, as above.
2. **Removed `\renewcommand{\qedsymbol}{$\blacksquare$}`.**
   `latex/latex-template.sty` already sets exactly this and stays the one place
   that does — it was moved out of that file's `\if@skipenv` block so that it
   still runs when the environments come from `latex/theorem-envs.sty`, which
   ends a proof with `\qed` and sets no symbol of its own. (Leaving it inside
   would have quietly swapped the book's filled square for amsthm's hollow
   `\openbox` in all 671 proofs.) (`\qedadd` is kept.)
3. **Removed the language strings** (`\lgthm` … `\lgtip`) and
   `\DeclareOption{french}`. The book's environment names come from its config.
4. **Removed sleek's own environments** — `thm`, `lem`, `prop`, `defn`, `hyp`,
   `meth`, `quest`, `answ`, `expl`, `rmk`, `note`, `tip` and the `framedthm` …
   `framedexpl` declarations built on them. They are the wrong set (this book's
   environments are whatever `config/config.json` lists) and the wrong colours
   (hard-coded `red!20`, `blue!20`, `black!20!green!20`, `orange!20`).
   `latex/theorem-envs.sty` declares the book's environments instead, with the
   same `\declaretheorem[mdthm={style=..., ...}]` pattern, taking each colour
   from the config.
5. **Kept `mdfdefinestyle{thicc}` unchanged.** `latex/theorem-envs.sty` derives
   from it with `style=thicc` and overrides two things: `nobreak=false` (a
   single proof in this book can be longer than a page, and a book this long of
   unbreakable boxes is not an option) and the book's split skips.

## Things that are *not* changes but are worth knowing

- `latex/chapter-style.sty` sets `\parskip=0.8em` and `\parindent=0pt` after
  `parskip` has been loaded, so the paragraph style is the book's, as before.
- `sleek.sty` sets `\pagestyle{fancy}` with `\fancyhead[R]{\leftmark}`;
  `latex/chapter-style.sty` is loaded later and replaces the whole `\fancyhf`,
  so the book keeps its own running heads (CHAPTER <part>.<chapter>, the title
  in capitals, over the 0.4pt rule) and its `\part` banners.
- `\arraystretch` becomes 1.1 (sleek's value; the book had the 1.0 default).
- `latex/theorem-envs.sty` differs from the sibling book's in one place: its
  `\booksmallenv` only defines the plain `\begin{<env>}` form when the name is
  free or is already an environment. `config/config.json` lists a `check`
  environment, and taking that name would have destroyed LaTeX's `\check` math
  accent, which two sections use (`\check{\X}`, `\check{\bbeta}`). Nothing is
  lost: `filters/theorems.lua` writes the head into the content and emits
  `checkinline`, so the plain form is a fallback this book never uses.

## Effect on the book

Measured on `_build/tmp/latex-book/book.log`, a full `./build.py all` before and
after, in a copy of the tree:

| | before | after |
|---|---|---|
| overfull `\hbox` | 947 | **877** |
| pages carrying one | 338 | **272** |
| worst | 157.3pt | **154.4pt** |
| total pages | 2078 | **2033** |
| `!` errors | 98 | **84** |

Across the 351 section PDFs the same way: 356 overfull boxes before, **286**
after; 449 `!` errors before, **435** after.

The errors that went are the ten from `\begin{algorithm}`, which no edition of
this book ever defined (`latex/latex-template.sty` had an `algorithmbox` but no
`algorithm`, and `filters/theorems.lua` has no `algorithm` entry in
`SMALL_ENV_INLINE`, so the plain name was emitted and was undefined) — the five
algorithm blocks are typeset for the first time — and four `\mathbf allowed only
in math mode` from "Yohai, Víctor" in `src/ch20-residuals-leverage-influence`,
which T1's precomposed `í` no longer routes through an accent macro. Every other
error is unchanged: they are the book's own, they were there before this, and
none of them is new.

Section by section the overfull boxes move both ways; the largest single
regression is `src/ch06-projections/04-least-squares.md`, 2 boxes to 6, and the
largest single improvement is `src/ch01-matrix-algebra/03-inverses.md`, 5 to 0.
Nothing was changed in `src/` either way.
