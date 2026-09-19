"""
Build System Pandoc Module
==========================
Builds Pandoc commands for a page and runs them.
"""

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Optional
from urllib.parse import quote, urlencode

from .book import Book, Page, ENGINE_ROOT
from .utils import print_error, print_warning, run_with_crash_retry

MARKDOWN_FORMAT = "markdown+tex_math_single_backslash"


# =============================================================================
# Metadata
# =============================================================================

def page_metadata(book: Book, page: Page, output_format: str, extra: Optional[dict] = None) -> dict:
    """Metadata passed to the filters and templates for one page."""
    metadata = {
        "chapter-num": page.chapter_number,
        "section-num": page.section,
        "is-preface": page.is_preface,
        "environment_settings": book.environment_settings,
        "section-number": page.number,
        "asset-prefix": page.asset_prefix if output_format == "html" else "./",
        "title": book.title,
        "title-text": page.title,
        "pagetitle": f"{page.number} {page.title}" if page.number else page.title,
        "tikz-cache-dir": str(book.tikz_cache_dir),
        "page-shard": page.shard,
        "page-path": page.html_path,
        "book-root": str(book.root),
        "engine-root": str(ENGINE_ROOT),
        "source-path": str(page.source.relative_to(book.root)),
    }
    if book.config.get("deploy-domain"):
        metadata["page-url"] = f"https://{book.config['deploy-domain'].strip()}/{page.html_path}"
    if book.labels_file.exists():
        metadata["crossref-labels-file"] = str(book.labels_file)

    if output_format == "html":
        prefix = page.asset_prefix
        if page.chapter:
            metadata["breadcrumb-chapter-title"] = page.chapter.title
            metadata["breadcrumb-chapter-url"] = f"{prefix}{page.chapter.slug}/index.html"
        else:
            metadata["breadcrumb-chapter-title"] = page.title
            metadata["breadcrumb-chapter-url"] = f"{prefix}index.html"
        metadata["breadcrumb-page-title"] = page.title
        metadata["breadcrumb-page-url"] = f"{prefix}{page.html_path}"
        metadata["pdf-url"] = f"{prefix}pdf/{Path(page.html_path).with_suffix('.pdf').as_posix()}"
        metadata["book-pdf-url"] = f"{prefix}book/book.pdf"
        if book.config.get("repo-url"):
            metadata["repo-url"] = book.config["repo-url"]
        if book.config.get("issues-url"):
            # A new issue prefilled with the section and its address
            label = f"{page.number} {page.title}" if page.number else page.title
            address = metadata.get("page-url") or page.html_path
            query = urlencode({
                "title": f"Typo or mistake in {label}",
                "body": f"Page: {address}\n\nWhat is wrong (quote the sentence or formula):\n\n\nWhat it should say:\n",
            }, quote_via=quote)
            metadata["report-url"] = f"{book.config['issues-url']}?{query}"
        for direction in ("prev", "next"):
            neighbour = page.neighbours.get(direction)
            if neighbour:
                metadata[f"{direction}-page-url"] = f"{prefix}{neighbour.html_path}"
                metadata[f"{direction}-page-title"] = neighbour.title
                # For the aria-label attribute: rendered math markup would break the HTML there
                metadata[f"{direction}-page-label"] = plain_label(neighbour.title)
                metadata[f"{direction}-page-number"] = neighbour.number

    metadata.update(extra or {})
    return metadata


def plain_label(title: str) -> str:
    """A title as plain text for attributes: TeX math reduced to readable symbols."""
    from .extras import plain_title  # local import: extras imports this module
    return plain_title(title)


def write_metadata_file(book: Book, page: Page, output_format: str, metadata: dict) -> Path:
    key = hashlib.md5(f"{page.source}:{output_format}".encode()).hexdigest()
    meta_file = book.build_dir / "tmp" / "meta" / f"{key}.json"
    meta_file.parent.mkdir(parents=True, exist_ok=True)
    # Top-level numbers as strings, the form the filters and LaTeX template expect
    metadata = {k: str(v) if isinstance(v, int) and not isinstance(v, bool) else v for k, v in metadata.items()}
    meta_file.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return meta_file


# =============================================================================
# Pandoc Command Building
# =============================================================================

def build_pandoc_command(book: Book, page: Page, output_file: Path, output_format: str, metadata: dict) -> list[str]:
    """Pandoc command converting one page to HTML or LaTeX (output_format "html" or "pdf")."""
    filters = book.filters_dir
    cmd = [
        "pandoc", str(page.source),
        "-o", str(output_file),
        "--standalone",
        "-f", MARKDOWN_FORMAT,
        "--lua-filter", str(filters / "format-visibility.lua"),
        "--lua-filter", str(filters / "components.lua"),
        "--lua-filter", str(filters / "tikz.lua"),
        "--lua-filter", str(filters / "enumerate.lua"),
        "--lua-filter", str(filters / "theorems.lua"),
        "--metadata-file", str(write_metadata_file(book, page, output_format, metadata)),
    ]

    if output_format == "html":
        cmd += [
            "--template", str(book.html_template),
            "--lua-filter", str(filters / "numbering.lua"),
            "--mathjax",
        ]
    elif output_format == "pdf":
        # # -> chapter (0.4), ## -> section (0.4.1), ### -> subsection (0.4.1.1)
        cmd += [
            "--template", str(book.latex_template),
            "--top-level-division", "chapter",
            "--resource-path", ":".join([str(page.source.parent), str(book.src_dir)]),
            "--include-in-header", str(book.engine_file("latex/preamble.tex")),
        ]
        if book.macros_file.exists():
            cmd += ["--include-in-header", str(book.macros_file)]
    return cmd


# =============================================================================
# Pandoc Execution
# =============================================================================

def run_pandoc(cmd: list[str], *, cwd: Optional[Path] = None, env: Optional[dict] = None) -> bool:
    """Run a pandoc command and return success. Logs errors and warnings."""
    try:
        result = run_with_crash_retry(cmd, cwd=cwd, env=env)
    except subprocess.CalledProcessError as e:
        print_error(f"Pandoc failed for {cmd[1]}")
        if e.stdout and e.stdout.strip():
            print(f"  stdout:\n{e.stdout.strip()}")
        if e.stderr and e.stderr.strip():
            print_error(f"  stderr:\n{e.stderr.strip()}")
        return False
    if result.stderr and result.stderr.strip():
        print_warning(f"Pandoc stderr for {Path(cmd[1]).name}:")
        for line in result.stderr.strip().splitlines():
            print(f"  {line}")
    return True
