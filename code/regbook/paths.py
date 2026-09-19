import atexit
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
GENERATED = ROOT / "code" / "_generated"
_FIGURES: list[Path] = []


def chapter_dir(chapter: str) -> Path:
    """src/chNN-slug/ for a chapter tag such as "ch06"."""
    matches = sorted(SRC.glob(f"{chapter}-*"))
    if not matches:
        raise FileNotFoundError(f"no chapter directory src/{chapter}-*")
    return matches[0]


def figure_path(chapter: str, name: str) -> Path:
    """src/chNN-slug/<name>.pdf, next to the chapter's Markdown.

    The PDF is used by the print edition; an SVG copy for the website is made from it
    when the script exits (pdftocairo), so scripts only ever save the PDF.
    """
    path = chapter_dir(chapter) / f"{name}.pdf"
    _FIGURES.append(path)
    return path


def generated_path(chapter: str, name: str) -> Path:
    """code/_generated/<chapter>/<name>.json, read by tools/check_numbers.py."""
    d = GENERATED / chapter
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{name}.json"


@atexit.register
def _make_svgs():
    if not _FIGURES or not shutil.which("pdftocairo"):
        return
    for pdf in _FIGURES:
        if pdf.exists():
            subprocess.run(["pdftocairo", "-svg", str(pdf), str(pdf.with_suffix(".svg"))], check=True)
