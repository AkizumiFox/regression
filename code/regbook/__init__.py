"""Shared helpers for the book's code: plotting style, output paths, generated numbers."""
from .paths import figure_path, generated_path
from .style import use_book_style, COLORS
from .texout import Generated

__all__ = ["figure_path", "generated_path", "use_book_style", "COLORS", "Generated"]
