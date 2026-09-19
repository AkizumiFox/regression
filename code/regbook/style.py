import matplotlib as mpl

# The book's figure palette (also used for the colours in TikZ figures).
COLORS = {
    "ink": "#1B1F27",
    "accent": "#2F5F96",
    "second": "#A45C25",
    "third": "#14705F",
    "thread": "#74509A",
    "grid": "#D2D9E0",
    "muted": "#7A8492",
}


def use_book_style():
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Libertinus Serif", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.edgecolor": COLORS["ink"],
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.2,
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
    })
