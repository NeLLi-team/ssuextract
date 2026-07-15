"""Publication figure defaults adapted from the beautiful-data-viz skill asset."""

from __future__ import annotations

from typing import Optional

import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.figure import Figure


TUFTE_COLORS = {
    "light_bg": "#ffffff",
    "light_text": "#111111",
    "light_secondary": "#666666",
    "light_axis": "#cccccc",
}


def set_beautiful_style(
    *, medium: str = "paper", background: str = "light", dpi: int = 180
) -> None:
    """Set restrained Matplotlib defaults for a paper or notebook figure."""
    if background != "light":
        raise ValueError("SSUextract documentation figures use a light background")
    sizes = {
        "paper": {"label": 9, "tick": 8, "legend": 8},
        "notebook": {"label": 11, "tick": 10, "legend": 10},
    }.get(medium)
    if sizes is None:
        raise ValueError(f"unsupported medium: {medium}")
    mpl.rcParams.update(
        {
            "figure.dpi": dpi,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "figure.facecolor": TUFTE_COLORS["light_bg"],
            "axes.facecolor": TUFTE_COLORS["light_bg"],
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "text.color": TUFTE_COLORS["light_text"],
            "axes.labelcolor": TUFTE_COLORS["light_secondary"],
            "axes.labelsize": sizes["label"],
            "xtick.labelsize": sizes["tick"],
            "ytick.labelsize": sizes["tick"],
            "legend.fontsize": sizes["legend"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": TUFTE_COLORS["light_axis"],
            "axes.linewidth": 0.6,
            "axes.grid": False,
            "lines.linewidth": 1.4,
            "lines.markersize": 4,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
        }
    )


def finalize_axes(
    ax: Axes,
    *,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
) -> Axes:
    """Remove non-data ink and apply axis labels."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_alpha(0.8)
    ax.spines["bottom"].set_alpha(0.8)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.grid(False)
    figure = ax.figure
    if isinstance(figure, Figure):
        figure.canvas.draw_idle()
    return ax
