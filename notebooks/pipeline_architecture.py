"""Render the SSUextract pipeline diagram as SVG.

The diagram reads top to bottom: what the user supplies on the left, the
numbered pipeline steps in the middle, and the files each step writes on the
right. Marks are grey. The one accent colour marks the optional tree step,
which also has dashed lines and an "optional" label.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

plt.rcParams["svg.hashsalt"] = "ssuextract-docs"

INK = "#1F2328"
MUTED = "#57606A"
LINE = "#8C959F"
CARD = "#F3F4F6"
CARD_EDGE = "#D0D7DE"
BADGE = "#32383F"
ACCENT = "#D55E00"
WHITE = "#FFFFFF"
SANS = "DejaVu Sans"
MONO = "DejaVu Sans Mono"
DASHES = (0, (4, 3))

# The axes fill the figure and use inches, so one data unit is one inch.
FIGURE_WIDTH_IN = 8.9
FIGURE_HEIGHT_IN = 6.25
INPUT_X = 0.1
INPUT_WIDTH = 1.85
STEP_X = 2.4
STEP_WIDTH = 4.1
OUTPUT_X = 6.9
CARD_HEIGHT = 0.78
HEADER_Y = 6.03
FIRST_ROW_CENTER = 5.45
ROW_PITCH = 0.98
BADGE_OFFSET = 0.33
TEXT_OFFSET = 0.66


@dataclass(frozen=True)
class Stroke:
    """Colour and dash pattern of a line or a border."""

    color: str
    linestyle: str | tuple[int, tuple[int, int]] = "-"


@dataclass(frozen=True)
class Output:
    """One file a step writes, with a short description."""

    filename: str
    caption: str


@dataclass(frozen=True)
class Step:
    """One pipeline step: its badge, title, description and outputs."""

    badge: str
    title: str
    description: str
    outputs: tuple[Output, ...] = ()
    is_optional: bool = False


@dataclass(frozen=True)
class Source:
    """Something the user supplies, and the row of the step that reads it."""

    title: str
    detail: str
    row: int
    detail_font: str = SANS
    is_optional: bool = False


FLOW = Stroke(LINE)
BORDER = Stroke(CARD_EDGE)
OPTIONAL = Stroke(ACCENT, DASHES)

STEPS = (
    Step(
        "1",
        "Detect rRNA genes",
        "Infernal searches each assembly with the 16S\n"
        "and 18S rRNA gene models (RF00177, RF01960).",
    ),
    Step(
        "2",
        "Resolve overlaps",
        "Where both models hit one locus, the hit with\nthe lower E-value is kept.",
    ),
    Step(
        "3",
        "Extract genes",
        "Each kept hit is cut out on its strand. Partial\n"
        "hits in model order are joined.",
        (Output("extracted/*.fna", "gene sequences"),),
    ),
    Step(
        "4",
        "Search references",
        "BLAST compares each gene with the 16S or 18S\nindex of the database profile.",
        (Output("blast_top_hits.tsv", "ranked reference hits"),),
    ),
    Step(
        "5",
        "Assign taxonomy",
        "A gene identical to a reference gets its taxonomy.\n"
        "Other genes need a calibrated database profile.",
        (
            Output("cmsearch_summary.tsv", "one row per gene"),
            Output("cmsearch_summary.tab", "contig counts per sample"),
        ),
    ),
    Step(
        "+",
        "Tree classification",
        "Builds a tree of the gene and 100 references. The\n"
        "nearest neighbors must also support the call.",
        (Output("tree_nearest_neighbors.tsv", "references by tree distance"),),
        is_optional=True,
    ),
)
SOURCES = (
    Source("Assemblies", ".fna, .fa or .fasta files", row=0),
    Source("Database profile", "curated or img", row=3),
    Source(
        "Optional flag",
        "--tree_classification",
        row=5,
        detail_font=MONO,
        is_optional=True,
    ),
)


def build() -> Figure:
    """Draw the pipeline diagram and return the figure."""
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), dpi=150)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, FIGURE_WIDTH_IN)
    ax.set_ylim(0, FIGURE_HEIGHT_IN)
    ax.axis("off")
    fig.patch.set_facecolor(WHITE)

    _draw_headers(ax)
    for row, step in enumerate(STEPS):
        _draw_step(ax, step, _row_center(row))
    _draw_step_links(ax)
    for source in SOURCES:
        _draw_source(ax, source)
    return fig


def _row_center(row: int) -> float:
    return FIRST_ROW_CENTER - row * ROW_PITCH


def _draw_headers(ax: Axes) -> None:
    for x, label in (
        (INPUT_X, "INPUT"),
        (STEP_X, "PIPELINE STEPS"),
        (OUTPUT_X, "OUTPUT FILES"),
    ):
        ax.text(
            x,
            HEADER_Y,
            label,
            ha="left",
            va="center",
            color=MUTED,
            family=SANS,
            fontsize=8,
            fontweight="bold",
        )


def _draw_step(ax: Axes, step: Step, center_y: float) -> None:
    _draw_card(
        ax,
        STEP_X,
        center_y,
        STEP_WIDTH,
        WHITE if step.is_optional else CARD,
        OPTIONAL if step.is_optional else BORDER,
    )
    _draw_badge(ax, step, center_y)
    text_x = STEP_X + TEXT_OFFSET
    ax.text(
        text_x,
        center_y + 0.2,
        step.title,
        ha="left",
        va="center",
        color=INK,
        family=SANS,
        fontsize=11.5,
        fontweight="bold",
    )
    ax.text(
        text_x,
        center_y - 0.13,
        step.description,
        ha="left",
        va="center",
        color=MUTED,
        family=SANS,
        fontsize=8.8,
        linespacing=1.3,
    )
    if step.is_optional:
        ax.text(
            STEP_X + STEP_WIDTH - 0.16,
            center_y + 0.2,
            "OPTIONAL",
            ha="right",
            va="center",
            color=MUTED,
            family=SANS,
            fontsize=8,
            fontweight="bold",
        )
    _draw_outputs(ax, step, center_y)


def _draw_card(
    ax: Axes,
    x: float,
    center_y: float,
    width: float,
    facecolor: str,
    stroke: Stroke,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, center_y - CARD_HEIGHT / 2),
            width,
            CARD_HEIGHT,
            boxstyle="round,pad=0,rounding_size=0.1",
            linewidth=1.1,
            edgecolor=stroke.color,
            facecolor=facecolor,
            linestyle=stroke.linestyle,
            zorder=2,
        )
    )


def _draw_badge(ax: Axes, step: Step, center_y: float) -> None:
    badge_x = STEP_X + BADGE_OFFSET
    ax.add_patch(
        Circle(
            (badge_x, center_y),
            0.17,
            facecolor=WHITE if step.is_optional else BADGE,
            edgecolor=ACCENT if step.is_optional else BADGE,
            linewidth=1.2,
            zorder=3,
        )
    )
    ax.text(
        badge_x,
        center_y - 0.005,
        step.badge,
        ha="center",
        va="center",
        color=ACCENT if step.is_optional else WHITE,
        family=SANS,
        fontsize=10,
        fontweight="bold",
        zorder=4,
    )


def _draw_outputs(ax: Axes, step: Step, center_y: float) -> None:
    if not step.outputs:
        return
    _draw_arrow(
        ax,
        (STEP_X + STEP_WIDTH, center_y),
        (OUTPUT_X - 0.08, center_y),
        OPTIONAL if step.is_optional else FLOW,
    )
    block_height = 0.4
    top = center_y + block_height * (len(step.outputs) - 1) / 2
    for index, output in enumerate(step.outputs):
        y = top - index * block_height
        ax.text(
            OUTPUT_X,
            y + 0.075,
            output.filename,
            ha="left",
            va="center",
            color=INK,
            family=MONO,
            fontsize=8.6,
        )
        ax.text(
            OUTPUT_X,
            y - 0.085,
            output.caption,
            ha="left",
            va="center",
            color=MUTED,
            family=SANS,
            fontsize=8.2,
        )


def _draw_step_links(ax: Axes) -> None:
    badge_x = STEP_X + BADGE_OFFSET
    for (upper_row, _), (lower_row, lower) in pairwise(enumerate(STEPS)):
        _draw_arrow(
            ax,
            (badge_x, _row_center(upper_row) - CARD_HEIGHT / 2),
            (badge_x, _row_center(lower_row) + CARD_HEIGHT / 2),
            OPTIONAL if lower.is_optional else FLOW,
        )


def _draw_source(ax: Axes, source: Source) -> None:
    center_y = _row_center(source.row)
    stroke = OPTIONAL if source.is_optional else FLOW
    _draw_card(ax, INPUT_X, center_y, INPUT_WIDTH, WHITE, stroke)
    ax.text(
        INPUT_X + 0.16,
        center_y + 0.14,
        source.title,
        ha="left",
        va="center",
        color=INK,
        family=SANS,
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        INPUT_X + 0.16,
        center_y - 0.15,
        source.detail,
        ha="left",
        va="center",
        color=MUTED,
        family=source.detail_font,
        fontsize=8.4,
    )
    _draw_arrow(ax, (INPUT_X + INPUT_WIDTH, center_y), (STEP_X, center_y), stroke)


def _draw_arrow(
    ax: Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    stroke: Stroke,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=1.2,
            color=stroke.color,
            linestyle=stroke.linestyle,
            shrinkA=1.5,
            shrinkB=1.5,
            zorder=3,
        )
    )


def main() -> int:
    """Write the diagram to docs/assets/figures/pipeline-architecture.svg."""
    output = Path(__file__).resolve().parents[1] / "docs" / "assets" / "figures"
    figure = build()
    svg = output / "pipeline-architecture.svg"
    figure.savefig(
        svg,
        bbox_inches="tight",
        pad_inches=0.08,
        facecolor=WHITE,
        metadata={"Date": None, "Creator": "SSUextract"},
    )
    svg_text = svg.read_text(encoding="utf-8")
    svg.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(figure)
    svg.chmod(0o644)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
