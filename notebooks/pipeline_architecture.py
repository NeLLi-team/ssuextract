"""Render the SSUextract pipeline diagram as SVG and PNG.

The layer-band layout and orthogonal routing follow the visual pattern used by
``fmschulz/memd`` in ``docs/figures/architecture.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


plt.rcParams["svg.hashsalt"] = "ssuextract-docs"


INDIGO = "#3949AB"
INDIGO_SOFT = "#E8EAF6"
TEAL = "#00695C"
TEAL_SOFT = "#E0F2F1"
ORANGE = "#BF360C"
ORANGE_SOFT = "#FFF3E0"
SLATE = "#37474F"
SLATE_SOFT = "#ECEFF1"
GRID = "#CFD8DC"
TEXT = "#212121"
SUBTLE = "#546E7A"
WHITE = "#FFFFFF"


@dataclass(frozen=True)
class Box:
    x: float
    y: float
    width: float
    height: float
    title: str
    body: str
    edge: str
    fill: str
    mono: bool = False

    @property
    def cx(self) -> float:
        return self.x + self.width / 2

    @property
    def cy(self) -> float:
        return self.y + self.height / 2


def draw_band(ax, y: float, height: float, label: str) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (0.4, y),
            13.2,
            height,
            boxstyle="round,pad=0.02,rounding_size=0.22",
            linewidth=1.0,
            edgecolor=GRID,
            facecolor="#FAFAFA",
            zorder=1,
        )
    )
    ax.text(
        0.5,
        y + height + 0.04,
        label,
        ha="left",
        va="bottom",
        color=SLATE,
        family="DejaVu Sans",
        fontsize=9.5,
        fontweight="bold",
    )


def draw_box(ax, box: Box) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (box.x, box.y),
            box.width,
            box.height,
            boxstyle="round,pad=0.02,rounding_size=0.16",
            linewidth=1.4,
            edgecolor=box.edge,
            facecolor=box.fill,
            zorder=2,
        )
    )
    ax.text(
        box.cx,
        box.y + box.height - 0.20,
        box.title,
        ha="center",
        va="top",
        color=box.edge,
        family="DejaVu Sans",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        box.cx,
        box.y + box.height * 0.39,
        box.body,
        ha="center",
        va="center",
        color=TEXT,
        family="DejaVu Sans Mono" if box.mono else "DejaVu Sans",
        fontsize=8.2 if box.mono else 8.7,
    )


def line(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    linestyle: str = "-",
) -> None:
    ax.plot(
        [start[0], end[0]],
        [start[1], end[1]],
        color=color,
        linewidth=1.25,
        linestyle=linestyle,
        zorder=3,
    )


def arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    linestyle: str = "-",
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=11,
            linewidth=1.25,
            color=color,
            linestyle=linestyle,
            shrinkA=2,
            shrinkB=2,
            zorder=4,
        )
    )


def build() -> plt.Figure:
    fig, ax = plt.subplots(figsize=(11.0, 9.0), dpi=150)
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 11.5)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.patch.set_facecolor(WHITE)

    bands = {
        "input": (9.05, 1.25, "Inputs"),
        "detect": (6.75, 1.85, "Detection and extraction"),
        "search": (4.55, 1.75, "Marker-specific search"),
        "taxonomy": (2.75, 1.35, "Taxonomy selection"),
        "output": (0.15, 2.15, "Results"),
    }
    for y, height, label in bands.values():
        draw_band(ax, y, height, label)

    inputs = [
        Box(1.1, 9.20, 3.35, 0.95, "Assembly FASTA", ".fna  ·  .fa  ·  .fasta", INDIGO, INDIGO_SOFT, True),
        Box(5.30, 9.20, 3.35, 0.95, "Covariance models", "RF00177  ·  RF01960", INDIGO, INDIGO_SOFT, True),
        Box(9.50, 9.20, 3.35, 0.95, "Database profile", "curated  or  img", INDIGO, INDIGO_SOFT, True),
    ]
    detection = [
        Box(1.05, 6.95, 3.55, 1.45, "Infernal search", "cmsearch per\nassembly and model", TEAL, TEAL_SOFT),
        Box(
            5.23,
            6.95,
            3.55,
            1.45,
            "Model competition",
            "CL00111 same-strand overlap\nE-value → bit score",
            TEAL,
            TEAL_SOFT,
        ),
        Box(
            9.40,
            6.95,
            3.55,
            1.45,
            "Interval extraction",
            "complete interval · strand\nreverse complement\nextracted/*.fna",
            TEAL,
            TEAL_SOFT,
        ),
    ]
    searches = [
        Box(2.05, 4.75, 4.15, 1.35, "RF00177", "16S rRNA gene\nBLAST index", TEAL, TEAL_SOFT),
        Box(7.80, 4.75, 4.15, 1.35, "RF01960", "18S rRNA gene\nBLAST index", TEAL, TEAL_SOFT),
    ]
    taxonomy = [
        Box(
            1.20,
            2.92,
            5.20,
            1.00,
            "BLAST taxonomy",
            "best-score LCA  ·  ties and conflicts retained",
            ORANGE,
            ORANGE_SOFT,
        ),
        Box(
            7.10,
            2.92,
            5.70,
            1.00,
            "Tree neighbors (optional)",
            "top 100 · cmalign · trim · IQ-TREE 3 · neighbor LCA",
            ORANGE,
            ORANGE_SOFT,
        ),
    ]
    outputs = [
        Box(0.80, 0.38, 3.75, 1.65, "Per-hit taxonomy", "cmsearch_summary.tsv", SLATE, SLATE_SOFT, True),
        Box(5.12, 0.38, 3.75, 1.65, "Reference evidence", "blast_top_hits.tsv\ntree_nearest_neighbors.tsv", SLATE, SLATE_SOFT, True),
        Box(9.45, 0.38, 3.75, 1.65, "Category counts", "cmsearch_summary.tab", SLATE, SLATE_SOFT, True),
    ]

    for box in [*inputs, *detection, *searches, *taxonomy, *outputs]:
        draw_box(ax, box)

    sequence_rail = 8.95
    for source in inputs[:2]:
        line(ax, (source.cx, source.y), (source.cx, sequence_rail), TEAL)
    line(ax, (inputs[0].cx, sequence_rail), (inputs[1].cx, sequence_rail), TEAL)
    arrow(ax, (detection[0].cx, sequence_rail), (detection[0].cx, detection[0].y + detection[0].height), TEAL)

    arrow(ax, (detection[0].x + detection[0].width, detection[0].cy), (detection[1].x, detection[1].cy), TEAL)
    arrow(ax, (detection[1].x + detection[1].width, detection[1].cy), (detection[2].x, detection[2].cy), TEAL)

    search_rail = 6.52
    line(ax, (detection[2].cx, detection[2].y), (detection[2].cx, search_rail), TEAL)
    sequence_entries = [box.cx - 0.55 for box in searches]
    line(ax, (sequence_entries[0], search_rail), (detection[2].cx, search_rail), TEAL)
    for box in searches:
        entry_x = box.cx - 0.55
        arrow(ax, (entry_x, search_rail), (entry_x, box.y + box.height), TEAL)

    database_x = 13.25
    database_rail = 6.28
    reference_style = "--"
    line(ax, (inputs[2].cx, inputs[2].y), (database_x, inputs[2].y), ORANGE, reference_style)
    line(ax, (database_x, inputs[2].y), (database_x, database_rail), ORANGE, reference_style)
    database_entries = [box.cx + 0.55 for box in searches]
    line(ax, (database_entries[0], database_rail), (database_x, database_rail), ORANGE, reference_style)
    for box in searches:
        entry_x = box.cx + 0.55
        arrow(
            ax,
            (entry_x, database_rail),
            (entry_x, box.y + box.height),
            ORANGE,
            reference_style,
        )

    taxonomy_rail = 4.32
    for box in searches:
        line(ax, (box.cx, box.y), (box.cx, taxonomy_rail), TEAL)
    line(ax, (searches[0].cx, taxonomy_rail), (searches[1].cx, taxonomy_rail), TEAL)
    for box in taxonomy:
        arrow(ax, (box.cx, taxonomy_rail), (box.cx, box.y + box.height), TEAL)

    output_rail = 2.52
    for box in taxonomy:
        line(ax, (box.cx, box.y), (box.cx, output_rail), TEAL)
    line(ax, (outputs[0].cx, output_rail), (outputs[2].cx, output_rail), TEAL)
    for box in outputs:
        arrow(ax, (box.cx, output_rail), (box.cx, box.y + box.height), TEAL)

    ax.text(
        0.4,
        11.30,
        "SSUextract — marker-specific small-subunit rRNA extraction",
        ha="left",
        va="top",
        color=TEXT,
        family="DejaVu Sans",
        fontsize=12,
        fontweight="bold",
    )
    legend_y = 10.88
    ax.plot([0.4, 0.85], [legend_y, legend_y], color=TEAL, linewidth=2.4)
    ax.text(0.95, legend_y, "sequence flow", va="center", color=TEAL, fontsize=9, fontweight="bold")
    ax.plot(
        [2.65, 3.10],
        [legend_y, legend_y],
        color=ORANGE,
        linewidth=2.4,
        linestyle=reference_style,
    )
    ax.text(3.20, legend_y, "database profile", va="center", color=ORANGE, fontsize=9, fontweight="bold")

    return fig


def main() -> None:
    output = Path(__file__).resolve().parents[1] / "docs" / "assets" / "figures"
    figure = build()
    metadata = {"Date": None, "Creator": "SSUextract"}
    svg = output / "pipeline-architecture.svg"
    png = output / "pipeline-architecture.png"
    figure.savefig(svg, bbox_inches="tight", facecolor=WHITE, metadata=metadata)
    svg_text = svg.read_text(encoding="utf-8")
    svg.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    figure.savefig(png, bbox_inches="tight", facecolor=WHITE, dpi=200, metadata=metadata)
    plt.close(figure)
    svg.chmod(0o644)
    png.chmod(0o644)


if __name__ == "__main__":
    main()
