#!/usr/bin/env python3
"""Create the Lab 1 teaching figure with rounded, cited genome-wide values.

Yeast: Goffeau et al. (1996), https://pubmed.ncbi.nlm.nih.gov/8849441/
Human: NIGMS, https://www.nigms.nih.gov/biobeat/2024/04/genetics-by-the-numbers
Counts and genome sizes are intentionally rounded, not release-specific estimates.
Run with a Python environment containing matplotlib, then rebuild the course site.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


OUTPUT = Path(__file__).resolve().parents[1] / "labs/week-01/assets"
INK, MUTED, RULE = "#1b2533", "#5d6876", "#d9dfe7"
YEAST, HUMAN = "#237c72", "#2e6aa6"
GENES = [6000, 20000]
SIZES_MB = [12.1, 3200]
DENSITY = [genes / size for genes, size in zip(GENES, SIZES_MB)]


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
    fig = plt.figure(figsize=(11, 8.6), facecolor="white")
    fig.text(0.06, 0.945, "Yeast and human: gene density", fontsize=23, color=INK, weight="bold")
    fig.text(0.06, 0.905, "Protein-coding genes per megabase of genomic DNA", fontsize=13, color=MUTED)

    for x, label in [(0.06, "SPECIES"), (0.33, "APPROX. GENES"), (0.55, "GENOME SIZE"), (0.78, "GENES / Mb")]:
        fig.text(x, 0.845, label, fontsize=10, color=MUTED, weight="bold")
    for y, species, count, size, density, color in [
        (0.795, "Yeast (S. cerevisiae)", "6,000", "12.1 Mb", "~500", YEAST),
        (0.735, "Human", "20,000", "3,200 Mb", "~6", HUMAN),
    ]:
        for x, value in [(0.06, species), (0.33, count), (0.55, size), (0.78, density)]:
            fig.text(x, y, value, fontsize=13, color=color, weight="bold" if x == 0.78 else "normal")

    ax = fig.add_axes([0.20, 0.49, 0.69, 0.17])
    ax.barh([1, 0], DENSITY, color=[YEAST, HUMAN], height=0.44)
    ax.set_yticks([1, 0], ["Yeast", "Human"], fontsize=12, color=INK)
    ax.set_xlim(0, 550)
    ax.set_xticks(range(0, 501, 100))
    ax.tick_params(axis="both", length=0, labelcolor=MUTED)
    ax.set_xlabel("Protein-coding genes per Mb (same linear scale)", color=MUTED, fontsize=11, labelpad=10)
    ax.grid(axis="x", color=RULE, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    for y, density, label in [(1, DENSITY[0], "~500"), (0, DENSITY[1], "~6")]:
        ax.text(density + 8, y, label, va="center", color=INK, fontsize=12, weight="bold")

    fig.text(0.06, 0.405, "Yeast has about 80 times as many protein-coding genes per Mb.", fontsize=13, color=INK, weight="bold")
    fig.text(0.06, 0.373, "Density = gene count / genome size. This is different from the fraction covered by exons.", fontsize=11, color=MUTED)
    fig.text(0.06, 0.316, "How genes can be packed", fontsize=14, color=INK, weight="bold")
    fig.text(0.06, 0.287, "Schematic only: positions and lengths below are not to scale.", fontsize=10.5, color=MUTED)

    cartoon = fig.add_axes([0.06, 0.12, 0.88, 0.145])
    cartoon.set_xlim(0, 1000)
    cartoon.set_ylim(0, 145)
    cartoon.axis("off")
    for y, name, color in [(112, "Yeast", YEAST), (42, "Human", HUMAN)]:
        cartoon.text(0, y, name, va="center", fontsize=12, color=color, weight="bold")
        cartoon.plot([140, 970], [y, y], color=RULE, lw=2, zorder=1)
    for start, end in [(150, 220), (247, 305), (330, 412), (440, 503), (530, 596), (625, 710), (738, 802), (830, 941)]:
        cartoon.add_patch(Rectangle((start, 102), end - start, 20, facecolor=YEAST, zorder=3))
    for exons in [[(155, 175), (225, 238), (310, 335)], [(515, 535), (596, 610), (665, 685)], [(825, 843), (890, 905), (935, 955)]]:
        cartoon.plot([exons[0][0], exons[-1][1]], [42, 42], color=HUMAN, lw=1.6, zorder=2)
        for start, end in exons:
            cartoon.add_patch(Rectangle((start, 32), end - start, 20, facecolor=HUMAN, zorder=3))
    cartoon.annotate("intron", xy=(275, 42), xytext=(275, 79), ha="center", fontsize=10, color=MUTED,
                     arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
    cartoon.annotate("intergenic DNA", xy=(425, 42), xytext=(425, 3), ha="center", fontsize=10, color=MUTED,
                     arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})

    fig.text(0.06, 0.087, "Filled boxes: exons     Thin colored lines: introns     Gray line: intergenic DNA", fontsize=10.5, color=MUTED)
    fig.text(0.06, 0.048, "Rounded teaching values; actual counts depend on the genome assembly and annotation release.", fontsize=9.5, color=MUTED)
    fig.text(0.06, 0.027, "Sources: Goffeau et al., Science (1996); NIGMS, Genetics by the Numbers (2024).", fontsize=9.5, color=MUTED)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "svg"]:
        fig.savefig(OUTPUT / f"yeast-human-gene-density.{extension}", dpi=200, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
