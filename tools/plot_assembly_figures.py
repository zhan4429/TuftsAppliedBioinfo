#!/usr/bin/env python3
"""Regenerate the Week 3 assembly figures with Python and matplotlib.

The workflow and N50 example are teaching diagrams. The size/GC chart uses the
instructor's SeqKit output supplied in the course editing conversation, not
simulated contigs. See the assets README for values, derivations and references.
This script is optional: the site builder uses the checked-in images.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUTPUT = Path(__file__).resolve().parents[1] / "labs/week-03/assets/genome-assembly"
INK, MUTED, RULE = "#203348", "#5e6d7b", "#dbe3e8"
TEAL, BLUE, AMBER = "#18796f", "#316caa", "#b66d17"
PALE_TEAL, PALE_BLUE, PALE_AMBER = "#e9f4f1", "#eef4fb", "#fff3df"


def canvas(height, title, subtitle):
    fig = plt.figure(figsize=(11, height), facecolor="white")
    fig.text(.045, .93, title, fontsize=22, weight="bold", color=INK)
    fig.text(.045, .865, subtitle, fontsize=12.5, color=MUTED)
    return fig


def save(fig, name):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for extension in ("svg", "png"):
        fig.savefig(OUTPUT / f"{name}.{extension}", dpi=170, facecolor="white")
    plt.close(fig)


def box(ax, x, y, width, height, title, detail, color, fill):
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0.012,rounding_size=0.025",
        edgecolor=color, facecolor=fill, linewidth=1.2,
    ))
    ax.text(x + width / 2, y + height * .68, title, ha="center", va="center",
            fontsize=14, weight="bold", color=color)
    ax.text(x + width / 2, y + height * .30, detail, ha="center", va="center",
            fontsize=11.5, color=INK, linespacing=1.4)


def arrow(ax, start, end, color=MUTED):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=15,
                                linewidth=1.5, color=color))


def workflow():
    fig = canvas(5.8, "One HiFi input, two ways to run the assembly",
                 "Submit both jobs in class; inspect and compare the results after class.")
    ax = fig.add_axes([.045, .12, .91, .68])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    box(ax, .015, .365, .23, .29, "Same HiFi reads", "hifi.fastq.gz\nselected in Part 1", INK, "#f3f6f8")
    box(ax, .34, .62, .28, .28, "Part 1: direct", "hifiasm module\nyour batch script", TEAL, PALE_TEAL)
    box(ax, .34, .09, .28, .32, "Part 2: pipeline", "fastplong → hifiasm\nNextflow + containers", BLUE, PALE_BLUE)
    box(ax, .72, .62, .26, .28, "Measure by hand", "SeqKit: contiguity\ncompleasm: completeness", TEAL, PALE_TEAL)
    box(ax, .72, .09, .26, .32, "Read the reports", "QUAST: contiguity\nBUSCO: completeness", BLUE, PALE_BLUE)
    arrow(ax, (.25, .55), (.33, .75), TEAL)
    arrow(ax, (.25, .46), (.33, .26), BLUE)
    arrow(ax, (.635, .76), (.71, .76), TEAL)
    arrow(ax, (.635, .25), (.71, .25), BLUE)
    fig.text(.045, .055, "Both routes use hifiasm. Similar output is reproducibility evidence, not proof of accuracy.",
             fontsize=12, color=INK)
    save(fig, "two-assembly-routes")


def n50():
    lengths = [5, 3, 2, 1, 1]
    assert sum(lengths) == 12 and sum(lengths[:1]) < 6 <= sum(lengths[:2])
    fig = canvas(5.1, "N50: add contig lengths from longest to shortest",
                 "Worked example only: five contigs, 12 Mb in total. Half of the assembly is 6 Mb.")
    ax = fig.add_axes([.13, .22, .44, .54])
    ax.barh(range(5), lengths, color=[TEAL, TEAL, "#bdcbd5", "#bdcbd5", "#bdcbd5"], height=.62)
    ax.set(yticks=range(5), yticklabels=[f"Contig {i}" for i in range(1, 6)],
           xlim=(0, 5.8), xticks=range(6))
    ax.invert_yaxis()
    ax.set_xlabel("Contig length (Mb)", color=MUTED, fontsize=12)
    ax.tick_params(axis="both", length=0, labelsize=11, labelcolor=MUTED)
    ax.grid(axis="x", color=RULE, linewidth=.7)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)
    for y, length in enumerate(lengths):
        ax.text(length + .12, y, f"{length} Mb", va="center", color=INK, fontsize=12)
    fig.text(.64, .70, "1. Start with the longest", fontsize=13, color=INK, weight="bold")
    fig.text(.64, .635, "5 Mb: below the 6 Mb halfway point", fontsize=12, color=MUTED)
    fig.text(.64, .52, "2. Add the next contig", fontsize=13, color=INK, weight="bold")
    fig.text(.64, .455, "5 + 3 = 8 Mb: reaches at least half", fontsize=12, color=MUTED)
    fig.text(.64, .32, "N50 = 3 Mb", fontsize=22, color=TEAL, weight="bold")
    fig.text(.64, .245, "L50 = 2 contigs", fontsize=15, color=TEAL)
    fig.text(.045, .065, "N50 describes contiguity. Long contigs can still contain errors or omit genes.",
             fontsize=12, color=INK)
    save(fig, "n50-worked-example")


def measured_size_gc():
    total_bp, long_bp = 17_499_613, 13_298_638
    total_gc, long_gc = 35.57, 38.84
    short_bp = total_bp - long_bp
    short_gc = (total_bp * total_gc - long_bp * long_gc) / short_bp
    assert short_bp == 4_200_975
    fig = canvas(5.6, "Where is the extra assembled sequence?",
                 "Instructor's hifiasm 0.25.0 run: 134 contigs, 17.50 Mb, overall GC 35.57%.")
    ax = fig.add_axes([.07, .45, .86, .32])
    ax.barh([0], [long_bp / 1e6], color=TEAL, height=.6)
    ax.barh([0], [short_bp / 1e6], left=[long_bp / 1e6], color=AMBER, height=.6)
    ax.set(xlim=(0, 18.2), ylim=(-.7, .7), yticks=[], xticks=[0, 3, 6, 9, 12, 15, 18])
    ax.set_xlabel("Assembled sequence (Mb)", color=MUTED, fontsize=12)
    ax.tick_params(axis="x", length=0, labelsize=11, labelcolor=MUTED)
    ax.grid(axis="x", color=RULE, linewidth=.7)
    ax.set_axisbelow(True)
    for side in ax.spines.values():
        side.set_visible(False)
    ax.text(long_bp / 2e6, 0, "13.30 Mb", ha="center", va="center", color="white", fontsize=18, weight="bold")
    ax.text((long_bp + short_bp / 2) / 1e6, 0, "4.20 Mb", ha="center", va="center", color="white", fontsize=17, weight="bold")
    ax.axvline(12.1, color=INK, linestyle=(0, (3, 3)), linewidth=1.5)
    ax.text(12.1, .52, "Expected size: 12.1 Mb", ha="right", fontsize=11.5, color=INK)
    fig.text(.07, .30, "25 contigs ≥100 kb", color=TEAL, weight="bold", fontsize=16)
    fig.text(.07, .235, "Measured GC: 38.84%", color=INK, fontsize=13)
    fig.text(.61, .30, "109 contigs <100 kb", color=AMBER, weight="bold", fontsize=16)
    fig.text(.61, .235, f"Estimated GC: ~{short_gc:.1f}%", color=INK, fontsize=13)
    fig.text(.07, .12, "Length groups are diagnostic subsets, not assignments of nuclear or mitochondrial identity.",
             fontsize=11.5, color=INK)
    fig.text(.07, .065, "Short-contig GC is inferred from the rounded SeqKit summaries; no contig-level depths are plotted.",
             fontsize=10.5, color=MUTED)
    save(fig, "observed-assembly-size-gc")


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "font.size": 12})
    workflow()
    n50()
    measured_size_gc()


if __name__ == "__main__":
    main()
