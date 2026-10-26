#!/usr/bin/env python3
"""Build static course pages from the repository Markdown files.

This intentionally supports the Markdown subset used by the course labs and
tutorials, without external dependencies, so GitHub Actions can run it on the
default runner.
"""

from __future__ import annotations

import html
import os
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"

SUMMARY_OVERRIDES = {
    "labs/week-01/01-command-line-git.md": (
        "Use Tufts HPC, inspect a GFF3 file, practice shell pipelines, and create a homework Git repository."
    ),
    "labs/week-01/02-hpc-environments-fetchngs.md": (
        "Request compute resources, submit SLURM jobs, test containers, and start an nf-core/fetchngs download."
    ),
    "tutorials/guide-git-github-on-pax.md": (
        "Configure Git identity, add an SSH key to GitHub, and connect Pax repositories safely."
    ),
}


REPO_URL = "https://github.com/zhan4429/TuftsAppliedBioinfo"

COURSE_TITLE = "Applied Bioinformatics"
COURSE_LEDE = (
    "Hands-on bioinformatics with real datasets: reproducible workflows on the "
    "Tufts HPC cluster, critical evaluation of results, and clear scientific "
    "reporting. Weekly labs and setup guides are posted here."
)

# Shown under the course description on the home page. Fill in what applies;
# entries left as "" are not displayed.
COURSE_DETAILS = [
    ("Course", ""),          # course number, e.g. as listed in SIS
    ("Term", ""),
    ("Instructor", ""),
    ("Class meetings", ""),
    ("Computing", "Tufts HPC cluster (Pax)"),
]

# Week-by-week schedule for the home page. Lab pages are attached to a row
# automatically from the week number in their path (labs/week-NN/...).
SCHEDULE = [
    {
        "label": "0",
        "weeks": [0],
        "title": "Optional preparation",
        "topics": "Self-paced practice with shell navigation and basic Python or R syntax.",
    },
    {
        "label": "1",
        "weeks": [1],
        "title": "Command line, HPC, and reproducible environments",
        "topics": (
            "Linux command line, text processing, Git/GitHub basics, HPC job submission with "
            "SLURM, conda/mamba environments, containers, public databases, and dataset inspection."
        ),
    },
    {
        "label": "2",
        "weeks": [2],
        "title": "Sequencing design and quality control",
        "topics": (
            "Sequencing technologies, experimental design, replication, confounding, batch effects, "
            "FASTQ format, FastQC, fastp, MultiQC, and QC interpretation."
        ),
    },
    {
        "label": "3",
        "weeks": [3],
        "title": "Genome assembly and annotation",
        "topics": (
            "Assembly QC, contamination screening, GFF3/GTF structure, functional annotation, "
            "annotation version issues, and IGV inspection."
        ),
    },
    {
        "label": "4",
        "weeks": [4],
        "title": "Read alignment and variant calling",
        "topics": (
            "SAM/BAM processing, alignment QC, variant calling with GATK HaplotypeCaller and "
            "DeepVariant, filtering, benchmarking, annotation, and VCF interpretation."
        ),
    },
    {
        "label": "5",
        "weeks": [5],
        "title": "Bulk RNA-seq quantification",
        "topics": (
            "RNA-seq design, pseudoalignment, alignment-based quantification, library type and "
            "strandedness, nf-core/rnaseq, count matrices, metadata, and exploratory analysis."
        ),
    },
    {
        "label": "6",
        "weeks": [6],
        "title": "Differential expression",
        "topics": (
            "DESeq2, design formulas, covariates and batch effects, multiple testing, "
            "visualization, enrichment analysis, and notebook/report preparation."
        ),
    },
    {
        "label": "7",
        "weeks": [7],
        "title": "Chromatin accessibility and binding",
        "topics": (
            "ATAC-seq, peak calling with MACS3, nf-core/atacseq, ChIP-seq and CUT&RUN/CUT&Tag "
            "concepts, QC metrics, peak annotation, motif enrichment, and integration with "
            "expression results."
        ),
    },
    {
        "label": "8",
        "weeks": [8],
        "title": "Single-cell RNA-seq preprocessing",
        "topics": (
            "10x barcodes and UMIs, Cell Ranger, STARsolo, alevin-fry, nf-core/scrnaseq, AnnData "
            "and Seurat objects, QC, filtering, normalization, and feature selection."
        ),
    },
    {
        "label": "9",
        "weeks": [9],
        "title": "Single-cell analysis and spatial transcriptomics",
        "topics": (
            "Dimensionality reduction, clustering, marker genes, cell type annotation, batch "
            "integration, pseudobulk differential expression, and spatial transcriptomics."
        ),
    },
    {
        "label": "10",
        "weeks": [10],
        "title": "Microbiome and metagenomics",
        "topics": (
            "16S amplicon workflows, QIIME2/DADA2, compositional data, shotgun taxonomic and "
            "functional profiling, metagenome-assembled genomes, diversity metrics, and "
            "differential abundance."
        ),
    },
    {
        "label": "11",
        "weeks": [11],
        "title": "Mini-project preparation",
        "topics": "Dataset acquisition and QC as assigned.",
    },
    {
        "label": "12\u201314",
        "weeks": [12, 13, 14],
        "title": "Mini-project",
        "topics": (
            "Reproduce key findings from a published paper, maintain a reproducible GitHub "
            "repository, prepare notebooks and reports, and present results."
        ),
    },
]

REPRODUCIBILITY_POINTS = [
    "Keep commands, parameters, sample IDs, and software versions visible.",
    "Use relative paths when possible.",
    "Separate raw data, intermediate files, final outputs, and reports.",
    "Document where input data came from.",
    "Preserve small configuration files and sample sheets.",
    "Avoid committing large generated files unless instructed.",
    "Record failed or unexpected results when they affect interpretation.",
]

AI_USE_POINTS = [
    "Disclose AI use in submitted notebooks or reports, including the tool used and what it helped with.",
    "Verify generated commands and code against documentation, error messages, and expected outputs.",
    "You are responsible for every line you submit.",
    "Written interpretation, discussion, and limitations sections should be your own scientific reasoning.",
    "Watch for hallucinated arguments, deprecated APIs, plausible but inappropriate statistical "
    "choices, and code that runs while producing the wrong analysis.",
]


def source_key(source: Path) -> str:
    return source.relative_to(ROOT).as_posix()


def natural_sort_key(path: Path) -> list[int | str]:
    parts: list[int | str] = []
    for piece in re.split(r"(\d+)", path.relative_to(ROOT).as_posix().lower()):
        parts.append(int(piece) if piece.isdigit() else piece)
    return parts


def strip_inline_markdown(text: str) -> str:
    text = text.strip()
    text = re.sub(r"^>\s*", "", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"_([^_]+)_", r"\1", text)
    text = re.sub(r"\*([^*]+)\*", r"\1", text)
    return html.unescape(text).strip()


def is_metadata_line(text: str) -> bool:
    normalized = text.strip("*_ ")
    return normalized.startswith("Applied Bioinformatics") or normalized.startswith("Time:")


def shorten(text: str, limit: int = 180) -> str:
    text = re.sub(r"\s+", " ", text).strip()
    if len(text) <= limit:
        return text
    return text[:limit].rsplit(" ", 1)[0].rstrip(".,;:") + "..."


def extract_summary(markdown: str, fallback: str) -> str:
    in_code = False
    seen_h1 = False
    paragraph: list[str] = []

    for raw in markdown.splitlines():
        line = raw.strip()

        if line.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        if not line:
            if paragraph:
                return shorten(" ".join(paragraph))
            continue
        if line.startswith("# "):
            seen_h1 = True
            continue
        if not seen_h1:
            continue
        if re.fullmatch(r"[-*_]{3,}", line):
            continue
        if line.startswith("<") or line.startswith("|"):
            continue

        plain = strip_inline_markdown(line)
        if not plain or is_metadata_line(plain):
            continue
        paragraph.append(plain)

    if paragraph:
        return shorten(" ".join(paragraph))
    return fallback


def lab_week_session(source: Path, title: str) -> tuple[int | None, int | None]:
    week_match = re.search(r"week-(\d+)", source.as_posix(), re.IGNORECASE)
    session_match = re.search(r"session[-_ ]?(\d+)", source.stem, re.IGNORECASE)
    if not session_match:
        session_match = re.search(r"week[-_ ]?\d+[-_ ]?(\d+)", source.stem, re.IGNORECASE)
    if not session_match:
        session_match = re.search(r"week\s+\d+\s+session\s+(\d+)", title, re.IGNORECASE)
    week = int(week_match.group(1)) if week_match else None
    session = int(session_match.group(1)) if session_match else None
    return week, session


def lab_label(week: int | None, session: int | None) -> str:
    if week is not None and session is not None:
        return f"Week {week}, Session {session}"
    if week is not None:
        return f"Week {week}"
    return "Lab"


def display_title(title: str) -> str:
    title = re.sub(r"^Week\s+\d+\s+Session\s+\d+:\s*", "", title, flags=re.IGNORECASE)
    title = re.sub(r"^Lab\s+\d+[A-Za-z]?\s+[-\u2013\u2014]\s*", "", title, flags=re.IGNORECASE)
    title = re.sub(r"^Tutorial:\s*", "", title, flags=re.IGNORECASE)
    return title


def discover_lab_pages() -> list[dict[str, str | Path]]:
    sources = sorted(
        (
            path
            for path in (ROOT / "labs").glob("**/*.md")
            if path.name.lower() != "readme.md"
        ),
        key=natural_sort_key,
    )
    return [
        {
            "source": source,
            "output": DOCS / source.relative_to(ROOT).with_suffix(".html"),
            "back_label": "Back to schedule",
            "back_fragment": "labs",
            "section": "labs",
        }
        for source in sources
    ]


def discover_tutorial_pages() -> list[dict[str, str | Path]]:
    sources = sorted((ROOT / "tutorials").glob("*.md"), key=natural_sort_key)
    return [
        {
            "source": source,
            "output": DOCS / source.relative_to(ROOT).with_suffix(".html"),
            "label": "Tutorial",
            "back_label": "Back to setup and tutorials",
            "back_fragment": "setup",
            "section": "tutorials",
        }
        for source in sources
    ]


def slugify(text: str) -> str:
    text = re.sub(r"<[^>]+>", "", text)
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "-", text)
    return text.strip("-") or "section"


def inline_markdown(text: str) -> str:
    placeholders: list[str] = []

    def stash(value: str) -> str:
        placeholders.append(value)
        return f"\x00{len(placeholders) - 1}\x00"

    def code_repl(match: re.Match[str]) -> str:
        return stash(f"<code>{html.escape(match.group(1))}</code>")

    def link_repl(match: re.Match[str]) -> str:
        label = inline_markdown(match.group(1))
        href = html.escape(match.group(2), quote=True)
        return stash(f'<a href="{href}">{label}</a>')

    escaped = html.escape(text)
    escaped = re.sub(r"`([^`]+)`", code_repl, escaped)
    escaped = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link_repl, escaped)
    escaped = re.sub(
        r"&lt;(https?://[^&]+)&gt;",
        lambda match: stash(
            f'<a href="{html.escape(match.group(1), quote=True)}">{html.escape(match.group(1))}</a>'
        ),
        escaped,
    )
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    escaped = re.sub(r"_([^_]+)_", r"<em>\1</em>", escaped)
    escaped = re.sub(r"\*([^*]+)\*", r"<em>\1</em>", escaped)

    for index, value in enumerate(placeholders):
        escaped = escaped.replace(f"\x00{index}\x00", value)
    return escaped


def render_table(rows: list[str]) -> str:
    split_rows = []
    for row in rows:
        cells = [cell.strip().replace(r"\|", "|") for cell in row.strip().strip("|").split("|")]
        split_rows.append(cells)
    if len(split_rows) < 2:
        return "\n".join(rows)

    header = split_rows[0]
    body = split_rows[2:]
    out = ["<div class=\"table-wrap\"><table>", "<thead><tr>"]
    out.extend(f"<th>{inline_markdown(cell)}</th>" for cell in header)
    out.append("</tr></thead>")
    if body:
        out.append("<tbody>")
        for row in body:
            out.append("<tr>")
            out.extend(f"<td>{inline_markdown(cell)}</td>" for cell in row)
            out.append("</tr>")
        out.append("</tbody>")
    out.append("</table></div>")
    return "\n".join(out)


def render_markdown(markdown: str) -> tuple[str, str, str, str, list[tuple[int, str, str]]]:
    lines = markdown.splitlines()
    html_parts: list[str] = []
    toc: list[tuple[int, str, str]] = []
    title = "Lab"
    meta = ""
    time_note = ""
    paragraph: list[str] = []
    list_stack: list[str] = []
    code_lines: list[str] | None = None
    code_lang = ""
    table_lines: list[str] = []
    seen_h1 = False

    def close_paragraph() -> None:
        nonlocal paragraph
        if paragraph:
            html_parts.append(f"<p>{inline_markdown(' '.join(paragraph))}</p>")
            paragraph = []

    current_item: list[str] = []

    def flush_item() -> None:
        nonlocal current_item
        if current_item:
            html_parts.append(f"<li>{inline_markdown(' '.join(current_item))}</li>")
            current_item = []

    def close_lists(target_depth: int = 0) -> None:
        flush_item()
        while len(list_stack) > target_depth:
            html_parts.append(f"</{list_stack.pop()}>")

    def close_table() -> None:
        nonlocal table_lines
        if table_lines:
            html_parts.append(render_table(table_lines))
            table_lines = []

    def heading_anchor(text: str) -> str:
        base = slugify(text)
        used = {anchor for _level, _text, anchor in toc}
        anchor = base
        counter = 2
        while anchor in used:
            anchor = f"{base}-{counter}"
            counter += 1
        return anchor

    for raw in lines:
        line = raw.rstrip()

        if code_lines is not None:
            if line.startswith("```"):
                lang_class = f" language-{html.escape(code_lang)}" if code_lang else ""
                code = html.escape("\n".join(code_lines))
                html_parts.append(f'<pre class="lab-code"><code class="{lang_class}">{code}</code></pre>')
                code_lines = None
                code_lang = ""
            else:
                code_lines.append(raw)
            continue

        if line.startswith("```"):
            close_paragraph()
            close_lists()
            close_table()
            code_lines = []
            code_lang = line[3:].strip()
            continue

        if line.startswith("|") and line.endswith("|"):
            close_paragraph()
            close_lists()
            table_lines.append(line)
            continue
        else:
            close_table()

        if not line.strip():
            close_paragraph()
            close_lists()
            continue

        if re.fullmatch(r"[-*_]{3,}", line.strip()):
            close_paragraph()
            close_lists()
            close_table()
            html_parts.append("<hr>")
            continue

        if line.startswith("<details>") or line.startswith("</details>"):
            close_paragraph()
            close_lists()
            html_parts.append(line)
            continue

        if line.startswith("<summary>") or line.startswith("</summary>"):
            close_paragraph()
            close_lists()
            html_parts.append(line)
            continue

        heading = re.match(r"^(#{1,6})\s+(.+)$", line)
        if heading:
            close_paragraph()
            close_lists()
            level = len(heading.group(1))
            text = heading.group(2).strip()
            if level == 1 and not seen_h1:
                title = re.sub(r"<[^>]+>", "", text)
                seen_h1 = True
                continue
            anchor = heading_anchor(text)
            if level <= 3:
                toc.append((level, re.sub(r"`", "", text), anchor))
            html_parts.append(f'<h{level} id="{anchor}">{inline_markdown(text)}</h{level}>')
            continue

        plain_line = strip_inline_markdown(line)
        if seen_h1 and not meta and plain_line.startswith("Applied Bioinformatics"):
            meta = plain_line
            continue
        if seen_h1 and not time_note and not html_parts and not paragraph and plain_line.startswith("Time:"):
            time_note = plain_line
            continue

        quote = re.match(r"^>\s*(.+)$", line)
        if quote:
            close_paragraph()
            close_lists()
            html_parts.append(f"<blockquote>{inline_markdown(quote.group(1))}</blockquote>")
            continue

        ordered = re.match(r"^(\d+)\.\s+(.+)$", line)
        unordered = re.match(r"^-\s+(.+)$", line)
        if ordered or unordered:
            close_paragraph()
            kind = "ol" if ordered else "ul"
            text = ordered.group(2) if ordered else unordered.group(1)
            if not list_stack or list_stack[-1] != kind:
                close_lists()
                html_parts.append(f"<{kind}>")
                list_stack.append(kind)
            flush_item()
            current_item = [text]
            continue

        if current_item:
            # Wrapped continuation line of the current list item.
            current_item.append(line.strip())
            continue

        paragraph.append(line.strip())

    if code_lines is not None:
        code = html.escape("\n".join(code_lines))
        html_parts.append(f'<pre class="lab-code"><code>{code}</code></pre>')
    close_paragraph()
    close_lists()
    close_table()
    return title, meta, time_note, "\n".join(html_parts), toc


def relative_prefix(output: Path) -> str:
    prefix = Path(*([".."] * len(output.parent.relative_to(DOCS).parts)))
    return "" if str(prefix) == "." else prefix.as_posix()


def prefixed(prefix: str, path: str) -> str:
    return f"{prefix}/{path}" if prefix else path


def esc(text: str) -> str:
    return html.escape(text, quote=True)


def head_markup(title: str, description: str, prefix: str, with_code_copy: bool) -> str:
    favicon_href = prefixed(prefix, "assets/favicon.svg")
    styles_href = prefixed(prefix, "styles.css")
    font_href = prefixed(prefix, "assets/fonts/source-serif-4-latin-wght-normal.woff2")
    script = ""
    if with_code_copy:
        script = f'\n    <script src="{prefixed(prefix, "assets/lab-code-copy.js")}" defer></script>'
    return f"""  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="{esc(description)}">
    <title>{esc(title)}</title>
    <link rel="icon" href="{favicon_href}" type="image/svg+xml">
    <link rel="preload" href="{font_href}" as="font" type="font/woff2" crossorigin>
    <link rel="stylesheet" href="{styles_href}">{script}
  </head>"""


def institution_bar() -> str:
    return """    <div class="institution-bar">
      <div class="wrap">
        <a href="https://www.tufts.edu">Tufts University</a>
        <span>Department of Biology</span>
      </div>
    </div>"""


def site_footer(prefix: str) -> str:
    return f"""    <footer class="site-footer">
      <div class="wrap">
        <p>
          {esc(COURSE_TITLE)}, Department of Biology, Tufts University.
          Course policies, grading, and due dates are set by the official syllabus.
        </p>
        <a href="{REPO_URL}">Course repository on GitHub</a>
      </div>
    </footer>"""


def page_template(page: dict) -> str:
    prefix = relative_prefix(Path(page["output"]))
    index_href = prefixed(prefix, "index.html")
    back_href = prefixed(prefix, f'index.html#{page["back_fragment"]}')
    is_lab_page = page["section"] == "labs"
    title = str(page["title"])
    heading = display_title(title)
    toc_items = "\n".join(
        f'            <a class="toc-level-{level}" href="#{anchor}">{esc(text)}</a>'
        for level, text, anchor in page["toc"]
        if level <= 3
    )
    time_note = f'\n            <p class="lab-time">{esc(page["time_note"])}</p>' if page["time_note"] else ""

    pager = ""
    if is_lab_page and (page.get("prev") or page.get("next")):
        links = []
        for key, label, cls in (("prev", "Previous lab", "pager-prev"), ("next", "Next lab", "pager-next")):
            other = page.get(key)
            if not other:
                links.append(f'<span class="{cls}"></span>')
                continue
            href = Path(os.path.relpath(other["output"], Path(page["output"]).parent)).as_posix()
            links.append(
                f'<a class="{cls}" href="{esc(href)}"><span>{label}: {esc(str(other["label"]))}</span>'
                f'{esc(display_title(str(other["title"])))}</a>'
            )
        pager = '\n          <nav class="lab-pager" aria-label="Lab sequence">\n            ' + "\n            ".join(links) + "\n          </nav>"

    return f"""<!doctype html>
<html lang="en">
{head_markup(f"{title} | {COURSE_TITLE}", f"{title} for {COURSE_TITLE}, Tufts University.", prefix, True)}
  <body class="doc-page">
    <a class="skip-link" href="#main">Skip to main content</a>
{institution_bar()}
    <header class="site-header">
      <div class="wrap site-header-inner">
        <a class="site-title" href="{index_href}">{esc(COURSE_TITLE)}</a>
        <nav class="site-nav" aria-label="Main navigation">
          <a href="{prefixed(prefix, 'index.html#labs')}">Schedule</a>
          <a href="{prefixed(prefix, 'index.html#setup')}">Setup</a>
          <a href="{REPO_URL}">GitHub</a>
        </nav>
      </div>
    </header>
    <main class="lab-main wrap" id="main">
      <div class="lab-shell">
        <aside class="lab-toc" aria-label="Table of contents">
          <a class="back-link" href="{back_href}">{esc(page["back_label"])}</a>
          <details class="toc-details" open>
            <summary>On this page</summary>
            <nav>
{toc_items}
            </nav>
          </details>
        </aside>
        <article class="lab-article">
          <header class="lab-hero">
            <p class="lab-label">{esc(str(page["label"]))}</p>
            <h1>{esc(heading)}</h1>{time_note}
          </header>
          {page["body"]}{pager}
        </article>
      </div>
    </main>
{site_footer(prefix)}
    <script>
      (function () {{
        var toc = document.querySelector(".toc-details");
        if (toc && window.matchMedia("(max-width: 860px)").matches) toc.removeAttribute("open");
      }})();
    </script>
  </body>
</html>
"""


def course_details_markup() -> str:
    items = [(label, value) for label, value in COURSE_DETAILS if value]
    if not items:
        return ""
    rows = "\n".join(
        f"          <div><dt>{esc(label)}</dt><dd>{esc(value)}</dd></div>" for label, value in items
    )
    return f'\n        <dl class="course-details">\n{rows}\n        </dl>'


def schedule_row(label: str, title: str, topics: str, labs: list[dict]) -> str:
    status = "has-labs" if labs else "upcoming"
    lab_items = "\n".join(
        f'                <li><span class="session">'
        f'{esc(f"Session {lab["session"]}" if lab["session"] is not None else "Lab")}</span>'
        f'<a href="{esc(Path(lab["output"]).relative_to(DOCS).as_posix())}">'
        f'{esc(display_title(str(lab["title"])))}</a></li>'
        for lab in labs
    )
    labs_cell = f'\n              <ul class="session-list">\n{lab_items}\n              </ul>\n            ' if labs else ""
    topics_markup = f'\n              <p class="topic-detail">{esc(topics)}</p>' if topics else ""
    return f"""          <tr class="{status}">
            <th scope="row" class="col-week"><span class="week-num">{esc(label)}</span></th>
            <td class="col-topic">
              <span class="topic-title">{esc(title)}</span>{topics_markup}
            </td>
            <td class="col-labs">{labs_cell}</td>
          </tr>"""


def schedule_markup(lab_pages: list[dict]) -> str:
    by_week: dict[int | None, list[dict]] = {}
    for lab in lab_pages:
        by_week.setdefault(lab["week"], []).append(lab)

    rows = []
    covered: set[int] = set()
    for entry in SCHEDULE:
        covered.update(entry["weeks"])
        labs = [lab for week in entry["weeks"] for lab in by_week.get(week, [])]
        rows.append(schedule_row(entry["label"], entry["title"], entry["topics"], labs))

    # Labs for weeks not listed in SCHEDULE still get a row, so nothing is hidden.
    for week in sorted(w for w in by_week if w is not None and w not in covered):
        rows.append(schedule_row(str(week), f"Week {week}", "", by_week[week]))
    if None in by_week:
        rows.append(schedule_row("", "Other labs", "", by_week[None]))
    return "\n".join(rows)


def tutorial_markup(tutorial_pages: list[dict]) -> str:
    if not tutorial_pages:
        return '            <p class="muted">Setup guides will be posted here.</p>'
    items = "\n".join(
        f'              <li><a href="{esc(Path(t["output"]).relative_to(DOCS).as_posix())}">'
        f'{esc(display_title(str(t["title"])))}</a><p>{esc(str(t["summary"]))}</p></li>'
        for t in tutorial_pages
    )
    return f'            <ul class="resource-list">\n{items}\n            </ul>'


def bullet_list(points: list[str]) -> str:
    return "\n".join(f"              <li>{esc(point)}</li>" for point in points)


def index_template(lab_pages: list[dict], tutorial_pages: list[dict]) -> str:
    return f"""<!doctype html>
<html lang="en">
{head_markup(f"{COURSE_TITLE} | Tufts University", f"{COURSE_TITLE}: course schedule, labs, and setup guides. Department of Biology, Tufts University.", "", True)}
  <body class="home-page">
    <a class="skip-link" href="#main">Skip to main content</a>
{institution_bar()}
    <header class="course-header">
      <div class="wrap course-header-main">
        <h1 class="course-title">{esc(COURSE_TITLE)}</h1>
        <p class="course-lede">{esc(COURSE_LEDE)}</p>{course_details_markup()}
      </div>
    </header>
    <nav class="course-nav" aria-label="Main navigation">
      <div class="wrap">
        <a href="#labs">Schedule</a>
        <a href="#setup">Setup</a>
        <a href="#expectations">Expectations</a>
        <a href="{REPO_URL}">GitHub</a>
      </div>
    </nav>

    <main id="main" class="wrap home-main">
      <section class="home-section" id="labs" aria-labelledby="schedule-title">
        <h2 id="schedule-title">Schedule and labs</h2>
        <p class="section-note">
          Lab links appear in each week as they are posted. Due dates, grading,
          and course policies are set by the official syllabus.
        </p>
        <div class="table-scroll">
          <table class="schedule">
            <thead>
              <tr>
                <th scope="col" class="col-week">Week</th>
                <th scope="col" class="col-topic">Topics</th>
                <th scope="col" class="col-labs">Labs</th>
              </tr>
            </thead>
            <tbody>
{schedule_markup(lab_pages)}
            </tbody>
          </table>
        </div>
      </section>

      <section class="home-section" id="setup" aria-labelledby="setup-title">
        <h2 id="setup-title">Setup and tutorials</h2>
        <div class="two-col setup-cols">
          <div>
            <h3>Get the course repository</h3>
            <p>Clone the repository once on Pax:</p>
            <pre class="lab-code"><code>git clone {REPO_URL}.git</code></pre>
            <p>Then pull before each class to get the latest labs:</p>
            <pre class="lab-code"><code>cd TuftsAppliedBioinfo
git pull</code></pre>
            <p class="muted">
              Keep large data and generated results on Pax or course storage, not in GitHub.
            </p>
          </div>
          <div>
            <h3 id="tutorials">Tutorials</h3>
{tutorial_markup(tutorial_pages)}
          </div>
        </div>
      </section>

      <section class="home-section" id="expectations" aria-labelledby="expectations-title">
        <h2 id="expectations-title">Course expectations</h2>
        <div class="two-col">
          <div>
            <h3>Reproducible work</h3>
            <p>
              Bioinformatics work should be understandable and rerunnable by someone else,
              including your future self. A pipeline result is not automatically correct
              because it ran successfully.
            </p>
            <ul class="plain-list">
{bullet_list(REPRODUCIBILITY_POINTS)}
            </ul>
          </div>
          <div>
            <h3>Using AI assistants</h3>
            <p>
              AI coding assistants are permitted and expected as part of the course,
              but they do not replace scientific judgment.
            </p>
            <ul class="plain-list">
{bullet_list(AI_USE_POINTS)}
            </ul>
          </div>
        </div>
      </section>
    </main>

{site_footer("")}
  </body>
</html>
"""


def main() -> None:
    pages = [*discover_lab_pages(), *discover_tutorial_pages()]

    for page in pages:
        source = Path(page["source"])
        markdown = source.read_text(encoding="utf-8")
        title, _meta, time_note, body, toc = render_markdown(markdown)
        page.update(title=title, time_note=time_note, body=body, toc=toc)
        if page["section"] == "labs":
            week, session = lab_week_session(source, title)
            page.update(week=week, session=session, label=lab_label(week, session))
        page["summary"] = SUMMARY_OVERRIDES.get(
            source_key(source),
            extract_summary(markdown, f"Open {display_title(title)} for this course session."),
        )

    lab_pages = [page for page in pages if page["section"] == "labs"]
    tutorial_pages = [page for page in pages if page["section"] == "tutorials"]
    for i, page in enumerate(lab_pages):
        page["prev"] = lab_pages[i - 1] if i > 0 else None
        page["next"] = lab_pages[i + 1] if i + 1 < len(lab_pages) else None

    for page in pages:
        output = Path(page["output"])
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(page_template(page), encoding="utf-8")
        print(f"Built {output.relative_to(ROOT)}")

    index_output = DOCS / "index.html"
    index_output.write_text(index_template(lab_pages, tutorial_pages), encoding="utf-8")
    print(f"Built {index_output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
