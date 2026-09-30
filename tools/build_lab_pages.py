#!/usr/bin/env python3
"""Build static lab pages from the repository Markdown labs.

This intentionally supports the Markdown subset used by the course labs, without
external dependencies, so GitHub Actions can run it on the default runner.
"""

from __future__ import annotations

import html
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LABS = [
    {
        "source": ROOT / "labs/week-01/01-command-line-git.md",
        "output": ROOT / "docs/labs/week-01/01-command-line-git.html",
        "session": "Session 1",
        "summary": "Use Pax, inspect a GFF3 file, practice shell pipelines, and create a homework Git repository.",
    },
    {
        "source": ROOT / "labs/week-01/02-hpc-environments-fetchngs.md",
        "output": ROOT / "docs/labs/week-01/02-hpc-environments-fetchngs.html",
        "session": "Session 2",
        "summary": "Request compute resources, submit SLURM jobs, test containers, and start an nf-core/fetchngs download.",
    },
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


def render_markdown(markdown: str) -> tuple[str, str, str, list[tuple[int, str, str]]]:
    lines = markdown.splitlines()
    html_parts: list[str] = []
    toc: list[tuple[int, str, str]] = []
    title = "Lab"
    meta = ""
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

    def close_lists(target_depth: int = 0) -> None:
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

        if seen_h1 and not meta and line.startswith("Applied Bioinformatics"):
            meta = line
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
            html_parts.append(f"<li>{inline_markdown(text)}</li>")
            continue

        paragraph.append(line.strip())

    if code_lines is not None:
        code = html.escape("\n".join(code_lines))
        html_parts.append(f'<pre class="lab-code"><code>{code}</code></pre>')
    close_paragraph()
    close_lists()
    close_table()
    return title, meta, "\n".join(html_parts), toc


def page_template(title: str, meta: str, body: str, toc: list[tuple[int, str, str]], lab: dict[str, str]) -> str:
    toc_items = "\n".join(
        f'<a class="toc-level-{level}" href="#{anchor}">{html.escape(text)}</a>'
        for level, text, anchor in toc
        if level <= 3
    )
    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="{html.escape(title)} for Applied Bioinformatics.">
    <title>{html.escape(title)} | Applied Bioinformatics</title>
    <link rel="icon" href="../../assets/favicon.svg" type="image/svg+xml">
    <link rel="stylesheet" href="../../styles.css">
  </head>
  <body class="lab-page">
    <a class="skip-link" href="#main">Skip to main content</a>
    <header class="site-header" id="top">
      <nav class="nav" aria-label="Main navigation">
        <a class="brand" href="../../index.html#labs" aria-label="Applied Bioinformatics labs">
          <span class="brand-mark" aria-hidden="true">AB</span>
          <span>Applied Bioinformatics</span>
        </a>
        <div class="nav-links">
          <a href="../../index.html#labs">Labs</a>
          <a href="https://github.com/zhan4429/TuftsAppliedBioinfo">GitHub</a>
        </div>
      </nav>
    </header>
    <main class="lab-main" id="main">
      <div class="lab-shell">
        <aside class="lab-toc" aria-label="Lab table of contents">
          <a class="back-link" href="../../index.html#labs">Back to labs</a>
          <p class="section-kicker">{html.escape(lab["session"])}</p>
          <nav>{toc_items}</nav>
        </aside>
        <article class="lab-article">
          <header class="lab-hero">
            <p class="section-kicker">{html.escape(lab["session"])}</p>
            <h1>{html.escape(title)}</h1>
            <p>{html.escape(lab["summary"])}</p>
            <span>{html.escape(meta)}</span>
          </header>
          {body}
        </article>
      </div>
    </main>
    <footer class="site-footer">
      <p>Applied Bioinformatics - Tufts University Department of Biology</p>
      <a href="../../index.html#labs">Course labs</a>
    </footer>
  </body>
</html>
"""


def main() -> None:
    for lab in LABS:
        source = Path(lab["source"])
        output = Path(lab["output"])
        title, meta, body, toc = render_markdown(source.read_text(encoding="utf-8"))
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(page_template(title, meta, body, toc, lab), encoding="utf-8")
        print(f"Built {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
