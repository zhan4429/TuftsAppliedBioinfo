#!/usr/bin/env python3
"""Build static course pages from the repository Markdown files.

This intentionally supports the Markdown subset used by the course labs and
tutorials, without external dependencies, so GitHub Actions can run it on the
default runner.
"""

from __future__ import annotations

import html
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


def lab_label(source: Path, title: str) -> str:
    week_match = re.search(r"week-(\d+)", source.as_posix(), re.IGNORECASE)
    session_match = re.search(r"session[-_ ]?(\d+)", source.stem, re.IGNORECASE)
    if not session_match:
        session_match = re.search(r"week[-_ ]?\d+[-_ ]?(\d+)", source.stem, re.IGNORECASE)
    if not session_match:
        title_match = re.search(r"week\s+\d+\s+session\s+(\d+)", title, re.IGNORECASE)
        session_match = title_match

    if week_match and session_match:
        return f"Week {int(week_match.group(1))} Session {int(session_match.group(1))}"
    if week_match:
        return f"Week {int(week_match.group(1))}"
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
            "back_label": "Back to labs",
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
            "back_label": "Back to tutorials",
            "back_fragment": "tutorials",
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


def relative_prefix(output: Path) -> str:
    prefix = Path(*([".."] * len(output.parent.relative_to(DOCS).parts)))
    return "" if str(prefix) == "." else prefix.as_posix()


def prefixed(prefix: str, path: str) -> str:
    return f"{prefix}/{path}" if prefix else path


def page_template(title: str, meta: str, body: str, toc: list[tuple[int, str, str]], page: dict[str, str]) -> str:
    prefix = relative_prefix(Path(page["output"]))
    index_href = prefixed(prefix, f'index.html#{page["back_fragment"]}')
    labs_href = prefixed(prefix, "index.html#labs")
    tutorials_href = prefixed(prefix, "index.html#tutorials")
    favicon_href = prefixed(prefix, "assets/favicon.svg")
    styles_href = prefixed(prefix, "styles.css")
    code_copy_href = prefixed(prefix, "assets/lab-code-copy.js")
    output_path = Path(page["output"])
    is_lab_page = output_path.relative_to(DOCS).parts[0] == "labs"
    code_copy_script = f'\n    <script src="{code_copy_href}" defer></script>' if is_lab_page else ""
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
    <link rel="icon" href="{favicon_href}" type="image/svg+xml">
    <link rel="stylesheet" href="{styles_href}">{code_copy_script}
  </head>
  <body class="lab-page">
    <a class="skip-link" href="#main">Skip to main content</a>
    <header class="site-header" id="top">
      <nav class="nav" aria-label="Main navigation">
        <a class="brand" href="{labs_href}" aria-label="Applied Bioinformatics">
          <span class="brand-mark" aria-hidden="true">AB</span>
          <span>Applied Bioinformatics</span>
        </a>
        <div class="nav-links">
          <a href="{labs_href}">Labs</a>
          <a href="{tutorials_href}">Tutorials</a>
          <a href="https://github.com/zhan4429/TuftsAppliedBioinfo">GitHub</a>
        </div>
      </nav>
    </header>
    <main class="lab-main" id="main">
      <div class="lab-shell">
        <aside class="lab-toc" aria-label="Page table of contents">
          <a class="back-link" href="{index_href}">{html.escape(page["back_label"])}</a>
          <p class="section-kicker">{html.escape(page["label"])}</p>
          <nav>{toc_items}</nav>
        </aside>
        <article class="lab-article">
          <header class="lab-hero">
            <p class="section-kicker">{html.escape(page["label"])}</p>
            <h1>{html.escape(title)}</h1>
            <p>{html.escape(page["summary"])}</p>
            <span>{html.escape(meta)}</span>
          </header>
          {body}
        </article>
      </div>
    </main>
    <footer class="site-footer">
      <p>Applied Bioinformatics - Tufts University Department of Biology</p>
      <a href="{labs_href}">Course labs</a>
    </footer>
  </body>
</html>
"""


def card_markup(cards: list[dict[str, str]]) -> str:
    return "\n".join(
        f"""          <a
            class="lab-card"
            href="{html.escape(card["href"], quote=True)}"
          >
            <span class="lab-session">{html.escape(card["label"])}</span>
            <strong>{html.escape(card["title"])}</strong>
            <span>{html.escape(card["summary"])}</span>
          </a>"""
        for card in cards
    )


def index_template(lab_cards: list[dict[str, str]], tutorial_cards: list[dict[str, str]]) -> str:
    return f"""<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta
      name="description"
      content="Applied Bioinformatics course lab hub for Tufts University Department of Biology."
    >
    <title>Applied Bioinformatics | Tufts Biology</title>
    <link rel="icon" href="assets/favicon.svg" type="image/svg+xml">
    <link rel="stylesheet" href="styles.css">
  </head>
  <body>
    <a class="skip-link" href="#main">Skip to main content</a>

    <header class="site-header" id="top">
      <nav class="nav" aria-label="Main navigation">
        <a class="brand" href="#top" aria-label="Applied Bioinformatics home">
          <span class="brand-mark" aria-hidden="true">AB</span>
          <span>Applied Bioinformatics</span>
        </a>
        <div class="nav-links">
          <a href="#labs">Labs</a>
          <a href="#tutorials">Tutorials</a>
          <a href="https://github.com/zhan4429/TuftsAppliedBioinfo">GitHub</a>
        </div>
      </nav>
    </header>

    <main id="main">
      <section class="masthead" aria-labelledby="page-title">
        <div class="masthead-inner">
          <p class="eyebrow">Tufts University - Department of Biology</p>
          <h1 id="page-title">Applied Bioinformatics</h1>
          <p class="lede">
            A lab hub for hands-on bioinformatics work on Pax, with copyable
            commands, GitHub-friendly Markdown, and short notes students can
            use during class.
          </p>
          <div class="hero-actions" aria-label="Primary actions">
            <a class="button primary" href="#labs">Open labs</a>
            <a class="button secondary" href="#tutorials">Open tutorials</a>
            <a class="button secondary" href="https://github.com/zhan4429/TuftsAppliedBioinfo">GitHub repository</a>
          </div>
        </div>
      </section>

      <section class="quick-strip" aria-label="Course workflow highlights">
        <div class="quick-item">
          <span class="quick-number">01</span>
          <strong>Clone once</strong>
          <span>Keep the lab files local for class.</span>
        </div>
        <div class="quick-item">
          <span class="quick-number">02</span>
          <strong>Pull often</strong>
          <span>Refresh before each lab session.</span>
        </div>
        <div class="quick-item">
          <span class="quick-number">03</span>
          <strong>Open a lab</strong>
          <span>Use the course pages during class.</span>
        </div>
        <div class="quick-item">
          <span class="quick-number">04</span>
          <strong>Run carefully</strong>
          <span>Use Tufts HPC and course storage as directed.</span>
        </div>
      </section>

      <section class="lab-hub" id="labs" aria-labelledby="labs-title">
        <div class="section-copy wide">
          <p class="section-kicker">Available Labs</p>
          <h2 id="labs-title">Course labs</h2>
          <p>
            Open these labs during class. They are organized
            for reading, copying commands, checking expected outputs, and
            revisiting answers after trying each question.
          </p>
        </div>

        <div class="lab-card-grid">
{card_markup(lab_cards)}
        </div>

        <div class="lab-setup">
          <div class="section-copy">
            <p class="section-kicker">Using the Repository</p>
            <h2>Keep your copy current</h2>
            <p>
              Clone once, then pull before each class. Large data and generated
              results should stay on Pax or course storage, not in GitHub.
            </p>
          </div>
          <div class="command-panel" aria-label="Basic repository commands">
            <div class="panel-bar">
              <span></span>
              <span></span>
              <span></span>
            </div>
            <pre><code>git clone &lt;repository-url&gt;
cd &lt;repository-name&gt;
git pull</code></pre>
          </div>
        </div>
      </section>

      <section class="lab-hub" id="tutorials" aria-labelledby="tutorials-title">
        <div class="section-copy wide">
          <p class="section-kicker">Tutorials</p>
          <h2 id="tutorials-title">Setup guides</h2>
          <p>
            Use these short guides when you need a one-time setup step or a
            refresher for course tools.
          </p>
        </div>

        <div class="lab-card-grid tutorial-card-grid">
{card_markup(tutorial_cards)}
        </div>
      </section>
    </main>

    <footer class="site-footer">
      <p>Applied Bioinformatics - Tufts University Department of Biology</p>
      <a href="https://github.com/zhan4429/TuftsAppliedBioinfo">GitHub repository</a>
    </footer>
  </body>
</html>
"""


def build_card(page: dict[str, str | Path], title: str) -> dict[str, str]:
    output = Path(page["output"])
    return {
        "href": output.relative_to(DOCS).as_posix(),
        "label": str(page["label"]),
        "title": display_title(title),
        "summary": str(page["summary"]),
        "section": str(page["section"]),
    }


def main() -> None:
    cards: list[dict[str, str]] = []

    for page in [*discover_lab_pages(), *discover_tutorial_pages()]:
        source = Path(page["source"])
        output = Path(page["output"])
        markdown = source.read_text(encoding="utf-8")
        title, meta, body, toc = render_markdown(markdown)
        page["label"] = str(page.get("label") or lab_label(source, title))
        page["summary"] = SUMMARY_OVERRIDES.get(
            source_key(source),
            extract_summary(markdown, f"Open {display_title(title)} for this course session."),
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(page_template(title, meta, body, toc, page), encoding="utf-8")
        print(f"Built {output.relative_to(ROOT)}")
        cards.append(build_card(page, title))

    lab_cards = [card for card in cards if card["section"] == "labs"]
    tutorial_cards = [card for card in cards if card["section"] == "tutorials"]
    index_output = DOCS / "index.html"
    index_output.write_text(index_template(lab_cards, tutorial_cards), encoding="utf-8")
    print(f"Built {index_output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
