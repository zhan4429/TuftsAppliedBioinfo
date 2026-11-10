# AWK Cheat Sheet

Applied Bioinformatics - Tufts University Department of Biology

Use `awk` to select rows, extract columns, count features, and calculate totals in
text files. This reference starts with a small practice file, then applies the same
ideas to genome annotations. See the [sed cheat sheet](sed-cheatsheet.html) for
replacing text and editing line prefixes.

## The Basic Pattern

```text
awk -F'\t' 'condition { action }' file.tsv
```

Read this as: split each line into tab-separated fields, check the condition, and
perform the action on matching lines. Replace `condition`, `action`, and `file.tsv`
with the rule and input you need.

| Piece | Meaning |
| --- | --- |
| `-F'\t'` | Fields are separated by tabs. Without `-F`, awk normally splits on whitespace. |
| `$1`, `$2`, `$3` | The first, second, and third fields on the current line. |
| `$0` | The whole current line. |
| `NR` | The number of input lines read so far. |
| `NF` | The number of fields on the current line. |
| `'...'` | Single quotes protect the awk program from expansion by the shell. |

With a condition but no action, awk prints the matching line. With an action but
no condition, it performs that action on every line.

## Make a Small Practice File

Run this in your own work directory. The values are separated by **tabs** and there
is no header row. The columns are sample name, group, and read count.

```bash
cat > example.tsv <<'EOF'
sample_A	control	1200
sample_B	treated	1800
sample_C	control	900
EOF
```

`cat` writes the lines between the two `EOF` markers into `example.tsv`.
`>` replaces that output file if it already exists, so use a new filename for practice.

## Select Columns

Print just the sample names:

```bash
awk -F'\t' '{ print $1 }' example.tsv
```

Print sample names and read counts, separated by a tab:

```bash
awk -F'\t' '{ print $1 "\t" $3 }' example.tsv
```

The `"\t"` inserts a tab between the two fields. By comparison, `print $1, $3`
uses a space between them by default.

## Select Rows

Keep the control samples:

```bash
awk -F'\t' '$2 == "control"' example.tsv
```

This prints the complete lines for `sample_A` and `sample_C`. `==` compares values;
the text value `"control"` goes in double quotes inside the single-quoted program.

Keep samples with at least 1,500 reads:

```bash
awk -F'\t' '$3 >= 1500' example.tsv
```

Only `sample_B` passes. Combine conditions with `&&`, meaning **and**:

```bash
awk -F'\t' '$2 == "control" && $3 >= 1000' example.tsv
```

Only `sample_A` passes both conditions.

## Count Rows and Add Values

Count the selected control rows:

```bash
awk -F'\t' '$2 == "control"' example.tsv | wc -l
```

A pipe (`|`) sends awk's output to the next command. `wc -l` counts those lines,
so the answer is **2**.

Add the read counts from all three samples:

```bash
awk -F'\t' '{ total += $3 } END { print total }' example.tsv
```

`total += $3` adds the current read count to a running total. Awk treats the initial
numeric value of `total` as zero. `END` runs after the last input line, so the total
is printed once: **3900**. To sum only control samples, add `$2 == "control"`
before the first pair of braces.

## Count Categories with a Pipeline

```bash
awk -F'\t' '{ print $2 }' example.tsv |
  sort | uniq -c | sort -rn
```

| Step | Purpose |
| --- | --- |
| `print $2` | Extract the group names. |
| `sort` | Put identical names next to each other. |
| `uniq -c` | Count each group of adjacent identical lines. |
| `sort -rn` | Sort the counts numerically, largest first. |

The counts are **2 control** and **1 treated**. Keep the first `sort`: `uniq`
groups only adjacent matches. A line ending in `|` continues onto the next line;
you can also type the pipeline on one line.

## Work with GFF3 Annotations

In a GFF3 data row, column 1 is the sequence name, column 3 is the feature type,
and columns 4 and 5 are the start and end. Use tabs as the field separator because
the attributes column can contain spaces. Lines starting with `#` are comments.

Run the following examples in a directory containing `anno.gff3`, such as your
Lab 1 folder.

### Inspect Data Rows

```bash
awk -F'\t' '!/^#/ { print $1 "\t" $3 }' anno.gff3 | head -n 5
```

`/^#/` matches lines starting with `#`; `!` means **not**. This skips comments,
extracts sequence names and feature types, and displays the first five results.

### Count a Feature Type

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 | wc -l
```

Match the feature-type field exactly. Searching for the word `gene` anywhere in a
line can also match other feature types whose attributes mention a gene.

### Rank Sequences by Gene Count

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 |
  cut -f1 | sort | uniq -c | sort -rn | head -n 5
```

`cut -f1` keeps column 1. The next three commands count and rank the sequence
names; `head -n 5` keeps the five largest counts. Build the pipeline one step at a
time if you want to see what each command contributes.

### Sum Annotated Exon Lengths

```bash
awk -F'\t' '$3 == "exon" { total += $5 - $4 + 1 } END { print total }' anno.gff3
```

Use the lab's GFF3 length formula, `$5 - $4 + 1`, for each exon. This sums feature
lengths; it does not merge overlapping exons into distinct genomic positions.

## Save Output and Avoid Common Mistakes

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 > genes.gff3
head -n 5 genes.gff3
wc -l genes.gff3
```

- `>` writes to a file and replaces any existing contents. Keep the output filename different from the input filename.
- Keep the awk program in straight single quotes so the shell does not expand `$1`, `$3`, and other field references.
- Field numbers depend on the input file. Check the columns before reusing a command.
- If a TSV has a header, `NR > 1` skips its first line. Apply this only when a header is actually present.

## Further Reading

- [GNU awk manual on fields](https://www.gnu.org/software/gawk/manual/html_node/Fields.html)
- [GNU awk manual on output separators](https://www.gnu.org/software/gawk/manual/html_node/Output-Separators.html)
- [GNU awk user's guide](https://www.gnu.org/s/gawk/manual/gawk.html)
- [sed cheat sheet](sed-cheatsheet.html)
- [Command Line, Text Processing, and Git lab](../labs/week-01/01-command-line-git.html)
