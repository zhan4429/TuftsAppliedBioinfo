# sed Cheat Sheet

Applied Bioinformatics - Tufts University Department of Biology

Use `sed` to replace text, add or remove prefixes, and select lines from text files.
The examples below print the edited text or save it to a new file. See the
[AWK cheat sheet](awk-cheatsheet.html) for filtering by columns and doing arithmetic.

## The Basic Pattern

```text
sed 's/old/new/' input.txt
```

`s` means **substitute**: replace the first match of `old` with `new` on each line.
The slashes separate the parts of the instruction. The single quotes keep the shell
from interpreting the sed program.

By default, sed prints the result to the terminal. It does not change the input
file. To save the result, use `>` and a different output filename.

## Make a Small Practice File

Run this in your own work directory, using a new filename if `samples.txt` already exists:

```bash
cat > samples.txt <<'EOF'
control_A.fastq.gz
control_B.fastq.gz
treated_A.fastq.gz
EOF
```

`cat` writes the lines between the two `EOF` markers into `samples.txt`.

## Replace Text

```bash
sed 's/control/reference/' samples.txt
```

Output:

```text
reference_A.fastq.gz
reference_B.fastq.gz
treated_A.fastq.gz
```

Lines without a match are printed unchanged. Without the `g` flag, only the first
match on each line is replaced. Compare:

```bash
echo 'read read read' | sed 's/read/READ/'
echo 'read read read' | sed 's/read/READ/g'
```

The first command prints `READ read read`; the second prints `READ READ READ`.
`g` means replace **all matches on each line**. The pipe (`|`) sends the text from
`echo` to sed.

## Match a Prefix or Suffix

| Pattern | Meaning |
| --- | --- |
| `^` | Start of a line. |
| `$` | End of a line. |
| `.` | Any single character. |
| `\.` | A literal dot. |

Add a prefix to every line:

```bash
sed 's/^/raw_/' samples.txt
```

The first line becomes `raw_control_A.fastq.gz`. `^` matches the start of the line,
so nothing in the middle is changed.

Remove a leading prefix:

```bash
sed 's/^control_//' samples.txt
```

The first two lines become `A.fastq.gz` and `B.fastq.gz`. The replacement is empty,
which removes the matched text. `treated_A.fastq.gz` stays unchanged.

Remove an exact filename suffix:

```bash
sed 's/\.fastq\.gz$//' samples.txt
```

The first line becomes `control_A`. Escape the dots with `\.` so they mean actual
dots; `$` restricts the match to the end of a line.

## Select or Exclude Lines

Print only the first five lines:

```bash
sed -n '1,5p' samples.txt
```

`-n` suppresses automatic printing. `1,5` selects lines 1 through 5, and `p` prints
them. This practice file has only three lines, so all three appear.

Print only lines starting with `control`:

```bash
sed -n '/^control/p' samples.txt
```

Exclude lines starting with `#`, such as comments in an annotation file:

```bash
sed '/^#/d' anno.gff3
```

`d` removes matching lines from the output. The original file is unchanged.
Here `anno.gff3` means your annotation file, not the small practice file above.

## Change Chromosome Prefixes in GFF3

GFF3 data rows start with a sequence name, such as `I` or `chrI`. Lines beginning
with `#` are comments or directives. Run these examples in the directory containing
your Lab 1 annotation files.

Add `chr` only to data rows:

```bash
sed '/^#/!s/^/chr/' anno.gff3 > anno_chr.gff3
```

Read the program in two parts:

- `/^#/!` means apply the following edit only to lines that do **not** start with `#`.
- `s/^/chr/` adds `chr` at the start of those lines.

For example, a data row beginning with `IV` now begins with `chrIV`. The comment
lines stay unchanged. Before using a renamed annotation with other tools, also check
that reference sequence names and any `##sequence-region` directives are consistent.

Remove the leading `chr` again:

```bash
sed 's/^chr//' anno_chr.gff3 > anno_no_prefix.gff3
```

`^chr` matches only a leading `chr`, so it does not remove `chr` from text in the
attributes column. Comments beginning with `#` stay unchanged.

## Save and Check the Result

```bash
sed 's/^control_/reference_/' samples.txt > samples_renamed.txt
head -n 5 samples_renamed.txt
wc -l samples_renamed.txt
```

`head` previews the saved text, and `wc -l` reports its line count. A substitution
changes text without adding or removing lines; a filtering command such as `d`
can change the line count.

- Always use a different output filename: `>` can empty the input before sed reads it if both names are the same.
- These examples use a new output file, so you can inspect the edit before replacing anything.
- A sed search pattern is a regular expression. Use anchors for prefixes and suffixes and escape literal dots.
- Sed edits text patterns; for selecting a GFF3 feature type specifically in column 3, use [awk](awk-cheatsheet.html).

## Further Reading

- [GNU sed manual](https://www.gnu.org/software/sed/manual/sed.html)
- [GNU sed manual on selecting lines](https://www.gnu.org/software/sed/manual/html_node/Addresses-overview.html)
- [GNU sed common commands](https://www.gnu.org/software/sed/manual/html_node/Common-Commands.html)
- [AWK cheat sheet](awk-cheatsheet.html)
- [Command Line, Text Processing, and Git lab](../labs/week-01/01-command-line-git.html)
