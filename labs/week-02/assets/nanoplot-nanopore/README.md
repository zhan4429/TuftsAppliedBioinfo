# Nanopore NanoPlot worked example

Source: instructor-provided `NanoPlot-report.html`. The report does not name an input FASTQ. Its 701,791 reads, 9,633,074,833 bases, length N50 of 20,948 bp, and maximum length of 276,396 bp match the Nanopore dataset used in the lab; that match is the basis for the dataset label.

`NanoPlot-report.html` is an unchanged copy of the supplied report. It loads Plotly 3.7.0 from the Plotly CDN and needs an internet connection for interactive plots.

The three PNGs are static exports of the report's original Plotly figures, rendered with its own Plotly library at 1200 by 750 CSS pixels, scale 1.5 (1800 by 1125 PNG pixels). Original data, axes, binning, titles, colors, N50 markers and display ranges were retained. No smoothing, filtering, rebinning or image edits were applied.

Selected original plot IDs:

- Read-count histogram: `69cda055-e944-4e5e-baf1-7551f021ac1f`.
- Base-weighted histogram: `88c403ef-3899-4b70-99fd-454ddee140f8`.
- Length-versus-quality dot plot: `88020048-734d-43d6-9587-f422ed0183f0` (10,000 sampled reads, as embedded in the supplied report; original displayed axis range retained).

The figures summarize reported basecall qualities; they do not measure accuracy against a reference. Per-read quality cutoffs differ from per-base Q20/Q30 percentages.

Source SHA-256: `e953bc117ef96da8db17abfd2ecfdcaea6bf3ccf0013a39aa9a6b03f15ae0a74`.
