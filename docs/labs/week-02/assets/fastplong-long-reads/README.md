# Long-read filtering examples

Source: instructor-provided fastplong 0.7.1 HTML reports, generated 2026-11-12:

- `SRX13548509_SRR17374240.fastplong.html` (Nanopore)
- `SRX14356846_SRR18210286.fastplong.html` (PacBio HiFi)

Both use the lab settings: `--mean_qual 10 --length_required 1000 --unqualified_percent_limit 100`.

The PNGs export each report's native **Before filtering: Read mean quality distribution** data using the report's Plotly 1.2.0 library. The original grouped bars, percentages, and trace colors are retained. Display changes: descriptive platform titles, shared x-axis (0–100) and y-axis (0–14%), explicit arithmetic-mean axis wording, a dashed Q10 cutoff, and larger labels/export dimensions. No smoothing, resampling, or filtering of the plotted data was applied. Each mean-quality histogram bin represents a rounded whole-number score; the exact filtering result cannot be recovered by simply summing bins below the cutoff.

Blue bars show percentages of reads; purple bars show percentages of bases contributed by reads in the same mean-quality bin. The purple bars do not show the percentage of individual bases at that Phred score. These arithmetic mean scores must not be equated with SeqKit's error-probability-based `AvgQual`, nor with measured alignment accuracy.

`figure-data.json` preserves both before- and after-filtering native chart arrays, source plot IDs, and SHA-256 hashes for the two source HTML files. Only the before-filtering charts are displayed: the lab's filters remove six Nanopore reads and no HiFi reads, so the after-filtering distributions are visually unchanged.

The additional default-settings Nanopore report was reviewed but is not included in the lab figures because the default-run exercise was removed from the class workflow.

Metric definitions: [fastplong filter implementation](https://github.com/OpenGene/fastplong/blob/main/src/filter.cpp), [fastplong histogram implementation](https://github.com/OpenGene/fastplong/blob/main/src/stats.cpp), and [SeqKit statistics documentation](https://bioinf.shenwei.me/seqkit/usage/#stats).
