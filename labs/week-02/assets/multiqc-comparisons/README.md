# MultiQC class-report comparisons

Source: instructor-provided `Week-2-QC-across-four-platforms_multiqc_report.html`.

The two PNG figures are static Plotly exports from that report. Adapter content combines four short-read inputs. The quality plot combines the two NanoStat summaries and uses the report's own Percentages control so each library sums to 100%.

Exports retain all plotted data and the report's native colors. Display-only changes: the adapter legend is visible, sample and adapter labels are shortened, the long-read rows are labeled PacBio HiFi and Nanopore instead of their output-directory names, the percentage axis has an explicit title, and margins, background, and export dimensions are set for the lab page. No data filtering, smoothing, or resampling was applied.

Native source plot IDs: `fastqc_adapter_content_plot` and `nanostat_quality_dist`. The selected original dataset definitions, including the count and percentage values, are retained in `figure-data.json`. No input filesystem paths or report metadata are included in these assets.

The quality bars describe per-read mean qualities, not percentages of individual Q20 or Q30 bases. Differences among these separately prepared libraries do not isolate a platform effect.

Source report SHA-256: `7a66c3305713920d71198811dc94894412496a315d2ca883622f08d65a1fe63b`.
