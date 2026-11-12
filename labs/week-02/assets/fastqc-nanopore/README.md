# Nanopore FastQC teaching example

Source: instructor-provided `SRX13548509_SRR17374240_fastqc.html` for `SRX13548509_SRR17374240.fastq.gz`.

The two PNGs are the original embedded per-base quality and read-length distribution plots, decoded without resizing, cropping, or other image changes. `nanopore-fastqc.html` is an unchanged copy of the supplied report.

Basic Statistics reports 701,791 reads, lengths 2,000–276,396 bp, and 38% GC. The selected modules are Per base sequence quality (FAIL) and Sequence Length Distribution (WARNING). Variable read length is expected for Nanopore data. These flags do not by themselves establish that the reads are unusable, nor do they establish that all quality concerns are artifacts.
