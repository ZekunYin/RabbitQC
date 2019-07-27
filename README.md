# RabbitQC
A tool designed to provide high-speed scalable quality control for sequencing data which can take full advantage of modern hardware.
It includes a variety of function modules and supports different sequencing technologies (Illumina, Oxford Nanopore, PacBio). RabbitQC achieves speedups between one and two orders-of-magnitude compared to other state-of-the-art tools.

# simple usage
## for short read
* for single end data (not compressed)
```
RabbitQC -w nthreads -i in.fq -o out.fq
```
* for paired end data (gzip compressed)
```
RabbitQC -w nthreads -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```
## for long read
```
RabbitQC -w nthreads -D in.fq
```

By default, the HTML report is saved to `RabbitQC.html` (can be specified with `-h` option), and the JSON report is saved to `RabbitQC.json` (can be specified with `-j` option).

# Options
RabbitQC suports all fastp options for short read quality control and all NanoQC optiions for long read quality control. For details please to [fastp](https://github.com/OpenGene/fastp) and [NanoQC](https://github.com/wdecoster/nanoQC).

# examples of report
`RabbitQC` creates reports in both HTML and JSON format.

# build
cd rabbit\_qc && make

# citation
RabbitQC paper is under review now.

## If you use RabbitQC for short read quality control please cite:

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

## If you use RabbitQC for long read quality control please cite:

De Coster W, D’Hert S, Schultz D T, et al. NanoPack: visualizing and processing long-read sequencing data[J]. Bioinformatics, 2018, 34(15): 2666-2669.

