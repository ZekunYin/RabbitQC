![RabbitQC](qc.png)

# RabbitQC
A tool designed to provide high-speed scalable quality control for sequencing data which can take full advantage of modern hardware.
It includes a variety of function modules and supports different sequencing technologies (Illumina, Oxford Nanopore, PacBio). RabbitQC achieves speedups between one and two orders-of-magnitude compared to other state-of-the-art tools.

# Simple usage
## For short read
* For single end data (not compressed)
```
rabbit_qc -w nthreads -i in.fq -o out.fq
```
* For paired end data (gzip compressed)
```
rabbit_qc -w nthreads -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```
## For long read
```
rabbit_qc -w nthreads -D -i in.fq
```

## For large gz files
A more efficient strategy to process large gzip compressed FASTQ files is to decompress files using pugz and then process them using RabbitQC. Pugz has been integrated into RabbitQC project.

```
cd RabbitQC/pugz && make asserts=0
./gunzip -t nthreads in.fq.gz
```

# Options
For more help information, please refer to `rabbit_qc -h`.

If `-w` opition is not specified, RabbitQC will set working thread number to total CPU cores - 2.
By default, the HTML report is saved to `RabbitQC.html` (can be specified with `-h` option), and the JSON report is saved to `RabbitQC.json` (can be specified with `-j` option).

RabbitQC suports all fastp options for short read quality control and all NanoQC optiions for long read quality control. For details please refer to [fastp](https://github.com/OpenGene/fastp) and [NanoQC](https://github.com/wdecoster/nanoQC).

# Examples of report
`RabbitQC` creates reports in both HTML and JSON format.

# Build

**For Linux and OSX:**

```bash
cd RabbitQC && make
```
**For Windows:**

We provide a prebuild binary for x64 windows (tested on 64bit Windows 10) [here](https://github.com/ZekunYin/RabbitQC/releases). Or you can build RabbitQC using MYSY2.

```bash
cd RabbitQC && make
```



# Citation

Zekun Yin, Hao Zhang, Meiyang Liu, Wen Zhang, Honglei Song, Haidong Lan, Yanjie Wei, Beifang Niu, Bertil Schmidt, Weiguo Liu, RabbitQC: High-speed scalable quality control for sequencing data, Bioinformatics, , btaa719, https://doi.org/10.1093/bioinformatics/btaa719
<!--

## If you use RabbitQC for short read quality control please cite:

Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

## If you use RabbitQC for long read quality control please cite:

De Coster W, D’Hert S, Schultz D T, et al. NanoPack: visualizing and processing long-read sequencing data[J]. Bioinformatics, 2018, 34(15): 2666-2669.
-->
