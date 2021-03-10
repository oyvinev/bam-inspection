Simple Python script for parsing and inspecting BAM and BAI-files, closely following the [BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)

## Dependencies
Python >= 3.8

No thirdparty dependencies

## Usage
*python inspect_bam.py BAM_OR_BAI [--fast]*

## About

Intended to debug issues with BAM/BAI-files (with your own necessary modifications).

Can also be used as a tool to understand the [BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)

Note that it is incomplete w.r.t. parsing `seq`, `qual` and auxiliary data.