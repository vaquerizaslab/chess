# CHESS

CHESS is a tool for the assessment of similarity between Hi-C matrices using structural similarity.

<!-- TOC depthFrom:1 depthTo:8 withLinks:1 updateOnSave:1 orderedList:0 -->

- [CHESS](#chess)
    - [Requirements]
    - [Installation](#installation)
    - [Usage](#usage)
        - [`plot`](#plot)
            - [Launching from the command line](#launching-from-the-command-line)
            - [Interactive plotting window](#interactive-plotting-window)
        - [`tads`](#tads)
        - [`subset`](#subset)

<!-- /TOC -->

## Requirements

CHESS was written and tested using Python 3.6.1 on Linux and macOS.

## Installation

You can install CHESS from the command line using PyPI:

```
pip install chess-hic
```

or download the code from our [GitHub repo](https://github.com/vaquerizaslab/chess/) and install manually

```bash
python setup.py install
```

## Usage

There two basic modes if running CHESS: 
  - for inter-species comparisons (default)
  - for intra-species comparisons using the -genome_scan flag.

CHESS takes the same mandatory, positional arguments for both modes:

* A Hi-C matrix file in sparse matrix format for the reference sample. This should be tab delimited, where each line has three columns:
      \<row index\> \<column index\> \<value\>)

> NOTE: Using CHESS with large matrices at high resolution can require a lot of memory,
         especially when using many threads. Typically, every thread will load one full
         chromosome Hi-C map into memory (the number of threads can be controlled with the -p flag).

* A BED file with region information for the reference Hi-C matrix.
  This should be a tab-delimited file where each row contains chromosome name, start, and end coordinates (exclusive) of the region. This file must 
  not contain any headers. If a fourth column is present, it is assumed to be a unique identifier
  (index/name), which is then used to refer to that region in sparse matrix format
  (see above).

* A Hi-C matrix file in sparse matrix format for the query sample.

* A BED file with region information for the query Hi-C matrix.

* A BEDPE file (pairs) file that specifies which regions in the reference matrix should
  be compared to which regions in the query matrix. This should be a tab-delimited file
  with columns:
  \<ref chromosome>\> \<ref start>\> \<ref end>\> \<qry chromosome>\> \<qry start>\> \<qry end>\> \<comparison name/id>\> \<anything, not used>\> \<ref strand>\> \<qry strand>\>
  The 8th column is required in order to match the BEDPE standart format (see http://bedtools.readthedocs.io/en/latest/content/general-usage.html).
  This file must not contain any headers. The end coordinates are exclusive.

* The path to the output file.

When called with only these arguments, CHESS will compare the specified regions between the reference and query matrices, using the query matrix as a background for the computation of p- and z-values. The output will consist of three files:

* OUT, a tab-delimited file with columns: 
  \<comparison name/ID>\> \<p-value>\> \<structural similarity score (ssim)>\> \<z-score>\>
  The comparison name/ID is the same that is specified in the BEDPE input file.
  The p-values and z-scores can be used to rank the regions according to their similarity.
  The raw structural similarity score is sensitive to sizes and size differences between the compared regions and should be handled with care.

* OUT.FULL_RAW.json, a Python dictionary of format:
```python
  {
      name/ID: {
        region_1: score,
        ...,
        region_n: score
      }  
  }
```
  This dictionary has an entry for every pair in the input BEDPE that was tested, and contains every region from the background that the reference was compared to and the corresponding structural similarity score.

* OUT.FULL_ROUNDED_QUERIES.json, a Python dictionary of format:
```python
  {
      name/ID: region
  }
```
  This dictionary contains the rounded positions of the input queries that were used in the actual comparisons. The specified start, end positions usually have to be rounded to the nearest bin.
  These positions can be used to discriminate the 'true' query region from the rest in the FULL_RAW.json file.

When the -genome_scan flag is set, OUT will be the only output file, as no background calculation is done in that case.

Optional arguments give you more control:

* `--reference_ID <string>` Species / Sample identifier for the reference Hi-C data. Will be used as a prefix for the intermediate file names corresponding to the reference Hi-C. Default: REF

* `--query_ID <string>` Species / Sample identifier for the query Hi-C data. Will be used as a prefix for the intermediate file names corresponding to the query Hi-C. Default: QRY

* `-d`, `--set-wd <string>` lets you set the working directory of CHESS. All intermediate files will be written to this directory.

* `--no-clean` lets you keep all intermediate files (split chromosome files and helpers).

* `--fast-input` lets you use intermediate files generated in previous runs located in the working directory. Reduces the runtime, as the input sparse matrix won't have to be split.

* `--genome-scan` will run without any background comparisons and report only the raw ssim scores for the pairs defined in the pairs file. All regions in the pairs file are expected to be of similar size to ensure comparability of the scores. Use this if you want to find differences between two Hi-C datasets mapped to the same genome.

* `--limit-background` reduces the computation of the background score distribution to the chromosome the query region (defined in the pairs file) is located on. Reduces runtime.



#### Interactive plotting window

After starting TADtool with the `plot` option as instructed above, the interactive plotting window will open.