# CHESS

CHESS is a tool for the assessment of similarity between Hi-C matrices using structural similarity.

<!-- TOC depthFrom:1 depthTo:8 withLinks:1 updateOnSave:1 orderedList:0 -->

- [CHESS](#chess)
    - [Requirements](#requirements)
    - [Installation](#installation)
    - [Quick start](#quick-start)
    - [Usage](#usage)

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

## Quick start

After installation, you can use the example data to familiarize yourself with the applications of CHESS implemented in the current version.

Use `chess -h` for quick help and orientation.

### Finding changing regions between biological conditions in the same species

This will run a comparison of 250 kb submatrices of chromosome X in wildtype _Drosophila melanogaster_ to a _zld_ knockdown.
The output directory is set to `examples/Dmel_genome_scan` and should, after the successful run, contain a comparison_results.tsv file with the raw similarity scores that can be used to rank the regions according to their similarity.
(Hi-C data from Hug et al. 2017)
```bash
chess sim \
examples/Dmel_genome_scan/zld_X.sparse.gz \
examples/Dmel_genome_scan/zld_X.regions.gz \
examples/Dmel_genome_scan/wt_X.sparse.gz \
examples/Dmel_genome_scan/wt_X.regions.gz \
examples/Dmel_genome_scan/Dmel_zld_kd_wt_nc14_chrm_X_250kwindow_25kbstep.pairs.gz \
examples/Dmel_genome_scan/comparison_results.tsv \
--set-wd examples/Dmel_genome_scan --genome-scan
```
> NOTE: Running this example might take some time. To speed it up, you can use the `-p <int>` flag to specify the number of cores to use (default: 1) and split the workload.


### Comparing regions accross species

This will run a comparison of 4 regions of varying sizes located on chromosome 19 in mouse and on chromosome 10 in human.
The output directory is set to `examples/Mmus_Hsap_syntenic` and should, after the successful run, contain a comparison_results.tsv file with the raw similarity scores along with p-values and z-scores which can be used to rank the regions according to their similarity.
(Hi-C data from Rao et al. 2014)
```bash
chess sim \
examples/Mmus_Hsap_syntenic/rao2014_hg19_chr10_25kb.sparse.gz \
examples/Mmus_Hsap_syntenic/rao2014_hg19_chr10_25kb.regions.gz \
examples/Mmus_Hsap_syntenic/rao2014_mm10_chr19_25kb.sparse.gz \
examples/Mmus_Hsap_syntenic/rao2014_mm10_chr19_25kb.regions.gz \
examples/Mmus_Hsap_syntenic/hg19_mm10_syntenic_examples.pairs.gz \
examples/Mmus_Hsap_syntenic/comparison_results.tsv \
--set-wd examples/Mmus_Hsap_syntenic
```

## Usage

CHESS has two basic commands: `sim`, which wraps the matrix comparison features and `oe`, which you can use to transform normalized Hi-C matrices into observed / expected matrices.

### sim

The main purpose of CHESS is the assessment of similarity between two Hi-C matrices, one of which is called the 'reference' and the other the 'query'. For this purpose, it has two basic modes:
  - inter-species comparisons (default)
  - intra-species comparisons using the `--genome-scan` flag.

The default mode will use the the queries whole genome Hi-C to check how the similarity between the reference and query compares to the similarity between the reference and any other region of the same size as the query. Using that information, the comparison between the reference and the query will be assigned a z-score and a p-value which can be used to rank multiple reference - query pairs according to their similarity.

The intra-species mode will only allow comparisons between matrices of equal size and is intended for finding changes in Hi-C data between biological conditions. The ouput will only contain a similarity score. The regions with the lowest scores are the least similar between both conditions, as assessed by CHESS.

CHESS takes the same mandatory, positional arguments for both modes:

* A normalized Hi-C matrix file in sparse matrix format for the reference sample. This should be tab delimited, where each line has three columns:
      \<row index\> \<column index\> \<value\>)

> NOTE: Using CHESS with large matrices at high resolution can require a lot of memory,
         especially when using many threads. Typically, every thread will load one full
         chromosome Hi-C map into memory (the number of threads can be controlled with the `-p` flag).

> NOTE: By default, the input matrices will be converted to observed / expected matrices. In case you already have observed / expected matrices, you can use the `--converted-input` flag to skip the observed / expected conversion.

* A BED file with region information for the reference Hi-C matrix.
  This should be a tab-delimited file where each row contains chromosome name, start, and end coordinates (exclusive) of the region. This file must 
  not contain any headers. If a fourth column is present, it is assumed to be a unique identifier
  (index/name), which is then used to refer to that region in sparse matrix format
  (see above).

* A normalized Hi-C matrix file in sparse matrix format for the query sample.

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

  You can read the .json output in python using the json package:
```python
  import json
  with open('path_to_json_file', 'r') as f:
    my_dict = json.load(f)
```

When the -genome_scan flag is set, OUT will be the only output file, as no background calculation is done in that case.

Optional arguments give you more control:

* `--reference_ID <string>` Species / Sample identifier for the reference Hi-C data. Will be used as a prefix for the intermediate file names corresponding to the reference Hi-C. Default: REF

* `--query_ID <string>` Species / Sample identifier for the query Hi-C data. Will be used as a prefix for the intermediate file names corresponding to the query Hi-C. Default: QRY

* `-d`, `--set-wd <string>` lets you set the working directory of CHESS. All intermediate files will be written to this directory (default: './'.

* `--no-clean` lets you keep all intermediate files (split chromosome files, helpers and full input matrices converted to observed / expected matrices).

* `--fast-input` lets you use intermediate files generated in previous runs located in the working directory. Reduces the runtime, as the input sparse matrix won't have to be split.

* `--genome-scan` will run without any background comparisons and report only the raw ssim scores for the pairs defined in the pairs file. All regions in the pairs file are expected to be of similar size to ensure comparability of the scores. Use this if you want to find differences between two Hi-C datasets mapped to the same genome.

* `--limit-background` reduces the computation of the background score distribution to the chromosome the query region (defined in the pairs file) is located on. Reduces runtime. However, it is in general advisable to use the whole genome as a background. Otherwise, some bias might be introduced due to differing chromosome sizes.

* `--converted-input` will skip the observed / expected conversion of the input matrices. Use if you already have observed / expected matrices or want to compare matrices in different format.

* `--no-raw` will only write only the standard OUT file. Set if you have no use for the OUT.FULL_RAW.json and OUT.FULL_ROUNDED_QUERIES.json files.

* `-p <int>` lets you choose the number of cores that CHESS will use (default: 1).

* `--keep-unmappable-bins` disables the deletion of deletion of unmappable bins from matrices before comparison. By default, bins that are marked as unmappable (have not contacts) in either of two matrices in a comparison are deleted from both matrices. Disabling this might give high similarity scores for matrices that happen to have unmappable bins at the same positions.

* `--mappability-cutoff <float>` maximum allowed content of unmappable bins in a matrix. In case a matrix has more unmappable bins, it will not be compared. (default: 0.1)

* `-r <float>, --relative-windowsize <float>` Fraction of the matrix size that will be used for the window size parameter of the structural similarity function. (default: 1)

* `-a <float>, --absolute-windowsize <float>` Absolute value for the window size parameter of the structural similarity function. Overwrites `-r`.

### oe

`chess oe` lets you convert your input matrix to observed / expected format. Calling `chess sim` this is done autmatically, but if you want to convert your matrix for other reasons, you can use this.

`oe` takes three positional arguments (in that order):

* A Hi-C matrix file in sparse matrix format. This should be tab delimited, where each line has three columns:
      \<row index\> \<column index\> \<value\>)

* A BED file with region information for the Hi-C matrix.
  This should be a tab-delimited file where each row contains chromosome name, start, and end coordinates (exclusive) of the region. This file must 
  not contain any headers. If a fourth column is present, it is assumed to be a unique identifier
  (index/name), which is then used to refer to that region in sparse matrix format
  (see above).

* A path to the output matrix file. Will in the same format as the input matrix.