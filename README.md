# CHESS

CHESS is a tool for the comparison and automatic feature extraction for chromatin contact data.

<!-- TOC depthFrom:1 depthTo:8 withLinks:1 updateOnSave:1 orderedList:0 -->

- [CHESS](#chess)
    - [Requirements](#requirements)
    - [Installation](#installation)
    - [Quick start](#quick-start)
    - [Usage](#usage)
      - [sim](#sim)
      - [oe](#oe)
      - [pairs](#pairs)
      - [filter](#filter)
      - [extract](#extract)
      - [crosscorrelate](#crosscorrelate)
    - [References](#references)

<!-- /TOC -->

## Requirements

CHESS was written and tested using Python 3.6 - 3.8.2 on Scientific Linux (6.10) and OS X El Capitan (10.11.6) and Ubuntu Bionic Beaver (18.04) and Ubuntu Focal Fossa (20.04).

It requires the packages `cython` (0.29.16), `scipy` (1.0.0), `numpy` (1.14.0), `scikit-image` (0.13.1), `pandas` (0.22.0),
`pathos` (0.2.1), `kneed` (0.6.0), `tqdm` (4.43.0), `intervaltree` (3.0.2), `pybedtools` (0.8.1), `future` (0.16.0) and `fanc` (0.8.28).

## Installation

**When installing with a python version < 3.8, Please make sure to have `cython` (0.29.16) installed before installling CHESS via `pip` or `python install setup.py`, or use a `pip` version >= 20.0.2.**

### Via `pip` from the Python Package Index (PyPi)

CHESS can be installed within just a few minutes from the command line using PyPi:

```
pip install chess-hic
```

### From source

You can also download the code from our [GitHub repo](https://github.com/vaquerizaslab/chess)
and install CHESS manually. Make sure to have `cython` installed. The other dependencies should be downloaded autmatically; if you encounter problems, try to install the other [requirements](#requirements) manually.

```bash
git clone https://github.com/vaquerizaslab/chess  # or download manually
pip install chess
```

or

```bash
git clone https://github.com/vaquerizaslab/chess  # or download manually
cd chess
python setup.py install
```


## Quick start

After installation, you can use the example data to familiarize yourself with the applications
of CHESS implemented in the current version. If you installed via `pip`, you can download the
examples from our [GitHub repo](https://github.com/vaquerizaslab/chess).

Use `chess -h` for quick help and orientation.

### Finding changing regions between biological conditions in the same species

The following will run a comparison of 250 kb submatrices of chromosome X in wildtype 
_Drosophila melanogaster_ to a _zld_ knockdown (Hi-C data from Hug et al. 2017):

```bash
chess sim \
examples/Dmel_genome_scan/zld_X.sparse.gz \
examples/Dmel_genome_scan/zld_X.regions.gz \
examples/Dmel_genome_scan/wt_X.sparse.gz \
examples/Dmel_genome_scan/wt_X.regions.gz \
examples/Dmel_genome_scan/Dmel_zld_kd_wt_nc14_chrm_X_250kwindow_25kbstep.pairs.gz \
examples/Dmel_genome_scan/comparison_results.tsv
```

> NOTE: This example run should finish within 2 minutes on a single core. 
> However, to speed it up, you can use the `-p <int>` flag to specify the 
> number of cores to use (default: 1) and split the workload.

The output file `examples/Dmel_genome_scan/comparison_results.tsv` 
contains four columns:

- ID: the ID of the original pair, as in the supplied region pairs file
- SN: the signal-to-noise ratio of the matrices in the comparison, 
  which can be used for filtering. Larger SN means more signal in comparison
  to noise.
- ssim: The raw similarity score, useful to rank matrices by similarity.
- z_ssim: the similarity z-score, in the context of all calculated similarity
  scores. Useful to assess if the similarity score is exceptionally high or low
  compared to other scores in the same run. For z-scores in the context of a 
  proper background model, see below.

When running `chess sim` to compare samples of different conditions, the goal often is
to identify regions that with particularly striking changes. As the 3D structure might 
change slightly across most of the genome, the most interesting regions might be the ones
with the largest deviation from the average change: this is reflected in large negative
values in the *z-score* column. 

Z-scores cannot be compared between runs of CHESS, as they are a type of 
within sample normalisation. To compare between runs of CHESS, the ssim values
should be used. For instance, a change of experimental conditions might not lead to marked
changes in particular regions of the genome, but instead induce an uniform genome wide effect
(e.g. global decompaction of TADs). Then, one might be interested in comparing the average
ssim from `chess sim` runs on replicates of one
condition to the average ssim of a run between conditions.

### Comparing regions across species

The following will run a comparison of 4 regions of varying sizes located on chromosome 
19 in mouse and on chromosome 10 in human (Hi-C data from Rao et al. 2014):

```bash
chess sim \
examples/Mmus_Hsap_syntenic/rao2014_hg19_chr10_25kb.sparse.gz \
examples/Mmus_Hsap_syntenic/rao2014_hg19_chr10_25kb.regions.gz \
examples/Mmus_Hsap_syntenic/rao2014_mm10_chr19_25kb.sparse.gz \
examples/Mmus_Hsap_syntenic/rao2014_mm10_chr19_25kb.regions.gz \
examples/Mmus_Hsap_syntenic/hg19_mm10_syntenic_examples.pairs.gz \
examples/Mmus_Hsap_syntenic/comparison_results.tsv --background-query
```

Note the `--background-query`, which enables the query-genome based background 
model.
The output file `examples/Mmus_Hsap_syntenic/comparison_results.tsv` has six columns:

* ID: the ID of the original pair, as in the supplied region pairs file
* SN: the signal-to-noise ratio of the matrices in the comparison, 
  which can be used for filtering
* ssim: The raw similarity score, useful to rank matrices by similarity
* z_ssim: the similarity z-score, in the context of all calculated similarity
  scores. Useful to assess if the similarity score is exceptionally high or low
  compared to other scores in the same run.
* z_bg: the similarity z-score based on the background matrix comparisons.
* p_bg: The similarity score p-value. This assesses how likely it is to obtain 
  a similarity of ssim or higher given the comparisons in the background model
  - in this case, all other regions in the query genome of the same window size

In contrast to the intra-species comparison described above, where we tried to find 
strong differences between samples, we are here concerned with assessing whether
two regions are significantly similar. To answer that questions, we first need to
determine the distribution of similarities between our reference and random regions
(here called 'background'), and then evaluate in this context the similarity between
the reference and query pair. CHESS reports on the significance and effect size of 
the similarity of a reference and query pair with the `p_bg` p-value and `z_bg` z-score
described above.
There are many conceivable ways of assembling the background set of matrices. The default way,
chosen here with the `--background-query` flag, is to segment the genome into windows of
query size, and compare the reference to all of those.

## Usage

CHESS has these basic commands: 
* `sim` contains the matrix comparison features 
* `oe` can be used to transform normalized Hi-C matrices into 
  observed / expected matrices 
* `pairs` helps you generate the pairs input file for comparing 
  Hi-C data mapped to the same genome between biological conditions
* `background` can be used to generate simple background model BED files
* `filter` can be used to filter results obtained by `sim`, for example 
  by signal-to-noise (SN) ratio
* `extract` allows you to the specific features from significant different 
  regions identified by `chess sim`
* `crosscorrelate` allows to classify the extracted features by `extract` 
  in order to obtain the main structural clusters.


### sim

The main purpose of CHESS is the assessment of similarity between two Hi-C matrices, 
one of which is called the 'reference' and the other the 'query'. The following arguments 
are mandatory to run `sim`:

* A normalized Hi-C matrix file in sparse matrix format for the reference sample. 
  This should be tab-delimited, where each line has three columns:
  \<row index\> \<column index\> \<value\>)

> NOTE: Using CHESS with large matrices at high resolution can require a lot of memory,
> especially when using many threads. Typically, every thread will load one full
> chromosome Hi-C map into memory (the number of threads can be controlled with the `-p` flag).

> NOTE: By default, the input matrices will be converted to observed / expected matrices. 
> In case you already have observed / expected matrices, you can use the `--converted-input` 
> flag to skip the observed / expected conversion. You can pre-compute observed/expected 
> matrices with `chess oe`

* A BED file with region information for the reference Hi-C matrix.
  This should be a tab-delimited file where each row contains chromosome name, start, 
  and end coordinates (exclusive) of the region. This file must 
  not contain any headers. If a fourth column is present, it is assumed to be a unique 
  identifier (index/name), which is then used to refer to that region in sparse matrix 
  format (see above).

> NOTE: The required formats for the Hi-C matrix (sparse) and regions (BED) are compatible 
> with the standard output of HiC-Pro (Servant et al. 2015).

* A normalized Hi-C matrix file in sparse matrix format for the query sample.

* A BED file with region information for the query Hi-C matrix.

* A BEDPE file (pairs) file that specifies which regions in the reference matrix should
  be compared to which regions in the query matrix. This should be a tab-delimited file
  with columns:
  \<ref chromosome\> \<ref start\> \<ref end\> \<qry chromosome\> \<qry start\> 
  \<qry end\> \<comparison name/id\> \<anything, not used\> \<ref strand\> \<qry strand\>
  The 8th column is required in order to match the BEDPE standard format 
  (see http://bedtools.readthedocs.io/en/latest/content/general-usage.html).
  This file must not contain any headers. The end coordinates are exclusive.

* The path to the output file.

The default mode with only the mandatory arguments will compare regions in the reference 
to the query using the region definitions in the pairs file. The output then contains the 
following columns:

- ID: the ID of the original pair, as in the supplied region pairs file
- SN: the signal-to-noise ratio of the matrices in the comparison, 
  which can be used for filtering
- ssim: The raw similarity score, useful to rank matrices by similarity
- z_ssim: the similarity z-score, in the context of all calculated similarity
  scores. Useful to assess if the similarity score is exceptionally high or low
  compared to other scores in the same run. For z-scores in the context of a 
  proper background model, see below.

#### Background models

To assess the statistical significance of your comparisons, you can enable a *background model*.

For inter-species comparisons, for example, it is useful to compare the original similarity 
score to the similarities obtained  by comparing the reference matrix to all other matrices in the
query. To enable this background model, use the `--background-query` parameter.

For a custom background model, you can include your own BED file with background regions using
`--background-regions`. You can generate these from a genome, for example, using `chess background`.

With an active background model, there are more columns in the output:

* ID: the ID of the original pair, as in the supplied region pairs file
* SN: the signal-to-noise ratio of the matrices in the comparison, 
  which can be used for filtering
* ssim: The raw similarity score, useful to rank matrices by similarity
* z_ssim: the similarity z-score, in the context of all calculated similarity
  scores. Useful to assess if the similarity score is exceptionally high or low
  compared to other scores in the same run.
* z_bg: the similarity z-score based on the background matrix comparisons.
* p_bg: The similarity score p-value. This assesses how likely it is to obtain 
  a similarity of ssim or higher given the comparisons in the background model
  - in this case, all other regions in the query genome of the same window size


#### Optional arguments

Optional arguments give you more control:

* `--converted-input` will skip the observed / expected conversion of the input matrices. 
  Use if you already have observed / expected matrices or want to compare matrices in different format.

* `-p <int>` lets you choose the number of cores that CHESS will use (default: 1).

* `--keep-unmappable-bins` disables the deletion of deletion of unmappable bins from matrices 
  before comparison. By default, bins that are marked as unmappable (have not contacts) 
  in either of two matrices in a comparison are deleted from both matrices. Disabling this 
  might give high similarity scores for matrices that happen to have unmappable bins at the 
  same positions.

* `--mappability-cutoff <float>` maximum allowed content of unmappable bins in a matrix. 
  In case a matrix has more unmappable bins, it will not be compared. (default: 0.1)

* `-r <float>, --relative-windowsize <float>` Fraction of the matrix size that will be used 
  for the window size parameter of the structural similarity function. (default: 1)

* `-a <float>, --absolute-windowsize <float>` Absolute value for the window size parameter 
  of the structural similarity function. Overwrites `-r`.


### oe

`chess oe` lets you convert your input matrix to observed / expected format. 
 Calling `chess sim` this is done automatically, but if you want to convert 
 your matrix for other reasons, you can use this.

`oe` takes three positional arguments (in that order):

* A Hi-C matrix file in sparse matrix format. This should be tab delimited, 
  where each line has three columns: \<row index\> \<column index\> \<value\>)

* A BED file with region information for the Hi-C matrix.
  This should be a tab-delimited file where each row contains chromosome name, 
  start, and end coordinates (exclusive) of the region. This file must 
  not contain any headers. If a fourth column is present, it is assumed to be a unique identifier
  (index/name), which is then used to refer to that region in sparse matrix format
  (see above).

* A path to the output matrix file. Will in the same format as the input matrix.


### pairs

`chess pairs` helps you generate the pairs input file to `chess sim`. There is no need 
to use this function for generating the pairs file, it is just here to make your life 
easier. `chess pairs` will perform a window slide with the specified input parameters. 
All positions of the window will be written to the specified output file in the format 
required by `chess sim`.

`pairs` takes four positional arguments (in that order):

* A genome id as recognized by pybedtools 
 (see [here](https://daler.github.io/pybedtools/autodocs/pybedtools.helpers.chromsizes.html)) 
  or a path to a tab-delimited file indicating the chromosome names and sizes of your genome 
  (no header, two columns: \<chromosome name\> \<chromosome size\>).

* The size of the window.

* The step of the window slide.

* The path to the output file.

You can use the following optional arguments to tweak `chess pairs`' behaviour:

* `--file-input` will force `pairs` to read from file. Won't attempt to use pybedtools 
  in that case, which is by default tried first. You can use this if you gave your 
  chromosome.sizes file a name that is also used as a UCSC identifier.

* `--chromosome` lets you restrict the pair generation to a specific chromosome, in 
  case you don't want all of them.


### filter

`chess filter` helps you filter the results file that is generated by `chess sim` 
and converts it to BED format, where the name (4th) column contains the comparison 
pair IDs. It outputs one BED file for all references and one for all queries.
You can then use `chess filter ... --genome-scan` to create a single BED file
if you run it on results from a `chess sim` without any background.

`filter` takes three positional arguments (in that order):

* A results file generated in a `chess sim` run.

* The pairs file used in the same `chess sim` run.

* The path to the output file. (REF and QUERY or genome_scan will be appended to that filename)

You can use the following parameters to specify how to filter the results. `mode` has to 
be `geq` (greated equal), `leq` (less equal), `l` (less) or `g` (greater). `value` 
defines the threshold value for the filter. You can combine any number of filters, 
as long as the corresponding columns are present in the results file. If you do not 
specify any filters, the results will simply be converted to BED format.

* `-p mode value`: filter by p-value (p_bg column)

* `-z mode value`: filter by z-score based on background comparisons (z_bg column)

* `-zs mode value`: filter by z-score based on all ssim values (z_ssim column)

* `-s mode value`: filter by ssim column

* `-n mode value`: filter by SN column

The following parameters give you some control over the output:

* `--score <string>` column in the results file that should be used for the score column in the output BED file(s). 

* `--genome-scan`: write only a single output BED. You can use this to avoid generating redundant files when the results you want to filter were generated in a `chess sim` run without background comparisons.

### extract
`chess extract` allows you to extract the specific features from significant different regions identified by `chess sim`. It outputs the coordinates of the gained and lost features, and its submatrix.

`chess extract` takes following mandatory arguments:

* A BEDPE file (pairs) file that specifies which regions in the reference matrix should
  be compared to which regions in the query matrix that have been observed to be significantly different by `chess sim`. This should be a tab-delimited file
  with columns:
  \<ref chromosome\> \<ref start\> \<ref end\> \<qry chromosome\> \<qry start\> \<qry end\> \<comparison name/id\> \<anything, not used\> \<ref strand\> \<qry strand\>
  The 8th column is required in order to match the BEDPE standard format (see http://bedtools.readthedocs.io/en/latest/content/general-usage.html).
  This file must not contain any headers. The end coordinates are exclusive.
  
* A normalized Hi-C matrix file in sparse matrix format for the reference sample.

* A BED file with region information for the reference Hi-C matrix.

* A normalized Hi-C matrix file in sparse matrix format for the query sample.

* A BED file with region information for the query Hi-C matrix.

* The path to the output directory.

### crosscorrelate
`chess crosscorrelate` allows to classify the extracted features by `chess extract` in order to obtain the main structural clusters. It will compute a 2D-crosscorrelation between all the features, and will calculate the optimal number of clusters according to the Elbow-index. It outputs the number of clusters and the classification of the specific features according to their structural pattern similarity.

`chess crosscorrelate` takes following mandatory arguments:

* The output file from `chess extract` that contains the submatrices of the features.

* A BEDPE file (pairs) file that specifies which regions in the reference matrix should
  be compared to which regions in the query matrix that have been observed to be significantly different by `chess sim`. This should be a tab-delimited file
  with columns:
  \<ref chromosome\> \<ref start\> \<ref end\> \<qry chromosome\> \<qry start\> \<qry end\> \<comparison name/id\> \<anything, not used\> \<ref strand\> \<qry strand\>
  The 8th column is required in order to match the BEDPE standard format (see http://bedtools.readthedocs.io/en/latest/content/general-usage.html).
  This file must not contain any headers. The end coordinates are exclusive.
  
 * Path to store the outputs.


## References

* Hug, Clemens B., Alexis G. Grimaldi, Kai Kruse, and Juan M. Vaquerizas. 2017. “Chromatin Architecture Emerges during Zygotic Genome Activation Independent of Transcription.” Cell 169 (2). Elsevier: 216–228.e19. doi:10.1016/j.cell.2017.03.024.
* Rao, Suhas S.P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. 2014. “A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping.” Cell 159 (7). Elsevier: 1665–80. doi:10.1016/j.cell.2014.11.021.
* Servant, Nicolas, Nelle Varoquaux, Bryan R. Lajoie, Eric Viara, Chong-Jian Chen, Jean-Philippe Vert, Edith Heard, Job Dekker, and Emmanuel Barillot. 2015. “HiC-Pro: An Optimized and Flexible Pipeline for Hi-C Data Processing.” Genome Biology 16 (1). BioMed Central: 259. doi:10.1186/s13059-015-0831-x.
