##########################
Frequently Asked Questions
##########################


*****************************
Why is ``chess sim`` so slow?
*****************************
A common cause for long runtimes of ``chess sim`` is when the input files are
in `Cooler <https://github.com/mirnylab/cooler>`_ format. ``chess sim`` compares
observed/expected transformed matrices, and these are often not precomputed
and stored in this format. If your runs spend a long time in the 
"INFO Loading query contact data" step, consider converting your data to
`Juicer <https://github.com/aidenlab/juicer>`_, or `FANC <https://github.com/vaquerizaslab/fanc>`_ 
format, as these usually store observed/expected values
(see `here <https://fan-c.readthedocs.io/en/latest/fanc-executable/compatibility.html#cooler-cool-and-mcool>`_).

******************************************************************************
Why do I get "ERROR No valid region pairs found" when using Juicer .hic files?
******************************************************************************

Juicer `removes the "chr" prefix from chromosome names <https://github.com/aidenlab/juicer/issues/118#issuecomment-510181960>`_.
This can cause chromosome names in Juicer .hic files to not match chromosome names in pairs BEDPE files created with
``chess pairs``. Please check whether the chromosome names in your pairs file match those in your .hic file.

********************************************************************************
Can I use CHESS to compare the results of loop calling software between samples?
********************************************************************************

Some users have asked whether it's possible to use CHESS to compare loop calls between samples, by supplying the BEDPE
file output by a loop calling tool as the input pairs file for CHESS. This not recommended, for reasons explained
below.

BEDPE files used to specify loops typically specify a small (< 10x10 bins) region that is located away from the
diagonal of the matrix. In contrast, the BEDPE files produced by ``chess pairs`` and used as input to ``chess sim``
specify large (e.g., 100x100 bins) regions that are located at the diagonal of the matrix.

While it may be possible to use CHESS to examine off-diagonal regions, this is so far **untested** and may
not work as expected! Instead, we recommend creating on-diagonal windows that span the width of the whole loop, or using
``chess pairs`` to create windows that tile the genome and run a genome-wide comparison between your samples.

************************************************************************************************
I used the SSIM Z-score and SN thresholds from the CHESS paper. Why don't I see any differences?
************************************************************************************************

Some users have reported that they have selected regions with SSIM Z-score < -1.2 and SN > 0.6,
as was done in Galan, Machnik, et al. 2020, but they don't get any regions that pass both these
thresholds, or that regions that pass these thresholds don't have any visible differences, or
that the results are very different between different datasets or resolutions. This could be because
there are no differences between the datasets, but it could also happen because these thresholds are
not appropriate for your data. Both SSIM and SN values are influenced by the parameters you choose,
like matrix resolution and region size, as well as by sequencing depth as this influences how noisy
the Hi-C matrix is at a given resolution. Therefore, it is important to choose SSIM and SN thresholds
that make sense for your data.

In particular, it is crucial to filter based on an appropriate SN threshold as well as SSIM, as some
genomic regions may have low SSIM due to noise. This occurs more often in relatively unstructured
regions of the genome. As shown in Galan, Machnik, et al. 2020
(`Extended Data Figure 2 <https://www.nature.com/articles/s41588-020-00712-y/figures/8>`_), the SSIM
score is largely influenced by the structure and contrast of the matrices. Note that SSIM is
calculated using the observed/expected matrices. In these matrices, a region with little 3D genome
structure will have low contrast and structure (low standard deviation of values across the
observed/expected matrix) as the observed interactions are similar to the expected values. Therefore,
the SSIM will also be low. The SN will also be low when two regions with little structure are compared,
so filtering for regions above a certain SN threshold will remove these regions.

***************************************************************
How do I choose an appropriate SSIM / SN threshold for my data?
***************************************************************

As discussed above, since the absolute values of both SSIM and SN vary between experiments due to different Hi-C
resolutions, window sizes, sequencing depth, etc, it is impossible to choose a single threshold that works for all
datasets.

For SSIM, the initial recommendation was to use the Z-normalised SSIM score from the CHESS comparison of two datasets.
This provides an easy to interpret way of identifying regions with larger differences than the genome-wide average.
We still recommend this approach.

The same approach can be applied to SN, but may be overly stringent. For example, for the DLBCL-control comparison
shown in Galan, Machnik et al. 2020 (`Figure 5 <https://www.nature.com/articles/s41588-020-00712-y/figures/5>`_),
selecting regions with an SN Z-score > 1 and a SSIM Z-score < -1.2
returns only 16 regions, compared to 61 with the original thresholds used in the publication.

Ideally, one would identify the distributions of SSIM and SN values for comparisons where no biologically meaningful
differences are expected, and use these to choose thresholds. For example, analogous to a 10% FDR, one might choose
the 90th percentile of such a reference SN distribution as the SN threshold for a comparison of interest, and the
10th percentile of such a reference SSIM distribution (note that the combination of two thresholds becomes more
stringent than a straightforward 10% FDR; you may prefer to relax these thresholds and use, for example, the 80th
percentile of the reference SN distribution and the 20th percentile of SSIM).

To produce these reference distributions of SSIM and SN, one could compare biological replicates or compare to
similar published data. We have tested this approach using comparisons of biological replicates (using Hi-C data
from Rao et al. 2014 and Ing-Simmons et al. 2021), and it performed well. If more than two biological replicates are
available, and comparisons between these do not give highly similar thresholds, you can consider using the most or
least stringent pair of thresholds, depending on how conservative you wish to be for your downstream analysis.

Because SSIM and SN are impacted by the parameters used for CHESS (window size, resolution) and also by noise, which
is impacted by sequencing depth, these thresholds should only be applied to CHESS comparisons performed using the
same parameters and similar sequencing depths. Since sequencing depth can vary per-chromosome (e.g. in samples from
XY cells, X and Y chromosomes will have relatively less coverage than the autosomes), you may want to consider using
per-chromosome thresholds. In order to apply thresholds to datasets with different sequencing depths (e.g. merged
replicates), Z-normalisation of SSIM and SN reference distributions can be performed before identifying percentile
thresholds. We have not yet extensively tested this approach, so it should be used with caution. Note that
Z-normalisation is not appropriate for comparisons where there is an expectation of globally reduced similarity,
e.g. cohesin / CTCF depletions that abolish TAD structures.

*******************************************
What if I don't have biological replicates?
*******************************************

It's not always feasible to have biological replicates. In this case, you could try running a CHESS comparison to a
reference dataset where you don't expect any biological differences. This could be comparable control data from
another lab, technical replicates, or even synthetic mixtures of your control and treatment Hi-C datasets. You can
then subtract your reference SSIM profiles from your real treatment-control SSIM profiles. Any remaining local minima
should highlight regions where the differences between treatment and control are greater than the differences between
the reference datasets. We have used this approach in some of our publications (Galan et al. 2020 and Ing-Simmons et
al. 2021).