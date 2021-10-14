********
Overview
********

CHESS is a command line tool for the quantitative comparison and automatic feature extraction for chromatin contact
data.

Please note that the structural similarity index cannot distinguish between differences in Hi-C datasets due to random
noise and
differences due to biologically relevant changes in genome organisation. In order to address this, CHESS outputs a
signal-to-noise metric (SN), which measures the mean magnitude of values in a difference matrix of the input Hi-C
datasets, divided by the variance across the matrix. This was introduced in the original CHESS publication in order
to control for noise effects, based on the observation that in noisy regions, values of neighbouring pixels of a
difference matrix often have opposite signs, fluctuating randomly, which leads to high variance. In regions with
visually striking differences such as differential domain organisation, the values in the differential region of the
matrix have low variance and high mean absolute value.

The level of noise in the Hi-C matrices used as input to CHESS will depend on the window size and resolution used, as
well as the sequencing depth of the input samples. Therefore, SN and SSIM distributions are not comparable across
comparisons that use different parameters (window size, resolution, sequencing depth). For ease of interpretation of
the output, we advise using Hi-C datasets with similar sequencing depth, especially if you are performing multiple
pairwise comparisons. See `our preprint <https://www.biorxiv.org/>`_)
for more information, and the FAQ for further advice on selecting SSIM and SN thresholds.

