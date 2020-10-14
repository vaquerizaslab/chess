****************
Example analysis
****************

CHESS has two typical use-cases, which are showcased in the publication (ADD LINK):

1. Finding genomic regions with differences in their chromatin conformation between
   two samples. This is useful for example for the identification of
   disease driving changes.

2. Comparing the chromatin conformation between two distinct regions in the same
   or in different genomes. This is useful for example for quantifying the grade
   of conservation of chromatin conformation in evolution.

===========================================================
Finding differences between samples of different conditions
===========================================================

This example illustrates how CHESS can help to identify genomic regions
that differ in their chromatin conformation between two samples, by the
example of samples from a `DLBCL <https://en.wikipedia.org/wiki/Diffuse_large_B-cell_lymphoma>`
patient and a control person. The analysis of these data is also showcased in
the `chess paper <SOME LINK>`, Figure 5.

For this analysis, we will need only three input files for CHESS:

1. The Hi-C data for the first sample
2. The Hi-C data for the second sample
3. The coordinates of the regions that should be compared in bedpe format.

1. and 2. can be in `Juicer <https://github.com/aidenlab/juicer>`_,
`Cooler <https://github.com/mirnylab/cooler>`_ or `FANC <https://github.com/vaquerizaslab/fanc>`_ format.

3. should be supplied in `BEDPE <https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format>`_
format.

If we want to scan a whole chromosome or genome for differences,
we can use the :ref:`chess pairs <chess-pairs>` subcommand to generate the
3. input file.

To follow along with the example analysis, first move to the data location:

.. code:: bash

  cd examples/dlbcl

-----------------------------
Generating a pairs input file
-----------------------------

In this example analysis, we will search the entire chromosome 2 for differences.
For this, we first need to generate the pairs file with the
:ref:`chess pairs <chess-pairs>` subcommand.
We will compare regions of 2 Mb size with a step size of 100 kb.
In addition, ``chess pairs`` needs to know the sizes of the chromosomes for which
we want to generate the pairs. Here we supply 'hg38' as an genome identifier:

.. code:: bash

  chess pairs hg38 2000000 100000 \
  ./hg38_chr2_2mb_win_100kb_step.bed --chromosome chr2

------------------
Running the search
------------------

With all necessary input files prepared, we can run the search with
the :ref:`chess sim <chess-sim>` subcommand:

.. code:: bash

  chess sim \
  ukm_control_fixed_le_25kb_chr2.hic \
  ukm_patient_fixed_le_25kb_chr2.hic \
  hg38_chr2_2mb_win_100kb_step.bed \
  ukm_chr2_2mb_control_vs_patient_chess_results.tsv

The output data are stored in the
ukm_chr2_2mb_control_vs_patient_chess_results.tsv file.

----------------------
Inspecting the results
----------------------

To inspect which regions show differences between the two samples,
it is useful to first inspect the similarity index along the compared
chromosome. Aside from the calculated similarity index, the signal-to-noise
ratio is important. One way to display this information together is this
(please check out the jupyter notebook at `examples/dlbcl/example_analysis.ipynb`
for the code to generate these plots):

.. figure:: plots/chr2_2mb_results_track.png
   :name: result-track

Here we display only regions with a signal-to-noise ration > 0.6 as colored
dots, the rest is gray. The most interesting regions are highlighted by
deep dips with above threshold signal-to-noise ratios, for example the
region around 148 Mb:

.. figure:: plots/chr2_148mb_example_region.png
   :name: result-region

----------------------
Choosing a window size
----------------------

In this analysis, we compared windows of 2 Mb size between our samples.
In general, choosing a different window size should be correlated,
with large windows simply averaging over the effects observed in smaller
windows. If we repeat the analysis above for 1 Mb and 4 Mb windows, and compare
the results by the window midpoints, we get the following relationships:

.. figure:: plots/ws_corr.png
   :name: ws-corr

.. figure:: plots/ws_track_corr.png
   :name: ws-track-corr

Despite the correlation, different window sizes can yield different results
in some regions:

* Larger windows cover more and longer long-range interactions;
  - If you are interested in changes of large effects stretching over 
    long genomic distances, choose a larger window size.
  - However, long-range interactions tend to be more noisy.
    The larger the window size, the smaller the number of regions that will
    pass a given signal-to-noise threshold. If your analysis does not return
    any regions of strong dissimilarity above your signal-to-noise threshold,
    lower the threshold or try a smaller window size.
* The larger the window, the smaller the effect of small changes;
  - If you are interested in finding changes in single TAD boundaries, 
    choose a small window. Large windows will cover multiple boundaries 
    and the score of the window will reflect their combined change.

=========================================
Comparing regions between Mouse and Human
=========================================
