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
we can use the `chess pairs` subcommand to generate the 3. input file.

-----------------------------
Generating a pairs input file
-----------------------------

In this example analysis, we will search the entire X chromosome for differences.
For this, we first need to generate the pairs file with the ``chess pairs``
subcommand.
We will compare regions of 1 Mb and 2 Mb size with a step size of 100 kb.
In addition, `chess pairs` needs to know the sizes of the chromosomes for which
we want to generate the pairs. Here we supply these with the chrom.sizes.tsv
file in the examples folder.

.. code:: bash

  chess pairs examples/Dmel_genome_scan/dm6.chrom.sizes.tsv 2000000 100000 \
  ./dm6_2mb_win_100kb_step.bed --chromosome X --file-input

  chess pairs examples/Dmel_genome_scan/dm6.chrom.sizes.tsv 1000000 100000 \
  ./dm6_1mb_win_100kb_step.bed --chromosome X --file-input

------------------
Running the search
------------------

With all necessary input files prepared, we can run the search with ``chess sim``:

.. code:: bash

  chess sim \
  examples/Dmel_genome_scan/juicer/zld.hic@25000 \
  examples/Dmel_genome_scan/juicer/wt.hic@25000 \
  ./dm6_2mb_win_100kb_step.bed \
  ./dm6_2mb_win_100kb_step_chess_results.tsv

  chess sim \
  examples/Dmel_genome_scan/juicer/zld.hic@25000 \
  examples/Dmel_genome_scan/juicer/wt.hic@25000 \
  ./dm6_1mb_win_100kb_step.bed \
  ./dm6_1mb_win_100kb_step_chess_results.tsv

The output data are stored in the 
./dm6_1mb_win_100kb_step_chess_results.tsv
and 
./dm6_2mb_win_100kb_step_chess_results.tsv
files.

=========================================
Comparing regions between Mouse and Human
=========================================
