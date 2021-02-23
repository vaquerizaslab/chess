###############################
The chess-hic command line tool
###############################

.. contents::
   :depth: 1


.. _chess-pairs:

***************
``chess pairs``
***************

This command allows the automatic creation of a ``bedpe`` file as required
by ``chess sim`` if its intended use is to search a whole genome or a number
of chromosomes for regions with differences in their chromatin conformation
between two samples.

The window size should not be smaller than 20x the bin size of the data
that will be compared with :ref:`chess sim <chess-sim>`. We recommend 
to use window sizes not smaller than 100x the bin size of the data.

.. argparse::
   :module: chess.commands
   :func: pairs_parser
   :prog: chess
   :nodescription:

  
.. _chess-sim:

***************
``chess sim``
***************

This command runs the comparisons between two sets of chromatin contact data.

Please note that the size of the regions passed in the `pairs` file cannot
be smaller than 20x the bin size of the data. We recommend to use regions
spanning at least 100x the bin size of the data.

.. argparse::
   :module: chess.commands
   :func: sim_parser
   :prog: chess
   :nodescription:

   
.. _chess-extract:

*****************
``chess extract``
*****************

This command extracts features from a set of input regions.

Please note that the input parameters have to be fine tuned depending on the
size of the analyzed regions and the target features.
For now, some experimentation by the user is required, but we are planning to 
release a guide to this in the future.

This command will write two files into the output directory:
gained_features.tsv and lost_features.tsv,
for the gained and lost features in the query matrix compared to the reference,
respectively.
These files contain the information about the position, and the shape of the features.
The first value of each row correspond to the region ID (same as in the `chess sim` output)
the second to an ID of the feature. The following four values correspond to the
position of the corners of the rectangle that contain the feature
(xmin, xmax, ymin and ymax) in the region matrix.
The following columns contain the contact values for the portion of the region
matrix belonging to the feature in the query matrix (gained_features.tsv)
or reference matrix (lost_features.tsv).

.. argparse::
   :module: chess.commands
   :func: extract_parser
   :prog: chess
   :nodescription:


.. _chess-crosscorrelate:

************************
``chess crosscorrelate``
************************

This command clusters extracted features by their topology.

.. argparse::
   :module: chess.commands
   :func: crosscorrelate_parser
   :prog: chess
   :nodescription: