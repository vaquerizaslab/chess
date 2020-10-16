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

.. argparse::
   :module: chess.commands
   :func: sim_parser
   :prog: chess
   :nodescription:

   
.. _chess-extract

***************
``chess extract``
***************

This command extracts features from a set of input regions.

Please note that the input parameters have to be fine tuned depending on the
size of the analyzed regions and the target features.
For now, some experimentation by the user is required, but we are planning to 
release a guide to this in the future.

.. argparse::
   :module: chess.commands
   :func: extract_parser
   :prog: chess
   :nodescription: