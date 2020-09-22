###############
Getting started
###############

.. contents::
   :depth: 3

********
Overview
********

CHESS is a command line tool for the quantitative comparison and automatic feature extraction for chromatin contact data.


.. _chess_installation:

************
Installation
************

.. _prerequisites:

=============
Prerequisites
=============

CHESS was written and tested using Python 3.6 - 3.8.2 on Scientific Linux (6.10) and OS X El Capitan (10.11.6) and Ubuntu Bionic Beaver (18.04) and Ubuntu Focal Fossa (20.04).

It requires the packages `cython` (0.29.16), `scipy` (1.0.0), `numpy` (1.14.0), `scikit-image` (0.13.1), `pandas` (0.22.0),
`pathos` (0.2.1), `kneed` (0.6.0), `tqdm` (4.43.0), `intervaltree` (3.0.2), `pybedtools` (0.8.1), `future` (0.16.0) and `fanc` (0.8.28).

**When installing with a python version < 3.8.2, Please make sure to have `cython` (0.29.16) installed before installling CHESS via `pip` or `python install setup.py`, or use a `pip` version >= 20.0.2.**

=================================================
Installation from the Python Package Index (PyPi)
=================================================

CHESS can be installed quickly from the command line using PyPi:

.. code:: bash
  
  pip install chess-hic


========================
Installation from source
========================

You can also download the code from our `GitHub repo <https://github.com/vaquerizaslab/chess>`_
and install CHESS manually. Make sure to have `cython` installed. The other dependencies should be downloaded autmatically; if you encounter problems, try to install the other :ref:`prerequisites <prerequisites>` manually.

.. code:: bash

  git clone https://github.com/vaquerizaslab/chess  # or download manually
  pip install chess


or

.. code:: bash

  git clone https://github.com/vaquerizaslab/chess  # or download manually
  cd chess
  python setup.py install


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



=========================================
Comparing regions between Mouse and Human
=========================================
