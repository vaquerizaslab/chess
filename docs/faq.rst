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