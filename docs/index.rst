.. pyGenomeTracks documentation master file, created by
   sphinx-quickstart on Thu Jan 23 17:49:39 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyGenomeTracks's documentation!
==========================================

Standalone program and library to plot beautiful genome browser tracks
----------------------------------------------------------------------

pyGenomeTracks aims to produce high-quality genome browser tracks that
are highly customizable. Currently, it is possible to plot:

 * bigwig
 * bed/gtf (many options)
 * bedgraph
 * bedgraph matrices (like TAD-separation scores)
 * epilogos
 * narrow peaks
 * links (represented as arcs, triangles or squares)
 * Hi-C matrices (as triangle or squares)
 * fasta
 * maf (multiple alignment format)

Here is a scheme which describe how pyGenomeTracks is working (graphical abstract of `Lopez-Delisle et al. 2020 <https://doi.org/10.1093/bioinformatics/btaa692>`_):

.. image:: content/images/graphicalabstract.png

pyGenomeTracks can make plots with or without Hi-C data. The following is an example output of
pyGenomeTracks from `Ramírez et al. 2017 <https://www.nature.com/articles/s41467-017-02525-w>`_.

.. image:: content/images/hic_example_nat_comm_small.png


There are 3 ways for using pyGenomeTracks:

* **Galaxy usage** --  the public `European Galaxy server <http://usegalaxy.eu>`_ let's you use pyGenomeTracks within the familiar Galaxy framework without the need to master the command line
* **command line usage** -- simply download and install the tool (see :doc:`content/installation` and :doc:`content/usage`)



Table of content
----------------
.. toctree::
   :maxdepth: 1
   :caption: Contents:

   content/installation
   content/usage
   content/citation
   content/all_tracks
   content/examples
   content/possible-parameters
   content/adding-new-tracks
   content/releases
   content/faq
