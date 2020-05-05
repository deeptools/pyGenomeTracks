[![PyPI Version](https://img.shields.io/pypi/v/pyGenomeTracks.svg?style=plastic)](https://pypi.org/project/pyGenomeTracks/) [![bioconda-badge](https://img.shields.io/conda/vn/bioconda/pyGenomeTracks.svg?style=plastic)](https://anaconda.org/bioconda/pygenometracks)
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=plastic)](http://bioconda.github.io)
[![Build Status](https://travis-ci.org/deeptools/pyGenomeTracks.svg?branch=master)](https://travis-ci.org/deeptools/pyGenomeTracks)


pyGenomeTracks
==============

Standalone program and library to plot beautiful genome browser tracks
----------------------------------------------------------------------

pyGenomeTracks aims to produce high-quality genome browser tracks that
are highly customizable. Currently, it is possible to plot:

 * bigwig
 * bed/gtf (many options)
 * bedgraph
 * epilogos
 * narrow peaks
 * links (represented as arcs)
 * Hi-C matrices

pyGenomeTracks can make plots with or without Hi-C data. The following is an example output of
pyGenomeTracks from [Ramírez et al. 2017](https://www.nature.com/articles/s41467-017-02525-w)

![pyGenomeTracks example](./docs/content/images/hic_example_nat_comm_small.png)

Table of content
----------------
  * [Installation](#installation)
  * [Usage](#usage)
  * [Citation](#citation)
  * [Documentation](#documentation)
  * [External users](#external-users)



Installation
------------
pyGenomeTracks works with python >=3.6.

Currently, the best way to install pyGenomeTracks is with anaconda

```bash
$ conda install -c bioconda -c conda-forge pygenometracks
```

Also, pyGenomeTracks can be installed using pip

```bash
$ pip install pyGenomeTracks
```

If the latest version wants to be installed use:

```bash
$ pip install  git+https://github.com/deeptools/pyGenomeTracks.git
```


Usage
-----
To run pyGenomeTracks a configuration file describing the tracks is required. The easiest way to create this file is using the program `make_tracks_file` which creates a configuration file with
defaults that can be easily changed. The format is:

```bash
$ make_tracks_file --trackFiles <file1.bed> <file2.bw> ... -o tracks.ini
```

`make_tracks_file` uses the file ending to guess the file type.

Then, a region can be plotted using:

```bash
$ pyGenomeTracks --tracks tracks.ini --region chr2:10,000,000-11,000,000 --outFileName nice_image.pdf
```

The ending `--outFileName` defines the image format. If `.pdf` is used, then the resulting image is a pdf. The options are pdf, png and svg.

Description of other possible arguments:
<!--- Start of possible arguments of pgt -->
``` text
optional arguments:
  -h, --help            show this help message and exit
  --tracks TRACKS       File containing the instructions to plot the tracks.
                        The tracks.ini file can be genarated using the
                        `make_tracks_file` program.
  --region REGION       Region to plot, the format is chr:start-end
  --BED BED             Instead of a region, a file containing the regions to
                        plot, in BED format, can be given. If this is the
                        case, multiple files will be created using a prefix
                        the value of --outFileName
  --width WIDTH         figure width in centimeters (default is 40)
  --height HEIGHT       Figure height in centimeters. If not given, the figure
                        height is computed based on the heights of the tracks.
                        If given, the track height are proportionally scaled
                        to match the desired figure height.
  --title TITLE, -t TITLE
                        Plot title
  --outFileName OUTFILENAME, -out OUTFILENAME
                        File name to save the image, file prefix in case
                        multiple images are stored
  --fontSize FONTSIZE   Font size for the labels of the plot (default is 0.3 *
                        figure width)
  --dpi DPI             Resolution for the image in case the ouput is a raster
                        graphics image (e.g png, jpg) (default is 72)
  --trackLabelFraction TRACKLABELFRACTION
                        By default the space dedicated to the track labels is
                        0.05 of the plot width. This fraction can be changed
                        with this parameter if needed.
  --trackLabelHAlign {left,right,center}
                        By default, the horizontal alignment of the track
                        labels is left. This alignemnt can be changed to right
                        or center.
  --decreasingXAxis     By default, the x-axis is increasing. Use this option
                        if you want to see all tracks with a decreasing
                        x-axis.
  --version             show program's version number and exit
```
<!--- End of possible arguments of pgt -->

Citation
---------
If you use pyGenomeTracks in your analysis, you can cite the following paper :

Fidel Ramírez, Vivek Bhardwaj, Laura Arrigoni, Kin Chung Lam, Björn A. Grüning, José Villaveces, Bianca Habermann, Asifa Akhtar & Thomas Manke. High-resolution TADs reveal DNA sequences underlying genome organization in flies. Nature Communications (2018) [doi:10.1038/s41467-017-02525-w](https://www.nature.com/articles/s41467-017-02525-w)

Documentation
-------------

Our [documentation](http://pygenometracks.readthedocs.io/) provide [examples](http://pygenometracks.readthedocs.org/en/latest/content/examples.html), as well as the [full list of possible parameters](http://pygenometracks.readthedocs.org/en/latest/content/possible-parameters.html) and [guidelines for developers who would like to add a new track type](http://pygenometracks.readthedocs.org/en/latest/content/adding-new-tracks.html).

<!-- I do not know what to do with that, is it External users? 
pyGenomeTracks is used by [HiCExporer](https://hicexplorer.readthedocs.io/) and [HiCBrowser](https://github.com/maxplanck-ie/HiCBrowser) (See e.g. [Chorogenome navigator](http://chorogenome.ie-freiburg.mpg.de/) which is made with HiCBrowser)
 -->
External users
--------------

* [CoolBox](https://github.com/GangCaoLab/CoolBox) is an interactive genomic data explorer for Jupyter Notebooks
* [Galaxy](https://usegalaxy.eu/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/pygenometracks/pygenomeTracks) integration offers a graphical user-interface to create PGT plots. It is also possible to include PGT into workflows and automatic pipelines.
