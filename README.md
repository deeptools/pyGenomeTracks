[![PyPI Version](https://img.shields.io/pypi/v/pyGenomeTracks.svg?style=plastic)](https://pypi.org/project/deepTools/) [![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=plastic)](https://anaconda.org/bioconda/pygenometracks)

pyGenomeTracks
==============

Standalone program and library to plot beautiful genome browser tracks
----------------------------------------------------------------------

pyGenomeTracks aims to produce high-quality genome browser tracks that
are highly customizable. Currently, it is possible to plot:

 * bigwig 
 * bed (many options)
 * bedgraph
 * links (represented as arcs) 

![pyGenomeTracks example](https://github.com/maxplanck-ie/pyGenomeTracks/raw/master/pygenometracks/tests/test_data/master_plot.png)

The configuration file for this image is [here](https://github.com/maxplanck-ie/pyGenomeTracks/blob/master/pygenometracks/tests/test_data/browser_tracks.ini)

Installation
------------
pyGenomeTracks works with python 2.7 and python 3.6.

Currently, the best way to install pyGenomeTracks is

```bash
$ pip install pyGenomeTracks
```

If the latest version wants to be installed use:

```bash
$ pip install  git+https://github.com/maxplanck-ie/pyGenomeTracks.git
```


Usage
-----
To run pyGenomeTracks a configuration file describing the tracks is required. The easiest way to create this file is using `make_tracks_file` which setups up a file with defaults that can be easily changed. The format is:

```bash
$ make_tracks_file --trackFiles <file1.bed> <file2.bw> ... -o tracks.ini
```

`make_tracks_file` uses the file ending to guess the file type. 

Then, a region can be plotted using:

```bash
$ pyGenomeTracks --tracks tracks.ini --region chr2:10,000,000-11,000,000 -o nice_image.pdf
```

pyGenomeTracks is used by [HiCExporer](https://hicexplorer.readthedocs.io/) and [HiCBrowser](https://github.com/maxplanck-ie/HiCBrowser) (See e.g. [Chorogenome navigator](http://chorogenome.ie-freiburg.mpg.de/) which is made with HiCBrowser)
