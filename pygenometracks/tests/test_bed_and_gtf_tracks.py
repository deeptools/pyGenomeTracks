# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

browser_tracks = """
[x-axis]
where = top
title = where=top

[spacer]
height = 0.05

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style=flybase;fontsize=10
style=UCSC
fontsize = 10

[spacer]
height = 1

[test bed6]
file = dm3_genes.bed6.gz
height = 7
title = bed6 border color=black; gene rows=10; fontsize=8; color=Reds (when a color map is used for the color (e.g. coolwarm, Reds) the bed score column mapped to a color)
fontsize = 7
file_type = bed
color = Reds
border color = black
gene rows = 10

[spacer]
height = 1

[test bed4]
file = dm3_genes.bed4.gz
height = 10
title = bed4 fontsize=10; line width=1.5; global max row=yes (global max row sets the number of genes per row as the maximum found anywhere in the genome, hence the white space at the bottom)
fontsize = 10
file_type = bed
global max row = yes
interval_height = 200
line width = 1.5


[spacer]
height = 1

[test gtf]
file = dm3_subset_BDGP5.78.gtf.gz
height = 10
title = gtf from ensembl
fontsize = 12
file_type = bed


[spacer]
height = 1

[test bed]
file = dm3_subset_BDGP5.78_asbed_sorted.bed.gz
height = 10
title = gtf from ensembl in bed12
fontsize = 12
file_type = bed


[spacer]
height = 1

[x-axis]
fontsize=30
title = fontsize=30

"""
with open(ROOT + "bed_and_gtf_tracks.ini", 'w') as fh:
    fh.write(browser_tracks)


tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_bed_and_gtf():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0}/bed_and_gtf_tracks.ini --region X:3000000-3300000 --trackLabelFraction 0.2 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_bed_and_gtf.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
