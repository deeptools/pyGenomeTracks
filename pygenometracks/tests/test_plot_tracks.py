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

[tads]
file = tad_classification.bed
title = TADs color=bed_rgb; border color = black
file_type = domains
border color = black
color = bed_rgb
height = 5

[tads 2]
file = tad_classification.bed
title = TADs orientation=inverted; color=#cccccc; border color=red
file_type = domains
border color = red
color = #cccccc
orientation = inverted
height = 3

[spacer]
height = 0.5

[tad state]
file = chromatinStates_kc.bed.gz
height = 1.2
title = bed display=interleaved; labels=off
display = interleaved
labels = off

[spacer]
height = 0.5

[tad state]
file = chromatinStates_kc.bed.gz
height = 0.5
title = bed display=collapsed; color=bed_rgb
display = color from bed rgb
labels = off
color = bed_rgb
display = collapsed

[spacer]
height = 0.5

[test bedgraph]
file = bedgraph_chrx_2e6_5e6.bg
color = blue
height = 3
title = bedgraph color=blue
max_value = 100

[test arcs]
file = test.arcs
title = links orientation=inverted
orientation = inverted
line style = dashed
height = 2

[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
height = 3
title = bigwig number of bins=2000
number of bins = 2000

[test bigwig overlay]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
title = color:red; max_value=50; number of bins=100 (next track: overlay previous=yes; max_value=50; show data range=no; color=#0000FF80 (blue, with alpha 0.5))
min_value =0
max_value = 50
height = 4
number of bins = 100

[test bigwig overlay]
file = bigwig_chrx_2e6_5e6.bw
color = #0000FF80
title =
min_value =0
max_value = 50
show data range = no
overlay previous = yes
number of bins = 100

[test arcs]
file = test.arcs
color = red
line width = 3
title = links line width=3
height = 3

[spacer]
height = 0.5
title = height=0.5

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style=flybase;fontsize=10
style=flybase
fontsize = 10

[spacer]
height = 1

[test gene rows]
file = dm3_genes.bed.gz
height = 3
title = max num rows=3; style=UCSC
fontsize = 8
style = UCSC
gene rows = 3

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

[test bed6]
file = dm3_genes.bed6.gz
height = 10
title = bed6 fontsize=10; line width=1.5; global max row=yes (global max row sets the number of genes per row as the maximum found anywhere in the genome, hence the white space at the bottom)
fontsize = 10
file_type = bed
global max row = yes
interval_height = 200
line width = 1.5

[x-axis]
fontsize=30
title = fontsize=30

[vlines]
file = tad_classification.bed
type = vlines

"""
with open(ROOT + "browser_tracks.ini", 'w') as fh:
    fh.write(browser_tracks)

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth=200000; transform=log1p;min_value=5
depth = 200000
min_value=5
transform = log1p
file_type = hic_matrix
show_masked_bins = no

[hic matrix]
file = Li_et_al_2015.h5
title = depth=250000; orientation=inverted; colormap=PrRd; min_value=5; max_value=70
min_value = 5
max_value = 70
depth = 250000
colormap = PuRd
file_type = hic_matrix
show_masked_bins = no
orientation = inverted

[spacer]
height = 0.5

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth=300000; transform=log1p; colormap Blues (TADs: overlay previous=yes; line width=1.5)
colormap = Blues
min_value = 10
max_value = 150
depth = 300000
transform = log1p
file_type = hic_matrix

[tads]
file = tad_classification.bed
#title = TADs color=bed_rgb; border color = black
file_type = domains
border color = black
color = none
height = 5
line width = 1.5
overlay previous = share-y

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth=250000; transform=log1p; colormap=bone_r (links: overlay previous=share-y; links type=triangles; color=darkred; line style=dashed, bigwig: color=red)
colormap = bone_r
min_value = 15
max_value = 200
depth = 250000
transform = log1p
file_type = hic_matrix
show_masked_bins = no

[test arcs]
file = links2.links
title =
links type = triangles
line style = dashed
overlay previous = share-y
line width = 0.8
color = darkred

[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title =
overlay previous = yes
min_value = 0
max_value = 50

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth=200000; show_masked_bins=yes;
depth = 200000
file_type = hic_matrix
show_masked_bins = yes

[spacer]
height = 0.1

[x-axis]

"""

with open(ROOT + "browser_tracks_hic.ini", 'w') as fh:
    fh.write(browser_tracks_with_hic)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0}/browser_tracks.ini --region X:3000000-3500000 --trackLabelFraction 0.2 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_plot.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0}/browser_tracks_hic.ini --region X:2500000-3500000 --trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_plot_hic.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
