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
title = TADs color=bed_rgb; border_color = black
file_type = domains
border_color = black
color = bed_rgb
height = 5

[tads 2]
file = tad_classification.bed
title = TADs orientation=inverted; color=#cccccc; border_color=red
file_type = domains
border_color = red
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
height = 1.5
title = bedgraph color=blue
max_value = 100

[test arcs]
file = test.arcs
title = links orientation=inverted
orientation = inverted
line_style = dashed
height = 2

[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
height = 1.5
title = bigwig number of bins=2000
number_of_bins = 2000

[spacer]

[test bigwig overlay]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
title = color:red; max_value=50; number_of_bins=100 (next_track: overlay_previous=yes;
        max_value=50; show_data_range=no; color=#0000FF80 (blue, with alpha 0.5))
min_value =0
max_value = 50
height = 2
number_of_bins = 100

[test bigwig overlay]
file = bigwig_chrx_2e6_5e6.bw
color = #0000FF80
title =
min_value =0
max_value = 50
show_data_range = no
overlay_previous = yes
number_of_bins = 100

[spacer]
height = 1

[tads 3]
file = tad_classification.bed
title = TADs color=#cccccc; border_color=red (next track: overlay_previous=share-y links; type = loops)
file_type = domains
border_color = red
color = #cccccc
height = 3

[test arcs overlay]
file = test.arcs
color = red
line_width = 10
links_type = loops
overlay_previous = share-y

[test arcs]
file = test.arcs
line width = 3
color = RdYlGn
title = links line width=3 color RdYlGn
height = 3

[spacer]
height = 0.5
title = height=0.5

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style=flybase;fontsize=10
style =flybase
fontsize = 10

[spacer]
height = 1

[test gene rows]
file = dm3_genes.bed.gz
height = 3
title = max num rows=3; style=UCSC
fontsize = 8
style = UCSC
gene_rows = 3

[spacer]
height = 1

[test bed6]
file = dm3_genes.bed6.gz
height = 7
title = bed6 border_color=black; gene_rows=10; fontsize=8; color=Reds (when a
        color map is used for the color (e.g. coolwarm, Reds) the bed score
        column mapped to a color)
fontsize = 7
file_type = bed
color = Reds
border color = black
gene rows = 10

[test bed6]
file = dm3_genes.bed6.gz
height = 10
title = bed6 fontsize=10; line_width=1.5; global_max_row=yes
        (global max row sets the number of genes per row as the maximum found
        anywhere in the genome, hence the white space at the bottom)
fontsize = 10
file_type = bed
global_max_row = yes
interval_height = 200
line_width = 1.5

[x-axis]
fontsize=30
title = fontsize=30

[vlines]
file = tad_classification.bed
type = vlines

"""
with open(ROOT + "browser_tracks.ini", 'w') as fh:
    fh.write(browser_tracks)


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


def test_plot_tracks_existing_chr_empty_tracks():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0}/browser_tracks.ini --region X:0-1000000 --trackLabelFraction 0.2 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_plot_2.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_missing_chr():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0}/browser_tracks.ini --region Y:0-1000000 --trackLabelFraction 0.2 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_plot_3.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
