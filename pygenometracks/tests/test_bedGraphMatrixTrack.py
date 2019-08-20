# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

tracks = """

[test bedgraph tabix]
file = bedgraph_chrx_2e6_5e6.bg.bgz
color = blue
height = 2
title = tabix color=blue; type=fill
min_value = -5
type = fill

[spacer]
width = 0.01

[test bedgraph tabix]
file = bedgraph_chrx_2e6_5e6.bg.bgz
color = brown
height = 2
title = tabix color=brown; type=line:0.5;nans to zeros=True; orientation=inverted
nans to zeros = True
min_value = -5
max_value = 50
type = line:0.5
orientation = inverted

[spacer]
width = 0.01

[test bedgraph]
file = bedgraph_chrx_2e6_5e6.bg
color = red
height = 2
title = color=red;type=points:0.5
min_value = 0
max_value = 50
type = points:0.5

[spacer]

[test bedgraph matrix]
file = tad_separation_score.bm.gz
title = bedgraph matrix (file type is tabix)
height = 5
file_type = bedgraph_matrix

[spacer]

[test bedgraph matrix lines]
file = tad_separation_score.bm.gz
title = type=lines
height = 5
file_type = bedgraph_matrix
type = lines
pos score in bin = block

[spacer]
width = 0.01

[test bedgraph matrix lines]
file = tad_separation_score_with_gap.bm.bgz
title = type=lines; show data range=no (file type is tabix with a gap)
height = 5
file_type = bedgraph_matrix
type = lines
pos score in bin = center
show data range = no

[test bedgraph matrix lines tabix]
file = tad_separation_score.bm.gz
title = type=lines;plot horizontal lines=yes
height = 6
file_type = bedgraph_matrix
type = lines
plot horizontal lines=yes
[spacer]

[x-axis]
"""

with open(ROOT + "bedgraph.ini", 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_bedgraphmatrix_track():
    region = "X:2850000-3150000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_', delete=False)
    args = "--tracks {root}/bedgraph.ini --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName  {outfile}".format(root=ROOT, outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_bedgraph.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
