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
height = 3
title = tabix color=blue; type=fill
min_value = -5
type = fill

[test bedgraph tabix]
file = bedgraph_chrx_2e6_5e6.bg.bgz
color = brown
height = 3
title = tabix color=brown; type=line:0.5;nans to zeros=True
nans to zeros = True
min_value = -5
type = line:0.5

[test bedgraph]
file = bedgraph_chrx_2e6_5e6.bg
color = red
height = 3
title = color=red;type=points:0.5
min_value = -5
type = points:0.5

# [test bedgraph]
# file = bedgraph_chrx_2e6_5e6.bg
# color = red
# height = 3
# title = color=blue
# min_value = 0
#
# [test bedgraph matrix]
# file = tad_score.gz
# title = bedgraph matrix
# height = 8
# file_type = bedgraph_matrix
#
# [test bedgraph matrix lines]
# file = tad_score.bm.bgz
# title = type=lines
# height = 8
# file_type = bedgraph_matrix
# type = lines
#
# [test bedgraph matrix lines tabix]
# file = tad_score.bm.bgz
# title = type=lines (file type is tabix)
# height = 8
# file_type = bedgraph_matrix
# type = lines
#
# [test bedgraph matrix lines tabix 2]
# file = tad_score.bm.bgz
# title = type=lines (file type is tabix)
# height = 8
# file_type = bedgraph_matrix
# type = lines
# plot horizontal lines=True

[spacer]

[x-axis]
"""

with open(ROOT + "bedgraph.ini", 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_narrow_track():
    region = "X:3000000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_', delete=False)
    args = "--tracks {root}/bedgraph.ini --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName  {outfile}".format(root=ROOT, outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_bedgraph.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
