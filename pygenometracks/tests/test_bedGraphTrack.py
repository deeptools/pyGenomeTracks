# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

tracks = """
[test bedgraph lines]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = gray
height = 2
type = line
title = orientation=inverted; show data range=no
orientation = inverted
show data range = no
max_value = 50

[test bedgraph lines:0.2]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = red
height = 2
type = line:0.2
title = type=line:0.2

[spacer]

[test bedgraph points]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = black
height = 2
min_value = -15
max_value = 100
type = points:0.5
title = type=point:0.5; min_value=0;max_value=100

[spacer]

[test bedgraph nans to zeros]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = red
height = 2
nans to zeros = True
title = nans to zeros =True

[spacer]


[spacer]

[test bedgraph negative color]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = red
negative color = blue
max_value = 15
min_value = -15
#type=line
height = 3
title = negative color = blue


[test bedgraph negative color line]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = red
negative color = blue
max_value = 15
min_value = -15
type=line
height = 3
title = negative color = blue type=line
type = line

[test bedgraph negative color points]
file = bedgraph_chrx_2e6_5e6_with_negative.bg
color = red
negative color = black
max_value = 15
min_value = -15
type=line
height = 3
title = negative color = blue type=line
type = points

[spacer]

[x-axis]
"""

with open(ROOT + "bedgraph.ini", 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixel difference tolerance


def test_narrow_track():
    region = "X:2700000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_', delete=False)
    args = "--tracks {root}/bedgraph.ini --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName  {outfile}".format(root=ROOT, outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_bedgraph.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
