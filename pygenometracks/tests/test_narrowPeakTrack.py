# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

tracks = """
[narrow]
file = test.narrowPeak
height = 4
max_value = 40
title = max_value=40

[narrow 2]
file = test.narrowPeak
height = 2
show labels = no
show data range =  no
color = #00FF0080
use summit = no
title = show labels=no; show data range=no; use summit=no;color=#00FF0080
[spacer]

[narrow 3]
file = test.narrowPeak
height = 2
show labels = no
color = #0000FF80
use summit = no
width adjust = 4
title = show labels=no;width adjust=3

[spacer]

[narrow 4]
file = test.narrowPeak
height = 3
type = box
color = blue
title = type=box;color=blue;

[x-axis]
"""

with open(ROOT + "narrow_peak.ini", 'w') as fh:
    fh.write(tracks)

with open(ROOT + "narrow_peak2.ini", 'w') as fh:
    fh.write(tracks.replace('test.narrowPeak', 'test2.narrowPeak'))

tolerance = 13  # default matplotlib pixed difference tolerance


def test_narrow_track():
    region = "X:2760000-2802000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='narrowTrack_test_', delete=False)
    args = "--tracks {root}/narrow_peak.ini --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName  {outfile}".format(root=ROOT, outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_narrowPeak.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_narrow_track_2():
    region = "X:2760000-2802000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='narrowTrack2_test_', delete=False)
    args = "--tracks {root}/narrow_peak2.ini --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName  {outfile}".format(root=ROOT, outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(ROOT + '/master_narrowPeak2.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
