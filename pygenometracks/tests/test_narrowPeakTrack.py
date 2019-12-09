# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[narrow]
file = test.narrowPeak
height = 4
max_value = 40
title = max_value = 40

[narrow 2]
file = test.narrowPeak
height = 2
show_labels = false
show_data_range =  false
color = #00FF0080
use_summit = false
title = show_labels = false; show_data_range = false; use_summit = false; color = #00FF0080

[spacer]

[narrow 3]
file = test.narrowPeak
height = 2
show_labels = false
color = #0000FF80
use_summit = false
width_adjust = 4
title = show_labels = false; use_summit = false; width_adjust = 4

[spacer]

[narrow 4]
file = test.narrowPeak
height = 3
type = box
color = blue
title = type = box; color = blue;

[x-axis]
"""

with open(os.path.join(ROOT, "narrow_peak.ini"), 'w') as fh:
    fh.write(tracks)

with open(os.path.join(ROOT, "narrow_peak2.ini"), 'w') as fh:
    fh.write(tracks.replace('test.narrowPeak', 'test2.narrowPeak'))

tolerance = 13  # default matplotlib pixed difference tolerance


def test_narrow_track():
    region = "X:2760000-2802000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='narrowTrack_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, 'narrow_peak.ini'),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_narrowPeak.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_narrow_track_2():
    region = "X:2760000-2802000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='narrowTrack2_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, 'narrow_peak2.ini'),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_narrowPeak2.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
