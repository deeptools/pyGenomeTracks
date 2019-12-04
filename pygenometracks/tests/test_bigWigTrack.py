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
[test bigwig lines]
file = bigwig2_X_2.5e6_3.5e6.bw
color = gray
height = 2
type = line
title = orientation = inverted; show_data_range = false
orientation = inverted
show_data_range = false
max_value = 50

[test bigwig lines:0.2]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 2
type = line:0.2
title = type = line:0.2

[spacer]

[test bigwig points]
file = bigwig_chrx_2e6_5e6.bw
color = black
height = 2
min_value = 0
max_value = 100
type = points:0.5
title = type = point:0.5; min_value = 0; max_value = 100

[spacer]

[test bigwig nans to zeros]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 2
nans_to_zeros = true
title = nans_to_zeros = true

[spacer]

[test bigwig mean]
file = bigwig2_X_2.5e6_3.5e6.bw
color = gray
height = 5
title = gray:summary_method = mean; blue:summary_method = max;
        red:summary_method = min
type = line
summary_method = mean
max_value = 150
min_value = -5
show_data_range = false
number_of_bins = 300

[test bigwig max]
file = bigwig2_X_2.5e6_3.5e6.bw
#title = test
color = blue
type = line
summary_method = max
max_value = 150
min_value = -15
show_data_range = false
overlay_previous = share-y
number_of_bins = 300

[test bigwig min]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
type = line
summary_method = min
max_value = 150
min_value = -25
overlay_previous = share-y
number_of_bins = 300

[spacer]

[x-axis]
"""

with open(os.path.join(ROOT, "bigwig.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_bigwig_track():
    region = "X:2700000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "bigwig.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_bigwig.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_alpha():
    region = "X:2700000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_alpha_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "alpha.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_alpha.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
