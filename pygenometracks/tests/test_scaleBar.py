# -*- coding: utf-8 -*-
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[sb1]
file_type = scalebar
title = scalebar; where = top
where = top

[spacer]

[sb1b]
file_type = scalebar
title = scalebar; height = 3; where = top
height = 3
where = top

[spacer]

[sb5]
file_type = scalebar
fontsize = 6
title = scalebar; fontsize = 6; where = top
where = top

[spacer]

[sb2]
file_type = scalebar
title = scalebar

[spacer]

[sb3]
file_type = scalebar
where = right
title = scalebar; where = right

[spacer]

[sb4]
file_type = scalebar
where = bottom
title = scalebar; where = bottom

[spacer]

[sb5]
file_type = scalebar
where = right
x_center = 3200000
size = 100002
fontsize = 8
line_width = 2
color = red
alpha = 0.5
title = scalebar; where = right;x_center = 3200000;size = 100002;fontsize = 8;line_width =2;color = red;alpha - 0.5
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "scale_bar.ini"), 'w') as fh:
    fh.write(browser_tracks)


tolerance = 13  # default matplotlib pixed difference tolerance


def test_scale_bar():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region X:3000000-3300000 --trackLabelFraction 0.2" \
           " --width 38 --dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, "scale_bar.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_scale_bar.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region X:3200000-3300000 --trackLabelFraction 0.2" \
           " --width 38 --dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, "scale_bar.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_scale_bar_zoom.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
