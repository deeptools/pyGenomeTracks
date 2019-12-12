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
[x-axis]
title = A very long title that does not fit into a single line whatever is the width of the label
"""
with open(os.path.join(ROOT, "title.ini"), 'w') as fh:
    fh.write(browser_tracks)


tolerance = 13  # default matplotlib pixed difference tolerance


def test_regular_width_label():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region X:3000000-3500000 --trackLabelFraction 0.2" \
           " --width 38 --dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, "title.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_title_0.2.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region X:3000000-3500000 --trackLabelFraction 0.5" \
           " --width 38 --dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, "title.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_title_0.5.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_big_font():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region X:3000000-3500000 --trackLabelFraction 0.5" \
           " --width 38 --dpi 130 --outFileName {1} --fontSize 20" \
           "".format(os.path.join(ROOT, "title.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_title_0.5_fs20.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
