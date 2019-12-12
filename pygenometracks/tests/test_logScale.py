# -*- coding: utf-8 -*-
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[test bigwig log]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = log1p
title = bigwig transform = log1p

[spacer]

[test bigwig log]
file = bigwig_chrx_2e6_5e6.bw
color = red
min_value = 0
height = 5
transform = log1p
title = bigwig transform = log1p min_value = 0

[x-axis]
"""

with open(os.path.join(ROOT, "log1p.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue
transform = log

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = red
height = 5
title = bedgraph color = red orientation = inverted
transform = log
orientation = inverted

[x-axis]
"""

with open(os.path.join(ROOT, "log.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_log1p_track():
    region = "X:2700000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='log1p_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "log1p.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_log1p.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)

def test_log_tracks():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_log_', delete=False)
    args = "--tracks {0} --region chr2:73,800,000-75,744,000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, "log.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_log.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
