import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[test with categories]
file = epilog.qcat.bgz
height = 5
title = height = 5; categories_file = epilog_cats.json
categories_file = epilog_cats.json

[spacer]

[test bedgraph tabix]
file = epilog.qcat.bgz
height = 5
title = height = 5

[test bedgraph tabix]
categories_file = epilog_cats.json
file = epilog.qcat.bgz
height = 5
title = height = 5
orientation = inverted

[spacer]
height = 0.5

[x-axis]
"""

with open(os.path.join(ROOT, "epilogos.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_epilogos_track():
    region = "X:3100000-3150000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "epilogos.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_epilogos.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
