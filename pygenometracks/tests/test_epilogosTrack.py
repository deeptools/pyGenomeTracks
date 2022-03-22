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
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "epilogos.ini")
    region = "X:3100000-3150000"
    expected_file = os.path.join(ROOT, 'master_epilogos.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_epilogos_track_overlap_chr_end():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "epilogos.ini")
    region = "chrX:3490000-24000000"
    expected_file = os.path.join(ROOT, 'master_epilogos_overlap_end.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_epilogos_track_over_chr_end():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "epilogos.ini")
    region = "X:25000000-30000000"
    expected_file = os.path.join(ROOT, 'master_epilogos_over_end.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
