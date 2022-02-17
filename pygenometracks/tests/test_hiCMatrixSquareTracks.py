import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks_with_hic = """
[hic]
file = Li_et_al_2015.h5

[hic_matrix_square]
file = Li_et_al_2015.h5
title = no region2 on h5 file
file_type = hic_matrix_square

[hic_matrix_square2]
file = Li_et_al_2015.cool
region2 = X:3000000-3500000
title = region2 = X:3000000-3500000 on cool file
file_type = hic_matrix_square

[hic_matrix_square3]
file = Li_et_al_2015.h5
region2 = chrX:3000000-3500000
title = region2 = chrX:3000000-3500000 on h5 file forced height
height = 2
file_type = hic_matrix_square

[hic_matrix_square4]
file = Li_et_al_2015.cool
transform = log1p
show_masked_bins = true
orientation = inverted
title = cool file log1p show_masked_bins inverted
file_type = hic_matrix_square

[x-axis]
"""

with open(os.path.join(ROOT, "hic_square.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_hic_square():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "hic_square.ini")
    region = "X:2500000-3500000"
    expected_file = os.path.join(ROOT, 'master_hic_square.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_hic_square_dec():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "hic_square.ini")
    region = "X:2500000-3500000"
    expected_file = os.path.join(ROOT, 'master_hic_square_dec.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 --dpi 130 "\
           "--decreasingXAxis "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
