import os.path
import shutil
from tempfile import NamedTemporaryFile

import matplotlib as mpl
from get_matplotlib_CI_version import get_CI_mpl_version
from matplotlib.testing.compare import compare_images

import pygenometracks.plotTracks

mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[gwas]
file = gwas_1.gwas
height = 4
transform = no
title = test_1 transform = no orientation = inverted color = black
orientation = inverted
color = black

[spacer]

[gwas]
file = gwas_1.gwas
height = 4
title = test_1 default values

[spacer]

[gwas_2]
file = gwas_2.gwas
file_has_header = True
height = 2
title = test_2 file_has_header = true color = #50E3C2 border_color = red line_width = 2 marker_size = 90 show_data_range = false
color = #50E3C2
border_color = red
line_width = 2
marker_size = 90
show_data_range = false

[spacer]

[gwas_2]
file = gwas_1.gwas
height = 4
title = test_1 default values min_value = 1 max_value = 1e-15
min_value = 1
max_value = 1e-15

[spacer]

[gwas1]
file = gwas_1.gwas
height = 4
title = test_1 orange overlayed with test_2 blue
color = orange

[gwas2]
file = gwas_2.gwas
file_has_header = True
color = blue
overlay_previous = share-y

[x-axis]
"""

with open(os.path.join(ROOT, "gwas.ini"), 'w') as fh:
    fh.write(tracks)


tracks = """
[gwas]
file = head_all_hg38_qcd_LE1_simgwas_quant1a.simgwas_quant1.glm.linear
height = 4
file_has_header = True
grid = true

title = glm.linear grid = true

[spacer]

[gwas raster]
file = head_all_hg38_qcd_LE1_simgwas_quant1a.simgwas_quant1.glm.linear
height = 4
file_has_header = True
grid = true
rasterize = true
title = glm.linear grid = true rasterize = true

[x-axis]
"""

with open(os.path.join(ROOT, "gwas2.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance
default_mpl_version = get_CI_mpl_version()


def test_gwas_track():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 26
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas.ini")
    region = "X:3000000-3200000"
    expected_file = os.path.join(ROOT, 'master_gwas.png')
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.3 --plotWidth 30 --dpi 130 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_gwas_track_chrX():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 15
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas.ini")
    region = "chrX:3000000-3200000"
    expected_file = os.path.join(ROOT, 'master_gwas.png')
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.3 --plotWidth 30 --dpi 130 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance + 14)  # 14 corresponds to the 'chr' on the x axis
    assert res is None, res

    os.remove(outfile.name)


def test_gwas_track_chrY():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = 26
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.png', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas.ini")
    region = "chrY:3000000-3200000"
    expected_file = os.path.join(ROOT, 'master_gwas_chrY.png')
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.2 --dpi 130 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_gwas_track_raster():

    if mpl.__version__ != default_mpl_version:
        my_tolerance = tolerance
    else:
        my_tolerance = tolerance

    outfile = NamedTemporaryFile(suffix='.pdf', prefix='gwas_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gwas2.ini")
    region = "1:0-1300000"
    expected_file = os.path.join(ROOT, 'master_gwas2.pdf')
    # matplotlib compare on pdf will create a png next to it.
    # To avoid issues related to write in test_data folder
    # We copy the expected file into a temporary place
    new_expected_file = NamedTemporaryFile(suffix='.pdf',
                                           prefix='pyGenomeTracks_test_',
                                           delete=False)
    shutil.copy(expected_file, new_expected_file.name)
    expected_file = new_expected_file.name
    args = f"--tracks {ini_file} --region {region} " \
           "--trackLabelFraction 0.2 --dpi 10 " \
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, my_tolerance)
    assert res is None, res

    os.remove(outfile.name)
