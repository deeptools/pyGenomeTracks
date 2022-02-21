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
[test bedgraph tabix]
file = bedgraph_chrx_2e6_5e6.bg.bgz
color = blue
height = 2
title = tabix color = blue; type = fill
min_value = -5
type = fill

[spacer]
height = 0.01

[test bedgraph tabix without tbi]
file = bedgraph_chrx_2e6_5e6_2.bg.bgz
color = brown
height = 2
title = tabix color = brown; type = line:0.5; nans_to_zeros = true; orientation = inverted
nans_to_zeros = true
min_value = -5
max_value = 50
type = line:0.5
orientation = inverted

[spacer]
height = 0.01

[test bedgraph]
file = bedgraph_chrx_2e6_5e6.bg
color = red
height = 2
title = color = red; type = points:0.5
min_value = 0
max_value = 50
type = points:0.5

[spacer]

[test bedgraph like bigwig]
file = bedgraph_chrx_2e6_5e6.bg
color = red
height = 2
title = color = red; type = points:0.5; summary_method = max
summary_method = max
min_value = 0
max_value = 50
type = points:0.5

[spacer]

[test bedgraph like bigwig_2]
file = bedgraph_chrx_2e6_5e6.bg
color = red
height = 2
title = color = red; type = points:0.5; summary_method = max; nbins = 3000
summary_method = max
number_of_bins = 3000
min_value = 0
max_value = 50
type = points:0.5

[spacer]

[test bedgraph matrix]
file = tad_separation_score.bm.gz
title = bedgraph matrix (file type is tabix) colormap = pink
height = 5
colormap = pink
file_type = bedgraph_matrix

[spacer]

[test bedgraph matrix]
file = tad_separation_score.bm.gz
title = bedgraph matrix (file type is tabix) colormap = red (this is not valid and viridis will be used)
height = 5
colormap = red
file_type = bedgraph_matrix

[spacer]

[test bedgraph matrix lines]
file = tad_separation_score.bm.gz
title = type = lines
height = 5
file_type = bedgraph_matrix
type = lines
pos_score_in_bin = block

[spacer]
height = 0.01

[test bedgraph matrix lines]
file = tad_separation_score_with_gap.bm.bgz
title = type = lines; show_data_range = false (file type is tabix with a gap)
height = 5
file_type = bedgraph_matrix
type = lines
pos_score_in_bin = center
show_data_range = false

[test bedgraph matrix lines tabix]
file = tad_separation_score.bm.gz
title = type = lines; plot_horizontal_lines = true
height = 6
file_type = bedgraph_matrix
type = lines
plot_horizontal_lines = true

[spacer]

[x-axis]
"""

with open(os.path.join(ROOT, "bedgraph.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bedgraph matrix]
file = tad_separation_score.bm.gz
type = lines
title = bedgraph matrix type = lines default colors
height = 5
file_type = bedgraph_matrix

[spacer]

[test bedgraph matrix]
file = tad_separation_score.bm.gz
type = lines
title = bedgraph matrix type = lines summary_color = red
summary_color = red
height = 5
file_type = bedgraph_matrix

[spacer]

[test bedgraph matrix]
file = tad_separation_score.bm.gz
type = lines
title = bedgraph matrix type = lines summary_color = red individual_color = (0, 1, 0)
summary_color = red
individual_color = (0, 1, 0)
height = 5
file_type = bedgraph_matrix

[spacer]

[test bedgraph matrix]
file = tad_separation_score.bm.gz
type = lines
title = bedgraph matrix type = lines individual_color = (.34, 1, 0)
individual_color = (.34, 1, 0)
height = 5
file_type = bedgraph_matrix

[x-axis]
"""

with open(os.path.join(ROOT, "bedgraph_matrix_line_colors.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_bedgraphmatrix_track():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bedgraph.ini")
    region = "X:2850000-3150000"
    expected_file = os.path.join(ROOT, 'master_bedgraph.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_bedgraphmatrix_track_chr1():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bedgraph.ini")
    region = "chr1:0-100000"
    expected_file = os.path.join(ROOT, 'master_bedgraph_chr1.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_bedgraph_matrix_line_colors():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bedgraph_matrix_line_colors.ini")
    region = "X:2850000-3150000"
    expected_file = os.path.join(ROOT, 'master_bedgraph_matrix_line_colors.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
