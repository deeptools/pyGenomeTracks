# -*- coding: utf-8 -*-
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

with open(os.path.join(ROOT, "bed_and_gtf_tracks.ini"), 'r') as f:
    with open(os.path.join(ROOT, "bed_and_gtf_tracks_dep.ini"), 'w') as fo:
        flag = False
        for line in f:
            if 'dm3_subset_BDGP5.78.gtf.gz' in line:
                flag = True
            if 'file_type' in line and flag:
                line = 'file_type = bed\n'
                flag = False
            fo.write(line)

browser_tracks = """
[test bigwig lines]
file = bigwig2_X_2.5e6_3.5e6.bw
color = gray
height = 2
type = line
#title = orientation=inverted; show data range=no
title = orientation = inverted; show_data_range = false
orientation = inverted
show data range = no
max_value = 50

[test bigwig lines:0.2]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 2
type = line:0.2
title = type=line:0.2

[spacer]

[test bigwig points]
file = bigwig_chrx_2e6_5e6.bw
color = black
height = 2
min_value = 0
max_value = 100
type = points:0.5
# title = type=point:0.5; min_value=0;max_value=100
title = type = point:0.5; min_value = 0; max_value = 100

[spacer]

[test bigwig nans to zeros]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 2
nans to zeros = True
# title = nans to zeros =True
title = nans_to_zeros = true

[spacer]

[test bigwig mean]
file = bigwig2_X_2.5e6_3.5e6.bw
color = gray
height = 5
# title = gray:summary method=mean; blue:summary method=max; red:summary method=min
title = gray:summary_method = mean; blue:summary_method = max;
        red:summary_method = min
type = line
summary method = mean
max_value = 150
min_value = -5
show data range = no
number of bins = 300

[test bigwig max]
file = bigwig2_X_2.5e6_3.5e6.bw
#title = test
color = blue
type = line
summary method = max
max_value = 150
min_value = -15
show data range = no
overlay previous = share-y
number of bins = 300

[test bigwig min]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
type=line
summary method = min
max_value = 150
min_value = -25
overlay previous = share-y
number of bins = 300

[spacer]

[x-axis]
"""
with open(os.path.join(ROOT, "bigwig_dep.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks_with_hic = """
[hic matrix]
file = small_test2.cool
title = cool with few interactions show_masked_bins = false (default)
depth = 200000
file_type = hic_matrix
boudaries_file = tad_classification.bed
height = 5

[spacer]

[hic matrix]
file = small_test2.cool
title = cool with few interactions show_masked_bins = true
depth = 200000
file_type = hic_matrix
height = 5
show_masked_bins = true

[x-axis]
"""
with open(os.path.join(ROOT, "browser_tracks_hic_small_test_boundaries_file.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)


tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_bed_and_gtf_dep():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'bed_and_gtf_tracks_dep.ini')
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_bed_and_gtf.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
    os.remove(ini_file)


def test_bigwig_track_dep():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bigwig_dep.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_bigwig.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
    os.remove(ini_file)


def test_plot_tracks_with_hic_small_file_boundaries():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    region = '1:0-200000'
    expected_file = os.path.join(ROOT, 'master_plot_hic_small_test.png')

    ini_file = os.path.join(ROOT, 'browser_tracks_hic_small_test_boundaries_file.ini')

    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.23 --width 38 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
    os.remove(ini_file)
