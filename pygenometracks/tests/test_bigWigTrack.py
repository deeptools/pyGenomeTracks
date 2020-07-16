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
min_value = auto
max_value = auto

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

tracks = """
[test hlines]
color = red
line_width = 2
line_style = dashed
y_values = 10, 200
min_value = 0
show_data_range = true
height = 5
title = hlines: color = red; line_width = 2; line_style = dashed; y_values = 10, 200
file_type = hlines

[spacer]

[test bigwig fill]
file = bigwig2_X_2.5e6_3.5e6.bw
color = gray
height = 2
type = fill
title = bigwig: gray fill overlayed with hlines at 10 and 200 blue dotted
max_value = 50

[test hlines ovelayed]
color = blue
line_style = dotted
y_values = 10, 200
overlay_previous = share-y
file_type = hlines

[spacer]

[x-axis]
"""

with open(os.path.join(ROOT, "hlines.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bigwig1]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
height = 4
title = first bw
min_value = 0
max_value = 30

[test bigwig2]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title = second bw
min_value = 0
max_value = 30
orientation = inverted

[spacer]
height = 0.5

[test bigwig dif]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
negative_color = red
height = 8
title = operation = file - second_file
operation = file - second_file
min_value = -30
max_value = 30
nans_to_zeros = true

[spacer]
height = 0.5

[test bigwig op]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
negative_color = red
height = 8
title = operation = log10((1 + file)/(1 + second_file))
operation = log10((1 + file)/(1 + second_file))
nans_to_zeros = true

[spacer]
height = 0.5

[test bigwig op2]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title = operation = 2 + second_file
operation = 2 + second_file
nans_to_zeros = true
max_value = 32
min_value = 0

[test bigwig op2bis]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
height = 4
title = operation = 1 + 2 * file
operation = 1 + 2 * file
nans_to_zeros = true
max_value = 32
min_value = 0

[spacer]
height = 0.5

[test bigwig op3]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = green
height = 4
title = operation = max(file, second_file) in green overlayed with (file + second_file) / 2 in lime overlayed with min(file, second_file) in yellow
operation = max(file, second_file)
nans_to_zeros = true
max_value = 30
min_value = 0

[test bigwig op4]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = lime
operation = (file + second_file) / 2
nans_to_zeros = true
overlay_previous = share-y

[test bigwig op5]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig2_X_2.5e6_3.5e6.bw
color = yellow
operation = min(file, second_file)
nans_to_zeros = true
overlay_previous = share-y

[spacer]
height = 0.5


[x-axis]
"""

with open(os.path.join(ROOT, "operation.ini"), 'w') as fh:
    fh.write(tracks)

with open(os.path.join(ROOT, "incorrect_operation.ini"), 'w') as fh:
    fh.write(tracks.replace('operation = min(file, second_file)\n', 'operation = min(file, second_file)\ny_axis_values = original\n'))

tracks = """
[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
height = 7
title = No alpha:
        (bigwig color=blue 2000 bins) overlaid with (bigwig color = (0.6, 0, 0) max over 300 bins) overlaid with (bigwig mean color = green 200 bins)
number_of_bins = 2000
min_value = 0
max_value = 30

[test bigwig max]
file = bigwig2_X_2.5e6_3.5e6.bw
color = (0.6, 0, 0)
summary_method = max
number_of_bins = 300
overlay_previous = share-y

[test bigwig mean]
file = bigwig2_X_2.5e6_3.5e6.bw
color = green
type = fill
number_of_bins = 200
overlay_previous = share-y

[spacer]

[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
height = 7
title = alpha
        (bigwig color = blue 2000 bins) overlaid with (bigwig color = (0.6, 0, 0) alpha = 0.5 max over 300 bins) overlaid with (bigwig mean color = green alpha = 0.5 200 bins)
number_of_bins = 2000
min_value = 0
max_value = 30

[test bigwig max]
file = bigwig2_X_2.5e6_3.5e6.bw
color = (0.6, 0, 0)
alpha = 0.5
summary_method = max
number_of_bins = 300
overlay_previous = share-y

[test bigwig mean]
file = bigwig2_X_2.5e6_3.5e6.bw
color = green
alpha = 0.5
type = fill
number_of_bins = 200
overlay_previous = share-y

[spacer]

[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
height = 7
title = alpha for lines/points:
        (bigwig color=(0.6, 0, 0) alpha = 0.5 max) overlaid with (bigwig mean color = green alpha = 0.5 line:2) overlaid with (bigwig min color = blue alpha = 0.5 points:2)
color = (0.6, 0, 0)
alpha = 0.5
summary_method = max
number_of_bins = 300
min_value = 0
max_value = 30

[test bigwig mean]
file = bigwig2_X_2.5e6_3.5e6.bw
color = green
type = line:2
alpha = 0.5
summary_method = mean
number_of_bins = 300
overlay_previous = share-y

[test bigwig min]
file = bigwig2_X_2.5e6_3.5e6.bw
color = blue
summary_method = min
number_of_bins = 1000
type = points:3
alpha = 0.5
overlay_previous = share-y

[x-axis]
"""

with open(os.path.join(ROOT, "alpha.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bigwig fill]
file = bigwig2_X_2.5e6_3.5e6.bw
color = black
height = 2
type = fill
title = bigwig: black fill (height = 2)

[test bigwig fill with grid]
file = bigwig2_X_2.5e6_3.5e6.bw
color = black
height = 2
type = fill
grid = true
title = bigwig: black fill with grid (height = 2)

[spacer]

[test bigwig fill with grid]
file = bigwig2_X_2.5e6_3.5e6.bw
color = black
height = 5
type = fill
grid = true
title = bigwig: black fill with grid (height = 5)

[spacer]

[test bigwig fill with grid]
file = bigwig2_X_2.5e6_3.5e6.bw
color = black
height = 5
type = fill
grid = true
max_value = 50
title = bigwig: black fill with grid (height = 5 max_value = 50)

[spacer]

[test bigwig fill with grid]
file = bigwig2_X_2.5e6_3.5e6.bw
color = black
height = 15
type = fill
grid = true
max_value = 50
title = bigwig: black fill with grid (height = 15 max_value = 50)

[x-axis]
"""

with open(os.path.join(ROOT, "grid.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[bigwig file test]
file = bigwig_chrx_2e6_5e6.bw
# height of the track in cm (optional value)
height = 4
title = bigwig
min_value = 0
max_value = 30
"""

with open(os.path.join(ROOT, "example_bigwig.ini"), 'w') as fh:
    fh.write(tracks)

with open(os.path.join(ROOT, "example_bigwig_invalid_custom_color.ini"), 'w') as fh:
    fh.write(tracks + 'color = (a')

with open(os.path.join(ROOT, "example_bigwig_invalid_custom_color2.ini"), 'w') as fh:
    fh.write(tracks + 'color = (1)')

with open(os.path.join(ROOT, "example_bigwig_invalid_custom_color3.ini"), 'w') as fh:
    fh.write(tracks + 'color = Reds')

with open(os.path.join(ROOT, "example_bigwig_invalid_transform.ini"), 'w') as fh:
    fh.write(tracks + 'transform = myfunction')

tracks = """
[bigwig op test]
file = bigwig3_X_2.5e6_3.5e6.bw
second_file = bigwig_chrx_2e6_5e6.bw
# height of the track in cm (optional value)
operation = file - second_file
height = 4
title = file - second_file
"""

with open(os.path.join(ROOT, "example_op.ini"), 'w') as fh:
    fh.write(tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_bigwig_track():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bigwig.ini")
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


def test_alpha():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_alpha_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "alpha.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_alpha.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_hlines():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_hlines_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "hlines.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_hlines.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_grid():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_grid_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "grid.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_grid.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_op():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_op_test_',
                                 delete=False)
    for pref in ['', 'incorrect_']:
        ini_file = os.path.join(ROOT, f"{pref}operation.ini")
        region = "X:2700000-3100000"
        expected_file = os.path.join(ROOT, 'master_operation.png')
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        if 'incorrect' in ini_file:
            os.remove(ini_file)


def test_op_fakeChr():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_op_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "operation.ini")
    region = "fakeChr:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_operation_fakeChr.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_defaults():
    region = "X:2,500,000-3,000,000"
    for suf in [''] + ['_invalid_custom_color' + s for s in ['', '2', '3']] + \
            ['_invalid_transform']:
        outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_test_',
                                     delete=False)
        ini_file = os.path.join(ROOT, f"example_bigwig{suf}.ini")
        expected_file = os.path.join(ROOT, 'master_example_bigwig.png')
        args = f"--tracks {ini_file} --region {region} "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        if 'invalid' in ini_file:
            os.remove(ini_file)


def test_op_chr_in_only_one_bw():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bigwig_op_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "example_op.ini")
    region = "2L:0-1000"
    expected_file = os.path.join(ROOT, 'master_operation_2L.png')
    args = f"--tracks {ini_file} --region {region} "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
