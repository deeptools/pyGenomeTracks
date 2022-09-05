# -*- coding: utf-8 -*-
import unittest
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

tracks = """
[test bigwig]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = no
title = bigwig transform = no

[spacer]

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
title = bigwig transform = log1p min_value = 0 y_axis_values = original
y_axis_values = original

[x-axis]
"""

with open(os.path.join(ROOT, "log1p.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bigwig]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = no
title = bigwig transform = no
[test bigwig]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = no
orientation = inverted
grid = true
title = bigwig transform = no orientation = inverted grid = true


[spacer]

[test bigwig log]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = log1p
title = bigwig transform = log1p
[test bigwig log]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = log1p
orientation = inverted
grid = true
title = bigwig transform = log1p orientation = inverted grid = true

[spacer]

[test bigwig log]
file = bigwig_chrx_2e6_5e6.bw
color = red
min_value = 0
height = 5
transform = log1p
title = bigwig transform = log1p min_value = 0 y_axis_values = original
y_axis_values = original

[test bigwig log]
file = bigwig_chrx_2e6_5e6.bw
color = red
min_value = 0
height = 5
transform = log1p
orientation = inverted
grid = true
title = bigwig transform = log1p min_value = 0 y_axis_values = original orientation = inverted grid = true
y_axis_values = original

[x-axis]
"""

with open(os.path.join(ROOT, "log1p_grid.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bigwig]
file = bigwig_chrx_2e6_5e6.bw
color = red
height = 5
transform = log
title = bigwig transform = log
"""

with open(os.path.join(ROOT, "log_neg.ini"), 'w') as fh:
    fh.write(tracks)
with open(os.path.join(ROOT, "mlog_neg.ini"), 'w') as fh:
    fh.write(tracks.replace('log', '-log'))


tracks = """
[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue transform = no
transform = no

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue transform = log
transform = log

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = red
height = 5
title = bedgraph color = red transform = log min_value = 1
min_value = 1
transform = log

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = green
height = 5
title = bedgraph color = green transform = log log_pseudocount = 2 min_value = 0
transform = log
log_pseudocount = 2
min_value = 0

[test bedgraph with operation]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = green
height = 5
title = bedgraph color = green operation = log(2+file) min_value = 0.7
operation = log(2+file)
min_value = 0.7

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = black
height = 5
title = bedgraph color = black transform = log2 log_pseudocount = 1 min_value = 0
transform = log2
log_pseudocount = 1
min_value = 0

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = black
height = 5
title = bedgraph color = black operation = log2(1+file) min_value = 0
operation = log2(1+file)
min_value = 0

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = black
height = 5
title = bedgraph color = black transform = log2 log_pseudocount = 1 min_value = 0 y_axis_values = original
transform = log2
log_pseudocount = 1
min_value = 0
y_axis_values = original

[x-axis]
"""

with open(os.path.join(ROOT, "log.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
grid = true
title = bedgraph color = blue transform = no grid=true
transform = no

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
grid = true
title = bedgraph color = blue transform = log grid=true
transform = log

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = red
height = 5
grid = true
title = bedgraph color = red transform = log min_value = 1 grid=true
min_value = 1
transform = log

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = green
height = 5
grid = true
title = bedgraph color = green transform = log log_pseudocount = 2 min_value = 0 grid=true
transform = log
log_pseudocount = 2
min_value = 0

[test bedgraph with operation]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = green
height = 5
grid = true
title = bedgraph color = green operation = log(2+file) min_value = 0.7 grid=true
operation = log(2+file)
min_value = 0.7

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = black
height = 5
grid = true
title = bedgraph color = black transform = log2 log_pseudocount = 1 min_value = 0 grid=true
transform = log2
log_pseudocount = 1
min_value = 0

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = black
height = 5
grid = true
title = bedgraph color = black operation = log2(1+file) min_value = 0 grid=true
operation = log2(1+file)
min_value = 0

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = black
height = 5
grid = true
title = bedgraph color = black transform = log2 log_pseudocount = 1 min_value = 0 y_axis_values = original grid=true
transform = log2
log_pseudocount = 1
min_value = 0
y_axis_values = original

[x-axis]
"""

with open(os.path.join(ROOT, "log_grid.ini"), 'w') as fh:
    fh.write(tracks)

tracks = """
[x-axis]

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue transform = log
transform = log

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue transform = log y_axis_values = original
transform = log
y_axis_values = original

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue transform = log10 y_axis_values = original
transform = log10
y_axis_values = original

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph color = blue transform = -log y_axis_values = original
transform = -log
y_axis_values = original
"""

with open(os.path.join(ROOT, "log_more.ini"), 'w') as fh:
    fh.write(tracks)

with open(os.path.join(ROOT, "log_more_incorrect.ini"), 'w') as fh:
    fh.write(tracks + 'type = invalid:a\n')

tolerance = 13  # default matplotlib pixed difference tolerance


def test_log1p_track():
    outfile = NamedTemporaryFile(suffix='.png', prefix='log1p_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "log1p.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_log1p.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_log1p_grid():
    outfile = NamedTemporaryFile(suffix='.png', prefix='log1p_grid_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "log1p_grid.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_log1p_grid.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_log_tracks():

    outfile = NamedTemporaryFile(suffix='.png',
                                 prefix='pyGenomeTracks_test_log_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "log.ini")
    region = "chr2:73,800,000-75,744,000"
    expected_file = os.path.join(ROOT, 'master_log.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_log_grid():

    outfile = NamedTemporaryFile(suffix='.png',
                                 prefix='pyGenomeTracks_test_loggrid_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "log_grid.ini")
    region = "chr2:73,800,000-75,744,000"
    expected_file = os.path.join(ROOT, 'master_log_grid.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


class TestLogNegMethods(unittest.TestCase):
    def test_log_tracks_with_0values(self):
        outfile = NamedTemporaryFile(suffix='.png', prefix='log_test_',
                                     delete=False)
        for pref in ['', 'm']:
            ini_file = os.path.join(ROOT, f"{pref}log_neg.ini")
            region = "X:2700000-3100000"
            args = f"--tracks {ini_file} --region {region} "\
                   "--trackLabelFraction 0.2 --dpi 130 "\
                   f"--outFileName {outfile.name}".split()
            with self.assertRaises(Exception) as context:
                pygenometracks.plotTracks.main(args)

            assert "coverage contains values smaller or equal to" in str(context.exception)
            os.remove(ini_file)


def test_log_more():

    outfile = NamedTemporaryFile(suffix='.png',
                                 prefix='pyGenomeTracks_test_log_',
                                 delete=False)
    for suf in ['', '_incorrect']:
        ini_file = os.path.join(ROOT, f"log_more{suf}.ini")
        region = "chr2:73,800,000-75,744,000"
        expected_file = os.path.join(ROOT, 'master_log_more.png')
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
    # I remove the incorrect ini file
    os.remove(ini_file)
