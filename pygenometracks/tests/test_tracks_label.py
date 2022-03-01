# -*- coding: utf-8 -*-
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile, TemporaryDirectory
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

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_ral():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_ral.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           "--trackLabelHAlign right "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_cal():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_cal.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           "--trackLabelHAlign center "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_cal_dpi250():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_cal_d250.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 250 "\
           "--trackLabelHAlign center "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_large_width_label_big_font():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.5_fs20.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
           "--fontSize 20 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fixed_height():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_force_height.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--height 10 --title force_height "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_non_existing_dir():

    outdir = TemporaryDirectory()
    output_file = os.path.join(outdir.name, "pGT_test", "test.png")
    ini_file = os.path.join(ROOT, "title.ini")
    region = "X:3000000-3500000"
    expected_file = os.path.join(ROOT, 'master_title_0.2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {output_file}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         output_file, tolerance)
    assert res is None, res

    outdir.cleanup()
