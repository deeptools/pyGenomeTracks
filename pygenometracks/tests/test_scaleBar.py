# -*- coding: utf-8 -*-
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[sb1]
file_type = scalebar
title = scalebar; where = top
where = top

[spacer]

[sb1b]
file_type = scalebar
title = scalebar; height = 3; where = top
height = 3
where = top

[spacer]

[sb5]
file_type = scalebar
fontsize = 6
title = scalebar; fontsize = 6; where = top
where = top

[spacer]

[sb2]
file_type = scalebar
title = scalebar

[spacer]

[sb3]
file_type = scalebar
where = right
title = scalebar; where = right

[spacer]

[sb4]
file_type = scalebar
where = bottom
title = scalebar; where = bottom

[spacer]

[sb5]
file_type = scalebar
where = right
x_center = 3200000
size = 100002
fontsize = 8
line_width = 2
color = red
alpha = 0.5
title = scalebar; where = right;x_center = 3200000;size = 100002;fontsize = 8;line_width =2;color = red;alpha - 0.5
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "scale_bar.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[sb1]
file_type = scalebar
title = scalebar where = right
where = right

[spacer]

[sb2]
file_type = scalebar
title = scalebar x_boundary1 = 3200000
x_boundary1 = 3200000

[spacer]

[sb3]
file_type = scalebar
title = scalebar x_boundary1 = 3200000 x_boundary2 = 3250000
x_boundary2 = 3250000
x_boundary1 = 3200000

[spacer]

[sb4]
file_type = scalebar
title = scalebar x_boundary1 = 3250000 x_boundary2 = 3200000
x_boundary1 = 3250000
x_boundary2 = 3200000

[spacer]

[sb5]
file_type = scalebar
title = scalebar x_boundary1 = 3200000 x_center = 3250000
x_center = 3250000
x_boundary1 = 3200000

[spacer]

[sb6]
file_type = scalebar
title = scalebar x_boundary1 = 3200000 size = 50000
size = 50000
x_boundary1 = 3200000

[spacer]

[sb7]
file_type = scalebar
title = scalebar x_boundary2 = 3200000 size = 50000
size = 50000
x_boundary2 = 3200000

[x-axis]
"""
with open(os.path.join(ROOT, "scale_bar_xboundary.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[sb]
file_type = scalebar
x_boundary1 = 2000000
x_boundary2 = 2500000
# x_center is not valid
x_center = 2000000
size = 500000
"""
with open(os.path.join(ROOT, "incompatible_param1.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[sb]
file_type = scalebar
x_boundary1 = 2000000
x_boundary2 = 2500000
x_center = 2250000
# size is not valid
size = 5000000
"""
with open(os.path.join(ROOT, "incompatible_param2.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[sb]
file_type = scalebar
x_boundary1 = 2000000
x_boundary2 = 2500000
# x_center is not valid
x_center = 2000000
"""
with open(os.path.join(ROOT, "incompatible_param3.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[sb]
file_type = scalebar
x_boundary1 = 2000000
x_boundary2 = 2500000
# size is not valid
size = 2000000
"""
with open(os.path.join(ROOT, "incompatible_param4.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[sb]
file_type = scalebar
x_center = 2000000
x_boundary2 = 2500000
# size is not valid
size = 2000000
"""
with open(os.path.join(ROOT, "incompatible_param5.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_scale_bar():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "scale_bar.ini")
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_scale_bar.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "scale_bar.ini")
    region = "X:3200000-3300000"
    expected_file = os.path.join(ROOT, 'master_scale_bar_zoom.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_xboundary():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "scale_bar_xboundary.ini")
    region = "X:3000000-3600000"
    expected_file = os.path.join(ROOT, 'master_scale_bar_xboundary.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_xboundary_dec():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "scale_bar_xboundary.ini")
    region = "X:3000000-3600000"
    expected_file = os.path.join(ROOT, 'master_scale_bar_xboundary_dec.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--decreasingXAxis "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_xboundary_outside():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "scale_bar_xboundary.ini")
    region = "X:2000000000-2500000000"
    expected_file = os.path.join(ROOT, 'master_scale_bar_xboundary_outside.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_xboundary_superzoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "scale_bar_xboundary.ini")
    region = "X:3199500-3201000"
    expected_file = os.path.join(ROOT, 'master_scale_bar_xboundary_superzoom.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_scale_bar_incompatible_param():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    region = "X:3000000-3300000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    for suf in ['1', '2', '3', '4', '5']:
        ini_file = os.path.join(ROOT, f"incompatible_param{suf}.ini")
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        try:
            pygenometracks.plotTracks.main(args)
        except Exception as e:
            assert 'is incompatible with' in str(e)
        else:
            raise Exception(f"incompatible_param{suf} should fail.")
        os.remove(ini_file)
