# -*- coding: utf-8 -*-
import matplotlib as mpl
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
from pygenometracks.utilities import InputError
mpl.use('agg')

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[arcs]
title = default
file = short_long.arcs
color = bwr
height = 5

[spacer]

[arcs]
file = short_long.arcs
color = bwr
height = 5
title = ylim = 6000000 (6Mb)
ylim = 6000000

[spacer]

[arcs]
file = short_long.arcs
color = bwr
height = 5
title = ylim = 200000 (200kb)
ylim = 200000

[spacer]

[arcs]
title = compact_arcs_level = 1
file = short_long.arcs
color = bwr
height = 5
compact_arcs_level = 1

[spacer]

[arcs]
title = compact_arcs_level = 1 ylim = 6000000 (6Mb)
ylim = 6000000
file = short_long.arcs
color = bwr
height = 5
compact_arcs_level = 1

[spacer]

[arcs]
title = compact_arcs_level = 2
file = short_long.arcs
color = bwr
height = 5
compact_arcs_level = 2

[spacer]

[arcs]
file = short_long.arcs
color = bwr
height = 5
title = ylim = 6000000 (6Mb) links_type = triangles
ylim = 6000000
links_type = triangles

[spacer]

[arcs]
file = short_long.arcs
color = bwr
height = 5
title = ylim = 6000000 (6Mb) links_type = triangles compact_arcs_level = 1
compact_arcs_level = 1
ylim = 6000000
links_type = triangles

[spacer]

[arcs]
file = short_long.arcs
color = bwr
height = 5
title = links_type = triangles compact_arcs_level = 2
compact_arcs_level = 2
links_type = triangles

[x-axis]
where = bottom
"""
with open(os.path.join(ROOT, "short_long_arcs.ini"), 'w') as fh:
    fh.write(browser_tracks)

with open(os.path.join(ROOT, "short_long_arcs_incorrect.ini"), 'w') as fh:
    fh.write(browser_tracks.replace('compact_arcs_level = 2', 'compact_arcs_level = 2\nylim=1000000'))

browser_tracks = """
[x-axis]
where = top

[spacer]
height = 0.05
[hic matrix]
file = Li_et_al_2015.cool
title = depth = 300000; transform = log1p; min_value = 5 (next track: overlay_previous = share-y links_type = loops)
depth = 300000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false

[test arcs overlay]
file = test_wide.arcs
color = red
line_width = 5
links_type = loops
overlay_previous = share-y

[test arcs]
file = test_wide.arcs
line_width = 3
color = RdYlGn
title = links line_width = 3 color RdYlGn
height = 3
orientation = inverted

[spacer]
height = 1

[hic matrix]
file = Li_et_al_2015.cool
title = depth = 300000; transform = log1p; min_value = 5 (next track: overlay_previous = share-y links_type = loops)
depth = 300000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false

[test arcs overlay]
file = test_wide.arcs
color = red
line_width = 5
links_type = loops
overlay_previous = share-y

[test arcs]
file = test_wide.arcs
line_width = 3
color = RdYlGn
title = links line_width = 3 color RdYlGn use_middle = true
use_middle = true
height = 3
orientation = inverted
"""
with open(os.path.join(ROOT, "arcs_use_middle.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[arcs]
title = loop with scores
file = test_high_score.arcs
color = cividis
height = 2
links_type = loops

[arcs]
title = arcs with scores
file = test_high_score.arcs
color = cividis
height = 2
min_value = 0
max_value = 80

[arcs]
title = arcs without scores
file = test_noscore.arcs
color = blue
line_width = 0.5
height = 2

"""
with open(os.path.join(ROOT, "arcs_no_score.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[arcs]
title = loop with scores
file = test_high_score.arcs
color = cividis
height = 2
links_type = loops

[arcs]
title = arcs with scores
file = test_high_score.arcs
color = cividis
height = 2
min_value = 0
max_value = 80

[arcs]
title = arcs without scores
file = test_noscore.arcs
color = cividis
height = 2

"""
with open(os.path.join(ROOT, "arcs_no_score_incorrect.ini"), 'w') as fh:
    fh.write(browser_tracks)

with open(os.path.join(ROOT, "arcs_no_score_invalid_score.ini"), 'w') as fh:
    fh.write(browser_tracks.replace('test_noscore.arcs',
                                    'arcs_invalid_score.arcs'))

with open(os.path.join(ROOT, "arcs_no_score_invalid_score2.ini"), 'w') as fh:
    fh.write(browser_tracks.replace('test_noscore.arcs',
                                    'arcs_invalid_score2.arcs'))

for suf in ['', '2']:
    browser_tracks = f"""
[arcs]
file = arcs_invalid{suf}.arcs
"""
    with open(os.path.join(ROOT, f"arcs_invalid{suf}.ini"), 'w') as fh:
        fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_short_long_arcs():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    for suf in ['', '_incorrect']:
        ini_file = os.path.join(ROOT, f"short_long_arcs{suf}.ini")
        region = "chr11:40000000-46000000"
        expected_file = os.path.join(ROOT, 'master_short_long_arcs.png')
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        
        # Remove incorrect ini file
        if 'incorrect' in ini_file:
            os.remove(ini_file)


def test_use_middle_arcs():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "arcs_use_middle.ini")
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_arcs_use_middle.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_arcs_no_score():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    for suf in ['', '_incorrect', '_invalid_score', '_invalid_score2']:
        ini_file = os.path.join(ROOT, f"arcs_no_score{suf}.ini")
        region = "X:3000000-3300000"
        expected_file = os.path.join(ROOT, 'master_arcs_no_score.png')
        args = f"--tracks {ini_file} --region {region} "\
            "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
            f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        # Remove the incorrect ini files or using invalid files
        if 'incorrect' in ini_file or 'invalid' in ini_file:
            os.remove(ini_file)


def test_arcs_invalid():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=True)
    ini_file = os.path.join(ROOT, "arcs_invalid.ini")
    region = "X:3000000-3300000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'not enough values to unpack (expected 6, got 5)' in str(e)
    else:
        raise Exception("The arcs_invalid should fail.")
    os.remove(ini_file)


def test_arcs_invalid2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=True)
    ini_file = os.path.join(ROOT, "arcs_invalid2.ini")
    region = "X:3000000-3300000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'One of the fields is not an integer.' in str(e)
    else:
        raise Exception("The arcs_invalid2 should fail.")
    os.remove(ini_file)
