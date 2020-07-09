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

tolerance = 13  # default matplotlib pixed difference tolerance


def test_short_long_arcs():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "short_long_arcs.ini")
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
