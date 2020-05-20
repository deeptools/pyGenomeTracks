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


tolerance = 13  # default matplotlib pixed difference tolerance


def test_short_long_arcs():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region chr11:40000000-46000000 --trackLabelFraction 0.2" \
           " --width 38 --dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, "short_long_arcs.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print(f"saving test to {outfile.name}")
    res = compare_images(os.path.join(ROOT, 'master_short_long_arcs.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
