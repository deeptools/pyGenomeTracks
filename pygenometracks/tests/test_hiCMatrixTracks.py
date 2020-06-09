import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log1p; min_value = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 250000; orientation = inverted; colormap = PuRd; min_value = 5;
        max_value = 70
min_value = 5
max_value = 70
depth = 250000
colormap = PuRd
file_type = hic_matrix
show_masked_bins = false
orientation = inverted

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 300000; transform = log1p; colormap Blues (TADs:
        overlay_previous = share-y; line_width = 1.5)
colormap = Blues
min_value = 10
max_value = 150
depth = 300000
transform = log1p
file_type = hic_matrix

[tads]
file = tad_classification.bed
#title = TADs color = none; border_color = black
file_type = domains
border_color = black
color = none
height = 5
line_width = 1.5
overlay_previous = share-y
show_data_range = false

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 250000; transform = log1p; colormap = bone_r (links: overlay_previous = share-y;
        links_type = triangles; color = darkred; line_style = dashed, bigwig: color = red)
colormap = bone_r
min_value = 15
max_value = 200
depth = 250000
transform = log1p
file_type = hic_matrix
show_masked_bins = false

[test arcs]
file = links2.links
title =
links_type = triangles
line_style = dashed
overlay_previous = share-y
line_width = 0.8
color = darkred
show_data_range = false


[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title =
overlay_previous = yes
min_value = 0
max_value = 50
show_data_range = false

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; show_masked_bins = true; colormap =
        ['blue', 'yellow', 'red']; max_value = 150
depth = 200000
colormap = ['blue', 'yellow', 'red']
max_value = 150
file_type = hic_matrix
show_masked_bins = true

[spacer]
height = 0.1

[x-axis]

"""

with open(os.path.join(ROOT, "browser_tracks_hic.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

with open(os.path.join(ROOT, "browser_tracks_cool.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic.replace(".h5", ".cool"))

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log1p; min_value = 5; height = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false
height = 5

[hic matrix 2]
file = Li_et_al_2015.h5
title = same but orientation=inverted; no height
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = false
orientation = inverted

[spacer]
height = 0.5

[hic matrix 3]
file = Li_et_al_2015.h5
title = same rasterize = false
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
rasterize = false
show_masked_bins = false

[x-axis]

"""

with open(os.path.join(ROOT, "browser_tracks_hic_rasterize_height.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

browser_tracks_with_hic = """
[hic matrix log]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log; show_masked_bins = true
depth = 200000
transform = log
min_value = -4
max_value = 4
file_type = hic_matrix
show_masked_bins = true

[spacer]

[hic matrix -log]
file = Li_et_al_2015.h5
title = depth = 200000; transform = -log; show_masked_bins = true
depth = 200000
transform = -log
min_value = -4
max_value = 4
file_type = hic_matrix
show_masked_bins = true

[x-axis]
"""

with open(os.path.join(ROOT, "browser_tracks_hic_log-log.ini"), 'w') as fh:
    fh.write(browser_tracks_with_hic)

browser_tracks_with_mcool = """
[x-axis]

[mcool1]
file = matrix.mcool::/0
depth = 1000000
file_type = hic_matrix

[mcool2]
file = matrix.mcool::/1
depth = 1000000
file_type = hic_matrix

[mcool3]
file = matrix.mcool::/2
depth = 1000000
file_type = hic_matrix

[mcool4]
file = matrix.mcool::/4
depth = 1000000
file_type = hic_matrix
"""

with open(os.path.join(ROOT, "mcool.ini"), 'w') as fh:
    fh.write(browser_tracks_with_mcool)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_with_hic():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, 'browser_tracks_hic.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_plot_hic.png'),
                         outfile.name, tolerance)
    print("saving test to {}".format(outfile.name))
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_dec():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName {1} --decreasingXAxis" \
           "".format(os.path.join(ROOT, 'browser_tracks_hic.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_plot_hic_dec.png'),
                         outfile.name, tolerance)
    print("saving test to {}".format(outfile.name))
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_cool_region():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, 'browser_tracks_cool.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_plot_hic.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_hic_logmlog():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT, 'browser_tracks_hic_log-log.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_plot_hic_log-log.png'),
                         outfile.name, tolerance)
    print("saving test to {}".format(outfile.name))
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_rasterize_height():

    outfile = NamedTemporaryFile(suffix='.pdf', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-2600000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 10 --outFileName {1}" \
           "".format(os.path.join(ROOT,
                                  'browser_tracks_hic_rasterize_height.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT,
                                      'master_plot_hic_rasterize_height.pdf'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_mcool():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName {1}" \
           "".format(os.path.join(ROOT,
                                  'mcool.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT,
                                      'master_mcool.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
