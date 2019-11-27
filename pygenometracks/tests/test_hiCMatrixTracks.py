import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log1p; min_value = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = no

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 250000; orientation = inverted; colormap = PrRd; min_value = 5;
        max_value = 70
min_value = 5
max_value = 70
depth = 250000
colormap = PuRd
file_type = hic_matrix
show_masked_bins = no
orientation = inverted

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 300000; transform = log1p; colormap Blues (TADs:
        overlay_previous = yes; line_width = 1.5)
colormap = Blues
min_value = 10
max_value = 150
depth = 300000
transform = log1p
file_type = hic_matrix

[tads]
file = tad_classification.bed
#title = TADs color = bed_rgb; border_color = black
file_type = domains
border_color = black
color = none
height = 5
line_width = 1.5
overlay_previous = share-y
show_data_range = no

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
show_masked_bins = no

[test arcs]
file = links2.links
title =
links_type = triangles
line_style = dashed
overlay_previous = share-y
line_width = 0.8
color = darkred
show_data_range = no


[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title =
overlay_previous = yes
min_value = 0
max_value = 50
show_data_range = no

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; show_masked_bins = yes; colormap =
        ['blue', 'yellow', 'red']; max_value = 150
depth = 200000
colormap = ['blue', 'yellow', 'red']
max_value = 150
file_type = hic_matrix
show_masked_bins = yes

[spacer]
height = 0.1

[x-axis]

"""

with open(ROOT + "browser_tracks_hic.ini", 'w') as fh:
    fh.write(browser_tracks_with_hic)


browser_tracks_with_cool = """
[hic matrix]
file = Li_et_al_2015.cool
title = depth = 200000; transform = log1p; min_value = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = no

[hic matrix]
file = Li_et_al_2015.cool
title = depth = 250000; orientation = inverted; colormap = PrRd; min_value = 5;
        max_value = 70
min_value = 5
max_value = 70
depth = 250000
colormap = PuRd
file_type = hic_matrix
show_masked_bins = no
orientation = inverted

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.cool
title = depth = 300000; transform = log1p; colormap Blues (TADs:
        overlay_previous = yes; line_width = 1.5)
colormap = Blues
min_value = 10
max_value = 150
depth = 300000
transform = log1p
file_type = hic_matrix

[tads]
file = tad_classification.bed
#title = TADs color = bed_rgb; border_color = black
file_type = domains
border_color = black
color = none
height = 5
line_width = 1.5
overlay_previous = share-y
show_data_range = no

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.cool
title = depth = 250000; transform = log1p; colormap = bone_r (links: overlay_previous = share-y;
        links_type = triangles; color = darkred; line_style = dashed, bigwig: color = red)
colormap = bone_r
min_value = 15
max_value = 200
depth = 250000
transform = log1p
file_type = hic_matrix
show_masked_bins = no

[test arcs]
file = links2.links
title =
links_type = triangles
line_style = dashed
overlay_previous = share-y
line_width = 0.8
color = darkred
show_data_range = no


[test bigwig]
file = bigwig2_X_2.5e6_3.5e6.bw
color = red
height = 4
title =
overlay_previous = yes
min_value = 0
max_value = 50
show_data_range = no

[spacer]
height = 0.5

[hic matrix]
file = Li_et_al_2015.cool
title = depth = 200000; show_masked_bins = yes; colormap =
        ['blue', 'yellow', 'red']; max_value = 150
depth = 200000
colormap = ['blue', 'yellow', 'red']
max_value = 150
file_type = hic_matrix
show_masked_bins = yes

[spacer]
height = 0.1

[x-axis]

"""

with open(ROOT + "browser_tracks_cool.ini", 'w') as fh:
    fh.write(browser_tracks_with_cool)

browser_tracks_with_hic = """
[hic matrix]
file = Li_et_al_2015.h5
title = depth = 200000; transform = log1p; min_value = 5; height = 5
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = no
height = 5

[hic matrix 2]
file = Li_et_al_2015.h5
title = same but orientation=inverted; no height
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
show_masked_bins = no
orientation = inverted

[spacer]
height = 0.5

[hic matrix 3]
file = Li_et_al_2015.h5
title = same rasterize = off
depth = 200000
min_value = 5
transform = log1p
file_type = hic_matrix
rasterize = no
show_masked_bins = no

[x-axis]

"""

with open(ROOT + "browser_tracks_hic_rasterize_height.ini", 'w') as fh:
    fh.write(browser_tracks_with_hic)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_with_hic():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0}/browser_tracks_hic.ini --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(ROOT + '/master_plot_hic.png', outfile.name, tolerance)
    print("saving test to {}".format(outfile.name))
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_cool_region():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0}/browser_tracks_cool.ini --region X:2500000-3500000 "\
           "--trackLabelFraction 0.23 --width 38 " \
           "--dpi 130 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(ROOT + '/master_plot_hic.png', outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_with_hic_rasterize_height():

    outfile = NamedTemporaryFile(suffix='.pdf', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0}/browser_tracks_hic_rasterize_height.ini --region "\
           "X:2500000-2600000 --trackLabelFraction 0.23 --width 38 " \
           "--dpi 10 --outFileName  {1}".format(ROOT, outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(ROOT + '/master_plot_hic_rasterize_height.pdf',
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
