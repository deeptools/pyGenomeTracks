# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[x-axis]
where = top

[spacer]
height = 0.05

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph rasterize = true
rasterize = true
max_value = 10

[test bedgraph]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph
max_value = 10

[test bedgraph use middle]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
height = 5
title = bedgraph with use_middle = true
max_value = 10
use_middle = true

[genes]
file = HoxD_cluster_regulatory_regions_mm10.bed
height = 3
title = HoxD genes and regulatory regions

"""
with open(os.path.join(ROOT, "bedgraph_useMid.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[test file]
file = bedgraph_chrx_2e6_5e6.bg.bgz
color = blue
height = 4
title = file summary_method = mean
summary_method = mean
min_value = 0
max_value = 30

[test file]
file = bedgraph2_X_2.5e6_3.5e6.bdg
color = red
height = 4
title = second_file summary_method = mean
summary_method = mean
min_value = 0
max_value = 30


[spacer]
height = 0.5

[test op0]
file = bedgraph_chrx_2e6_5e6.bg.bgz
color = blue
height = 4
title = operation = log1p(file) (no summary_method)
operation = log1p(file)


[spacer]
height = 0.5

[test op1]
file = bedgraph_chrx_2e6_5e6.bg.bgz
second_file = bedgraph2_X_2.5e6_3.5e6.bdg
color = blue
negative_color = red
height = 8
title = operation = file - second_file
operation = file - second_file
min_value = -30
max_value = 30

[spacer]
height = 0.5

[test op2]
file = bedgraph2_X_2.5e6_3.5e6.bdg
second_file = bedgraph_chrx_2e6_5e6.bg.bgz
color = red
negative_color = blue
height = 8
title = operation = file - second_file (but files were switched)
operation = file - second_file
min_value = -30
max_value = 30

[spacer]
height = 0.5

[test op2]
file = bedgraph2_X_2.5e6_3.5e6.bdg
second_file = bedgraph_chrx_2e6_5e6.bg.bgz
color = red
negative_color = blue
height = 8
title = idem but nans_to_zeros = true
operation = file - second_file
nans_to_zeros = true
min_value = -30
max_value = 30

[spacer]
height = 0.5

[x-axis]
"""
with open(os.path.join(ROOT, "operation_bdg.ini"), 'w') as fh:
    fh.write(browser_tracks)


tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_bedgraph_tracks():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region chr2:73,800,000-75,744,000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, "bedgraph_useMid.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_bedgraph_useMid.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_bedgraph_tracks_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region chr2:74,000,000-74,800,000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, "bedgraph_useMid.ini"),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_bedgraph_useMid_zoom.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_bedgraph_tracks_rasterize():

    outfile = NamedTemporaryFile(suffix='.pdf', prefix='pyGenomeTracks_test_', delete=False)
    args = "--tracks {0} --region chr2:73,800,000-75,744,000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bedgraph_useMid.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_bedgraph_useMid.pdf'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_op_bdg():
    region = "X:2700000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bdg_op_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "operation_bdg.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_operation_bdg.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
