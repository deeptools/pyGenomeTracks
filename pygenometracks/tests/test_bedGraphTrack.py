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
