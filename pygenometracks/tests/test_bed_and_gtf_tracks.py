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
title = where =top

[spacer]
height = 0.05

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = UCSC; fontsize = 10
style = UCSC
fontsize = 10

[genes 2bis]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = UCSC; arrow_interval=10; fontsize = 10
style = UCSC
arrow_interval = 10
fontsize = 10

[spacer]
height = 1

[test bed6]
file = dm3_genes.bed6.gz
height = 7
title = bed6 border_color = black; gene_rows=10; fontsize=7; color=Reds
        (when a color map is used for the color (e.g. coolwarm, Reds) the bed
        score column mapped to a color)
fontsize = 7
file_type = bed
color = Reds
border_color = black
gene_rows = 10

[spacer]
height = 1

[test bed4]
file = dm3_genes.bed4.gz
height = 10
title = bed4 fontsize = 10; line_width = 1.5; global_max_row = true
        (global_max_row sets the number of genes per row as the maximum found
        anywhere in the genome, hence the white space at the bottom)
fontsize = 10
file_type = bed
global_max_row = true
line_width = 1.5

[spacer]
height = 1

[test gtf]
file = dm3_subset_BDGP5.78.gtf.gz
height = 10
title = gtf from ensembl
fontsize = 12
file_type = bed

[spacer]
height = 1

[test bed]
file = dm3_subset_BDGP5.78_asbed_sorted.bed.gz
height = 10
title = gtf from ensembl in bed12
fontsize = 12
file_type = bed

[spacer]
height = 1

[test gtf collapsed]
file = dm3_subset_BDGP5.78.gtf.gz
height = 10
title = gtf from ensembl one entry per gene
merge_transcripts = true
prefered_name = gene_name
fontsize = 12
file_type = bed

[spacer]
height = 1

[x-axis]
fontsize = 30
title = fontsize = 30

"""
with open(os.path.join(ROOT, "bed_and_gtf_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[x-axis]
where = top
title = where =top

[spacer]
height = 0.05

[genes 2]
file = dm3_genes.bed.gz
height = 3
title = genes (bed12) style = UCSC; fontsize = 10
style = UCSC
fontsize = 10

[genes 2bis]
file = dm3_genes.bed.gz
height = 3
title = genes (bed12) style = UCSC; arrow_interval=10; fontsize = 10
style = UCSC
arrow_interval = 10
fontsize = 10

[spacer]
height = 1

[test bed6]
file = dm3_genes.bed6.gz
height = 3
title = bed6 border_color = black; fontsize=8; color=red
fontsize = 8
file_type = bed
color = red
border_color = black

[spacer]
height = 1

[test bed6 arrowhead_included]
file = dm3_genes.bed6.gz
height = 3
title = bed6 border_color = black; fontsize=8; color=red; arrowhead_included = true
fontsize = 8
file_type = bed
color = red
border_color = black
arrowhead_included = true

[spacer]
height = 1

[test bed4]
file = dm3_genes.bed4.gz
height = 3
title = bed4 fontsize = 10; line_width = 1.5
fontsize = 10
file_type = bed
line_width = 1.5

[spacer]
height = 1

[test bed]
file = dm3_subset_BDGP5.78_asbed_sorted.bed.gz
height = 8
title = gtf from ensembl in bed12
fontsize = 12
file_type = bed

[spacer]
height = 1

[test bed]
file = dm3_subset_BDGP5.78_asbed_sorted.bed.gz
height = 8
title = gtf from ensembl in bed12; arrowhead_included = true
fontsize = 12
file_type = bed
arrowhead_included = true

[spacer]
height = 1

[x-axis]
fontsize = 30
title = fontsize = 30

[vlines]
type = vlines
file = dm3_genes.bed4.gz

"""
with open(os.path.join(ROOT, "bed_arrow_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[x-axis]
where = top

[spacer]
height = 0.05

[genes 1]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = flybase; fontsize = 10
style = UCSC
fontsize = 10

[spacer]
height = 1

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = flybase; fontsize = 10; max_labels = 600
style = UCSC
fontsize = 10
max_labels = 600

"""
with open(os.path.join(ROOT, "bed_maxLab_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[x-axis]
where = top

[spacer]
height = 0.05

[genes 0]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = flybase; fontsize = 10
style = flybase
fontsize = 10

[spacer]
height = 1

[genes 1]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = flybase; fontsize = 10; color_utr = red
style = flybase
fontsize = 10
color_utr = red

[spacer]
height = 1

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = flybase; fontsize = 10; height_utr = 0.7
style = flybase
fontsize = 10
height_utr = 0.7

"""
with open(os.path.join(ROOT, "bed_flybase_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[x-axis]
where = top
title = where =top

[spacer]
height = 0.05

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = UCSC
style = UCSC
fontsize = 10


[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = UCSC all_labels_inside = true
style = UCSC
fontsize = 10
all_labels_inside = true


[genes 3]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = UCSC labels_in_margin = true
style = UCSC
fontsize = 10
labels_in_margin = true


[genes 4]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = UCSC all_labels_inside = true labels_in_margin = true
style = UCSC
fontsize = 10
labels_in_margin = true
all_labels_inside = true
"""
with open(os.path.join(ROOT, "bed_all_labels_inside.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[x-axis]
where = top

[spacer]
height = 0.05

[test bed4]
file = dm3_genes.bed4.gz
height = 7
title = bed4 style = tssarrow
style = tssarrow
fontsize = 10

[spacer]
height = 1

[test bed6]
file = dm3_genes.bed6.gz
height = 7
title = bed6 style = tssarrow
style = tssarrow
fontsize = 10

[spacer]
height = 1

[genes 0]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = tssarrow; fontsize = 10
style = tssarrow
fontsize = 10

[genes 0]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = tssarrow; fontsize = 10; arrow_length = 5000 (5kb)
style = tssarrow
arrow_length = 5000
fontsize = 10

[spacer]
height = 1

[genes 1]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = tssarrow; fontsize = 10; color_utr = red
style = tssarrow
fontsize = 10
color_utr = red

[spacer]
height = 1

[genes 2]
file = dm3_genes.bed.gz
height = 7
title = genes (bed12) style = tssarrow; fontsize = 10; height_utr = 0.7
style = tssarrow
fontsize = 10
height_utr = 0.7
"""

with open(os.path.join(ROOT, "bed_tssarrow_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[x-axis]
where = top

[spacer]
height = 0.5

[genes 0]
file = hoxd_genes_rgb.bed.gz
height = 7
title = genes (bed12) style = tssarrow; fontsize = 10; color = bed_rgb
style = tssarrow
fontsize = 10
color = bed_rgb

[spacer]

[genes 1]
file = hoxd_genes_rgb.bed.gz
height = 7
title = genes (bed12) style = tssarrow; fontsize = 10; color = bed_rgb; color_utr = bed_rgb
style = tssarrow
fontsize = 10
color = bed_rgb
color_utr = bed_rgb

[spacer]

[genes 2]
file = hoxd_genes_rgb.bed.gz
height = 7
title = genes (bed12) style = UCSC; fontsize = 10; color = bed_rgb; color_utr = bed_rgb
style = UCSC
fontsize = 10
color = bed_rgb
color_utr = bed_rgb

[spacer]

[genes 3]
file = hoxd_genes_rgb.bed.gz
height = 7
title = genes (bed12) style = flybase; fontsize = 10; color = bed_rgb; color_utr = bed_rgb
style = flybase
fontsize = 10
color = bed_rgb
color_utr = bed_rgb

[spacer]

[genes 4]
file = hoxd_genes_noGm_rgb.bed.gz
height = 3
title = genes except Gm (bed12) style = tssarrow; fontsize = 10; color = bed_rgb; color_utr = bed_rgb
style = tssarrow
fontsize = 10
color = bed_rgb
color_utr = bed_rgb

[spacer]

[genes 5]
file = hoxd_genes_noGm_rgb.bed.gz
height = 1
title = genes except Gm (bed12) style = tssarrow; fontsize = 10; color = bed_rgb; color_utr = bed_rgb; display = collapsed
style = tssarrow
fontsize = 10
color = bed_rgb
color_utr = bed_rgb
display = collapsed
"""
with open(os.path.join(ROOT, "bed_genes_rgb.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[scores]
file = dm3_genes_withrgbandscore.bed.gz
title = genes with scores in coding bed_rgb in utr
color = cool_r
color_utr = bed_rgb
border_color = black
height = 10.0

[scores in utr]
file = dm3_genes_withrgbandscore.bed.gz
title = genes with scores in utr, coding is blue, border is bed_rgb
color = blue
color_utr = viridis
border_color = bed_rgb
height = 10.0

[scores in border_color]
file = dm3_genes_withrgbandscore.bed.gz
title = genes with scores in border_color, coding and utr are blue
color = blue
color_utr = blue
border_color = ['red', 'orange']
height = 10.0

[scores in cod+utr]
file = dm3_genes_withrgbandscore.bed.gz
title = genes with scores both in coding and utr as ['blue', 'purple'] (this does not work for color_utr)
color = ['blue', 'purple']
color_utr = ['blue', 'purple']
height = 10.0

[scores in cod+utr 2]
file = dm3_genes_withrgbandscore.bed.gz
title = genes with scores both in coding and utr as Reds
color = Reds
color_utr = Reds
height = 10.0
"""
with open(os.path.join(ROOT, "bed_colormap_genes.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_bed_and_gtf():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3000000-3300000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_and_gtf_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_and_gtf.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_arrow():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3120000-3150000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_arrow_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_arrow.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_arrow_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3130000-3140000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_arrow_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_arrow_zoom.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_flybase():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3000000-3300000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_flybase_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_flybase.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_tssarrow():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3000000-3300000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_tssarrow_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_tssarrow.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_tssarrow_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3020000-3070000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_tssarrow_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_tssarrow_zoom.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_tssarrow_zoom2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3130000-3150000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_tssarrow_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_tssarrow_zoom2.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_with_maxLab():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:2000000-3500000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_maxLab_tracks.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_maxLab.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_genes_rgb():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region chr2:74,650,000-74,710,000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_genes_rgb.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_genes_rgb.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_all_label_inside():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3100000-3200000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 " \
           "--trackLabelHAlign right --outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_all_labels_inside.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_all_label_inside.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_scores():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    args = "--tracks {0} --region X:3000000-3300000 "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--outFileName {1}" \
           "".format(os.path.join(ROOT, 'bed_colormap_genes.ini'),
                     outfile.name).split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(os.path.join(ROOT, 'master_bed_colormap_genes.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
