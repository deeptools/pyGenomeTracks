# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import pygenometracks.plotTracks
from pygenometracks.utilities import InputError


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
title = bed6 border_color = black; gene_rows=10; fontsize=7; color=Reds (when a color map is used for the color (e.g. coolwarm, Reds) the bed score column mapped to a color)
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
title = bed4 fontsize = 10; line_width = 1.5; global_max_row = true (global_max_row sets the number of genes per row as the maximum found anywhere in the genome, hence the white space at the bottom)
fontsize = 10
file_type = bed
global_max_row = true
line_width = 1.5

[spacer]
height = 1

[test gtf]
file = dm3_subset_BDGP5.78_gtf.dat
height = 10
title = gtf from ensembl (with dat extension)
fontsize = 12
file_type = gtf

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
file_type = gtf

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

[genes 2bix]
file = hoxd_genes_rgb.bed.gz
height = 7
title = same but color_backbone = bed_rgb
style = UCSC
fontsize = 10
color = bed_rgb
color_utr = bed_rgb
color_backbone = bed_rgb

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

[genes 3bis]
file = hoxd_genes_rgb.bed.gz
height = 7
title = same but color_backbone = bed_rgb
style = flybase
fontsize = 10
color = bed_rgb
color_utr = bed_rgb
color_backbone = bed_rgb

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

with open(os.path.join(ROOT, "bed_genes_rgb_incorrect.ini"), 'w') as fh:
    fh.write(browser_tracks.replace('style = flybase', 'style = inexisting'))

wrong_track = """
[test gtf]
file = dm3_subset_BDGP5.78_gtf.dat
height = 10
title = gtf from ensembl (with dat extension)
fontsize = 12
file_type = bed
"""
with open(os.path.join(ROOT, "gtf_as_bed.ini"), 'w') as fh:
    fh.write(wrong_track)


wrong_track = """
[test bed]
file = dm3_genes.bed.gz
height = 10
title = bed file
fontsize = 12
file_type = gtf
"""
with open(os.path.join(ROOT, "bed_as_gtf.ini"), 'w') as fh:
    fh.write(wrong_track)

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
title = genes with scores both in coding and utr as Reds min_value = 0.2 max_value = 2
color = Reds
color_utr = Reds
height = 10.0
min_value = 0.2
max_value = 2
"""
with open(os.path.join(ROOT, "bed_colormap_genes.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[genes]
file = dm3_genes_withrgbandscore.bed.gz
title = bed color = Reds
color = Reds
height = 4

[spacer]

[genes]
file = dm3_genes_withrgbandscore_shuffled.bed.gz
title = bed color = Reds bed is not sorted
color = Reds
height = 4

[spacer]

[genes]
file = dm3_genes_withrgbandscore_shuffled.bed.gz
title = bed color = Reds global_max_row = true and bed is not sorted
color = Reds
global_max_row = true
height = 4

[x-axis]
"""
with open(os.path.join(ROOT, "bed_shuffle.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[genes]
file = dm3_genes_withrgbandscore.bed.gz
title = bed global_max_row = true labels = false
labels = false
global_max_row = true
height = 4

[spacer]

[x-axis]
title = centered title

[vlines]
file = tad_classification.bed
type = vlines
line_width = 3
"""
with open(os.path.join(ROOT, "bed_vlines.ini"), 'w') as fh:
    fh.write(browser_tracks)
with open(os.path.join(ROOT, "bed_vlines_incorrect.ini"), 'w') as fh:
    fh.write(browser_tracks + 'line_style = dashed\n')

browser_tracks = """
[genes1]
file = example.bed
title = bed style = flybase
style = flybase
height = 4

[spacer]

[genes2]
file = example.bed
title = bed style = UCSC
style = UCSC
height = 4

[spacer]

[genes1]
file = example.bed
title = bed style = tssarrow
style = tssarrow
height = 4

[spacer]

[x-axis]
"""
with open(os.path.join(ROOT, "bed_different_UTR.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[genes1]
file = no_exon.gtf
title = default

[spacer]

[genes2]
file = no_exon.gtf
title = merge_transcripts=true
merge_transcripts = true

[spacer]

[x-axis]
"""
with open(os.path.join(ROOT, "gtf_no_exon.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_tracks_bed_and_gtf():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'bed_and_gtf_tracks.ini')
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_bed_and_gtf.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_arrow():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'bed_arrow_tracks.ini')
    region = "X:3120000-3150000"
    expected_file = os.path.join(ROOT, 'master_bed_arrow.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_arrow_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)

    ini_file = os.path.join(ROOT, "bed_arrow_tracks.ini")
    region = "X:3130000-3140000"
    expected_file = os.path.join(ROOT, 'master_bed_arrow_zoom.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_flybase():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)

    ini_file = os.path.join(ROOT, "bed_flybase_tracks.ini")
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_bed_flybase.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_tssarrow():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_tssarrow_tracks.ini")
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_bed_tssarrow.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_tssarrow_zoom():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_tssarrow_tracks.ini")
    region = "X:3020000-3070000"
    expected_file = os.path.join(ROOT, 'master_bed_tssarrow_zoom.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_tssarrow_zoom2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_tssarrow_tracks.ini")
    region = "X:3130000-3150000"
    expected_file = os.path.join(ROOT, 'master_bed_tssarrow_zoom2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_with_maxLab():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_maxLab_tracks.ini")
    region = "X:2000000-3500000"
    expected_file = os.path.join(ROOT, 'master_maxLab.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_genes_rgb():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    for suf in ['', '_incorrect']:
        ini_file = os.path.join(ROOT, f"bed_genes_rgb{suf}.ini")
        region = "chr2:74,650,000-74,710,000"
        expected_file = os.path.join(ROOT, 'master_bed_genes_rgb.png')
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        # remove the incorrect ini file
        if 'incorrect' in ini_file:
            os.remove(ini_file)


def test_plot_tracks_bed_all_label_inside():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_all_labels_inside.ini")
    region = "X:3100000-3200000"
    expected_file = os.path.join(ROOT, 'master_bed_all_label_inside.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--trackLabelHAlign right "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_all_label_inside_Xdec():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_all_labels_inside.ini")
    region = "X:3215000-3240000"
    expected_file = os.path.join(ROOT, 'master_bed_all_label_inside_dec.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--trackLabelHAlign right --decreasingXAxis "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_gtf_as_bed():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "gtf_as_bed.ini")
    region = "X:3100000-3200000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'This is probably not a bed file.' in str(e)
    else:
        raise Exception("The gtf as bed should fail.")
    os.remove(ini_file)


def test_bed_as_gtf():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_as_gtf.ini")
    region = "X:3100000-3200000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'This is not a gtf file.' in str(e)
    else:
        raise Exception("The bed as gtf should fail.")
    os.remove(ini_file)


def test_plot_tracks_bed_scores():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_colormap_genes.ini")
    region = "X:3000000-3300000"
    expected_file = os.path.join(ROOT, 'master_bed_colormap_genes.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_bed_shuffle():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_shuffle.ini")
    bed_file = os.path.join(ROOT, 'regions_chr1XY.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['chr1:0-500000', 'chrX:2500000-2600000', 'chrY:0-1000000']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_bed_shuffle_'
                                     + region_str + extension)
        res = compare_images(expected_file,
                             output_file, tolerance)
        assert res is None, res

        os.remove(output_file)


def test_plot_tracks_bed_vlines():
    extension = '.png'
    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    bed_file = os.path.join(ROOT, 'regionsXfakeChr.bed')
    for suf in ['', '_incorrect']:
        ini_file = os.path.join(ROOT, f"bed_vlines{suf}.ini")
        args = f"--tracks {ini_file} --BED {bed_file} "\
               "--trackLabelFraction 0.5 --width 38 --dpi 130 "\
               "--trackLabelHAlign center "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        for region in ['X:3000000-3300000', 'fakeChr:0-100']:
            region_str = region.replace(':', '-')
            output_file = outfile.name[:-4] + '_' + region_str + extension
            expected_file = os.path.join(ROOT, 'master_bed_vlines_'
                                         + region_str + extension)
            res = compare_images(expected_file,
                                 output_file, tolerance)
            assert res is None, res

            os.remove(output_file)
        if 'incorrect' in ini_file:
            os.remove(ini_file)


def test_plot_tracks_bed_different_UTR():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_different_UTR.ini")
    region = "chr1:0-500"
    expected_file = os.path.join(ROOT, 'master_different_UTR.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_gtf_no_exon():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, 'gtf_no_exon.ini')
    region = "381:0-1000"
    expected_file = os.path.join(ROOT, 'master_gtf_no_exon.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
