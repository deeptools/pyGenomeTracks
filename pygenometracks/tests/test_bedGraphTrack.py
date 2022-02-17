# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os
import pygenometracks.plotTracks
from pygenometracks.utilities import InputError

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
[x-axis]
where = top

[spacer]
height = 0.05

[test bedgraph use middle]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
color = blue
alpha = 0.5
height = 5
title = bedgraph with use_middle = true 2 files red and blue and minimum black
min_value = 0
max_value = 10
use_middle = true

[test bedgraph2 use middle]
file = GSM3182415_E12PHL_WT_Hoxd11vp.bedgraph.gz
color = red
alpha = 0.5
use_middle = true
overlay_previous = share-y

[test bedgraph min use middle]
file = GSM3182415_E12PHL_WT_Hoxd11vp.bedgraph.gz
second_file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
operation = min(file, second_file)
color = black
alpha = 1
use_middle = true
overlay_previous = share-y

[diff use middle]
file = GSM3182416_E12DHL_WT_Hoxd11vp.bedgraph.gz
second_file = GSM3182415_E12PHL_WT_Hoxd11vp.bedgraph.gz
operation = file - second_file
color = blue
negative_color = red
alpha = 0.5
height = 5
title = bedgraph with use_middle = true difference
min_value = -5
max_value = 5
use_middle = true


[genes]
file = HoxD_cluster_regulatory_regions_mm10.bed
height = 3
title = HoxD genes and regulatory regions

"""
with open(os.path.join(ROOT, "bedgraph_useMid_op.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[x-axis]
where = top

[test bedgraph neg]
file = test_with_neg_values.bg.gz
color = blue
negative_color = red
title = color = blue; negative_color = red
height = 5

[test bedgraph neg]
file = test_with_neg_values.bg.gz
color = cyan
negative_color = darkred
type = line:2
title = color = cyan; negative_color = darkred; type = line:2
height = 5

[test bedgraph neg]
file = test_with_neg_values.bg.gz
color = black
negative_color = lime
type = points:2
title = color = black; negative_color = lime; type = points:2
height = 5
"""

with open(os.path.join(ROOT, "bedgraph_negative.ini"), 'w') as fh:
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
second_file = bedgraph_chrx_2e6_5e6_2.bg.bgz
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
with open(os.path.join(ROOT, "incorrect_operation_bdg.ini"), 'w') as fh:
    fh.write(browser_tracks.replace('operation = log1p(file)\n', 'operation = log1p(file)\ny_axis_values = original\n'))

bedgraph_withNA = """
[test bedgraph withNA]
file = bedgraph_withNA.bdg
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "bedgraph_withNA.ini"), 'w') as fh:
    fh.write(bedgraph_withNA)

for suff in ["", "2", "3"]:
    invalid_bedgraph = f"""
[invalid_bedgraph{suff}]
file = invalid_bedgraph{suff}.bdg
height = 3

[x-axis]
"""
    with open(os.path.join(ROOT, f"invalid_bedgraph{suff}.ini"), 'w') as fh:
        fh.write(invalid_bedgraph)

unsorted_bedgraph = """
[unsorted_bedgraph]
file = unsorted_bedgraph.bdg
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "unsorted_bedgraph.ini"), 'w') as fh:
    fh.write(unsorted_bedgraph)

# Write a bedgraph with negative values:
with open(os.path.join(ROOT, "bedgraph_chrx_2e6_5e6.bg"), 'r') as f:
    with open(os.path.join(ROOT, "bedgraph_chrx_2e6_5e6_m.bg"), 'w') as fo:
        for line in f:
            ls = line.split('\t')
            ls[3] = '-' + ls[3]
            fo.write('\t'.join(ls))

log1p_with_neg = """
[m_bedgraph]
file = bedgraph_chrx_2e6_5e6_m.bg
transform = log1p
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "log1pm_bedgraph.ini"), 'w') as fh:
    fh.write(log1p_with_neg)


bedgraph_end_not_covered = """
[bedgraph]
file = simple.bdg
height = 3
summary_method = max

[x-axis]
"""
with open(os.path.join(ROOT, "bedgraph_end_not_covered.ini"), 'w') as fh:
    fh.write(bedgraph_end_not_covered)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_plot_bedgraph_tracks_with_bed():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bedgraph_useMid.ini")
    bed_file = os.path.join(ROOT, 'regions_imbricated_chr2.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['chr2:73800000-75744000', 'chr2:74000000-74800000']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_bedgraph_useMid_'
                                     + region_str + extension)
        res = compare_images(expected_file,
                             output_file, tolerance)
        assert res is None, res

        os.remove(output_file)


def test_plot_bedgraph_tracks_individual():
    extension = '.png'

    for region in ['chr2:73800000-75744000', 'chr2:74000000-74800000']:
        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        ini_file = os.path.join(ROOT, "bedgraph_useMid.ini")
        region_str = region.replace(':', '-')
        expected_file = os.path.join(ROOT, 'master_bedgraph_useMid_'
                                     + region_str + extension)
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)


def test_plot_bedgraph_use_mid_op():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bdg_NA_', delete=False)
    ini_file = os.path.join(ROOT, "bedgraph_useMid_op.ini")
    region = "chr2:74000000-74800000"
    expected_file = os.path.join(ROOT, 'master_bedgraph_useMid_op.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_bedgraph_tracks_rasterize():

    outfile = NamedTemporaryFile(suffix='.pdf', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bedgraph_useMid.ini")
    region = "chr2:73,800,000-75,744,000"
    expected_file = os.path.join(ROOT, 'master_bedgraph_useMid.pdf')
    # matplotlib compare on pdf will create a png next to it.
    # To avoid issues related to write in test_data folder
    # We copy the expected file into a temporary place
    new_expected_file = NamedTemporaryFile(suffix='.pdf',
                                           prefix='pyGenomeTracks_test_',
                                           delete=False)
    os.system(f'cp {expected_file} {new_expected_file.name}')
    expected_file = new_expected_file.name
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
    os.remove(expected_file.replace('.pdf', '_pdf.png'))


def test_op_bdg():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bdg_op_test_',
                                 delete=False)
    for pref in ['', 'incorrect_']:
        ini_file = os.path.join(ROOT, f"{pref}operation_bdg.ini")
        region = "X:2700000-3100000"
        expected_file = os.path.join(ROOT, 'master_operation_bdg.png')
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res

        os.remove(outfile.name)
        if 'incorrect' in ini_file:
            os.remove(ini_file)


def test_bdg_withNA():
    outfile = NamedTemporaryFile(suffix='.png', prefix='bdg_NA_', delete=False)
    ini_file = os.path.join(ROOT, "bedgraph_withNA.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_bedgraph_withNA.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_bdg_unsorted():
    outfile = NamedTemporaryFile(suffix='.png', prefix='pgt_test_', delete=False)
    ini_file = os.path.join(ROOT, "unsorted_bedgraph.ini")
    region = "X:2700000-3100000"
    expected_file = os.path.join(ROOT, 'master_bedgraph_withNA.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_negative():
    region = "X:2700000-3100000"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_negative_test_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "bedgraph_negative.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_negative.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_invalid_bedgraph():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pgt_test_', delete=True)
    ini_file = os.path.join(ROOT, "invalid_bedgraph.ini")
    region = "X:2700000-3100000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'Only one field detected' in str(e)
    else:
        raise Exception("The invalid_bedgraph should fail.")

    os.remove(ini_file)


def test_invalid_bedgraph2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pgt_test_', delete=True)
    ini_file = os.path.join(ROOT, "invalid_bedgraph2.ini")
    region = "X:2700000-3100000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'The start field is not an integer.' in str(e)
    else:
        raise Exception("The invalid_bedgraph2 should fail.")

    os.remove(ini_file)


def test_invalid_bedgraph3():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pgt_test_', delete=True)
    ini_file = os.path.join(ROOT, "invalid_bedgraph3.ini")
    region = "X:2700000-3100000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except InputError as e:
        assert 'The end field is not an integer.' in str(e)
    else:
        raise Exception("The invalid_bedgraph3 should fail.")

    os.remove(ini_file)


def test_bedgraph_neg_log1p():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pgt_test_', delete=True)
    ini_file = os.path.join(ROOT, "log1pm_bedgraph.ini")
    region = "X:2700000-3100000"
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except Exception as e:
        assert 'coverage contains values below or equal to - 1' in str(e)
    else:
        raise Exception("The bedgraph_neg_log1p should fail.")

    os.remove(ini_file)
    os.remove(os.path.join(ROOT, "bedgraph_chrx_2e6_5e6_m.bg"))


def test_bedgraph_end_not_covered():
    region = "chr7:100-400"
    outfile = NamedTemporaryFile(suffix='.png', prefix='bedgraph_end_not_covered_', delete=False)
    args = "--tracks {ini} --region {region} --trackLabelFraction 0.2 " \
           "--dpi 130 --outFileName {outfile}" \
           "".format(ini=os.path.join(ROOT, "bedgraph_end_not_covered.ini"),
                     outfile=outfile.name, region=region).split()
    pygenometracks.plotTracks.main(args)
    print("saving test to {}".format(outfile.name))
    res = compare_images(os.path.join(ROOT, 'master_bedgraph_end_not_covered.png'),
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
