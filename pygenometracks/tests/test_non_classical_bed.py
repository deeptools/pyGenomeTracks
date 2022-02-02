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
[macs2 broadPeak]
file = broadPeak.broadPeak
title = broadPeak
file_type = bed

[spacer]

[macs2 gappedPeak]
file = gappedPeak.gappedPeak
title = gappedPeak
file_type = bed

[spacer]

[macs2 filteredbed]
file = filtered.results.bed
title = filtered.results.bed (strange format)
file_type = bed
"""
with open(os.path.join(ROOT, "bed_unusual_formats.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[strange_strand]
file = strange_strand.bed
"""
with open(os.path.join(ROOT, "strange_strand.ini"), 'w') as fh:
    fh.write(browser_tracks)


browser_tracks = """
[invalid_strand]
file = invalid_strand.bed
"""
with open(os.path.join(ROOT, "invalid_strand.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_rgb]
file = invalid_rgb.bed
title = first line valid but neg invalid
color = bed_rgb

[invalid_rgb2]
file = invalid_rgb2.bed
title = first line invalid
color = bed_rgb
"""
with open(os.path.join(ROOT, "invalid_rgb.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_blockCount]
file = invalid_blockCount.bed
title = invalid block count in first line
color = bed_rgb

[spacer]

[invalid_blockCount]
file = invalid_blockCount2.bed
title = invalid block count in not first line -> bedtools raised an error
style = UCSC
height = 2
"""
with open(os.path.join(ROOT, "invalid_blockCount.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_CDScoo]
file = invalid_CDScoo.bed
title = invalid CDS coordinate in first line rgb is ignored
color = bed_rgb

[spacer]

[invalid_CDScoo2]
file = invalid_CDScoo2.bed
title = invalid CDS coordinate in not first line rgb can be used
color = bed_rgb
height = 2

[spacer]

[invalid_CDScoo3]
file = invalid_CDScoo3.bed
title = invalid CDS coordinate in not first line bed12 can be used
color = bed_rgb
style = UCSC
height = 2
"""
with open(os.path.join(ROOT, "invalid_CDScoo.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_blocks]
file = invalid_blocks.bed
title = block_length is not integer in first line, block info is ignored
style = UCSC

[spacer]

[invalid_blocks2]
file = invalid_blocks2.bed
title = block_length is not integer in not first line blocks are used
style = UCSC
height = 4
"""
with open(os.path.join(ROOT, "invalid_blocks.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_blocks3]
file = invalid_blocks3.bed
title = block_count_inconsistent between block count and others
"""
with open(os.path.join(ROOT, "invalid_blocks3.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_blocks4]
file = invalid_blocks4.bed
title = block_count_inconsistent between lengths and others
"""
with open(os.path.join(ROOT, "invalid_blocks4.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[bed with 0 block]
file = example_zeroBlock.bed
title = example with number of block = 0 will fail
"""
with open(os.path.join(ROOT, "invalid_zero_block.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[invalid_score]
file = invalid_score.bed
title = invalid score in first line strand and rgb is ignored
color = bed_rgb

[spacer]

[invalid_score2]
file = invalid_score2.bed
title = invalid score in not first line strand and rgb can be used
color = bed_rgb
height = 2
"""
with open(os.path.join(ROOT, "invalid_score.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


# For this test I need to fix the tolerence at 19 to make it work
# On python 3.7 and 3.8
# The only difference between outputs is the font width...
def test_plot_tracks_bed_unusual_format():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "bed_unusual_formats.ini")
    region = "X:20000-40000"
    expected_file = os.path.join(ROOT, 'master_bed_unusual_formats.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance + 6)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_strange_strand():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "strange_strand.ini")
    region = "chr1:0-500"
    expected_file = os.path.join(ROOT, 'master_strange_strand.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_invalid_strand():

    tolerance = 5

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_strand.ini")
    region = "chr1:0-500"
    expected_file = os.path.join(ROOT, 'master_strange_strand.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
    os.remove(ini_file)


def test_plot_tracks_bed_invalid_rgb():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_rgb.ini")
    region = "chr1:0-500"
    expected_file = os.path.join(ROOT, 'master_invalid_rgb.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_invalid_blockCount():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_blockCount.ini")
    region = "chrX:15000-24000"
    expected_file = os.path.join(ROOT, 'master_invalid_blockCount.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_invalid_CDScoo():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_CDScoo.ini")
    region = "chrX:15000-24000"
    expected_file = os.path.join(ROOT, 'master_invalid_CDScoo.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_invalid_blocks():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_blocks.ini")
    region = "chrX:15000-24000"
    expected_file = os.path.join(ROOT, 'master_invalid_blocks.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_invalid_score():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_score.ini")
    region = "chrX:15000-24000"
    expected_file = os.path.join(ROOT, 'master_invalid_score.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_plot_tracks_bed_invalid_block_count():

    region = "chrX:15000-24000"

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    for suf in ['3', '4']:
        ini_file = os.path.join(ROOT, f"invalid_blocks{suf}.ini")
        args = f"--tracks {ini_file} --region {region} "\
               "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
               f"--outFileName {outfile.name}".split()
        try:
            pygenometracks.plotTracks.main(args)
        except Exception as e:
            assert 'The number of blocks:' in str(e)
        else:
            raise Exception(f"bed_invalid_block{suf} should fail.")
        os.remove(ini_file)


def test_plot_tracks_bed_invalid_zero_block():

    region = "chr1:0-500"

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "invalid_zero_block.ini")
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    try:
        pygenometracks.plotTracks.main(args)
    except Exception as e:
        assert 'File type detected is bed12 but line' in str(e)
    else:
        raise Exception("invalid_zero_block should fail.")
    os.remove(ini_file)
