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
[x-axis]

[spacer]
height = 0.5

[fasta_track]
file = fasta_track.fasta
title = Reference
height = 5

[annotation]
file = fasta_track.bed
title = Annotation
height = 2
color = darkblue
labels = false
fontsize = 10
file_type = bed
"""
with open(os.path.join(ROOT, "fasta_tracks.ini"), 'w') as fh:
    fh.write(browser_tracks)

browser_tracks = """
[x-axis]

[spacer]
height = 0.5

[fasta_track]
file = fasta_track_malformed.fasta
title = Reference
height = 5
"""
with open(os.path.join(ROOT, "fasta_track_malformed.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_fasta_zoomout():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "rDNA_unit_8919x2_bp:0-100"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_zoomout.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_zoomin():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "rDNA_unit_8919x2_bp:0-11"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_zoomin.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_inexistingchr():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "chr2:0-11"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_inexistingchr.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_zoomin_dec():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "rDNA_unit_8919x2_bp:0-11"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_zoomin_dec.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--decreasingXAxis "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_extrazoomin():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "rDNA_unit_8919x2_bp:0-3"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_extrazoomin.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_extrazoomin_height2():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "rDNA_unit_8919x2_bp:0-3"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_extrazoomin_height2.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--height 2 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_end_chr():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_tracks.ini")
    region = "rDNA_unit_8919x2_bp:17840-17850"
    expected_file = os.path.join(ROOT, 'master_fasta_tracks_end_chr.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)


def test_fasta_malformed():

    outfile = NamedTemporaryFile(suffix='.png', prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "fasta_track_malformed.ini")
    region = "test:48-60"
    expected_file = os.path.join(ROOT, 'master_fasta_track_malformed.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res

    os.remove(outfile.name)
