# -*- coding: utf-8 -*-
import matplotlib as mpl
mpl.use('agg')
from matplotlib.testing.compare import compare_images
from tempfile import NamedTemporaryFile
import os.path
import shutil
import pygenometracks.plotTracks


ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

browser_tracks = """
[maf]
file = first.maf
reference = mm10
title = default
height = 3

[spacer]

[maf]
file = first.maf
reference = mm10
title = choose order
species_order = hg19 rn5   sorAra1
species_labels = Human Rat  Shrew
height = 3

[spacer]

[maf]
file = first.maf
reference = mm10
title = specify only danRer11
species_order = danRer11
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "first_maf.ini"), 'w') as fh:
    fh.write(browser_tracks)

with open(os.path.join(ROOT, "first_maf_order_species_only.ini"), 'w') as fh:
    fh.write(browser_tracks.replace("species_order", "species_order_only = true\nspecies_order"))

browser_tracks = """
[maf]
file = first.maf
reference = mm10
title = display_ref_seq = true
display_ref_seq = true
height = 3

[x-axis]
"""
with open(os.path.join(ROOT, "first_maf_seq.ini"), 'w') as fh:
    fh.write(browser_tracks)

with open(os.path.join(ROOT, "first_maf_seq_invalid.ini"), 'w') as fh:
    fh.write(browser_tracks.replace("first.maf\nreference",
                                    "first.maf\nfile_index = first.maf.hg19.index\nreference"))

with open(os.path.join(ROOT, "first_maf_seq_invalid2.ini"), 'w') as fh:
    fh.write(browser_tracks.replace("first.maf\nreference",
                                    "first.maf\nspecies_labels = Mouse Rat\nreference"))

with open(os.path.join(ROOT, "first_maf_seq_invalid3.ini"), 'w') as fh:
    fh.write(browser_tracks.replace("first.maf\nreference",
                                    "first.maf\nspecies_labels = Mouse Rat\nspecies_order = rn5\nreference"))

browser_tracks = """
[maf]
file = first.maf
file_index = first.maf.hg19.index
reference = hg19
title = display_ref_seq = true rasterize = true
display_ref_seq = true
height = 3
rasterize = true

[x-axis]
"""
with open(os.path.join(ROOT, "first_maf_seq_hg19.ini"), 'w') as fh:
    fh.write(browser_tracks)

# To get the second_maf:
# wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/multiz60way/maf/chr2.maf.gz
# gunzip -k chr2.maf.gz
# maf_build_index.py chr2.maf -s mm10
# echo -e "74071500\t74077028" > isl2.txt
# echo -e "74060473\t74082287" > isl2.txt
# maf_extract_ranges_indexed.py chr2.maf -s mm10.chr2 < isl2.txt > mm10_chr2_isl2.maf
# maf_limit_to_species.py mm10,rn5,hetGla2,hg19,susScr3,felCat5,ailMel1,pteVam1,loxAfr3,triMan1,ornAna1,galGal4,xenTro3,latCha1,fr3,danRer7 < mm10_chr2_isl2.maf > mm10_chr2_isl2_lessspe.maf
# mouse rat naked_mole_rat human pig cat panda megabat elephant manatee platypus chicken x_tropicalis coelacanth fugu zebrafish
# maf_build_index.py pygenometracks/tests/test_data/first.maf pygenometracks/tests/test_data/first.maf.hg19.index -s hg19

browser_tracks = """
[maf]
file = mm10_chr2_isl2_lessspe.maf
reference = mm10
title = default
height = 3

[spacer]

[maf]
file = mm10_chr2_isl2_lessspe.maf
reference = mm10
title = choose order Platypus Elephant Megabat
species_order = ornAna1 loxAfr3 pteVam1
species_labels = Platypus Elephant Megatbat
height = 3

[spacer]

[maf]
file = mm10_chr2_isl2_lessspe.maf
reference = mm10
title = species_order_only Platypus Elephant Megabat and show seq
species_order = ornAna1 loxAfr3 pteVam1
species_labels = Platypus Elephant Megabat
species_order_only = true
display_ref_seq = true
height = 5

[x-axis]
"""
with open(os.path.join(ROOT, "maf_withe.ini"), 'w') as fh:
    fh.write(browser_tracks)

tolerance = 13  # default matplotlib pixed difference tolerance


def test_first_maf():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf.ini")
    bed_file = os.path.join(ROOT, 'regions_maf.bed')
    args = f"--tracks {ini_file} --BED {bed_file} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    for region in ['chr2:34704975-34705208', 'chr2:34705032-34707346']:
        region_str = region.replace(':', '-')
        output_file = outfile.name[:-4] + '_' + region_str + extension
        expected_file = os.path.join(ROOT, 'master_first_maf_'
                                     + region_str + extension)
        res = compare_images(expected_file,
                             output_file, tolerance)
        assert res is None, res

        os.remove(output_file)


def test_first_maf_empty_chr():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf.ini")
    region = "chr1:0-1000"
    expected_file = os.path.join(ROOT, 'master_first_maf_empty_chr.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_first_maf_order_species_only():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf_order_species_only.ini")
    region = "2:34705032-34707346"
    expected_file = os.path.join(ROOT, 'master_first_maf_order_species_only.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_first_maf_seq():
    extension = '.png'

    for suf in ['', '_invalid', '_invalid2', '_invalid3']:
        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        ini_file = os.path.join(ROOT, f"first_maf_seq{suf}.ini")
        region = "2:34704975-34705208"
        expected_file = os.path.join(ROOT, 'master_first_maf_seq.png')
        args = f"--tracks {ini_file} --region {region} "\
            "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
            f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)
        if 'invalid' in ini_file:
            os.remove(ini_file)


def test_first_maf_seq_zoom():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf_seq.ini")
    region = "2:34705120-34705150"
    expected_file = os.path.join(ROOT, 'master_first_maf_seq_zoom.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_first_maf_seq_zoom_dec():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf_seq.ini")
    region = "2:34705120-34705150"
    expected_file = os.path.join(ROOT, 'master_first_maf_seq_zoom_dec.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           "--decreasingXAxis "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_first_maf_seq_hg19():
    extension = '.png'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf_seq_hg19.ini")
    region = "9:128093930-128093970"
    expected_file = os.path.join(ROOT, 'master_first_maf_seq_hg19.png')
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)


def test_first_maf_seq_zoom_change_height():
    if mpl.__version__ == "3.1.1":
        my_tolerance = 17
    else:
        my_tolerance = tolerance
    extension = '.png'
    ini_file = os.path.join(ROOT, "first_maf_seq.ini")
    region = "2:34705120-34705150"

    for height in [2, 10]:
        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        expected_file = os.path.join(ROOT, f'master_first_maf_seq_zoom_h{height}.png')
        args = f"--tracks {ini_file} --region {region} "\
            "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
            f"--height {height} " \
            f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, my_tolerance)
        assert res is None, res
        os.remove(outfile.name)


def test_second_maf_withe():
    extension = '.png'
    ini_file = os.path.join(ROOT, "maf_withe.ini")

    for i, region in enumerate(['chr2:74,070,244-74,071,016',
                                'chr2:74,075,687-74,075,808'],
                               start=1):
        outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                     delete=False)
        expected_file = os.path.join(ROOT, f'master_maf_withe_region{i}.png')
        args = f"--tracks {ini_file} --region {region} "\
            "--trackLabelFraction 0.2 --width 38 --dpi 130 "\
            f"--outFileName {outfile.name}".split()
        pygenometracks.plotTracks.main(args)
        res = compare_images(expected_file,
                             outfile.name, tolerance)
        assert res is None, res
        os.remove(outfile.name)
    # Remove the index file:
    os.remove(os.path.join(ROOT, 'mm10_chr2_isl2_lessspe.maf.index'))


def test_first_maf_seq_hg19_pdf():
    extension = '.pdf'

    outfile = NamedTemporaryFile(suffix=extension, prefix='pyGenomeTracks_test_',
                                 delete=False)
    ini_file = os.path.join(ROOT, "first_maf_seq_hg19.ini")
    region = "9:128093930-128093970"
    expected_file = os.path.join(ROOT, 'master_first_maf_seq_hg19.pdf')
    # matplotlib compare on pdf will create a png next to it.
    # To avoid issues related to write in test_data folder
    # We copy the expected file into a temporary place
    new_expected_file = NamedTemporaryFile(suffix='.pdf',
                                           prefix='pyGenomeTracks_test_',
                                           delete=False)
    shutil.copy(expected_file, new_expected_file.name)
    expected_file = new_expected_file.name
    args = f"--tracks {ini_file} --region {region} "\
           "--trackLabelFraction 0.2 --width 38 --dpi 10 "\
           f"--outFileName {outfile.name}".split()
    pygenometracks.plotTracks.main(args)
    res = compare_images(expected_file,
                         outfile.name, tolerance)
    assert res is None, res
    os.remove(outfile.name)
    os.remove(expected_file.replace('.pdf', '_pdf.png'))
