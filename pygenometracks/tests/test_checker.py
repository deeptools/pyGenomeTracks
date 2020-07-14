# -*- coding: utf-8 -*-
"""
This test file will check each part of the
configuration file checker
which can raise an InputError
"""
import unittest
import pygenometracks.plotTracks
import os
from pygenometracks.utilities import InputError


ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "test_data")

test_tracks_1 = """
[vlines]
type = vlines
"""
with open(os.path.join(ROOT, "test_tracks_1.ini"), 'w') as fh:
    fh.write(test_tracks_1)

test_tracks_2 = """
[wrong file_type]
file_type = newFileTypeIJustInvented
"""
with open(os.path.join(ROOT, "test_tracks_2.ini"), 'w') as fh:
    fh.write(test_tracks_2)

test_tracks_3 = """
[undetermined section]
title = I do not know what
"""
with open(os.path.join(ROOT, "test_tracks_3.ini"), 'w') as fh:
    fh.write(test_tracks_3)

test_tracks_4 = """
[missing necessary option]
file_type = bed
"""
with open(os.path.join(ROOT, "test_tracks_4.ini"), 'w') as fh:
    fh.write(test_tracks_4)

test_tracks_5 = """
[problem with boolean]
file_type = bed
file = empty.bed
labels = 2
"""
with open(os.path.join(ROOT, "test_tracks_5.ini"), 'w') as fh:
    fh.write(test_tracks_5)

test_tracks_6 = """
[problem with float]
file_type = bed
file = empty.bed
max_value = a
"""
with open(os.path.join(ROOT, "test_tracks_6.ini"), 'w') as fh:
    fh.write(test_tracks_6)

test_tracks_7 = """
[problem with float range]
file_type = bed
file = empty.bed
fontsize = -2
"""
with open(os.path.join(ROOT, "test_tracks_7.ini"), 'w') as fh:
    fh.write(test_tracks_7)

test_tracks_8 = """
[problem with integer]
file_type = bed
file = empty.bed
gene_rows = 1a2
"""
with open(os.path.join(ROOT, "test_tracks_8.ini"), 'w') as fh:
    fh.write(test_tracks_8)

test_tracks_9 = """
[problem with integer range]
file_type = bed
file = empty.bed
gene_rows = -2
"""
with open(os.path.join(ROOT, "test_tracks_9.ini"), 'w') as fh:
    fh.write(test_tracks_9)

test_tracks_10 = """
[inexisting file]
file_type = bed
file = FileWhichDoesNotExists.txt
"""
with open(os.path.join(ROOT, "test_tracks_10.ini"), 'w') as fh:
    fh.write(test_tracks_10)

test_tracks_11 = """
[unguessable file]
file = FileWhichDoesNotExists.unknownextension
"""
with open(os.path.join(ROOT, "test_tracks_11.ini"), 'w') as fh:
    fh.write(test_tracks_11)

test_tracks_12 = """
[operation without second_file]
file = bigwig_chrx_2e6_5e6.bw
operation = file - second_file
"""
with open(os.path.join(ROOT, "test_tracks_12.ini"), 'w') as fh:
    fh.write(test_tracks_12)

test_tracks_12b = """
[operation without second_file]
file = bedgraph2_X_2.5e6_3.5e6.bdg
operation = file - second_file
"""
with open(os.path.join(ROOT, "test_tracks_12b.ini"), 'w') as fh:
    fh.write(test_tracks_12)

test_tracks_13 = """
[operation with transform]
file = bigwig_chrx_2e6_5e6.bw
operation = file + 1
transform = log
"""
with open(os.path.join(ROOT, "test_tracks_13.ini"), 'w') as fh:
    fh.write(test_tracks_13)

test_tracks_13b = """
[operation with transform]
file = bedgraph2_X_2.5e6_3.5e6.bdg
operation = file + 1
transform = log
"""
with open(os.path.join(ROOT, "test_tracks_13b.ini"), 'w') as fh:
    fh.write(test_tracks_13b)

test_tracks_14 = """
[invalid operation]
file = bigwig_chrx_2e6_5e6.bw
operation = file + a
"""
with open(os.path.join(ROOT, "test_tracks_14.ini"), 'w') as fh:
    fh.write(test_tracks_14)

test_tracks_14b = """
[invalid operation]
file = bedgraph2_X_2.5e6_3.5e6.bdg
operation = file + a
"""
with open(os.path.join(ROOT, "test_tracks_14b.ini"), 'w') as fh:
    fh.write(test_tracks_14b)

test_tracks_15 = """
[invalid operation with 2 files]
file = bigwig_chrx_2e6_5e6.bw
second_file = bigwig_chrx_2e6_5e6.bw
operation = file + a * second_file
"""
with open(os.path.join(ROOT, "test_tracks_15.ini"), 'w') as fh:
    fh.write(test_tracks_15)

test_tracks_15b = """
[invalid operation with 2 files]
file = bedgraph2_X_2.5e6_3.5e6.bdg
second_file = bedgraph2_X_2.5e6_3.5e6.bdg
operation = file + a * second_file
"""
with open(os.path.join(ROOT, "test_tracks_15b.ini"), 'w') as fh:
    fh.write(test_tracks_15b)


test_tracks_16 = """
[x-axis]

[bedgraph]
file = bedgraph2_X_2.5e6_3.5e6.bdg
overlay_previous = anything
"""
with open(os.path.join(ROOT, "test_tracks_16.ini"), 'w') as fh:
    fh.write(test_tracks_16)


test_tracks_17 = """
[x-axis]

[vlines]
file = tad_classification.bed
type = vlines
line_width = a
"""
with open(os.path.join(ROOT, "test_tracks_17.ini"), 'w') as fh:
    fh.write(test_tracks_17)


class TestCheckerMethods(unittest.TestCase):

    def test_vline_without_file(self):
        """
        This test check that if you provide a section with
        type = vlines but you do not provide any file
        you will have an error with a message containing
        there is no file
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_1.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("there is no file" in str(context.exception))
        os.remove(ini_file)

    def test_missing_file_type(self):
        """
        This test check that if you provide a section with
        a file_type which is not part of the tracks class
        you will have an error with a message containing
        file_type ... does not exists
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_2.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("file_type newFileTypeIJustInvented does not exists" in
               str(context.exception))
        os.remove(ini_file)

    def test_unguessable_file_type_no_file(self):
        """
        This test check that if you provide a section with
        no file_type and no file which is not a x-axis nor a spacer
        you will have an error with a message containing
         there is no file_type nor file
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_3.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert(" there is no file_type nor file " in
               str(context.exception))
        os.remove(ini_file)

    def test_missing_necessary_option(self):
        """
        This test check that if you provide a section with
        a track which has a necessary option
        (for example a bed track needs a file) but
        this necessary option is missing
        you will have an error with a message containing
        the necessary property
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_4.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("the necessary property" in str(context.exception))
        os.remove(ini_file)

    def test_boolean(self):
        """
        This test check that if you provide
        a number where you should have a boolean
        you will have an error with a message containing
        boolean
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_5.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("boolean" in str(context.exception))
        os.remove(ini_file)

    def test_float(self):
        """
        This test check that if you provide
        a letter where you should have a float
        you will have an error with a message containing
        float
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_6.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("float" in str(context.exception))
        os.remove(ini_file)

    def test_float_limit(self):
        """
        This test check that if you provide
        a float which is not in the expected range:
        for example a negative fontsize
        you will have an error with a message containing
        whereas it should be between
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_7.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas it should be between" in str(context.exception))
        os.remove(ini_file)

    def test_integer(self):
        """
        This test check that if you provide
        a mix of letter and numbers where you should have a integer
        you will have an error with a message containing
        integer
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_8.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("integer" in str(context.exception))
        os.remove(ini_file)

    def test_integer_limit(self):
        """
        This test check that if you provide
        a number which is not in the expected range:
        for example a negative gene_rows
        you will have an error with a message containing
        whereas it should be between
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_9.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas it should be between" in str(context.exception))
        os.remove(ini_file)

    def test_file_missing(self):
        """
        This test check that if you provide
        a file that does not exists
        you will have an error with a message containing
        File in section ... not found
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_10.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("File in section [1. [inexisting file]] not found" in
               str(context.exception))
        os.remove(ini_file)

    def test_unguessable_file_type(self):
        """
        This test check that if you provide
        a section with no file_type and that
        the extension of the file do not match
        the extensions of the tracks class
        you will have an error with a message containing
        can not identify file type
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_11.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("can not identify file type" in str(context.exception))
        os.remove(ini_file)

    def test_operation_bw_without_second_file(self):
        """
        This test check that if you provide
        a second_file if second_file is part of the
        operation.
        """
        outfile_name = "test.png"
        for suf in ['', 'b']:
            ini_file = os.path.join(ROOT, f"test_tracks_12{suf}.ini")
            region = "X:3000000-3300000"
            args = f"--tracks {ini_file} --region {region} "\
                   f"--outFileName {outfile_name}".split()
            with self.assertRaises(InputError) as context:
                pygenometracks.plotTracks.main(args)

            assert("requires to set the parameter second_file" in
                   str(context.exception))
            os.remove(ini_file)

    def test_operation_with_transform(self):
        """
        This test check that if you do not provide
        both an operation and a transform.
        """
        outfile_name = "test.png"
        for suf in ['', 'b']:
            ini_file = os.path.join(ROOT, f"test_tracks_13{suf}.ini")
            region = "X:3000000-3300000"
            args = f"--tracks {ini_file} --region {region} "\
                   f"--outFileName {outfile_name}".split()
            with self.assertRaises(InputError) as context:
                pygenometracks.plotTracks.main(args)

            assert("'operation' and 'transform' cannot be set at the same time."
                   in str(context.exception))
            os.remove(ini_file)

    def test_invalid_operation(self):
        """
        This test check that if you give an invalid operation
        it will fail.
        """
        outfile_name = "test.png"
        for suf in ['', 'b']:
            ini_file = os.path.join(ROOT, f"test_tracks_14{suf}.ini")
            region = "X:3000000-3300000"
            args = f"--tracks {ini_file} --region {region} "\
                   f"--outFileName {outfile_name}".split()
            with self.assertRaises(Exception) as context:
                pygenometracks.plotTracks.main(args)

            assert("could not be computed"
                   in str(context.exception))
            os.remove(ini_file)

    def test_invalid_operation2(self):
        """
        This test check that if you give an invalid operation
        it will fail.
        """
        outfile_name = "test.png"
        for suf in ['', 'b']:
            ini_file = os.path.join(ROOT, f"test_tracks_15{suf}.ini")
            region = "X:3000000-3300000"
            args = f"--tracks {ini_file} --region {region} "\
                   f"--outFileName {outfile_name}".split()
            with self.assertRaises(Exception) as context:
                pygenometracks.plotTracks.main(args)

            assert("could not be computed"
                   in str(context.exception))
            os.remove(ini_file)

    def test_wrong_overlay_previous(self):
        """
        This test check that if you provide
        an overlay_previous which is not no, yes, share-y
        you will have an error with a message containing
        Possible options are no, yes, share-y
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_16.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
               f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("Possible options are no, yes, share-y" in str(context.exception))
        os.remove(ini_file)

    def test_wrong_line_width_vlines(self):
        """
        This test check that if you provide
        an line_width which is not float
        you will have an error with a message containing
        whereas we should have a float
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "test_tracks_17.ini")
        region = "X:3000000-3300000"
        args = f"--tracks {ini_file} --region {region} "\
               f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas we should have a float" in str(context.exception))
        os.remove(ini_file)


class TestInputRegionMethods(unittest.TestCase):

    def test_region_format_onlychr(self):
        """
        This test check that you get an error
        if the region format is not good
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "empty.ini")
        region = "X"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("is not valid, it should be chr:start-end" in
               str(context.exception))

    def test_region_format_chr_only_start(self):
        """
        This test check that you get an error
        if the region format is not good
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "empty.ini")
        region = "X:12"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(Exception) as context:
            pygenometracks.plotTracks.main(args)

        assert("is not valid, it should be chr:start-end" in
               str(context.exception))

    def test_region_format_invalid_start(self):
        """
        This test check that you get an error
        if the region format is not good
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "empty.ini")
        for region in ['X:a-12', 'X:-3']:
            args = f"--tracks {ini_file} --region {region} "\
                   f"--outFileName {outfile_name}".split()
            with self.assertRaises(InputError) as context:
                pygenometracks.plotTracks.main(args)

            assert("is not valid, it should be chr:start-end" in
                   str(context.exception))

    def test_region_format_invalid_end(self):
        """
        This test check that you get an error
        if the region format is not good
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "empty.ini")
        region = "X:12-1d"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("is not valid, it should be chr:start-end" in
               str(context.exception))

    def test_region_format_start_smaller_end(self):
        """
        This test check that you get an error
        if the region format is not good
        """
        outfile_name = "test.png"
        ini_file = os.path.join(ROOT, "empty.ini")
        region = "X:12-1"
        args = f"--tracks {ini_file} --region {region} "\
            f"--outFileName {outfile_name}".split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("end is larger" in str(context.exception))
