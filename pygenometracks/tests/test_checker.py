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
                    "test_data", "uncorrect_ini_files")

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


class TestCheckerMethods(unittest.TestCase):

    def test_vline_without_file(self):
        """
        This test check that if you provide a section with
        type = vlines but you do not provide any file
        you will have an error with a message containing
        there is no file
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_1.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("there is no file" in str(context.exception))

    def test_missing_file_type(self):
        """
        This test check that if you provide a section with
        a file_type which is not part of the tracks class
        you will have an error with a message containing
        file_type ... does not exists
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_2.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("file_type newFileTypeIJustInvented does not exists" in
               str(context.exception))

    def test_unguessable_file_type_no_file(self):
        """
        This test check that if you provide a section with
        no file_type and no file which is not a x-axis nor a spacer
        you will have an error with a message containing
         there is no file_type nor file
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_3.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert(" there is no file_type nor file " in
               str(context.exception))

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
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_4.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("the necessary property" in str(context.exception))

    def test_boolean(self):
        """
        This test check that if you provide
        a number where you should have a boolean
        you will have an error with a message containing
        boolean
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_5.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("boolean" in str(context.exception))

    def test_float(self):
        """
        This test check that if you provide
        a letter where you should have a float
        you will have an error with a message containing
        float
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_6.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("float" in str(context.exception))

    def test_float_limit(self):
        """
        This test check that if you provide
        a float which is not in the expected range:
        for example a negative fontsize
        you will have an error with a message containing
        whereas it should be between
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_7.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas it should be between" in str(context.exception))

    def test_integer(self):
        """
        This test check that if you provide
        a mix of letter and numbers where you should have a integer
        you will have an error with a message containing
        integer
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_8.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("integer" in str(context.exception))

    def test_integer_limit(self):
        """
        This test check that if you provide
        a number which is not in the expected range:
        for example a negative gene_rows
        you will have an error with a message containing
        whereas it should be between
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_9.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas it should be between" in str(context.exception))

    def test_file_missing(self):
        """
        This test check that if you provide
        a file that does not exists
        you will have an error with a message containing
        File in section ... not found
        """
        outfile_name = "test.png"
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_10.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("File in section [1. [inexisting file]] not found" in
               str(context.exception))

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
        args = "--tracks {0} --region X:3000000-3300000 " \
            "--outFileName {1}" \
            "".format(os.path.join(ROOT, "test_tracks_11.ini"),
                      outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("can not identify file type" in str(context.exception))
