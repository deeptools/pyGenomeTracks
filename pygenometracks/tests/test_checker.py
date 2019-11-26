# -*- coding: utf-8 -*-
import unittest
import pygenometracks.plotTracks
import os
from pygenometracks.tracksClass import InputError


ROOT = os.path.dirname(os.path.abspath(__file__)) + \
    "/test_data/uncorrect_ini_files/"

test_tracks_1 = """
[vlines]
type = vlines
"""
with open(ROOT + "test_tracks_1.ini", 'w') as fh:
    fh.write(test_tracks_1)

test_tracks_2 = """
[wrong file_type]
file_type = newFileTypeIJustInvented
"""
with open(ROOT + "test_tracks_2.ini", 'w') as fh:
    fh.write(test_tracks_2)

test_tracks_3 = """
[undetermined section]
title = I do not know what
"""
with open(ROOT + "test_tracks_3.ini", 'w') as fh:
    fh.write(test_tracks_3)

test_tracks_4 = """
[missing necessary option]
file_type = bed
"""
with open(ROOT + "test_tracks_4.ini", 'w') as fh:
    fh.write(test_tracks_4)

test_tracks_5 = """
[problem with boolean]
file_type = bed
file = empty.bed
labels = 2
"""
with open(ROOT + "test_tracks_5.ini", 'w') as fh:
    fh.write(test_tracks_5)

test_tracks_6 = """
[problem with float]
file_type = bed
file = empty.bed
max_value = a
"""
with open(ROOT + "test_tracks_6.ini", 'w') as fh:
    fh.write(test_tracks_6)

test_tracks_7 = """
[problem with float range]
file_type = bed
file = empty.bed
fontsize = -2
"""
with open(ROOT + "test_tracks_7.ini", 'w') as fh:
    fh.write(test_tracks_7)

test_tracks_8 = """
[problem with integer]
file_type = bed
file = empty.bed
gene_rows = 1a2
"""
with open(ROOT + "test_tracks_8.ini", 'w') as fh:
    fh.write(test_tracks_8)

test_tracks_9 = """
[problem with integer range]
file_type = bed
file = empty.bed
gene_rows = -2
"""
with open(ROOT + "test_tracks_9.ini", 'w') as fh:
    fh.write(test_tracks_9)

test_tracks_10 = """
[inexisting file]
file_type = bed
file = FileWhichDoesNotExists.txt
"""
with open(ROOT + "test_tracks_10.ini", 'w') as fh:
    fh.write(test_tracks_10)

test_tracks_11 = """
[unguessable file]
file = FileWhichDoesNotExists.unknownextension
"""
with open(ROOT + "test_tracks_11.ini", 'w') as fh:
    fh.write(test_tracks_11)


class TestCheckerMethods(unittest.TestCase):

    def test_vline_without_file(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_1.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("there is no file" in str(context.exception))

    def test_missing_file_type(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_2.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("file_type newFileTypeIJustInvented does not exists" in
               str(context.exception))

    def test_unguessable_file_type_no_file(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_3.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert(" there is no file_type nor file " in
               str(context.exception))

    def test_missing_necessary_option(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_4.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("the necessary property" in str(context.exception))

    def test_boolean(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_5.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("boolean" in str(context.exception))

    def test_float(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_6.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("float" in str(context.exception))

    def test_float_limit(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_7.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas it should be between" in str(context.exception))

    def test_integer(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_8.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("integer" in str(context.exception))

    def test_integer_limit(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_9.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("whereas it should be between" in str(context.exception))

    def test_file_missing(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_10.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("File in section [1. [inexisting file]] not found" in
               str(context.exception))

    def test_unguessable_file_type(self):
        outfile_name = "test.png"
        args = "--tracks {0}/test_tracks_11.ini --region X:3000000-3300000 " \
            "--outFileName  {1}".format(ROOT, outfile_name).split()
        with self.assertRaises(InputError) as context:
            pygenometracks.plotTracks.main(args)

        assert("can not identify file type" in str(context.exception))
