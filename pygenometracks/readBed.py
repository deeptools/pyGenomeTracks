# -*- coding: utf-8 -*-
from __future__ import division
from past.builtins import map

import sys
import collections
from .utilities import to_string


class ReadBed(object):
    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Example:
    bed = readBed(open("file.bed", 'r'))
    for interval in bed:
        print bed['start']

    """

    def __init__(self, file_handle):
        """
        :param file_handle: file handle
        :return:
        """

        self.file_type = None
        self.file_handle = file_handle
        self.line_number = 0
        # guess file type
        fields = self.get_no_comment_line()
        fields = to_string(fields)
        fields = fields.split('\t')

        self.guess_file_type(fields)
        self.file_handle.seek(0)
        self.prev_chrom = None
        self.prev_start = -1
        self.prev_line = None

        # list of bed fields
        self.fields = ['chromosome', 'start', 'end',
                       'name', 'score', 'strand',
                       'thick_start', 'thick_end',
                       'rgb', 'block_count',
                       'block_sizes', 'block_starts']

        if self.file_type == 'bed12':
            self.BedInterval = collections.namedtuple('BedInterval', self.fields)
        elif self.file_type == 'bed9':
            self.BedInterval = collections.namedtuple('BedInterval', self.fields[:9])
        else:
            self.BedInterval = collections.namedtuple('BedInterval', self.fields[:6])

    def __iter__(self):
        return self

    def get_no_comment_line(self):
        """
        Skips comment lines starting with '#'
        "track" or "browser" in the bed files
        :return:
        """
        line = next(self.file_handle)
        line = to_string(line)
        if line.startswith("#") or line.startswith("track") or \
                line.startswith("browser") or line.strip() == '':
            line = self.get_no_comment_line()

        self.line_number += 1
        return line

    def guess_file_type(self, line_values):
        """try to guess type of bed file by counting the fields
        """
        if len(line_values) == 3:
            self.file_type = 'bed3'
        elif len(line_values) == 4:
            self.file_type = 'bedgraph'
        elif len(line_values) == 6:
            self.file_type = 'bed6'
        elif len(line_values) == 12:
            self.file_type = 'bed12'
        elif len(line_values) == 9:
            # this is a case where a specific color is encoded in the 10 field of the bed file
            self.file_type = 'bed9'
        elif len(line_values) > 6:
            # assume bed6
            self.file_type = 'bed6'
            sys.stderr.write("Number of fields in BED file is not standard. Assuming bed6\n")
        else:
            # assume bed3
            self.file_type = 'bed3'
            sys.stderr.write("Number of fields in BED file is not standard. Assuming bed3\n")
        return self.file_type

    def next(self):
        """
        :return: bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)
        if self.prev_chrom == bed.chromosome:
            assert self.prev_start <= bed.start, \
                "Bed file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self.file_handle.name, self.prev_line, line)

        self.prev_chrom = bed.chromosome
        self.prev_start = bed.start
        self.prev_line = line

        return bed

    def __next__(self):
        """
        :return: bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)
        if self.prev_chrom == bed.chromosome:
            assert self.prev_start <= bed.start, \
                "Bed file not sorted. Please use a sorted bed file.\n" \
                "File: {}\n" \
                "Previous line: {}\n Current line{} ".format(self.file_handle.name, self.prev_line, line)

        self.prev_chrom = bed.chromosome
        self.prev_start = bed.start
        self.prev_line = line

        return bed

    def get_bed_interval(self, bed_line):
        r"""
        Processes each bed line from a bed file, casts the values and returns
        a namedtuple object

        >>> bed_line="chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,100\t20,200,700"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed = bed_f.get_bed_interval(bed_line)
        >>> bed.chromosome
        'chr1'
        >>> bed.block_starts
        [20, 200, 700]

        >>> bed_line="chr2\t0\t1000\tgene_1\t0.5\t-\n"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed_f.get_bed_interval(bed_line)
        BedInterval(chromosome='chr2', start=0, end=1000, name='gene_1', score=0.5, strand='-')
        """

        line_data = bed_line.strip()
        line_data = to_string(line_data)
        line_data = line_data.split("\t")

        if self.file_type == 'bed12':
            assert len(line_data) == 12, "File type detected is bed12 but line {}: {} does " \
                                         "not have 12 fields.".format(self.line_number, bed_line)

        elif self.file_type == 'bed3':
            assert len(line_data) == 3, "File type detected is bed3 but line {}: {} does " \
                "not have 3 fields.".format(self.line_number, bed_line)

        elif self.file_type == 'bed6':
            assert len(line_data) == 6, "File type detected is bed6 but line {}: {} does " \
                "not have 6 fields.".format(self.line_number, bed_line)
        line_values = []
        for idx, r in enumerate(line_data):
            # first field is always chromosome/contig name
            # and should be cast as a string
            # same for field 3 (name)
            if idx in [0, 3]:
                line_values.append(r)
            # check field strand
            elif idx == 5:
                if r not in ['+', '-', '.']:
                    if r == '1':
                        r = '+'
                    elif r == '-1':
                        r = '-'
                    else:
                        sys.stderr.write("*Warning, invalid strand value found {} for line #{}:\n{}\n "
                                         "Setting strand to '.'\n".format(r, bed_line, self.line_number))
                        r = '.'
                line_values.append(r)

            elif idx in [1, 2, 6, 7, 9]:
                # start and end fields must be integers, same for thichStart(6),
                # and thickEnd(7) and blockCount(9) fields
                try:
                    line_values.append(int(r))
                except ValueError:
                    sys.stderr.write("Value: {} in field {} at line {} is not an integer\n".format(r, idx + 1,
                                                                                                   self.line_number))
                    return dict()
            # check item rgb
            elif idx == 8:
                r = to_string(r)
                rgb = r.split(",")
                if len(rgb) == 3:
                    try:
                        r = map(int, rgb)
                    except ValueError as detail:
                        sys.stderr.write("Error reading line: #{}. The rgb field {} is not "
                                         "valid.\nError message: {}\n".format(self.line_number, r, detail))
                line_values.append(r)

            elif idx in [10, 11]:
                # this are the block sizes and block start positions
                r = to_string(r)
                r_parts = r.split(',')
                try:
                    r = [int(x) for x in r_parts if x != '']
                except ValueError as detail:
                    sys.stderr.write("Error reading line #{}. The block field {} is not "
                                     "valid.\nError message: {}\n".format(self.line_number, r, detail))
                line_values.append(r)

            else:
                try:
                    tmp = float(r)
                except ValueError:
                    tmp = r
                except TypeError:
                    tmp = r
                line_values.append(tmp)

        assert line_values[2] > line_values[1], \
            "Start position larger or equal than end for line #{}:\n{}\n".format(self.line_number,
                                                                                 bed_line)

        if self.file_type == 'bed3':
            line_values = line_values[0:3]
            # in case of a bed3, the id, score and strand
            # values are added as ".", 0, "." respectively
            line_values.extend([".", 0, "."])
        elif self.file_type == 'bed6':
            line_values = line_values[0:6]

        return self.BedInterval._make(line_values)
