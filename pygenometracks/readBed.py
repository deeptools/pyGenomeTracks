# -*- coding: utf-8 -*-
import sys
import collections
from .utilities import to_string, InputError


class ReadBed(object):
    """
    Reads a bed file. Based on the number of fields
    it tries to guess the type of bed file used. Current options
    are bed3, bed6 and bed12

    Example:
    bed = ReadBed(open("file.bed", 'r'))
    for interval in bed:
        print interval.start

    """

    def __init__(self, file_handle):
        """
        :param file_handle: file handle
        :return:
        """

        # file_type can be bed6, bed8, bed9, or bed12
        self.file_type = None
        # The number of fields to read at each line
        # Can be 3 to 12
        self.fields_to_read = 12
        self.file_handle = file_handle
        self.line_number = 0
        # guess file type
        try:
            fields = self.get_no_comment_line()
        except StopIteration:
            self.file_type = 'bed6'
        else:
            self.get_bed_interval(fields, is_first_line=True)
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

        if self.fields_to_read <= 6:
            self.BedInterval = collections.namedtuple('BedInterval',
                                                      self.fields[:6])
        else:
            self.BedInterval = collections.namedtuple('BedInterval',
                                                      self.fields[:self.fields_to_read])

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

    def get_bed_interval(self, bed_line, is_first_line=False):
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

        if not is_first_line:
            if self.file_type != 'bed6':
                # When bed6 you can have less fields in one row
                # because there are default values
                assert len(line_data) >= self.fields_to_read, \
                    "File type detected is {} but line {}: {} does " \
                    "not have {} fields.".format(self.file_type,
                                                 self.line_number,
                                                 bed_line,
                                                 self.fields_to_read)
            line_data = line_data[:self.fields_to_read]

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
                        sys.stderr.write("*Warning, invalid strand value"
                                         "found {} for line #{}:\n{}\n "
                                         "Setting strand to '.'"
                                         "\n".format(r, bed_line,
                                                     self.line_number))
                        r = '.'
                line_values.append(r)

            elif idx in [1, 2]:
                # start and end fields must be integers
                try:
                    line_values.append(int(r))
                except ValueError:
                    raise InputError("Value: {} in field {} at line {}"
                                     " is not an integer. This is "
                                     "probably not a bed file."
                                     "\n".format(r, idx + 1,
                                                 self.line_number))
            elif idx in [6, 7, 9]:
                # thichStart(6), thickEnd(7) and blockCount(9) fields
                # Should be integer, if they are not we change the bed type
                try:
                    line_values.append(int(r))
                except ValueError:
                    if is_first_line:
                        if idx == 9:
                            self.file_type = 'bed8'
                            self.fields_to_read = 8
                        else:
                            self.file_type = 'bed6'
                            self.fields_to_read = 6
                        sys.stderr.write("Value: {} in field {}"
                                         " is not an integer"
                                         "\n Only the first {} fields will"
                                         " be used.\n".format(r, idx + 1,
                                                              self.fields_to_read))
                        break
                    else:
                        default_value = 0 if (idx == 9) else line_values[1]
                        sys.stderr.write("Value: {} in field {} at line {}"
                                         " is not an integer"
                                         "\n {} will be used.\n"
                                         "".format(r, idx + 1,
                                                   self.line_number,
                                                   default_value))
            # check item rgb
            elif idx == 8:
                passed = True
                try:
                    line_values.append(int(r))
                except ValueError:
                    r = to_string(r)
                    rgb = r.split(",")
                    if len(rgb) == 3:
                        try:
                            r = list(map(int, rgb))
                        except ValueError:
                            passed = False
                        else:
                            line_values.append(r)
                    else:
                        passed = False
                if not passed:
                    if is_first_line:
                        sys.stderr.write("Warning: "
                                         "The rgb field {} is not "
                                         "valid.\nOnly the first 8 fields"
                                         " will be used.\n"
                                         "".format(self.line_number, r))
                        self.file_type = 'bed8'
                        self.fields_to_read = 8
                        break
                    else:
                        sys.stderr.write("Warning: reading line: #{}. "
                                         "The rgb field {} is not "
                                         "valid.\n0"
                                         " will be used.\n"
                                         "".format(self.line_number, r))
                        line_values.append(0)

            elif idx in [10, 11]:
                # this are the block sizes and block start positions
                r = to_string(r)
                r_parts = r.split(',')
                try:
                    r = [int(x) for x in r_parts if x != '']
                except ValueError as detail:
                    if is_first_line:
                        sys.stderr.write("Warning: "
                                         "The block field {} is not "
                                         "valid.\nError message: {}"
                                         "\nOnly the first 9 fields"
                                         " will be used.\n"
                                         "".format(r,
                                                   detail))
                        self.file_type = 'bed9'
                        self.fields_to_read = 9
                        break
                    else:
                        sys.stderr.write("Warning: reading line #{}, "
                                         "the block field {} is not "
                                         "valid.\nError message: {}"
                                         "\nNo block will be used.\n"
                                         "".format(self.line_number, r,
                                                   detail))
                        line_values[9] = 1
                        line_values[10] = line_values[1]
                        line_values[11] = line_values[2] - line_values[1]
                        break
                else:
                    line_values.append(r)

            elif idx == 4:
                try:
                    tmp = float(r)
                except (ValueError, TypeError) as detail:
                    if is_first_line:
                        sys.stderr.write("Warning: "
                                         "The block field 5 (score) is not "
                                         "valid: {}.\nError message: {}"
                                         "\nOnly the first 4 fields "
                                         "will be used.\n".format(r,
                                                                  detail))
                        self.file_type = 'bed6'
                        self.fields_to_read = 4
                        break
                    else:
                        sys.stderr.write("Warning: reading line #{}, "
                                         "the block field 5 (score) is not "
                                         "valid: {}.\nError message: {}"
                                         "\n0 will be used.\n"
                                         "".format(self.line_number, r,
                                                   detail))
                        line_values.append(0.)
                else:
                    line_values.append(tmp)

        if is_first_line:
            if self.file_type is None:
                self.fields_to_read = max([i for i in [3, 4, 5, 6, 8, 9, 12] if i <= len(line_values)])
                if self.fields_to_read <= 6:
                    self.file_type = 'bed6'
                else:
                    self.file_type = 'bed{}'.format(self.fields_to_read)
            return()

        assert line_values[2] > line_values[1], \
            "Start position larger or equal than end" \
            " for line #{}:\n{}\n".format(self.line_number,
                                          bed_line)
        if len(line_values) == 12:
            check_bed12(line_values, self.line_number, bed_line)

        if len(line_values) < 6:
            assert len(line_values) > 2, \
                "The number of field is less than 3.\n" \
                "This is not a bed file.\n" \
                "File: {}\n" \
                "Current line: {}\n".format(self.file_handle.name,
                                            line_values)
            # If there is less than 6 fields, the default values will be added
            default = [".", 0, "."]
            line_values = [line_values[i] if i < len(line_values)
                           else default[i - 3]
                           for i in range(6)]

        return self.BedInterval._make(line_values)


def check_bed12(line_values, line_number, bed_line):
    block_counts = line_values[9]
    block_sizes = line_values[10]
    block_relative_starts = line_values[11]
    assert len(block_sizes) == block_counts, \
        "The number of blocks: {} does not correspond to" \
        "the number of blocks sizes: {}\nline #{}:\n" \
        "{}".format(block_counts, str(block_sizes),
                    line_number, bed_line)
    assert len(block_relative_starts) == block_counts, \
        "The number of blocks: {} does not correspond to" \
        "the number of blocks relative starts: {}\nline #{}:\n" \
        "{}".format(block_counts, str(block_relative_starts),
                    line_number, bed_line)
    for i in range(block_counts):
        block_start = line_values[1] + block_relative_starts[i]
        block_end = block_start + block_sizes[i]
        assert block_start <= line_values[2], \
            "The block number {} of line {} has a starting position " \
            "greater than the end of the feature:\n{}The" \
            " 12th field of a bed12 should contains relative start" \
            " positions of blocks.".format(i, line_number, bed_line)
        assert block_end <= line_values[2], \
            "The block number {} of line {} has an ending position " \
            "greater than the end of the feature:\n{}The" \
            " 12th field of a bed12 should contains relative start" \
            " positions of blocks and the 11th field should contains" \
            " the length of each block.".format(i, line_number, bed_line)
