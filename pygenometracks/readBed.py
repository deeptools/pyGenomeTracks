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

    def __next__(self):
        """
        :return: bedInterval object
        """
        line = self.get_no_comment_line()

        bed = self.get_bed_interval(line)

        return bed

    def get_bed_interval(self, bed_line, is_first_line=False):
        r"""
        Processes each bed line from a bed file, casts the values and returns
        a namedtuple object

        >>> bed_line="chr1\t0\t1000\tgene_1\t0.5\t-\t0\t1000\t0\t3\t10,20,300\t0,200,700"
        >>> with open('/tmp/test.bed', 'w') as fh:
        ...     foo = fh.write(bed_line)
        >>> bed_f = ReadBed(open('/tmp/test.bed','r'))
        >>> bed = bed_f.get_bed_interval(bed_line)
        >>> bed.chromosome
        'chr1'
        >>> bed.block_starts
        [0, 200, 700]

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
        if is_first_line:
            if len(line_data) == 1:
                if line_data[0].startswith("{\\rtf"):
                    raise InputError("The file is a rtf file."
                                     " Please save it as plain text.")
                else:
                    raise InputError("Only one field detected, you may use"
                                     " a bed delimited by space. This format "
                                     "is not supported by pyGenomeTracks.")
        else:
            if self.file_type != 'bed6':
                # When bed6 you can have less fields in one row
                # because there are default values
                assert len(line_data) >= self.fields_to_read, \
                    f"File type detected is {self.file_type} but line" \
                    f" {self.line_number}: {bed_line} does " \
                    f"not have {self.fields_to_read} fields."
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
                                         f"found {r} for line #{bed_line}:\n{self.line_number}\n "
                                         "Setting strand to '.'"
                                         "\n")
                        r = '.'
                line_values.append(r)

            elif idx in [1, 2]:
                # start and end fields must be integers
                try:
                    line_values.append(int(r))
                except ValueError:
                    raise InputError(f"Value: {r} in field {idx + 1} at line"
                                     f" {self.line_number}"
                                     " is not an integer. This is "
                                     "probably not a bed file."
                                     "\n")
            elif idx in [6, 7]:
                # thichStart(6) and thickEnd(7) fields
                # Should be integer, if they are not we change the bed type
                try:
                    line_values.append(int(r))
                except ValueError:
                    if is_first_line:
                        self.file_type = 'bed6'
                        self.fields_to_read = 6
                        sys.stderr.write(f"Value: {r} in field {idx + 1}"
                                         " is not an integer"
                                         "\n Only the first "
                                         f"{self.fields_to_read} fields will"
                                         " be used.\n")
                        break
                    else:
                        if idx == 6:
                            default_value = line_values[1]
                        else:  # idx is 7
                            default_value = line_values[6]
                        line_values.append(default_value)
                        sys.stderr.write(f"Value: {r} in field {idx + 1} at"
                                         f" line {self.line_number}"
                                         " is not an integer"
                                         f"\n {default_value} will be used.\n")
            elif idx == 9:
                # blockCount(9) field
                # Should be integer, if they are not we change the bed type
                try:
                    line_values.append(int(r))
                except ValueError:
                    if is_first_line:
                        self.file_type = 'bed8'
                        self.fields_to_read = 8
                        sys.stderr.write(f"Value: {r} in field {idx + 1}"
                                         " is not an integer"
                                         "\n Only the first "
                                         f"{self.fields_to_read} fields will"
                                         " be used.\n")
                        break
                    else:
                        sys.stderr.write("Warning: reading line "
                                         f"#{self.line_number}, "
                                         f"the block number {r} is not "
                                         "valid.\nNo block will be used.\n")
                        line_values.append(1)
                        line_values.append([line_values[2] - line_values[1]])
                        line_values.append([0])
                        break
            # check item rgb
            elif idx == 8:
                passed = True
                try:
                    # This is what happens in UCSC browser:
                    line_values.append([0, 0, int(r)])
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
                                         f"The rgb field {r} is not "
                                         "valid.\nOnly the first 8 fields"
                                         " will be used.\n")
                        self.file_type = 'bed8'
                        self.fields_to_read = 8
                        break
                    else:
                        sys.stderr.write("Warning: reading line: "
                                         f"#{self.line_number}. "
                                         f"The rgb field {r} is not "
                                         "valid.\n0,0,0"
                                         " will be used.\n")
                        line_values.append([0, 0, 0])

            elif idx in [10, 11]:
                # this are the block sizes and block start positions
                r = to_string(r)
                r_parts = r.split(',')
                try:
                    r = [int(x) for x in r_parts if x != '']
                except ValueError as detail:
                    if is_first_line:
                        sys.stderr.write("Warning: "
                                         f"The block field {r} is not "
                                         f"valid.\nError message: {detail}"
                                         "\nOnly the first 9 fields"
                                         " will be used.\n")
                        self.file_type = 'bed9'
                        self.fields_to_read = 9
                        break
                    else:
                        sys.stderr.write("Warning: reading line "
                                         f"#{self.line_number}, "
                                         f"the block field {r} is not "
                                         f"valid.\nError message: {detail}"
                                         "\nOne block will be used.\n")
                        line_values = line_values[:9]
                        line_values.append(1)
                        line_values.append([line_values[2] - line_values[1]])
                        line_values.append([0])
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
                                         f"valid: {r}.\nError message: {detail}"
                                         "\nOnly the first 4 fields "
                                         "will be used.\n")
                        self.file_type = 'bed6'
                        self.fields_to_read = 4
                        break
                    else:
                        sys.stderr.write("Warning: reading line "
                                         f"#{self.line_number}, "
                                         "the block field 5 (score) is not "
                                         f"valid: {r}.\nError message: {detail}"
                                         "\n0 will be used.\n")
                        line_values.append(0.)
                else:
                    line_values.append(tmp)

        if is_first_line:
            if self.file_type is None:
                self.fields_to_read = max([i for i in [3, 4, 5, 6, 8, 9, 12] if i <= len(line_values)])
                if self.fields_to_read <= 6:
                    self.file_type = 'bed6'
                else:
                    self.file_type = f'bed{self.fields_to_read}'
            return

        assert line_values[2] > line_values[1], \
            "Start position larger or equal than end" \
            f" for line #{self.line_number}:\n{bed_line}\n"

        if len(line_values) == 12:
            check_bed12(line_values, self.line_number, bed_line)

        if len(line_values) < 6:
            assert len(line_values) > 2, \
                "The number of field is less than 3.\n" \
                "This is not a bed file.\n" \
                f"File: {self.file_handle.name}\n" \
                f"Current line: {line_values}\n"
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
        f"The number of blocks: {block_counts} does not correspond to" \
        f"the number of blocks sizes: {str(block_sizes)}\nline #{line_number}:" \
        f"\n{bed_line}"
    assert len(block_relative_starts) == block_counts, \
        f"The number of blocks: {block_counts} does not correspond to" \
        f"the number of blocks relative starts: {str(block_relative_starts)}" \
        f"\nline #{line_number}:\n{bed_line}"
    for i in range(block_counts):
        block_start = line_values[1] + block_relative_starts[i]
        block_end = block_start + block_sizes[i]
        assert block_start <= line_values[2], \
            f"The block number {i} of line {line_number} has a starting" \
            f" position greater than the end of the feature:\n{bed_line}The" \
            " 12th field of a bed12 should contains relative start" \
            " positions of blocks."
        assert block_end <= line_values[2], \
            f"The block number {i} of line {line_number} has an ending position " \
            f"greater than the end of the feature:\n{bed_line}The" \
            " 12th field of a bed12 should contains relative start" \
            " positions of blocks and the 11th field should contains" \
            " the length of each block."
    assert min(block_relative_starts) == 0, \
        f"The blocks relative_starts of line\n{bed_line}\n" \
        "does not contain 0. BED blocks must span chromStart" \
        " to chromEnd."
    assert max([bstart + bsize for bstart, bsize in zip(block_relative_starts, block_sizes)]) == line_values[2] - line_values[1], \
        f"The blocks described in line\n{bed_line}\n" \
        "does not cover chromEnd. BED blocks must span chromStart" \
        " to chromEnd."
