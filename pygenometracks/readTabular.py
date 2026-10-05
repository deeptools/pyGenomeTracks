# -*- coding: utf-8 -*-
import collections

from .utilities import opener, to_string


class ReadTabular():
    """
    Reads a file with multiple named columns with
    different type (string/integers..., not a bedgraph-like)
    and generate an iterator
    """
    # List of fields
    fields = None
    comments_keywords = ["#"]

    def __init__(self, file_path):
        """
        :param file_path: file path to open and generate a iterator
        :type file_path: str
        """
        self.file_handle = opener(file_path)
        self.line_number = 0
        self.length = self.get_nb_records()
        self.adjust_fields()
        self.Record = collections.namedtuple('Record', self.used_fields)

    def adjust_fields(self):
        self.used_fields = self.fields

    def __iter__(self):
        return self

    def is_comment(self, line_string):
        return any([line_string.startswith(c) for c in self.comments_keywords]) or line_string.strip() == ''

    def get_no_comment_line(self):
        """
        Skips comment lines starting with a
        comment keyword
        :return:
        """
        line = next(self.file_handle)
        self.line_number += 1
        line = to_string(line)
        if self.is_comment(line):
            line = self.get_no_comment_line()
        return line

    def get_nb_records(self):
        n = 0
        for line in self.file_handle:
            if not self.is_comment(to_string(line)):
                n += 1
        # Put back the handler to start
        self.file_handle.seek(0)
        return n

    def __next__(self):
        """
        :return: Record object
        """
        line = self.get_no_comment_line()
        return self.get_record(line)

    def get_line_data(self, line):
        line_data = line.strip()
        line_data = to_string(line_data)
        return line_data.split("\t")

    def get_record(self, line):
        return self.Record._make(self.get_line_data(line))
