import sys
import gzip
import numpy as np
from intervaltree import IntervalTree, Interval


class InputError(Exception):
    """Exception raised for errors in the input."""
    pass


def to_string(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        assert(sys.version_info[0] != 2)
#        if sys.version_info[0] == 2:
#            return str(s)
        return s.decode('ascii')
    if isinstance(s, list):
        return [to_string(x) for x in s]
    return s


def to_bytes(s):
    """
    Like toString, but for functions requiring bytes in python3
    """
    assert(sys.version_info[0] != 2)
#    if sys.version_info[0] == 2:
#        return s
    if isinstance(s, bytes):
        return s
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [to_bytes(x) for x in s]
    return s


def opener(filename):
    """
    Determines if a file is compressed or not
    """
    f = open(filename, 'rb')
    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f


def file_to_intervaltree(file_name):
    """
    converts a BED like file into a bx python interval tree
    :param file_name: string file name
    :return: interval tree dictionary. They key is the chromosome/contig name and the
    value is an IntervalTree. Each of the intervals have as 'value' the fields[3:] if any.
    """
    # iterate over a BED like file
    # saving the data into an interval tree
    # for quick retrieval
    file_h = opener(file_name)
    line_number = 0
    valid_intervals = 0
    prev_chrom = None
    prev_start = -1
    prev_line = None
    interval_tree = {}
    min_value = float('Inf')
    max_value = -float('Inf')

    for line in file_h.readlines():
        line_number += 1
        line = to_string(line)
        if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        try:
            chrom, start, end = fields[0:3]
        except Exception as detail:
            msg = "Error reading line: {}\nError message: {}".format(line_number, detail)
            raise InputError(msg)

        try:
            start = int(start)
        except ValueError as detail:
            msg = "Error reading line: {}. The start field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            raise InputError(msg)

        try:
            end = int(end)
        except ValueError as detail:
            msg = "Error reading line: {}. The end field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            raise InputError(msg)

        if prev_chrom == chrom:
            assert prev_start <= start, \
                "Bed file not sorted. Please use a sorted bed file.\n{}{} ".format(prev_line, line)

        if chrom not in interval_tree:
            interval_tree[chrom] = IntervalTree()

        value = None

        if len(fields) > 3:
            value = fields[3:]
            try:
                line_min = min(map(float, value))
                if line_min < min_value:
                    min_value = line_min

                line_max = max(map(float, value))
                if line_max > max_value:
                    max_value = line_max
            except ValueError:
                pass

        assert end > start, "Start position larger or equal than end for line\n{} ".format(line)

        interval_tree[chrom].add(Interval(start, end, value))
        valid_intervals += 1

    if valid_intervals == 0:
        sys.stderr.write("No valid intervals were found in file {}".format(file_name))
    file_h.close()

    return interval_tree, min_value, max_value


def plot_coverage(ax, x_values, score_list, plot_type, size, color,
                  negative_color, alpha):
    if plot_type == 'line':
        if color == negative_color:
            ax.plot(x_values, score_list, '-', linewidth=size, color=color,
                    alpha=alpha)
        else:
            import warnings
            warnings.warn('Line plots with a different negative color might not look pretty')
            pos_x_values = x_values.copy()
            pos_x_values[score_list < 0] = np.nan
            ax.plot(pos_x_values, score_list, '-', linewidth=size, color=color,
                    alpha=alpha)

            neg_x_values = x_values.copy()
            neg_x_values[score_list >= 0] = np.nan
            ax.plot(neg_x_values, score_list, '-', linewidth=size,
                    color=negative_color, alpha=alpha)

    elif plot_type == 'points':
        if color == negative_color:
            ax.plot(x_values, score_list, '.', markersize=size,
                    color=color,
                    alpha=alpha)
        else:
            pos_x_values = x_values.copy()
            pos_x_values[score_list < 0] = np.nan
            ax.plot(pos_x_values, score_list, '.',
                    markersize=size,
                    color=color,
                    alpha=alpha)
            neg_x_values = x_values.copy()
            neg_x_values[score_list >= 0] = np.nan
            ax.plot(neg_x_values, score_list, '.',
                    markersize=size,
                    color=negative_color,
                    alpha=alpha)
    else:
        if plot_type != 'fill':
            import warnings
            warnings.warn('The plot type was not part of known types '
                          '(fill, line, points) will be fill.')
        if color == negative_color:
            ax.fill_between(x_values, score_list, linewidth=0.1,
                            color=color,
                            facecolor=color,
                            alpha=alpha)
        else:
            pos_x_values = x_values.copy()
            pos_x_values[score_list < 0] = np.nan
            ax.fill_between(pos_x_values, score_list, linewidth=0.1,
                            color=color,
                            facecolor=color,
                            alpha=alpha)
            neg_x_values = x_values.copy()
            neg_x_values[score_list >= 0] = np.nan
            ax.fill_between(neg_x_values, score_list, linewidth=0.1,
                            color=negative_color,
                            facecolor=negative_color,
                            alpha=alpha)
