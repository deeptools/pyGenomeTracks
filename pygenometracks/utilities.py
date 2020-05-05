import sys
import gzip
import numpy as np
from tqdm import tqdm
from intervaltree import IntervalTree, Interval
import warnings


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

    for line in tqdm(file_h.readlines()):
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
            neg_x_values[score_list > 0] = np.nan
            ax.fill_between(neg_x_values, score_list, linewidth=0.1,
                            color=negative_color,
                            facecolor=negative_color,
                            alpha=alpha)


def transform(score_list, transform, log_pseudocount, file):
    if transform == 'no':
        return(score_list)
    elif transform in ['log', 'log2', 'log10']:
        if np.nanmin(score_list) <= - log_pseudocount:
            msg = ("\n*ERROR*\ncoverage contains values smaller or equal to"
                   " - {0}.\n"
                   "{1}({0} + <values>) transformation can not be applied to "
                   "values in file: {2}".format(log_pseudocount, transform,
                                                file))
            raise Exception(msg)
        else:
            return(eval('np.' + transform + '(log_pseudocount + score_list)'))
    elif transform == 'log1p':
        if np.nanmin(score_list) <= - 1:
            msg = ("\n*ERROR*\ncoverage contains values below or equal to - 1.\n"
                   "log1p(<values>) transformation can not be applied to "
                   "values in file: {}".format(file))
            raise Exception(msg)
        else:
            return(np.log1p(score_list))
    elif transform == '-log':
        if np.nanmax(score_list.max) <= - log_pseudocount:
            msg = ("\n*ERROR*\ncoverage contains values smaller or equal to"
                   " - {0}.\n"
                   "- log( {0} + <values>) transformation can not be applied"
                   " to values in file: {1}".format(log_pseudocount, file))
            raise Exception(msg)
        else:
            return(- np.log(log_pseudocount + score_list))
    else:
        warnings.warn('The transform: {} for file {} is not valid.'
                      'will not use any transformation'.format(transform,
                                                               file))
        return(score_list)


def get_length_w(fig_width, region_start, region_end, fontsize):
    """
    to improve the visualization of the labels
    it is good to have an estimation of their length
    in base pairs. In the following code I try to get the
    length of a 'W' in base pairs.
    """
    # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
    inches_per_pt = 1.0 / 72.27
    font_in_inches = fontsize * inches_per_pt
    region_len = region_end - region_start
    bp_per_inch = region_len / fig_width
    font_in_bp = font_in_inches * bp_per_inch
    return font_in_bp


def count_lines(file_h, asBed=False):
    n = 0
    for line in file_h:
        if asBed:
            line = to_string(line)
            if line.startswith("#") or line.startswith("track") or \
               line.startswith("browser") or line.strip() == '':
                continue
        n += 1
    file_h.close()
    return(n)
