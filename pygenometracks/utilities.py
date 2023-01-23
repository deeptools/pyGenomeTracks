import sys
import os
import gzip
import numpy as np
from tqdm import tqdm
from intervaltree import IntervalTree, Interval
import pybedtools
import tempfile
import warnings
import logging
from matplotlib.ticker import Formatter
import math


FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
logging.basicConfig(format=FORMAT)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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
        assert sys.version_info[0] != 2
#        if sys.version_info[0] == 2:
#            return str(s)
        return s.decode('ascii')
    if isinstance(s, list):
        return [to_string(x) for x in s]
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


def temp_file_from_intersect(file_name, plot_regions=None, around_region=0):
    """
    intersect file_name with the plot_regions +/- around_region
    :param file_name: string file name
    :param plot_regions:a list of tuple [(chrom1, start1, end1), (chrom2, start2, end2)]
                        with the region to restrict the data to.
    :param around_region: integer with the bp to extend to plot_regions
    :return: temporary file with the intersection
    """
    file_to_open = file_name
    # Check if we can restrict the interval tree to a region:
    if plot_regions is not None:
        # We use pybedtools to overlap:
        original_file = pybedtools.BedTool(file_name)
        # We extend the start and end:
        plot_regions_ext = [(chrom, max(0, start - around_region), end + around_region) for chrom, start, end in plot_regions]
        # We will overlap with both version of chromosome name:
        plot_regions_as_bed = '\n'.join([f'{chrom}\t{start}\t{end}\n{change_chrom_names(chrom)}\t{start}\t{end}' for chrom, start, end in plot_regions_ext])
        regions = pybedtools.BedTool(plot_regions_as_bed, from_string=True)
        # Bedtools will put a warning because we are using inconsistent
        # nomenclature (with and without chr)
        temporary_file = tempfile.NamedTemporaryFile(delete=False)
        sys.stderr = open(temporary_file.name, 'w')
        try:
            file_to_open = original_file.intersect(regions, wa=True, u=True).fn
        except pybedtools.helpers.BEDToolsError as e:
            log.warning(f"BEDTools intersect raised: {e}"
                        "\nWill not subset the file.")
            file_to_open = file_name
        except NotImplementedError:
            log.warning("BEDTools is not installed pygenometracks"
                        " will be slower.")
            file_to_open = file_name
        except Exception as e:
            log.warning(f"BEDTools intersect raised: {e}"
                        "\nWill not subset the file.")
            file_to_open = file_name
        sys.stderr.close()
        sys.stderr = sys.__stderr__
        with open(temporary_file.name, 'r') as f:
            temp_std_error = f.readlines()
        os.remove(temporary_file.name)
        error_lines = [line for line in temp_std_error if 'error' in line.lower()]
        if len(error_lines) > 0:
            error_lines_printable = '\n'.join(error_lines)
            log.warning("BEDTools intersect raised an error:\n"
                        f"{error_lines_printable}\n"
                        "Will not use BEDTools.\n")
            file_to_open = file_name
    return file_to_open


def file_to_intervaltree(file_name, plot_regions=None):
    """
    converts a BED like file into a bx python interval tree
    :param file_name: string file name
    :param plot_regions:a list of tuple [(chrom1, start1, end1), (chrom2, start2, end2)]
                        with the region to restrict the data to.
    :return: interval tree dictionary. They key is the chromosome/contig name and the
    value is an IntervalTree. Each of the intervals have as 'value' the fields[3:] if any.
    """
    file_to_open = temp_file_from_intersect(file_name, plot_regions, 0)
    # iterate over a BED like file
    # saving the data into an interval tree
    # for quick retrieval
    file_h = opener(file_to_open)
    line_number = 0
    valid_intervals = 0
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
            msg = f"Error reading line: {line_number}\nError message: {detail}"
            if len(fields) == 1:
                if fields[0].startswith("{\\rtf"):
                    raise InputError(f"The file {file_name} is a rtf file."
                                     " Please save it as plain text.")
                else:
                    raise InputError(f"Only one field detected in {file_name}"
                                     ", you may use"
                                     " a bed-like delimited by space. This format "
                                     "is not supported by pyGenomeTracks.")
            raise InputError(msg)

        try:
            start = int(start)
        except ValueError as detail:
            msg = f"Error reading line: {line_number}. The start field is not " \
                  f"an integer.\nError message: {detail}"
            raise InputError(msg)

        try:
            end = int(end)
        except ValueError as detail:
            msg = f"Error reading line: {line_number}. The end field is not " \
                  f"an integer.\nError message: {detail}"
            raise InputError(msg)

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

        assert end > start, f"Start position larger or equal than end for line\n{line} "

        interval_tree[chrom].add(Interval(start, end, value))
        valid_intervals += 1

    if valid_intervals == 0:
        if file_to_open != file_name:
            suffix = " after intersection with the plotted region"
        else:
            suffix = ""
        log.warning(f"No valid intervals were found in file {file_name}{suffix}")
    file_h.close()

    return interval_tree, min_value, max_value


def plot_coverage(ax, x_values, score_list, plot_type, size, color,
                  negative_color, alpha, grid):
    if grid:
        ax.grid(axis='y', zorder=0)
    if plot_type == 'line':
        if color == negative_color:
            ax.plot(x_values, score_list, '-', linewidth=size, color=color,
                    alpha=alpha)
        else:
            warnings.warn('Line plots with a different negative color might not look pretty.\n')
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
                          '(fill, line, points) will be fill.\n')
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
        return score_list
    elif transform in ['log', 'log2', 'log10']:
        if np.nanmin(score_list) <= - log_pseudocount:
            msg = ("\n*ERROR*\ncoverage contains values smaller or equal to"
                   f" - {log_pseudocount}.\n"
                   f"{transform}({log_pseudocount} + <values>) transformation "
                   "can not be applied to "
                   f"values in file: {file}")
            raise Exception(msg)
        else:
            return eval('np.' + transform + '(log_pseudocount + score_list)')
    elif transform == 'log1p':
        if np.nanmin(score_list) <= - 1:
            msg = ("\n*ERROR*\ncoverage contains values below or equal to - 1.\n"
                   "log1p(<values>) transformation can not be applied to "
                   f"values in file: {file}")
            raise Exception(msg)
        else:
            return np.log1p(score_list)
    elif transform == '-log':
        if np.nanmin(score_list) <= - log_pseudocount:
            msg = ("\n*ERROR*\ncoverage contains values smaller or equal to"
                   f" - {log_pseudocount}.\n"
                   f"- log( {log_pseudocount} + <values>) transformation can "
                   f"not be applied to values in file: {file}")
            raise Exception(msg)
        else:
            return - np.log(log_pseudocount + score_list)
    else:
        warnings.warn(f"The transform: {transform} for file {file} is not "
                      "valid. Will not use any transformation.\n")
        return score_list


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


def get_optimal_fontsize(fig_width, region_start, region_end):
    """
    to improve the visualization of the letters (one per base)
    it is good to have an estimation of the fontsize.
    """
    # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
    inches_per_pt = 1.0 / 72.27
    region_len = region_end - region_start
    bp_per_inch = region_len / fig_width
    fontsize = 1 / (inches_per_pt * bp_per_inch)
    return fontsize


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
    return n


def change_chrom_names(chrom):
    """
    Changes UCSC chromosome names to ensembl chromosome names
    and vice versa.
    """
    # TODO: mapping from chromosome names like mithocondria is missing
    if chrom.startswith('chr'):
        # remove the chr part from chromosome name
        chrom = chrom[3:]
    else:
        # prefix with 'chr' the chromosome name
        chrom = 'chr' + chrom

    return chrom


def get_region(region_string):
    """
    splits a region string into
    a chrom, start_region, end_region tuple
    The region_string format is chr:start-end
    """
    if region_string:
        # separate the chromosome name and the location using the ':' character
        try:
            chrom, position = region_string.strip().split(":")
        except ValueError:
            raise InputError(f"The region provided ({region_string})"
                             " is not valid, it should be chr:start-end.\n")

        # clean up the position
        for char in ",.;|!{}()":
            position = position.replace(char, '')

        position_list = position.split("-")
        assert len(position_list) == 2, \
            f"The region provided ({region_string})" \
            " is not valid, it should be chr:start-end.\n"

        try:
            region_start = int(position_list[0])
        except ValueError:
            raise InputError(f"The start value ({position_list[0]}) in the"
                             " region provided"
                             " is not valid, it should be chr:start-end.\n")
        try:
            region_end = int(position_list[1])
        except ValueError:
            raise InputError(f"The start value ({position_list[0]}) in the"
                             " region provided"
                             " is not valid, it should be chr:start-end.\n")

        if region_end <= region_start:
            raise InputError("Please check that the region end is larger "
                             "than the region start.\n"
                             f"Values given:\nstart: {region_start}\n"
                             f"end: {region_end}\n"
                             "To plot tracks with a decreasing axis "
                             "consider using `--decreasingXAxis`.")

        return chrom, region_start, region_end


class MyBasePairFormatter(Formatter):
    """
    Format tick values as pretty numbers and add as offset the unit

    The units are "b", "Kb", "Mb"
    The choice is made based on distance between extreme visible locs

    """

    def __init__(self):
        self.format = ''
        self.exponent = 0
        self.unit = ""

    def __call__(self, x, pos=None):
        """
        Return the format for tick value *x* at position *pos*.
        """
        if len(self.locs) == 0:
            # This should never happen...
            return ''
        else:
            xp = (x) / (10. ** self.exponent)
            if abs(xp) < 1e-8:
                xp = 0
            if len(self.locs) < 2 or x == self.locs[-2]:
                return self.format.format(xp) + ' ' + self.unit
            else:
                return self.format.format(xp)

    def set_locs(self, locs):
        # docstring inherited
        self.locs = locs
        if len(self.locs) > 0:
            self._set_unit()
            self._set_format()

    def _set_unit(self):
        # restrict to visible ticks
        vmin, vmax = sorted(self.axis.get_view_interval())
        locs = np.asarray(self.locs)
        locs = locs[(vmin <= locs) & (locs <= vmax)]
        locs = np.abs(locs)
        if not len(locs):
            # I don't understand how this can happen
            self.exponent = 0
            self.unit = ""
            return
        else:
            if np.abs(locs[-1] - locs[0]) <= 1e3:
                self.exponent = 0
                self.unit = "b"
            elif np.abs(locs[-1] - locs[0]) <= 7e5:
                self.exponent = 3
                self.unit = "Kb"
            else:
                self.exponent = 6
                self.unit = "Mb"

    # This is adapted from ScalarFormatter
    def _set_format(self):
        # set the format string to format all the ticklabels
        if len(self.locs) < 2:
            # Temporarily augment the locations with the axis end points.
            _locs = [*self.locs, *self.axis.get_view_interval()]
        else:
            _locs = self.locs
        locs = np.asarray(_locs) / 10. ** self.exponent
        loc_range = np.ptp(locs)
        # Curvilinear coordinates can yield two identical points.
        if loc_range == 0:
            loc_range = np.max(np.abs(locs))
        # Both points might be zero.
        if loc_range == 0:
            loc_range = 1
        if len(self.locs) < 2:
            # We needed the end points only for the loc_range calculation.
            locs = locs[:-2]
        loc_range_oom = int(math.floor(math.log10(loc_range)))
        # first estimate:
        sigfigs = max(0, 3 - loc_range_oom)
        # refined estimate:
        thresh = 1e-3 * 10 ** loc_range_oom
        while sigfigs >= 0:
            if np.abs(locs - np.round(locs, decimals=sigfigs)).max() < thresh:
                sigfigs -= 1
            else:
                break
        sigfigs += 1
        self.format = '{:,.' + str(sigfigs) + 'f}'
