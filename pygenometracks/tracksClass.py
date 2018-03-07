# -*- coding: utf-8 -*-
from __future__ import division
from past.builtins import map
from past.builtins import zip

import sys
import numpy as np
import logging
from configparser import ConfigParser

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.textpath
import matplotlib.colors
import matplotlib.gridspec
import matplotlib.cm
import mpl_toolkits.axisartist as axisartist
from matplotlib.patches import Rectangle
import textwrap
import os.path
from .readBed import ReadBed
from .utilities import to_string, to_bytes

import hicexplorer.HiCMatrix as HiCMatrix
import hicexplorer.utilities
import scipy.sparse
import copy

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ndarray size changed")
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
warnings.simplefilter(action='ignore', category=ImportWarning)

# import warnings
# warnings.filterwarnings('error')

from collections import OrderedDict
from intervaltree import IntervalTree, Interval


FORMAT = "[%(levelname)s:%(filename)s:%(lineno)s - %(funcName)20s()] %(message)s"
logging.basicConfig(format=FORMAT)
log = logging.getLogger(__name__)
# logging.basicConfig()
# log = logging.getLogger("tracksClass")
log.setLevel(logging.DEBUG)

DEFAULT_BED_COLOR = '#1f78b4'
DEFAULT_BIGWIG_COLOR = '#33a02c'
DEFAULT_BEDGRAPH_COLOR = '#a6cee3'
DEFAULT_MATRIX_COLORMAP = 'RdYlBu_r'
DEFAULT_TRACK_HEIGHT = 0.5  # in centimeters
DEFAULT_FIGURE_WIDTH = 40  # in centimeters
# proportion of width dedicated to (figure, legends)
DEFAULT_WIDTH_RATIOS = (0.90, 0.1)
DEFAULT_MARGINS = {'left': 0.04, 'right': 0.92, 'bottom': 0.03, 'top': 0.97}


class MultiDict(OrderedDict):
    """
    Class to allow identically named
    sections in configuration file
    by appending the section number as
    for example:
    1. section name
    """
    _unique = 0

    def __setitem__(self, key, val):
        if isinstance(val, OrderedDict):
            self._unique += 1
            key = "{}. [{}]".format(str(self._unique), key)
        OrderedDict.__setitem__(self, key, val)


class PlotTracks(object):

    def __init__(self, tracks_file, fig_width=DEFAULT_FIGURE_WIDTH,
                 fig_height=None, fontsize=None, dpi=None, track_label_width=None):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.dpi = dpi
        self.vlines_intval_tree = None
        self.vlines_properties = None
        self.track_list = None
        start = self.print_elapsed(None)
        self.available_tracks = self.get_available_tracks()
        self.parse_tracks(tracks_file)
        if fontsize:
            fontsize = fontsize
        else:
            fontsize = float(fig_width) * 0.3
        # the track label width is the fraction of the figure width that is used
        # for the track 'title' or label.
        if track_label_width is None:
            self.width_ratios = DEFAULT_WIDTH_RATIOS
        else:
            self.width_ratios = (1 - track_label_width, track_label_width)

        font = {'size': fontsize}
        matplotlib.rc('font', **font)
        # initialize each track
        self.track_obj_list = []
        for idx, properties in enumerate(self.track_list):
            if 'spacer' in properties:
                self.track_obj_list.append(SpacerTrack(properties))
                continue
            elif 'x-axis' in properties:
                self.track_obj_list.append(XAxisTrack(properties))
                continue
            else:
                # for all other tracks that are not axis or spacer
                # the track_class is obtained from the available tracks
                track_class = self.available_tracks[properties['file_type']]
                self.track_obj_list.append(track_class(properties))

            if 'title' in properties:
                # adjust titles that are too long
                # if the track label space is small
                if track_label_width < 0.1:
                    if sys.version_info[0] == 2:
                        properties['title'] = textwrap.fill(properties['title'].encode("UTF-8"), 12)
                    else:
                        properties['title'] = textwrap.fill(properties['title'], 12)
                else:
                    if sys.version_info[0] == 2:
                        properties['title'] = textwrap.fill(properties['title'].encode("UTF-8"), 30)
                    else:
                        properties['title'] = textwrap.fill(properties['title'], 30)

        log.info("time initializing track(s):")
        self.print_elapsed(start)

    @staticmethod
    def get_available_tracks():
        """
        Finds all the children of the class GenomeTrack (as well as the children of GenomeTrack subclasses
         (eg PlotTADs, which inherits PlotBed that inherits GenomeTrack)
        :return:
        Dict of available GenomeTracks. key = track type (eg. bedgraph, bed) and value = tracks class
        """
        avail_tracks = {}
        work = [GenomeTrack]
        while work:
            parent = work.pop()
            for child in parent.__subclasses__():
                if child not in avail_tracks:
                    track_type = child.TRACK_TYPE
                    avail_tracks[track_type] = child
                    work.append(child)

        return avail_tracks

    def get_tracks_height(self, start_region=None, end_region=None):
        """
        The main purpose of the following loop is
        to get the height of each of the tracks
        because for the Hi-C the height is variable with respect
        to the range being plotted, the function is called
        when each plot is going to be printed.

        Args:
            start_region: start of the region to plot. Only used in case the plot is a Hi-C matrix
            end_region: end of the region to plot. Only used in case the plot is a Hi-C matrix

        Returns:

        """
        track_height = []
        for track_dict in self.track_list:
            # if overlay previous is set to a value other than no
            # then, skip this track height
            if track_dict['overlay previous'] != 'no':
                continue
            elif 'x-axis' in track_dict and track_dict['x-axis'] is True:
                height = track_dict['fontsize'] / 10
            elif 'height' in track_dict:
                height = track_dict['height']
            # compute the height of a Hi-C track
            # based on the depth such that the
            # resulting plot appears proportional
            #
            #      /|\
            #     / | \
            #    /  |d \   d is the depth that we want to be proportional
            #   /   |   \  when plotted in the figure
            # ------------------
            #   region len
            #
            # d (in cm) =  depth (in bp) * width (in cm) / region len (in bp)

            elif 'depth' in track_dict and track_dict['file_type'] == 'hic_matrix':
                # to compute the actual width of the figure the margins and the region
                # set for the legends have to be considered
                # DEFAULT_MARGINS[1] - DEFAULT_MARGINS[0] is the proportion of plotting area

                hic_width = \
                    self.fig_width * (DEFAULT_MARGINS['right'] - DEFAULT_MARGINS['left']) * self.width_ratios[0]
                scale_factor = 0.6  # the scale factor is to obtain a 'pleasing' result.
                depth = min(track_dict['depth'], (end_region - start_region))

                height = scale_factor * depth * hic_width / (end_region - start_region)
            else:
                height = DEFAULT_TRACK_HEIGHT

            track_height.append(height)

        return track_height

    def plot(self, file_name, chrom, start, end, title=None):
        track_height = self.get_tracks_height(start_region=start, end_region=end)

        if self.fig_height:
            fig_height = self.fig_height
        else:
            fig_height = sum(track_height)

        log.debug("Figure size in cm is {} x {}. Dpi is set to {}\n".format(self.fig_width,
                                                                            fig_height, self.dpi))
        fig = plt.figure(figsize=self.cm2inch(self.fig_width, fig_height))
        if title:
            fig.suptitle(title)

        grids = matplotlib.gridspec.GridSpec(len(track_height), 2,
                                             height_ratios=track_height,
                                             width_ratios=self.width_ratios)
        axis_list = []
        # skipped_tracks is the count of tracks that have the
        # 'overlay previous' parameter and should be skipped
        skipped_tracks = 0
        axis = None
        for idx, track in enumerate(self.track_obj_list):
            if idx == 0 and track.properties['overlay previous'] != 'no':
                log.warn("First track can not have the `overlay previous` option")
                track.properties['overlay previous'] = 'no'

            if track.properties['overlay previous'] in ['yes', 'share-y']:
                overlay = True
                skipped_tracks += 1
            else:
                overlay = False

            if track.properties['overlay previous'] == 'share-y':
                ylim = axis.get_ylim()
            else:
                idx -= skipped_tracks
                axis = axisartist.Subplot(fig, grids[idx, 0])
                fig.add_subplot(axis)
                # turns off the lines around the tracks
                axis.axis[:].set_visible(False)
                # to make the background transparent
                axis.patch.set_visible(False)
                label_axis = plt.subplot(grids[idx, 1])
                label_axis.set_axis_off()

            axis.set_xlim(start, end)
            track.plot(axis, label_axis, chrom, start, end)

            if track.properties['overlay previous'] == 'share-y':
                axis.set_ylim(ylim)

            if not overlay:
                axis_list.append(axis)

        if self.vlines_intval_tree:
            self.plot_vlines(axis_list, chrom, start, end)

        fig.subplots_adjust(wspace=0, hspace=0.0,
                            left=DEFAULT_MARGINS['left'],
                            right=DEFAULT_MARGINS['right'],
                            bottom=DEFAULT_MARGINS['bottom'],
                            top=DEFAULT_MARGINS['top'])

        fig.savefig(file_name, dpi=self.dpi, transparent=False)
        return fig.get_size_inches()

    def plot_vlines(self, axis_list, chrom_region, start_region, end_region):
        """
        Plots dotted lines from the top of the first plot to the bottom
        of the last plot at the specified positions.

        :param axis_list: list of plotted axis
        :param chrom_region chromosome name
        :param start_region start position
        :param end_region end position

        :return: None
        """
        vlines_list = []
        if 'line width' in self.vlines_properties:
            line_width = self.vlines_properties['line width']
        else:
            line_width = 0.5

        chrom_region = check_chrom_str_bytes(self.vlines_intval_tree, chrom_region)

        if chrom_region not in list(self.vlines_intval_tree):
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.vlines_intval_tree, chrom_region)

        for region in sorted(self.vlines_intval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            vlines_list.append(region.begin)

        for ax in axis_list:
            ymin, ymax = ax.get_ylim()

            ax.vlines(vlines_list, ymin, ymax, linestyle='dashed', zorder=10, linewidth=line_width,
                      color=(0, 0, 0, 0.7), alpha=0.5)

        return

    def parse_tracks(self, tracks_file):
        """
        Parses a configuration file

        :param tracks_file: file path containing the track configuration
        :return: array of dictionaries and vlines_file. One dictionary per track
        """
        from ast import literal_eval
        parser = ConfigParser(dict_type=MultiDict, strict=False)
        parser.read_file(open(tracks_file, 'r'))

        tracks_file_path = os.path.dirname(tracks_file)

        track_list = []
        for section_name in parser.sections():
            track_options = dict({"section_name": section_name})
            if section_name.endswith('[spacer]'):
                track_options['spacer'] = True
            elif section_name.endswith('[x-axis]'):
                track_options['x-axis'] = True
            for name, value in parser.items(section_name):
                if name in ['max_value', 'min_value', 'depth', 'height', 'line width',
                            'fontsize', 'scale factor', 'number of bins'] and value != 'auto':
                    track_options[name] = literal_eval(value)
                else:
                    track_options[name] = value

            if 'type' in track_options and track_options['type'] == 'vlines':
                self.vlines_properties = self.check_file_exists(track_options, tracks_file_path)
            elif 'skip' in track_options and track_options['skip'] != 'no':
                pass
            else:
                track_list.append(track_options)

        updated_track_list = []
        for track_dict in track_list:
            warn = None
            if 'file' in track_dict and track_dict['file'] != '':
                track_dict = self.check_file_exists(track_dict, tracks_file_path)
                if 'file_type' not in track_dict:
                    track_dict['file_type'] = self.guess_filetype(track_dict, self.available_tracks)

            if 'overlay previous' not in track_dict:
                track_dict['overlay previous'] = 'no'
            #  set some default values
            if 'title' not in track_dict:
                track_dict['title'] = ''
                if track_dict['overlay previous'] != 'no' or track_dict['section_name'].endswith('[x-axis]') \
                        or track_dict['section_name'].endswith('[spacer]'):
                    pass
                else:
                    warn = "\ntitle not set for 'section {}'\n".format(track_dict['section_name'])
            if warn:
                sys.stderr.write(warn)
            updated_track_list.append(track_dict)
        self.track_list = updated_track_list
        if self.vlines_properties:
            self.vlines_intval_tree, __, __ = file_to_intervaltree(self.vlines_properties['file'])

    @staticmethod
    def check_file_exists(track_dict, tracks_path):
        """
        Checks if a file or list of files exists. If the file does not exists
        tries to check if the file may be relative to the track_file path, in such case
        the path is updated.
        :param track_dict: dictionary of track values. Should contain
                            a 'file' key containing the path of the file
                            or files to be checked separated by space
                            For example: file1 file2 file3
        :param tracks_path: path of the tracks file
        :return: None
        """
        for file_type in ['boundaries_file', 'file']:
            if file_type in track_dict:
                file_names = [x for x in track_dict[file_type].split(" ") if x != '']
                full_path_file_names = []
                for file_name in file_names:
                    try:
                        open(file_name, 'r').close()
                        full_path_file_names.append(file_name)
                    except IOError:
                        try:
                            # try to find the file in the same path as the
                            # the path of the
                            name_with_tracks_path = tracks_path + "/" + file_name
                            open(name_with_tracks_path, 'r').close()
                            full_path_file_names.append(name_with_tracks_path)
                        except IOError:
                            sys.stderr.write("\n*ERROR*\nFile in section [{}] "
                                             "not found:\n{}\n\n".format(track_dict['section_name'],
                                                                         file_name))
                            sys.exit(1)

                track_dict[file_type] = " ".join(full_path_file_names)
        return track_dict

    @staticmethod
    def guess_filetype(track_dict, available_tracks):
        """

        :param track_dict: dictionary of track values with the 'file' key
                    containing a string path of the file or files. Only the ending
                     of the last file is used in case when there are more files
        :param: available_tracks: list of available tracks

        :return: string file type detected
        """
        file_ = track_dict['file'].strip()
        file_type = None
        for track_type, track_class in available_tracks.items():
            for ending in track_class.SUPPORTED_ENDINGS:
                if file_.endswith(ending):
                    if file_type == track_class.TRACK_TYPE:
                        log.error("file_type already defined in other GenomeTrack")
                        exit()
                    else:
                        file_type = track_class.TRACK_TYPE

        if file_type is None:
            sys.exit("Section {}: can not identify file type. Please specify "
                     "the file_type for {}".format(track_dict['section_name'], file))

        return file_type

    @staticmethod
    def cm2inch(*tupl):
        inch = 2.54
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    @staticmethod
    def print_elapsed(start):
        import time
        if start:
            log.info(time.time() - start)
        return time.time()


def opener(filename):
    """
    Determines if a file is compressed or not
    """
    import gzip
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
    min_value = np.inf
    max_value = -np.inf

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
            sys.exit(msg)

        try:
            start = int(start)
        except ValueError as detail:
            msg = "Error reading line: {}. The start field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            sys.exit(msg)

        try:
            end = int(end)
        except ValueError as detail:
            msg = "Error reading line: {}. The end field is not " \
                  "an integer.\nError message: {}".format(line_number, detail)
            sys.exit(msg)

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


class GenomeTrack(object):
    """
    The TrackPlot object is a holder for all tracks that are to be plotted.
    For example, to plot a bedgraph file a new class that extends TrackPlot
    should be created.

    It is expected that all GenomeTrack objects have a plot method.

    """
    SUPORTED_ENDINGS = []
    TRACK_TYPE = None
    OPTIONS_TXT = """
# title of track (plotted on the right side)
title =
# height of track in cm (ignored if the track is overlay on top the previous track)
height = 2
# if the track wants to be plotted upside-down:
# orientation = inverted
# if the track wants to be plotted on top of the previous track. Options are 'yes' or 'share-y'. For the 'share-y'
# option the y axis values is shared between this plot and the overlay plot. Otherwise, each plot use its own scale
#overlay previous = yes
"""

    def __init__(self, properties_dict):
        self.properties = properties_dict
        self.file_type = 'test'


class SpacerTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        pass


class XAxisTrack(GenomeTrack):
    SUPPORTED_ENDINGS = []
    TRACK_TYPE = None

    def __init__(self, *args, **kwargs):
        super(XAxisTrack, self).__init__(*args, **kwargs)
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 15

    def plot(self, ax, label_axis, chrom_region, region_start, region_end):
        ticks = ax.get_xticks()
        if ticks[-1] - ticks[1] <= 1e5:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"

        elif 1e5 < ticks[-1] - ticks[1] < 4e6:
            labels = ["{:,.0f}".format((x / 1e3))
                      for x in ticks]
            labels[-2] += " Kb"
        else:
            labels = ["{:,.1f} ".format((x / 1e6))
                      for x in ticks]
            labels[-2] += " Mbp"

        if 'where' in self.properties and self.properties['where'] == 'top':
            ax.axis["x"] = ax.new_floating_axis(0, 0.2)
            ax.axis["x"].set_axis_direction("top")
            label_y_pos = 0.99
            vert_align = 'top'
        else:
            ax.axis["x"] = ax.new_floating_axis(0, 0.9)
            label_y_pos = 0.01
            vert_align = 'bottom'
        ax.text(0.5, label_y_pos, chrom_region, horizontalalignment='center',
                fontsize=int(self.properties['fontsize']), verticalalignment=vert_align, transform=ax.transAxes)

        ax.axis["x"].axis.set_ticklabels(labels)
        ax.axis['x'].axis.set_tick_params(which='minor', bottom='on')

        ax.axis["x"].major_ticklabels.set(size=int(self.properties['fontsize']))
        label_axis.text(0.15, 0.5, self.properties['title'],
                        horizontalalignment='left', size='large',
                        verticalalignment='center')


class BedGraphTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bg', '.bg.gz']
    TRACK_TYPE = 'bedgraph'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
color = green
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, properties_dict):
        # super(self.__class__, self).__init__(*args, **kwargs)
        self.properties = properties_dict
        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BEDGRAPH_COLOR
        self.interval_tree, ymin, ymax = file_to_intervaltree(self.properties['file'])

        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = ymax

        if 'min_value' not in self.properties or self.properties['min_value'] == 'auto':
            self.properties['min_value'] = ymin

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        score_list = []
        pos_list = []

        chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)

        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)

        for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            score_list.append(float(region.data[0]))
            pos_list.append(region.begin + (region.end - region.begin) / 2)

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BEDGRAPH_COLOR

        if 'extra' in self.properties and self.properties['extra'][0] == '4C':
            # draw a vertical line for each fragment region center
            self.ax.fill_between(pos_list, score_list,
                                 facecolor=self.properties['color'],
                                 edgecolor='none')
            self.ax.vlines(pos_list, [0], score_list, color='olive', linewidth=0.5)
            self.ax.plot(pos_list, score_list, '-', color='slateblue', linewidth=0.7)
        else:
            try:
                self.ax.fill_between(pos_list, score_list, facecolor=self.properties['color'])
            except ValueError:
                sys.stderr.write("Invalid color {} for {}. "
                                 "Using gray instead.".format(self.properties['color'], self.properties['file']))
                self.ax.fill_between(pos_list, score_list, facecolor='gray')

        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)

        ymax = self.properties['max_value']
        ymin = self.properties['min_value']

        if float(ymax) % 1 == 0:
            ymax_print = int(ymax)
        else:
            ymax_print = "{:.1f}".format(ymax)
        self.ax.set_ylim(ymin, ymax)
        ydelta = ymax - ymin
        small_x = 0.01 * (end_region - start_region)

        if 'show data range' in self.properties and \
                self.properties['show data range'] == 'no':
            pass
        else:
            # by default show the data range
            self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                         "[{}-{}]".format(ymin, ymax_print),
                         horizontalalignment='left', size='small',
                         verticalalignment='bottom')

        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)


class BigWigTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.bw', '.bigwig']
    TRACK_TYPE = 'bigwig'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
color = #666666
# the default for min_value and max_value is 'auto' which means that the scale will go
# from the minimum value found in the region plotted to the maximum value found.
min_value = 0
#max_value = auto
# The number of bins takes the region to be plotted and divides it into the number of bins specified
# Then, at each bin the bigwig mean value is computed and plotted.
# A lower number of bins produces a coarser tracks
number of bins = 500
# to convert missing data (NaNs) into zeros. Otherwise, missing data is not plotted.
nans to zeros = True
# for type, the options are: line, points, fill. Default is fill
# to add the preferred line width or point size use:
# type = line:lw where lw (linewidth) is float
# similarly points:ms sets the point size (markersize (ms) to the given float
# type = line:0.5
# type = points:0.5
# set show data range to no to hide the text on the upper-left showing the data range
show data range = yes
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        import pyBigWig
        self.bw = pyBigWig.open(self.properties['file'])
        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BIGWIG_COLOR

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        formated_region = "{}:{}-{}".format(chrom_region, start_region, end_region)
        log.debug("plotting {}".format(self.properties['file']))
        # compute the score in bins of 10000 SLOW
    #    score = np.array([self.bw.query(region[0], x, x+10000,1)[0]['mean']
    #                      for x in range(region[1], region[2], 10000)])

        num_bins = 700
        if 'number of bins' in self.properties:
            try:
                num_bins = int(self.properties['number of bins'])
            except TypeError:
                num_bins = 700
                sys.stderr.write("'number of bins' value: {} for bigwig file {} "
                                 "is not valid. Using default value (700)".format(self.properties['number of bins'],
                                                                                  self.properties['file']))
        chrom_region = check_chrom_str_bytes(self.bw.chroms().keys(), chrom_region)

        if chrom_region not in self.bw.chroms().keys():
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.bw.chroms().keys(), chrom_region)

        if chrom_region not in self.bw.chroms().keys():
            sys.stderr.write("Can not read region {} from bigwig file:\n\n"
                             "{}\n\nPlease check that the chromosome name is part of the bigwig file "
                             "and that the region is valid".format(formated_region, self.properties['file']))

        # on rare occasions pyBigWig may throw an error, apparently caused by a corruption
        # of the memory. This only occurs when calling trackPlot from different
        # processors. Reloading the file solves the problem.
        num_tries = 0
        while num_tries < 5:
            num_tries += 1
            try:
                scores_per_bin = np.array(self.bw.stats(chrom_region, start_region,
                                                        end_region, nBins=num_bins)).astype(float)
            except Exception as e:
                import pyBigWig
                self.bw = pyBigWig.open(self.properties['file'])

                log.warning("error found while reading bigwig scores ({}).\nTrying again. Iter num: {}".
                            format(e, num_tries))
                pass
            else:
                if num_tries > 1:
                    log.warning("After {} the scores could be computed".format(num_tries))
                break

        x_values = np.linspace(start_region, end_region, num_bins)

        if 'type' in self.properties and self.properties != 'fill':
            if self.properties['type'].find(":") > 0:
                plot_type, size = self.properties['type'].split(":")
                try:
                    size = float(size)
                except ValueError:
                    exit("Invalid value: 'type = {}' in section: {}\n"
                         "A number was expected and found '{}'".format(self.properties['type'],
                                                                       self.properties['section_name'],
                                                                       size))
            else:
                plot_type = self.properties['type']
                size = None

            if plot_type == 'line':
                self.ax.plot(x_values, scores_per_bin, '-', linewidth=size, color=self.properties['color'])

            elif plot_type == 'points':
                self.ax.plot(x_values, scores_per_bin, '.', markersize=size, color=self.properties['color'])

            else:
                exit("Invalid: 'type = {}' in section: {}\n".format(self.properties['type'],
                                                                    self.properties['section_name'],
                                                                    size))
        else:
            self.ax.fill_between(x_values, scores_per_bin, linewidth=0.1,
                                 color=self.properties['color'],
                                 facecolor=self.properties['color'])

        ymin, ymax = self.ax.get_ylim()
        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            ymax = self.properties['max_value']
        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            ymin = self.properties['min_value']

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.ax.set_ylim(ymax, ymin)
        else:
            self.ax.set_ylim(ymin, ymax)

    #    self.ax.set_yticks([ymax])
        ydelta = ymax - ymin

    #    self.ax.set_yticklabels(["{}-{}".format(int(ymin), int(ymax))], size='large')
        # set min max
        if float(ymax) % 1 == 0:
            ymax_print = int(ymax)
        else:
            ymax_print = "{:.1f}".format(ymax)
        small_x = 0.01 * (end_region - start_region)
        if 'show data range' in self.properties and self.properties['show data range'] == 'no':
            pass
        else:
            # by default show the data range
            self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                         "[{}-{}]".format(int(ymin), ymax_print),
                         horizontalalignment='left',
                         verticalalignment='top')

        """
        self.ax.text(region_end, ymax - ydelta * 0.2, self.properties['title'],
                horizontalalignment='right', size='large',
                verticalalignment='bottom')

        """
        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center')

        return self.ax


class BedTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['bed', 'bed3', 'bed6', 'bed12', 'bed.gz', 'bed3.gz', 'bed6.gz', 'bed12.gz']
    TRACK_TYPE = 'bed'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# if the type=genes is given
# the the file is interpreted as gene
# file. If the bed file contains the exon
# structure (bed 12) then this is plotted. Otherwise
# a region **with direction** is plotted.
# if the bed file contains a column for color (column 9), then this color can be used by
# setting:
# color = bed_rgb
color = darkblue
#if color is a valid colormap name (like RbBlGn), then the score is mapped
# to the colormap. If the color is simply a color name, then this color is used and the score is not considered.
# For the colormap option, the the min_value and max_value for the score can be provided, otherwise
# the maximum score and minimum score found are used.
#color = RdYlBu
#min_value=0
#max_value=100
# height of track in cm
height = 5
# to turn off/on printing of labels
labels = off
# optional: font size can be given to override the default size
fontsize = 10
# optional: line width
#line width = 0.5
# the display parameter defines how the bed file is plotted.
# The options are ['colapsed', 'interleaved', 'triangles'] This options asume that the regions do not overlap.
# `collapsed`: The bed regions are plotted one after the other in one line.
# `interleaved`: The bed regions are plotted in two lines, first up, then down, then up etc.
# if display is not given, then each region is plotted using the gene style
#optional, default is black. To remove the background color, simply set 'color' and 'background color' to the
# same value
#border color = black
# style to plot the genes when they have exon information
#style = UCSC
#style = flybase
# maximum number of gene rows to be plotted. This
# field is useful to limit large number of close genes
# to be printed over many rows. When several images want
# to be combined this must be set to get equal size, otherwise, on each image the height of each gene changes
#gene rows = 10
# if the track wants to be plotted on top of the previous track. Options are 'yes' or 'share-y'. For the 'share-y'
# option the y axis values is shared between this plot and the overlay plot. Otherwise, each plot use its own scale
#overlay previous = yes
# by default the ymax is the number of
# rows occupied by the genes in the region plotted. However,
# by setting this option, the global maximum is used instead.
# This is useful to combine images that are all consistent and
# have the same number of rows.
#global max row = yes
# optional. If not given is guessed from the file ending.
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwarg):
        super(BedTrack, self).__init__(*args, **kwarg)
        self.bed_type = None  # once the bed file is read, this is bed3, bed6 or bed12
        self.len_w = None  # this is the length of the letter 'w' given the font size
        self.interval_tree = {}  # interval tree of the bed regions

        from matplotlib import font_manager
        if 'fontsize' not in self.properties:
            self.properties['fontsize'] = 12
        else:
            self.properties['fontsize'] = float(self.properties['fontsize'])

        self.fp = font_manager.FontProperties(size=self.properties['fontsize'])

        if 'color' not in self.properties:
            self.properties['color'] = DEFAULT_BED_COLOR
        if 'border color' not in self.properties:
            self.properties['border color'] = 'black'
        if 'labels' not in self.properties:
            self.properties['labels'] = 'on'
        if 'style' not in self.properties:
            self.properties['style'] = 'flybase'
        if 'display' not in self.properties:
            self.properties['display'] = 'stacked'
        if 'interval height' not in self.properties:
            self.properties['interval_height'] = 100
        if 'line width' not in self.properties:
            self.properties['line width'] = 0.5

        self.colormap = None

        # check if the color given is a color map
        if not matplotlib.colors.is_color_like(self.properties['color']) and self.properties['color'] != 'bed_rgb':
            # check if the color is a valid colormap name
            if self.properties['color'] not in matplotlib.cm.datad:
                log.warning("*WARNING* color: '{}' for section {} is not valid. Color has "
                            "been set to {}".format(self.properties['color'], self.properties['section_name'],
                                                    DEFAULT_BED_COLOR))
                self.properties['color'] = DEFAULT_BED_COLOR
            else:
                self.colormap = self.properties['color']

        # to set the distance between rows
        self.row_scale = self.properties['interval_height'] * 2.3

        self.interval_tree, min_score, max_score = self.process_bed()
        if self.colormap is not None:
            if 'min_value' in self.properties:
                min_score = self.properties['min_value']
            if 'max_value' in self.properties:
                max_score = self.properties['max_value']

            norm = matplotlib.colors.Normalize(vmin=min_score,
                                               vmax=max_score)

            cmap = matplotlib.cm.get_cmap(self.properties['color'])
            self.colormap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

    def get_length_w(self, fig_width, region_start, region_end):
        '''
        to improve the visualization of the genes it is good to have an estimation of the label
        length. In the following code I try to get the length of a 'W' in base pairs.
        '''
        if self.properties['labels'] == 'on':
            # from http://scipy-cookbook.readthedocs.org/items/Matplotlib_LaTeX_Examples.html
            inches_per_pt = 1.0 / 72.27
            font_in_inches = self.properties['fontsize'] * inches_per_pt
            region_len = region_end - region_start
            bp_per_inch = region_len / fig_width
            font_in_bp = font_in_inches * bp_per_inch
            self.len_w = font_in_bp
            log.debug("len of w set to: {} bp".format(self.len_w))
        else:
            self.len_w = 1

        return self.len_w

    def process_bed(self):

        bed_file_h = ReadBed(opener(self.properties['file']))
        self.bed_type = bed_file_h.file_type

        if 'color' in self.properties and self.properties['color'] == 'bed_rgb' and \
           self.bed_type not in ['bed12', 'bed9']:
            log.warning("*WARNING* Color set to 'bed_rgb', but bed file does not have the rgb field. The color has "
                        "been set to {}".format(DEFAULT_BED_COLOR))
            self.properties['color'] = DEFAULT_BED_COLOR

        valid_intervals = 0
        interval_tree = {}

        max_score = float('-inf')
        min_score = float('inf')
        for bed in bed_file_h:
            if bed.score < min_score:
                min_score = bed.score
            if bed.score > max_score:
                max_score = bed.score

            if bed.chromosome not in interval_tree:
                interval_tree[bed.chromosome] = IntervalTree()

            interval_tree[bed.chromosome].add(Interval(bed.start, bed.end, bed))
            valid_intervals += 1

        if valid_intervals == 0:
            log.warning("No valid intervals were found in file {}".format(self.properties['file_name']))

        return interval_tree, min_score, max_score

    def get_max_num_row(self, len_w, small_relative):
        ''' Process the whole bed regions at the given figure length and font size to
        determine the maximum number of rows required.
        :return:
        '''

        self.max_num_row = {}
        for chrom in self.interval_tree:
            row_last_position = []  # each entry in this list contains the end position
            self.max_num_row[chrom] = 0
            for region in sorted(self.interval_tree[chrom][0:500000000]):
                bed = region.data
                if self.properties['labels'] == 'on':
                    bed_extended_end = int(bed.end + (len(bed.name) * len_w))
                else:
                    bed_extended_end = (bed.end + 2 * small_relative)

                # get smallest free row
                if len(row_last_position) == 0:
                    free_row = 0
                    row_last_position.append(bed_extended_end)
                else:
                    # get list of rows that are less than bed.start, then take the min
                    idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                    if len(idx_list):
                        free_row = min(idx_list)
                        row_last_position[free_row] = bed_extended_end
                    else:
                        free_row = len(row_last_position)
                        row_last_position.append(bed_extended_end)

                if free_row > self.max_num_row[bed.chromosome]:
                    self.max_num_row[bed.chromosome] = free_row

        log.debug("max number of rows set to {}".format(self.max_num_row))
        return self.max_num_row

    def get_y_pos(self, free_row):
        """
        The y_pos is set such that regions to be plotted do not overlap (stacked). To override this
        the properties['collapsed'] needs to be set.

        The algorithm uses a interval tree (self.region_interval) to check the overlaps
        and a sort of coverage vector 'rows used' to identify the row in which to plot
        :return: int y position
        """

        # if the domain directive is given, ypos simply oscilates between 0 and 100
        if self.properties['display'] == 'interlaced':
            ypos = self.properties['interval_height'] if self.counter % 2 == 0 else 1

        elif self.properties['display'] == 'collapsed':
            ypos = 0

        else:
            ypos = free_row * self.row_scale
        return ypos

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.counter = 0
        self.small_relative = 0.004 * (end_region - start_region)
        self.get_length_w(ax.get_figure().get_figwidth(), start_region, end_region)
        if 'global max row' in self.properties and self.properties['global max row'] == 'yes':
            self.get_max_num_row(self.len_w, self.small_relative)

        chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)
        if chrom_region not in self.interval_tree.keys():
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)

        genes_overlap = sorted(self.interval_tree[chrom_region][start_region:end_region])

        # turn labels off when too many intervals are visible.
        if self.properties['labels'] != 'off' and len(genes_overlap) > 60:
            self.properties['labels'] = 'off'

        linewidth = self.properties['line width']
        max_num_row_local = 1
        max_ypos = 0
        # check for the number of other intervals that overlap
        #    with the given interval
        #            1         2
        #  012345678901234567890123456
        #  1=========       4=========
        #       2=========
        #         3============
        #
        # for 1 row_last_position = [9]
        # for 2 row_last_position = [9, 14]
        # for 3 row_last_position = [9, 14, 19]
        # for 4 row_last_position = [26, 14, 19]

        row_last_position = []  # each entry in this list contains the end position
        # of genomic interval. The list index is the row
        # in which the genomic interval was plotted.
        # Any new genomic interval that wants to be plotted,
        # knows the row to use by finding the list index that
        # is larger than its start

        # check for overlapping genes including
        # label size (if plotted)

        for region in genes_overlap:
            """
            BED12 gene format with exon locations at the end
            chrX    20850   23076   CG17636-RA      0       -       20850   23017   0       3       946,765,64,     0,1031,2162,

            BED9
            bed with rgb at end
            chr2L   0       70000   ID_5    0.26864549832   .       0       70000   51,160,44

            BED6
            bed without rgb
            chr2L   0       70000   ID_5    0.26864549832   .
            """
            self.counter += 1
            bed = region.data

            if self.properties['labels'] == 'on':
                num_name_characters = len(bed.name) + 2  # +2 to account for an space before and after the name
                bed_extended_end = int(bed.end + (num_name_characters * self.len_w))
            else:
                bed_extended_end = (bed.end + 2 * self.small_relative)

            # get smallest free row
            if len(row_last_position) == 0:
                free_row = 0
                row_last_position.append(bed_extended_end)
            else:
                # get list of rows that are less than bed.start, then take the min
                idx_list = [idx for idx, value in enumerate(row_last_position) if value < bed.start]
                if len(idx_list):
                    free_row = min(idx_list)
                    row_last_position[free_row] = bed_extended_end
                else:
                    free_row = len(row_last_position)
                    row_last_position.append(bed_extended_end)

            rgb, edgecolor = self.get_rgb_and_edge_color(bed)

            ypos = self.get_y_pos(free_row)

            # do not plot if the maximum interval rows to plot is reached
            if 'gene rows' in self.properties and free_row >= int(self.properties['gene rows']):
                continue

            if free_row > max_num_row_local:
                max_num_row_local = free_row
            if ypos > max_ypos:
                max_ypos = ypos

            if self.bed_type == 'bed12':
                if self.properties['style'] == 'flybase':
                    self.draw_gene_with_introns_flybase_style(ax, bed, ypos, rgb, edgecolor, linewidth)
                else:
                    self.draw_gene_with_introns(ax, bed, ypos, rgb, edgecolor, linewidth)
            else:
                self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)

            if self.properties['labels'] == 'off':
                pass
            elif bed.start > start_region and bed.end < end_region:
                ax.text(bed.end + self.small_relative, ypos + (float(self.properties['interval_height']) / 2),
                        bed.name, horizontalalignment='left',
                        verticalalignment='center', fontproperties=self.fp)

        if self.counter == 0:
            log.warning("*Warning* No intervals were found for file {} "
                        "in section '{}' for the interval plotted ({}:{}-{}).\n".
                        format(self.properties['file'], self.properties['section_name'], chrom_region, start_region, end_region))

        ymax = 0

        if 'global max row' in self.properties and self.properties['global max row'] == 'yes':
            ymin = self.max_num_row[chrom_region] * self.row_scale

        elif 'gene rows' in self.properties:
            ymin = int(self.properties['gene rows']) * self.row_scale
        else:
            ymin = max_ypos + self.properties['interval_height']

        log.debug("ylim {},{}".format(ymin, ymax))
        # the axis is inverted (thus, ymax < ymin)
        ax.set_ylim(ymin, ymax)

        if 'display' in self.properties:
            if self.properties['display'] == 'domain':
                ax.set_ylim(-5, 205)
            elif self.properties['display'] == 'collapsed':
                ax.set_ylim(-5, 105)

        label_ax.text(0.15, 1, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='top', transform=label_ax.transAxes)

    def get_rgb_and_edge_color(self, bed):
        rgb = self.properties['color']
        edgecolor = self.properties['border color']

        if self.colormap:
            # translate value field (in the example above is 0 or 0.2686...) into a color
            rgb = self.colormap.to_rgba(bed.score)

        if self.properties['color'] == 'bed_rgb':
            # if rgb is set in the bed line, this overrides the previously
            # defined colormap
            if self.bed_type in ['bed9', 'bed12'] and len(bed.rgb) == 3:
                try:
                    rgb = [float(x) / 255 for x in bed.rgb]
                    if 'border color' in self.properties:
                        edgecolor = self.properties['border color']
                    else:
                        edgecolor = self.properties['color']
                except IndexError:
                    rgb = DEFAULT_BED_COLOR
            else:
                rgb = DEFAULT_BED_COLOR
        return rgb, edgecolor

    def draw_gene_simple(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws an interval with direction (if given)
        """
        from matplotlib.patches import Polygon

        if bed.strand not in ['+', '-']:
            ax.add_patch(Rectangle((bed.start, ypos), bed.end - bed.start, self.properties['interval_height'],
                                   edgecolor=edgecolor, facecolor=rgb, linewidth=linewidth))
        else:
            vertices = self._draw_arrow(ax, bed.start, bed.end, bed.strand, ypos)
            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=rgb,
                                 linewidth=linewidth))

    def draw_gene_with_introns_flybase_style(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws a gene using different styles
        """
        from matplotlib.patches import Polygon
        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor)
            return
        half_height = float(self.properties['interval_height']) / 2
        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black',
                linewidth=linewidth, zorder=-1)

        # get start, end of all the blocks
        positions = []
        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x0 < bed.thick_start < x1:
                positions.append((x0, bed.thick_start, 'UTR'))
                positions.append((bed.thick_start, x1, 'coding'))

            elif x0 < bed.thick_end < x1:
                positions.append((x0, bed.thick_end, 'coding'))
                positions.append((bed.thick_end, x1, 'UTR'))

            else:
                if x1 < bed.thick_start or x0 > bed.thick_end:
                    type = 'UTR'
                else:
                    type = 'coding'

                positions.append((x0, x1, type))

        # plot all blocks as rectangles except the last if the strand is + or
        # the first is the strand is -, which are drawn as arrows.
        if bed.strand == '-':
            positions = positions[::-1]

        first_pos = positions.pop()
        if first_pos[2] == 'UTR':
            _rgb = 'grey'
        else:
            _rgb = rgb

        vertices = self._draw_arrow(ax, first_pos[0], first_pos[1], bed.strand, ypos)

        ax.add_patch(Polygon(vertices, closed=True, fill=True,
                             edgecolor=edgecolor,
                             facecolor=_rgb,
                             linewidth=linewidth))

        for start_pos, end_pos, _type in positions:
            if _type == 'UTR':
                _rgb = 'grey'
            else:
                _rgb = rgb
            vertices = [(start_pos, ypos), (start_pos, ypos + self.properties['interval_height']),
                        (end_pos, ypos + self.properties['interval_height']), (end_pos, ypos)]

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 edgecolor=edgecolor,
                                 facecolor=_rgb,
                                 linewidth=linewidth))

    def _draw_arrow(self, ax, start, end, strand, ypos):
        """
        Draws a filled arrow
        :param ax:
        :param start:
        :param end:
        :param strand:
        :param ypos:
        :param rgb:
        :return: None
        """
        half_height = float(self.properties['interval_height']) / 2
        if strand == '+':
            x0 = start
            x1 = end  # - self.small_relative
            y0 = ypos
            y1 = ypos + self.properties['interval_height']
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            -----------------\
            ---------------- /

            """

            vertices = [(x0, y0), (x0, y1), (x1, y1), (x1 + self.small_relative, y0 + half_height), (x1, y0)]

        else:
            x0 = start  # + self.small_relative
            x1 = end
            y0 = ypos
            y1 = ypos + self.properties['interval_height']
            """
            The vertices correspond to 5 points along the path of a form like the following,
            starting in the lower left corner and progressing in a clock wise manner.

            /-----------------
            \-----------------
            """
            vertices = [(x0, y0), (x0 - self.small_relative, y0 + half_height), (x0, y1), (x1, y1), (x1, y0)]

        return vertices

    def draw_gene_with_introns(self, ax, bed, ypos, rgb, edgecolor, linewidth):
        """
        draws a gene like in flybase gbrowse.
        """
        from matplotlib.patches import Polygon

        if bed.block_count == 0 and bed.thick_start == bed.start and bed.thick_end == bed.end:
            self.draw_gene_simple(ax, bed, ypos, rgb, edgecolor, linewidth)
            return
        half_height = float(self.properties['interval_height']) / 2
        quarter_height = float(self.properties['interval_height']) / 4
        three_quarter_height = quarter_height * 3

        # draw 'backbone', a line from the start until the end of the gene
        ax.plot([bed.start, bed.end], [ypos + half_height, ypos + half_height], 'black', linewidth=linewidth, zorder=-1)

        for idx in range(0, bed.block_count):
            x0 = bed.start + bed.block_starts[idx]
            x1 = x0 + bed.block_sizes[idx]
            if x1 < bed.thick_start or x0 > bed.thick_end:
                y0 = ypos + quarter_height
                y1 = ypos + three_quarter_height
            else:
                y0 = ypos
                y1 = ypos + self.properties['interval_height']

            if x0 < bed.thick_start < x1:
                vertices = ([(x0, ypos + quarter_height), (x0, ypos + three_quarter_height),
                             (bed.thick_start, ypos + three_quarter_height),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (bed.thick_start, ypos + self.properties['interval_height']),
                             (x1, ypos + self.properties['interval_height']), (x1, ypos),
                             (bed.thick_start, ypos), (bed.thick_start, ypos + quarter_height)])

            elif x0 < bed.thick_end < x1:
                vertices = ([(x0, ypos),
                             (x0, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + self.properties['interval_height']),
                             (bed.thick_end, ypos + three_quarter_height),
                             (x1, ypos + three_quarter_height),
                             (x1, ypos + quarter_height),
                             (bed.thick_end, ypos + quarter_height),
                             (bed.thick_end, ypos)])
            else:
                vertices = ([(x0, y0), (x0, y1), (x1, y1), (x1, y0)])

            ax.add_patch(Polygon(vertices, closed=True, fill=True,
                                 linewidth=linewidth,
                                 edgecolor='none',
                                 facecolor=rgb))

            if idx < bed.block_count - 1:
                # plot small arrows using the character '<' or '>' over the back bone
                intron_length = bed.block_starts[idx + 1] - (bed.block_starts[idx] + bed.block_sizes[idx])
                marker = 5 if bed.strand == '+' else 4
                if intron_length > 3 * self.small_relative:
                    pos = np.arange(x1 + 1 * self.small_relative,
                                    x1 + intron_length + self.small_relative, int(2 * self.small_relative))
                    ax.plot(pos, np.zeros(len(pos)) + ypos + half_height, '.', marker=marker,
                            fillstyle='none', color='blue', markersize=3)

                elif intron_length > self.small_relative:
                    intron_center = x1 + int(intron_length) / 2
                    ax.plot([intron_center], [ypos + half_height], '.', marker=5,
                            fillstyle='none', color='blue', markersize=3)


class LinksTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.arcs', '.arc' '.link', '.links']
    TRACK_TYPE = 'links'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# the file format for links is (tab separated)
#   chr1 start1 end1 chr2 start2 end2 score
# for example:
#   chr1 100 200 chr1 250 300 0.5
# depending on the links type either and arc or a 'triangle' can be plotted. If an arc,
# a line will be drawn from the center of the first region (chr1: 150, tot the center of the other region (chr1:275).
# if a triangle, the vertix of the triangle will be drawn at the center between the two points (also the center of
# each position is used)
# links whose start or end is not in the region plotted are not shown.
# color of the lines
color = red
# for the links type, the options are arcs and triangles, the triangles option is convenient to overlay over a
# Hi-C matrix to highlight the matrix pixel of the highlighted link
links type = arcs
# if line width is not given, the score is used to set the line width
# using the following formula (0.5 * square root(score)
# line width = 0.5
# options for line style are 'solid', 'dashed', 'dotted' etc. The full list of
# styles can be found here: https://matplotlib.org/gallery/lines_bars_and_markers/linestyles.html
line style = solid
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, *args, **kwarg):
        super(LinksTrack, self).__init__(*args, **kwarg)
        # the file format expected is similar to file format of links in
        # circos:
        # chr1 100 200 chr1 250 300 0.5
        # where the last value is a score.
        if 'line width' not in self.properties:
            self.properties['line width'] = 0.5
        if 'line style' not in self.properties:
            self.properties['line style'] = 'solid'
        if 'links type' not in self.properties:
            self.properties['links type'] = 'arcs'
        self.max_height = None
        valid_intervals = 0
        interval_tree = {}
        line_number = 0
        with open(self.properties['file'], 'r') as file_h:
            for line in file_h.readlines():
                line_number += 1
                if line.startswith('browser') or line.startswith('track') or line.startswith('#'):
                    continue
                try:
                    chrom1, start1, end1, chrom2, start2, end2, score = line.strip().split('\t')
                except Exception as detail:
                    msg = 'File not valid. The format is chrom1 start1, end1, ' \
                          'chrom2, start2, end2, score\nError: {}\n in line\n {}'.format(detail, line)
                    sys.exit(msg)

                try:
                    start1 = int(start1)
                    end1 = int(end1)
                    start2 = int(start2)
                    end2 = int(end2)
                except ValueError as detail:
                    msg = "Error reading line: {}. One of the fields is not " \
                          "an integer.\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                assert start1 <= end1, "Error in line #{}, end1 larger than start1 in {}".format(line_number, line)
                assert start2 <= end2, "Error in line #{}, end2 larger than start2 in {}".format(line_number, line)
                try:
                    score = float(score)
                except ValueError as detail:
                    msg = "Error reading line: {}. The score is not valid {}. " \
                          "\nError message: {}".format(line_number, detail)
                    sys.exit(msg)

                if chrom1 != chrom2:
                    sys.stderr.write("Only links in same chromosome are used. Skipping line\n{}\n".format(line))
                    continue

                if chrom1 not in interval_tree:
                    interval_tree[chrom1] = IntervalTree()

                if start2 < start1:
                    start1, start2 = start2, start1
                    end1, end2 = end2, end1

                # each interval spans from the smallest start to the largest end
                interval_tree[chrom1].add(Interval(start1, end2, score))
                valid_intervals += 1

        if valid_intervals == 0:
            sys.stderr.write("No valid intervals were found in file {}".format(self.properties['file']))

        file_h.close()
        self.interval_tree = interval_tree

        if 'color' not in self.properties:
            self.properties['color'] = 'blue'

        if 'alpha' not in self.properties:
            self.properties['alpha'] = 0.8

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        """
        Makes and arc connecting two points on a linear scale representing
        interactions between Hi-C bins.
        :param ax: matplotlib axis
        :param label_ax: matplotlib axis for labels
        """
        self.max_height = 0
        count = 0

        chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)
        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)

        arcs_in_region = sorted(self.interval_tree[chrom_region][region_start:region_end])

        for idx, interval in enumerate(arcs_in_region):
            # skip intervals whose start and end are outside the plotted region
            if interval.begin < region_start and interval.end > region_end:
                continue

            if 'line width' in self.properties:
                self.line_width = float(self.properties['line width'])
            else:
                self.line_width = 0.5 * np.sqrt(interval.data)

            if self.properties['links type'] == 'triangles':
                self.plot_triangles(ax, interval)
            else:
                self.plot_arcs(ax, interval)

            count += 1

        # the arc height is equal to the radius, the track height is the largest
        # radius plotted plus an small increase to avoid cropping of the arcs
        self.max_height += self.max_height * 0.1
        log.debug("{} were links plotted".format(count))
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(self.max_height, -1)
        else:
            ax.set_ylim(-1, self.max_height)

        log.debug('title is {}'.format(self.properties['title']))
        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')

    def plot_arcs(self, ax, interval):
        from matplotlib.patches import Arc

        diameter = (interval.end - interval.begin)
        radius = float(diameter) / 2
        center = interval.begin + float(diameter) / 2
        if radius > self.max_height:
            self.max_height = radius
        ax.plot([center], [diameter])
        ax.add_patch(Arc((center, 0), diameter,
                         diameter, 0, 0, 180, color=self.properties['color'],
                         linewidth=self.line_width, ls=self.properties['line style']))

    def plot_triangles(self, ax, interval):
        from matplotlib.patches import Polygon
        x1 = interval.begin
        x2 = x1 + float(interval.end - interval.begin) / 2
        x3 = interval.end
        y1 = 0
        y2 = (interval.end - interval.begin)

        triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]), closed=False,
                           facecolor='none', edgecolor=self.properties['color'],
                           linewidth=self.line_width,
                           ls=self.properties['line style'])
        ax.add_artist(triangle)
        if y2 > self.max_height:
            self.max_height = y2


class PlotTADs(BedTrack):
    SUPPORTED_ENDINGS = ['.domain', '.domains' '.tad', '.tads']
    TRACK_TYPE = 'domains'

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        """
        Plots the boundaries as triangles in the given ax.
        """
        from matplotlib.patches import Polygon
        ymax = 0.001
        valid_regions = 0
        chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)
        if chrom_region not in self.interval_tree:
            orig = chrom_region
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)
            log.info('Chromosome name: {} does not exists. Changing name to {}'.format(orig, chrom_region))

        for region in sorted(self.interval_tree[chrom_region][start_region:end_region]):
            """      ______ y2
                  /\
                 /  \
                /    \ _____ y1
            _____________________
               x1 x2 x3
            """
            x1 = region.begin
            x2 = x1 + float(region.end - region.begin) / 2
            x3 = region.end
            y1 = 0
            y2 = (region.end - region.begin)

            rgb, edgecolor = self.get_rgb_and_edge_color(region.data)

            triangle = Polygon(np.array([[x1, y1], [x2, y2], [x3, y1]]), closed=True,
                               facecolor=rgb, edgecolor=edgecolor, linewidth=self.properties['line width'])
            ax.add_artist(triangle)
            valid_regions += 1

            if y2 > ymax:
                ymax = y2

        if valid_regions == 0:
            log.warning("No regions found for section {}.".format(self.properties['section_name']))

        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            ax.set_ylim(ymax, 0)
        else:
            ax.set_ylim(0, ymax)

        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')


class HiCMatrixTrack(GenomeTrack):
    SUPPORTED_ENDINGS = ['.h5', '.cool' '.mcool']
    TRACK_TYPE = 'hic_matrix'
    OPTIONS_TXT = """
title =
# The different options for color maps can be found here: https://matplotlib.org/users/colormaps.html
# the default color map is RdYlBu_r (_r) stands for reverse
#colormap = RdYlBu_r
# depth is the maximum distance that should be plotted.
depth = 100000
# height of track (in cm) can be given. Otherwise, the height is computed such that the proportions of the
# hic matrix are kept (e.g. the image does not appear shrink or extended)
# height = 10
# min_value and max_value refer to the contacts in the matrix.
#min_value =2.8
#max_value = 3.0
# the matrix can be transformed using the log1 (or log, but zeros could be problematic)
transform = log1p
# show masked bins plots as white lines
# those bins that were not used during the correction
# the default is to extend neighboring bins to
# obtain an aesthetically pleasant output
show_masked_bins = no
# if the track wants to be plotted upside-down:
# orientation = inverted
# optional if the values in the matrix need to be scaled the
# following parameter can be used. This is useful to plot multiple hic-matrices on the same scale
# scale factor = 1
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, properties_dict):
        # to avoid the color bar to span all the
        # width of the axis I pass two axes
        # to plot_matrix
        self.properties = properties_dict

        log.debug("self.properties", self.properties)
        if self.properties['file'].endswith('.cool'):
            # just init the cooler matrix.
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'], color_only_init=True)
        else:
            self.hic_ma = HiCMatrix.hiCMatrix(self.properties['file'])

        if len(self.hic_ma.matrix.data) == 0:
            log.error("Matrix {} is empty".format(self.properties['file']))
            exit(1)
        if 'show_masked_bins' in self.properties and self.properties['show_masked_bins'] == 'yes':
            pass
        else:
            self.hic_ma.maskBins(self.hic_ma.nan_bins)

        # check that the matrix can be log transformed
        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                if self.hic_ma.matrix.data.min() + 1 < 0:
                    log.error("\n*ERROR*\nMatrix contains negative values.\n"
                              "log1p transformation can not be applied to \n"
                              "values in matrix: {}".format(self.properties['file']))
                    exit(1)

            elif self.properties['transform'] == '-log':
                if self.hic_ma.matrix.data.min() < 0:
                    log.error("\n*ERROR*\nMatrix contains negative values.\n"
                              "log(-1 * <values>) transformation can not be applied to \n"
                              "values in matrix: {}".format(self.properties['file']))
                    exit(1)

            elif self.properties['transform'] == 'log':
                if self.hic_ma.matrix.data.min() < 0:
                    log.error("\n*ERROR*\nMatrix contains negative values.\n"
                              "log transformation can not be applied to \n"
                              "values in matrix: {}".format(self.properties['file']))
                    exit(1)

        new_intervals = hicexplorer.utilities.enlarge_bins(self.hic_ma.cut_intervals)
        self.hic_ma.interval_trees, self.hic_ma.chrBinBoundaries = \
            self.hic_ma.intervalListToIntervalTree(new_intervals)

        self.hic_ma.cut_intervals = new_intervals
        binsize = self.hic_ma.getBinSize()
        max_depth_in_bins = int(self.properties['depth'] / binsize)

        # work only with the lower matrix
        # and remove all pixels that are beyond
        # 2 * max_depth_in_bis which are not required
        # (this is done by subtracting a second sparse matrix
        # that contains only the lower matrix that wants to be removed.
        limit = 2 * max_depth_in_bins
        self.hic_ma.matrix = scipy.sparse.triu(self.hic_ma.matrix, k=0, format='csr') - \
            scipy.sparse.triu(self.hic_ma.matrix, k=limit, format='csr')
        self.hic_ma.matrix.eliminate_zeros()

        # fill the main diagonal, otherwise it looks
        # not so good. The main diagonal is filled
        # with an array containing the max value found
        # in the matrix
        if sum(self.hic_ma.matrix.diagonal()) == 0:
            log.info("Filling main diagonal with max value "
                     "because it empty and looks bad...\n")
            max_value = self.hic_ma.matrix.data.max()
            main_diagonal = scipy.sparse.dia_matrix(([max_value] * self.hic_ma.matrix.shape[0], [0]),
                                                    shape=self.hic_ma.matrix.shape)
            self.hic_ma.matrix = self.hic_ma.matrix + main_diagonal

        self.plot_inverted = False
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            self.plot_inverted = True

        self.norm = None

        if 'colormap' not in self.properties:
            self.properties['colormap'] = DEFAULT_MATRIX_COLORMAP

        self.cmap = matplotlib.cm.get_cmap(self.properties['colormap'])
        self.cmap.set_bad('white')

        self.cmap.set_bad('black')

    def plot(self, ax, label_ax, chrom_region, region_start, region_end):
        self.cbar_ax = copy.copy(label_ax)
        self.label_ax = label_ax
        # self.label_ax = label_ax
        self.ax = ax

        chrom_sizes = self.hic_ma.get_chromosome_sizes()
        chrom_region = check_chrom_str_bytes(chrom_sizes, chrom_region)

        if chrom_region not in list(chrom_sizes):
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(chrom_sizes, chrom_region)

        if region_end > chrom_sizes[chrom_region]:
            log.error("*Error*\nThe region to plot extends beyond the chromosome size. Please check.\n")
            log.error("{} size: {}. Region to plot {}-{}\n".format(chrom_region, chrom_sizes[chrom_region],
                                                                   region_start, region_end))

        if self.properties['file'].endswith('.cool'):
            # load now the region to be plotted
            pass

        # expand region to plus depth on both sides
        # to avoid a 45 degree 'cut' on the edges

        # get bin id of start and end of region in given chromosome
        chr_start_id, chr_end_id = self.hic_ma.getChrBinRange(chrom_region)
        chr_start = self.hic_ma.cut_intervals[chr_start_id][1]
        chr_end = self.hic_ma.cut_intervals[chr_end_id - 1][1]
        start_bp = max(chr_start, region_start - self.properties['depth'])
        end_bp = min(chr_end, region_end + self.properties['depth'])

        idx, start_pos = zip(*[(idx, x[1]) for idx, x in
                               enumerate(self.hic_ma.cut_intervals)
                               if x[0] == chrom_region and x[1] >= start_bp and x[2] <= end_bp])

        idx = idx[0:-1]
        # select only relevant matrix part
        matrix = self.hic_ma.matrix[idx, :][:, idx]
        # limit the 'depth' based on the length of the region being viewed

        region_len = region_end - region_start
        depth = min(self.properties['depth'], int(region_len * 1.25))
        depth_in_bins = int(1.5 * region_len / self.hic_ma.getBinSize())

        if depth < self.properties['depth']:
            # remove from matrix all data points that are not visible.
            matrix = matrix - scipy.sparse.triu(matrix, k=depth_in_bins, format='csr')
        matrix = np.asarray(matrix.todense().astype(float))
        if 'scale factor' in self.properties:
            matrix = matrix * self.properties['scale factor']

        if 'transform' in self.properties:
            if self.properties['transform'] == 'log1p':
                matrix += 1
                self.norm = matplotlib.colors.LogNorm()

            elif self.properties['transform'] == '-log':
                mask = matrix == 0
                matrix[mask] = matrix[mask is False].min()
                matrix = -1 * np.log(matrix)

            elif self.properties['transform'] == 'log':
                mask = matrix == 0
                matrix[mask] = matrix[mask is False].min()
                matrix = np.log(matrix)

        if 'max_value' in self.properties and self.properties['max_value'] != 'auto':
            vmax = self.properties['max_value']

        else:
            # try to use a 'aesthetically pleasant' max value
            vmax = np.percentile(matrix.diagonal(1), 80)

        if 'min_value' in self.properties and self.properties['min_value'] != 'auto':
            vmin = self.properties['min_value']
        else:
            if depth_in_bins > matrix.shape[0]:
                depth_in_bins = matrix.shape[0] - 5

            # if the region length is large with respect to the chromosome length, the diagonal may have
            # very few values or none. Thus, the following lines reduce the number of bins until the
            # diagonal is at least length 5
            num_bins_from_diagonal = int(region_len / self.hic_ma.getBinSize())
            for num_bins in range(0, num_bins_from_diagonal)[::-1]:
                distant_diagonal_values = matrix.diagonal(num_bins)
                if len(distant_diagonal_values) > 5:
                    break

            vmin = np.median(distant_diagonal_values)

        log.info("setting min, max values for track {} to: {}, {}\n".format(self.properties['section_name'],
                                                                            vmin, vmax))
        img = self.pcolormesh_45deg(matrix, start_pos, vmax=vmax, vmin=vmin)
        img.set_rasterized(True)
        if self.plot_inverted:
            self.ax.set_ylim(depth, 0)
        else:
            self.ax.set_ylim(0, depth)

        self.cbar_ax.patch.set_alpha(0.0)
        try:
            if 'transform' in self.properties and \
                    self.properties['transform'] in ['log', 'log1p']:
                # get a useful log scale
                # that looks like [1, 2, 5, 10, 20, 50, 100, ... etc]

                # The following code is problematic with some versions of matplotlib.
                # Should be uncommented once the problem is clarified
                from matplotlib.ticker import LogFormatter
                formatter = LogFormatter(10, labelOnlyBase=False)
                aa = np.array([1, 2, 5])
                tick_values = np.concatenate([aa * 10 ** x for x in range(10)])
                cobar = plt.colorbar(img, ticks=tick_values, format=formatter, ax=self.cbar_ax, fraction=0.95)
            else:
                cobar = plt.colorbar(img, ax=self.cbar_ax, fraction=0.95)
            cobar.solids.set_edgecolor("face")
            cobar.ax.tick_params(labelsize='smaller')
            # cobar.ax.set_ylabel(self.properties['title'])

            # adjust the labels of the colorbar
            labels = cobar.ax.get_yticklabels()
            ticks = cobar.ax.get_yticks()
            if ticks[0] == 0:
                # if the label is at the start of the colobar
                # move it above avoid being cut or overlapping with other track
                labels[0].set_verticalalignment('bottom')
            if ticks[-1] == 1:
                # if the label is at the end of the colobar
                # move it a bit inside to avoid overlapping
                # with other labels
                labels[-1].set_verticalalignment('top')
            cobar.ax.set_yticklabels(labels)

        except ValueError:
            pass

        self.label_ax.text(0.30, 0.5, self.properties['title'],
                           size='large',
                           verticalalignment='center')

    def pcolormesh_45deg(self, matrix_c, start_pos_vector, vmin=None,
                         vmax=None):
        """
        Turns the matrix 45 degrees and adjusts the
        bins to match the actual start end positions.
        """
        import itertools
        # code for rotating the image 45 degrees
        n = matrix_c.shape[0]
        # create rotation/scaling matrix
        t = np.array([[1, 0.5], [-1, 0.5]])
        # create coordinate matrix and transform it
        matrix_a = np.dot(np.array([(i[1], i[0])
                                    for i in itertools.product(start_pos_vector[::-1],
                                                               start_pos_vector)]), t)
        # this is to convert the indices into bp ranges
        x = matrix_a[:, 1].reshape(n + 1, n + 1)
        y = matrix_a[:, 0].reshape(n + 1, n + 1)
        # plot
        im = self.ax.pcolormesh(x, y, np.flipud(matrix_c),
                                vmin=vmin, vmax=vmax, cmap=self.cmap, norm=self.norm)
        return im


class PlotBedGraphMatrix(BedGraphTrack):
    SUPPORTED_ENDINGS = ['.bm', '.bm.gz' '.bedgraphmatrix']
    TRACK_TYPE = 'bedgraph_matrix'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
# a bedgraph matrix file is like a bedgraph, except that per bin there
# are more than one value separated by tab: E.g.
# This file type is produced by HiCExplorer tool hicFindTads and contains
# the TAD-separation score at different window sizes
# chrX	18279	40131	0.399113	0.364118	0.320857	0.274307
# chrX	40132	54262	0.479340	0.425471	0.366541	0.324736
#min_value = 0.10
#max_value = 0.70
# if type is set as lines, then the TAD score lines are drawn instead
# of the matrix otherwise a heatmap is plotted
type = lines
#plot horizontal lines=False
file_type = {}
    """.format(TRACK_TYPE)

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax
        self.label_ax = label_ax
        """
        Plots a bedgraph matrix file, that instead of having
        a single value per bin, it has several values.
        """

        start_pos = []
        matrix_rows = []
        chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)

        # if type(next(iter(self.interval_tree))) is np.bytes_ or type(next(iter(self.interval_tree))) is bytes:
        #     chrom_region = toBytes(chrom_region)
        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)

            # if type(next(iter(self.interval_tree))) is np.bytes_ or type(next(iter(self.interval_tree))) is bytes:
            #     chrom_region = toBytes(chrom_region)

        for region in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            start_pos.append(region.begin)
            values = map(float, region.data)
            matrix_rows.append(values)

        matrix = np.vstack(matrix_rows).T
        if 'orientation' in self.properties and self.properties['orientation'] == 'inverted':
            matrix = np.flipud(matrix)

        if 'type' in self.properties and self.properties['type'] == 'lines':
            for row in matrix:
                self.ax.plot(start_pos, row, color='grey', linewidth=0.5)

            self.ax.plot(start_pos, matrix.mean(axis=0), linestyle="--", marker="|")
            ymax = self.properties['max_value']
            ymin = self.properties['min_value']
            self.ax.set_ylim(ymin, ymax)

            if 'show data range' in self.properties and self.properties['show data range'] == 'no':
                pass
            else:
                if float(ymax) % 1 == 0:
                    ymax_print = int(ymax)
                else:
                    ymax_print = "{:.1f}".format(ymax)

                if float(ymin) % 1 == 0:
                    ymin_print = int(ymin)
                else:
                    ymin_print = "{:.1f}".format(ymin)

                ydelta = ymax - ymin
                small_x = 0.01 * (end_region - start_region)
                # by default show the data range
                self.ax.text(start_region - small_x, ymax - ydelta * 0.2,
                             "[{}-{}]".format(ymin_print, ymax_print),
                             horizontalalignment='left',
                             verticalalignment='bottom')
            if 'plot horizontal lines' in self.properties and self.properties['plot horizontal lines']:
                # plot horizontal lines to compare values
                self.ax.hlines(np.arange(0, 1.1, 0.1), start_region, end_region, linestyle="--",
                               zorder=0, color='grey')

        else:
            x, y = np.meshgrid(start_pos, np.arange(matrix.shape[0]))
            shading = 'gouraud'
            vmax = self.properties['max_value']
            vmin = self.properties['min_value']

            img = self.ax.pcolormesh(x, y, matrix, vmin=vmin, vmax=vmax, shading=shading)
            img.set_rasterized(True)

        self.label_ax.text(0.15, 0.5, self.properties['title'],
                           horizontalalignment='left', size='large',
                           verticalalignment='center', transform=self.label_ax.transAxes)


class NarrowPeakTrack(BedGraphTrack):
    SUPPORTED_ENDINGS = ['.narrowPeak']
    TRACK_TYPE = 'narrow_peal'
    OPTIONS_TXT = GenomeTrack.OPTIONS_TXT + """
#max_value = 0.70
show data range = yes
show labels = yes
# the narrowPeak format provides the information of the
# peak summit. By default this information is used
# although some peaks may look crooked.
use summit = yes
# type of plot: either box or peak
# box will plot a rectangle of the peak width
# peak will plot the shape of the peak, whose height is the
# narrowPeak file signal value (usually peak coverage)
type = peak
# if the peaks look too thin, the can be adjusted
width adjust = 1.5
file_type = {}
    """.format(TRACK_TYPE)

    def __init__(self, properties_dict):
        self.properties = properties_dict
        if 'color' not in self.properties:
            self.properties['color'] = '#FF000080'  # red, alpha=0.55
        if 'show data range' not in self.properties:
            self.properties['show data range'] = 'yes'
        if 'show labels' not in self.properties:
            self.properties['show labels'] = 'yes'
        if 'use summits' not in self.properties:
            self.properties['use summits'] = 'yes'
        if 'width adjust' not in self.properties:
            self.properties['width adjust'] = 1.5
        else:
            self.properties['width adjust'] = float(self.properties['width adjust'])
        if 'type' not in self.properties:
            self.properties['type'] = 'peak'

        self.interval_tree, _, _ = file_to_intervaltree(self.properties['file'])

    def peak_plot(self, start, end, height, center=None, width_adjust=1.5):
        # uses bezier curves to plot a shape that
        # looks like a peak
        from matplotlib.path import Path
        import matplotlib.patches as patches
        peak_width = float(end - start)
        if center is None:
            center = peak_width / 2 + start
        if width_adjust != 1:
            start -= width_adjust * peak_width / 2
            end += width_adjust * peak_width / 2
            peak_width *= width_adjust

        path_data = [
            (Path.MOVETO, (start, 0)),
            (Path.CURVE4, (start + peak_width / 2, 0)),
            (Path.CURVE4, (start + peak_width * 0.4, height)),
            (Path.CURVE4, (center, height)),
            (Path.CURVE4, (end - peak_width * 0.4, height)),
            (Path.CURVE4, (end - peak_width / 2, 0)),
            (Path.CURVE4, (end, 0))]

        codes, verts = zip(*path_data)
        path = Path(verts, codes)
        patch = patches.PathPatch(path, facecolor=self.properties['color'])
        self.ax.add_patch(patch)

    def plot(self, ax, label_ax, chrom_region, start_region, end_region):
        self.ax = ax

        chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)
        if chrom_region not in list(self.interval_tree):
            chrom_region = change_chrom_names(chrom_region)
            chrom_region = check_chrom_str_bytes(self.interval_tree, chrom_region)
        max_signal = -1
        for peak in sorted(self.interval_tree[chrom_region][start_region - 10000:end_region + 10000]):
            name, score, strand, signal_value, p_value, q_value, summit = peak.data
            signal_value = float(signal_value)
            p_value = float(p_value)
            q_value = float(q_value)
            summit = int(summit)
            if summit > 0:
                summit = peak.begin + summit
            else:
                summit = None

            if self.properties['type'] == 'box':
                self.ax.add_patch(Rectangle((peak.begin, 0), peak.end - peak.begin,
                                            100, edgecolor='black',
                                            facecolor=self.properties['color']))
                max_signal = 110
            else:
                if signal_value > max_signal:
                    max_signal = signal_value
                self.peak_plot(peak.begin, peak.end, signal_value, center=summit,
                               width_adjust=self.properties['width adjust'])

            x_pos = peak.begin + float(peak.end - peak.begin) / 2
            y_pos = 0 - max_signal * 0.05
            if self.properties['show labels'] != 'no':
                self.ax.text(x_pos, y_pos, "{}\np-val:{:.1f}\nq-val:{:.1f}".format(name, p_value, q_value),
                             horizontalalignment='center', size='smaller',
                             verticalalignment='top')

        if 'max_value' not in self.properties or self.properties['max_value'] == 'auto':
            self.properties['max_value'] = max_signal

        ymax = self.properties['max_value']
        if self.properties['show labels'] != 'no':
            if self.properties['type'] == 'box':
                ymin = ymax * -3
            else:
                ymin = ymax * -0.8
        else:
            ymin = 0
        self.ax.set_ylim(ymin, ymax)
        if self.properties['show data range'] != 'no' and self.properties['type'] != 'box':
            if float(ymax) % 1 == 0:
                ymax_print = int(ymax)
            else:
                ymax_print = "{:.1f}".format(ymax)
            ymin_print = 0
            small_x = 0.01 * (end_region - start_region)
            # by default show the data range
            self.ax.text(start_region - small_x, ymax - ymax * 0.01,
                         "[{}-{}]".format(ymin_print, ymax_print),
                         horizontalalignment='left',
                         verticalalignment='top')

        label_ax.text(0.15, 0.5, self.properties['title'],
                      horizontalalignment='left', size='large',
                      verticalalignment='center')


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


def check_chrom_str_bytes(iteratable_obj, p_obj):
    # determine type
    if isinstance(p_obj, list) and len(p_obj) > 0:
        type_ = type(p_obj[0])
    else:
        type_ = type(p_obj)
    if not isinstance(type(next(iter(iteratable_obj))), type_):
        if type(next(iter(iteratable_obj))) is str:
            p_obj = to_string(p_obj)
        elif type(next(iter(iteratable_obj))) in [bytes, np.bytes_]:
            p_obj = to_bytes(p_obj)
    return p_obj
